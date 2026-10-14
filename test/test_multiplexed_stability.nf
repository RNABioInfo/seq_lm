#!/usr/bin/env nextflow

nextflow.enable.types = true

include { DifferentialExpressionResult ; Sample } from '../lib/sample.nf'
include { resolve_gene_contrasts } from '../lib/gene_contrasts.nf'
include { differential_stability ; read_stability_tsv ; stability_config } from '../lib/stability.nf'

params {
    test_scenario: String = 'shared'
}

def write_multiplexed_snapshot(root: Path, samples: List<Sample>, index: Integer, retry: Boolean) -> Path {
    def snapshot: Path = root.resolve("snapshot_${index}")
    java.nio.file.Files.createDirectories(snapshot)
    java.nio.file.Files.writeString(
        snapshot.resolve('sample_metadata.tsv'),
        'sample\tgroup\n' + samples.collect { sample -> "${sample.name}\t${sample.group}\n" }.join(''),
    )
    java.nio.file.Files.writeString(snapshot.resolve('feature_counts.tsv'), 'feature_id\ta1\ngene1\t100\ngene2\t100\n')
    ['stress_a', 'stress_b'].withIndex().each { group: String, contrast_index: Integer ->
        def contrast: Path = snapshot.resolve("contrast_${contrast_index + 1}")
        def log_fc: Number = group == 'stress_a' ? 2.0 : retry && index >= 5 ? 5.0 : 2.0 + Math.min(index, 2)
        java.nio.file.Files.createDirectories(contrast)
        java.nio.file.Files.writeString(
            contrast.resolve('edgeR_results.tsv'),
            "feature_id\tlogFC\tFDR\ngene1\t${log_fc}\t0.01\ngene2\t0\t0.5\n",
        )
    }
    java.nio.file.Files.writeString(snapshot.resolve('contrasts.tsv'),
        'contrast_id\ttarget_group\treference_group\ncontrast_1\tstress_a\tcontrol\ncontrast_2\tstress_b\tcontrol\n')
    return snapshot
}

workflow {
    def scenario: String = params.test_scenario
    assert scenario in ['shared', 'separate', 'excluded', 'retry', 'log', 'marker_failure']
    def root: Path = java.nio.file.Files.createTempDirectory('seq-lm-multiplexed-stability-')
    def excluded: Boolean = scenario == 'excluded'
    def retry: Boolean = scenario == 'retry'
    def behavior: String = scenario == 'log' ? 'log' : 'terminate'
    def samples: List<Sample> = [
        record(name: 'c1', group: 'control', order: null, bam_dir: root.resolve('c1'), is_live: false, protocol_run_id: 'shared-run'),
        record(name: 'c2', group: 'control', order: null, bam_dir: root.resolve('c2'), is_live: false, protocol_run_id: 'shared-run'),
        record(name: 'a1', group: 'stress_a', order: null, bam_dir: root.resolve('a1'), is_live: true, protocol_run_id: scenario == 'separate' ? 'run-a' : 'shared-run'),
        record(name: 'a2', group: 'stress_a', order: null, bam_dir: root.resolve('a2'), is_live: true, protocol_run_id: scenario == 'separate' ? 'run-a' : 'shared-run'),
        record(name: 'b1', group: 'stress_b', order: null, bam_dir: root.resolve('b1'), is_live: !excluded, protocol_run_id: scenario == 'separate' ? 'run-b' : 'shared-run'),
        record(name: 'b_restored', group: 'stress_b', order: null, bam_dir: root.resolve('b_restored'), is_live: true, protocol_run_id: 'shared-run'),
        record(name: 'b_closed', group: 'stress_b', order: null, bam_dir: root.resolve('b_closed'), is_live: true, protocol_run_id: 'shared-run'),
    ]
    samples.each { sample -> java.nio.file.Files.createDirectories(sample.bam_dir) }
    java.nio.file.Files.writeString(root.resolve('b_closed/STOP.manual'), '')
    if (scenario == 'marker_failure') {
        java.nio.file.Files.delete(root.resolve('a2'))
    }
    def active_samples: List<Sample> = samples.findAll { sample -> sample.name != 'b_restored' }
    def final_index: Integer = retry ? 8 : excluded ? 3 : 5
    def results: List<DifferentialExpressionResult> = (0..final_index).collect { index: Integer ->
        record(
            batch_index: index,
            analysis_index: index,
            report_sequence: index,
            results: write_multiplexed_snapshot(root, samples, index, retry),
        )
    }
    def settings: Map = [
        behavior: behavior,
        num_stable_batches: 2,
        max_feature_diff_fraction: 0.05,
        max_median_abs_lfc_delta: 0.1,
        min_jaccard_similarity: 0.9,
        max_fdr: 0.05,
        min_abs_lfc: 1.0,
        min_de_calls_for_fraction_metrics: 20,
        max_small_set_call_changes: 2,
    ]
    settings.contrasts = resolve_gene_contrasts(samples, null)
    settings.config = stability_config(settings, samples)
    def fake_log: Path = root.resolve('calls.tsv')
    def fake_state: Path = root.resolve('first-request')
    def fake_manager: Path = root.resolve('seq-run-manager')
    java.nio.file.Files.writeString(
        fake_manager,
        """#!/bin/sh
            run_id=''
            previous=''
            for argument in \"\$@\"; do
                if [ \"\$previous\" = '--run-id' ]; then run_id=\"\$argument\"; fi
                previous=\"\$argument\"
            done
            printf '%s\\n' \"\$run_id\" >> '${fake_log}'
            if [ '${scenario}' = 'retry' ] && [ ! -f '${fake_state}' ]; then
                : > '${fake_state}'
                exit 2
            fi
            exit 0
        """.stripIndent(),
    )
    assert fake_manager.toFile().setExecutable(true)
    stability_audits = differential_stability(
        channel.fromList(results),
        samples,
        active_samples,
        behavior,
        settings,
        [previous_results: file("${projectDir}/../data/OPTIONAL_FILE"), streaks: [:], eligible: [:]],
        [command: fake_manager.toString(), host: 'fake.local', port: 9501, client_certificate: null, client_private_key: null, ca_certificate: null],
    )
    stability_audits.collect().filter { _audits -> scenario != 'marker_failure' }.map { audits ->
        audits = audits.toSorted { left, right -> left.analysis_index <=> right.analysis_index }
        def rows_by_batch: Map = audits.collectEntries { audit ->
            [(audit.batch_index): read_stability_tsv(audit.sample_stability).collectEntries { row: Map -> [(row.sample): row] }]
        }
        assert rows_by_batch.size() == final_index + 1
        assert rows_by_batch[0].a1.action_result == 'none'
        assert rows_by_batch[1].a1.action_result == 'none'
        def calls: List<String> = java.nio.file.Files.exists(fake_log) ? java.nio.file.Files.readAllLines(fake_log) : []
        if (excluded || scenario == 'separate') {
            assert rows_by_batch[2].a1.action_result == 'stop_created'
            assert rows_by_batch[2].a2.action_result == 'stop_created'
            assert rows_by_batch[2].a1.run_eligible == 'true'
            assert rows_by_batch[2].a1.run_blocking_samples == ''
        }
        else {
            [2, 3].each { index: Integer ->
                assert rows_by_batch[index].a1.eligible == 'true'
                assert rows_by_batch[index].a1.run_eligible == 'false'
                assert rows_by_batch[index].a1.run_blocking_samples == 'stress_b/b1'
                assert rows_by_batch[index].a1.action_result == 'waiting_for_run_stability'
                assert rows_by_batch[index].b1.action_result == 'none'
            }
        }
        if (scenario == 'log') {
            assert calls.empty
            assert rows_by_batch[4].a1.action_result == 'logged'
            assert rows_by_batch[4].b1.action_result == 'logged'
            assert rows_by_batch[4].a1.run_eligible == 'true'
            assert samples.every { sample -> !java.nio.file.Files.exists(sample.bam_dir.resolve('STOP')) }
        }
        else {
            def stopped_names: List<String> = excluded ? ['a1', 'a2'] : ['a1', 'a2', 'b1']
            assert stopped_names.every { name: String -> java.nio.file.Files.exists(root.resolve("${name}/STOP")) }
            assert stopped_names.every { name: String -> rows_by_batch[final_index][name].action_result == 'stop_exists' }
            if (retry) {
                assert calls == ['shared-run', 'shared-run']
                assert stopped_names.every { name: String -> rows_by_batch[4][name].action_result == 'stop_failed' }
                assert rows_by_batch[5].a1.action_result == 'waiting_for_run_stability'
                assert rows_by_batch[6].a1.action_result == 'waiting_for_run_stability'
                assert stopped_names.every { name: String -> rows_by_batch[7][name].action_result == 'stop_created' }
            }
            else if (scenario == 'separate') {
                assert calls == ['run-a', 'run-b']
                assert rows_by_batch[2].b1.action_result == 'none'
                assert rows_by_batch[4].b1.action_result == 'stop_created'
            }
            else {
                assert calls == ['shared-run']
                if (!excluded) {
                    assert stopped_names.every { name: String -> rows_by_batch[4][name].action_result == 'stop_created' }
                }
            }
        }
        assert rows_by_batch[2].b_closed.action_result == 'stop_exists'
        assert rows_by_batch[2].b_restored.effectively_live == 'false'
        assert rows_by_batch[2].c1.effectively_live == 'false'
        assert rows_by_batch[2].b_restored.action_result == 'none'
        assert !java.nio.file.Files.exists(root.resolve('b_restored/STOP'))
        assert !java.nio.file.Files.exists(root.resolve('c1/STOP'))
        "Multiplexed stability scenario passed: ${scenario}"
    }.view()
}
