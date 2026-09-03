#!/usr/bin/env nextflow

nextflow.enable.types = true

include { DifferentialExpressionResult ; Sample } from '../lib/sample.nf'
include { differential_stability ; discover_stability_state ; stability_config } from '../lib/stability.nf'

params {
    test_behavior: String = 'log'
    test_reset: Boolean = false
    test_stop_mode: String = 'success'
}

def write_snapshot(root: Path, analysis_index: Integer, log_fc: Number) -> Path {
    def snapshot: Path = root.resolve("snapshot_${analysis_index}")
    def contrast: Path = snapshot.resolve('group_treated_vs_control')
    java.nio.file.Files.createDirectories(contrast)
    java.nio.file.Files.writeString(
        snapshot.resolve('sample_metadata.tsv'),
        'sample\tgroup\nc1\tcontrol\nc2\tcontrol\nt1\ttreated\nt2\ttreated\n',
    )
    java.nio.file.Files.writeString(
        snapshot.resolve('feature_counts.tsv'),
        'feature_id\tc1\nc1\t1\nc2\t2\n',
    )
    java.nio.file.Files.writeString(
        contrast.resolve('edgeR_results.tsv'),
        "feature_id\tlogFC\tFDR\nc1\t${log_fc}\t0.01\nc2\t0\t0.5\n",
    )
    return snapshot
}

workflow {
    def root: Path = java.nio.file.Files.createTempDirectory('seq-lm-stability-')
    def output_root: Path = Path.of(params.out_dir as String).toAbsolutePath().normalize()
    def include_run_ids: Boolean = params.test_stop_mode != 'missing'
    def samples: List<Sample> = [
        record(name: 'c1', group: 'control', order: null, bam_dir: root.resolve('c1'), is_live: true, protocol_run_id: include_run_ids ? 'run-c1' : null),
        record(name: 'c2', group: 'control', order: null, bam_dir: root.resolve('c2'), is_live: true, protocol_run_id: include_run_ids ? 'run-c2' : null),
        record(name: 't1', group: 'treated', order: null, bam_dir: root.resolve('t1'), is_live: true, protocol_run_id: include_run_ids ? 'run-t1' : null),
        record(name: 't2', group: 'treated', order: null, bam_dir: root.resolve('t2'), is_live: true, protocol_run_id: include_run_ids ? 'run-t2' : null),
    ]
    samples.each { sample -> java.nio.file.Files.createDirectories(sample.bam_dir) }
    def final_index: Integer = params.test_stop_mode == 'retry' ? 4 : params.test_reset ? 3 : 2
    def results: List<DifferentialExpressionResult> = (0..final_index).collect { analysis_index: Integer ->
        def log_fc: Number = params.test_reset && analysis_index >= 2 ? 3.0 : 2.0
        record(
            batch_index: analysis_index,
            analysis_index: analysis_index,
            report_sequence: analysis_index,
            results: write_snapshot(root, analysis_index, log_fc),
        )
    }
    def settings: Map = [
        behavior: params.test_behavior as String,
        num_stable_batches: 2,
        max_feature_diff_fraction: 0.05,
        max_median_abs_lfc_delta: 0.1,
        min_jaccard_similarity: 0.9,
        max_fdr: 0.05,
        min_abs_lfc: 1.0,
        min_de_calls_for_fraction_metrics: 20,
        max_small_set_call_changes: 2,
    ]
    settings.config = stability_config(settings, samples)
    def fake_log: Path = root.resolve('seq-run-manager-calls.tsv')
    def fake_state: Path = root.resolve('seq-run-manager-state')
    def fake_manager: Path = root.resolve('seq-run-manager')
    java.nio.file.Files.createDirectories(fake_state)
    java.nio.file.Files.writeString(
        fake_manager,
        """#!/bin/sh
            run_id=''
            previous=''
            for argument in \"\$@\"; do
                if [ \"\$previous\" = '--run-id' ]; then run_id=\"\$argument\"; fi
                previous=\"\$argument\"
            done
            printf '%s\\t%s\\n' \"\$run_id\" \"\$*\" >> '${fake_log}'
            if [ '${params.test_stop_mode}' = 'retry' ] && [ ! -f '${fake_state}/'\"\$run_id\" ]; then
                : > '${fake_state}/'\"\$run_id\"
                exit 2
            fi
            exit 0
        """.stripIndent(),
    )
    assert fake_manager.toFile().setExecutable(true)
    def credential_paths: List<Path> = ['client.pem', 'key.pem', 'ca.crt'].collect { name: String ->
        def path: Path = root.resolve(name)
        java.nio.file.Files.writeString(path, '')
        path
    }
    def minknow_connection: Map = [
        command: fake_manager.toString(),
        host: 'minknow.local',
        port: 9501,
        client_certificate: credential_paths[0],
        client_private_key: credential_paths[1],
        ca_certificate: credential_paths[2],
    ]
    differential_stability(
        channel.fromList(results),
        samples,
        samples,
        params.test_behavior,
        settings,
        [previous_results: file("${projectDir}/../data/OPTIONAL_FILE"), streaks: [:], eligible: [:]],
        minknow_connection,
    )
    differential_stability.out
        .collect()
        .map { audits ->
            assert audits.size() == final_index + 1
            def rows: List<String> = java.nio.file.Files.readAllLines(audits[-1].sample_stability)
            def expected_eligible: String = params.test_reset ? 'false' : 'true'
            assert rows.drop(1).every { row: String -> row.split('\t', -1)[8] == (params.test_reset ? '1' : params.test_stop_mode == 'retry' ? '4' : '2') }
            assert rows.drop(1).every { row: String -> row.split('\t', -1)[9] == expected_eligible }
            if (params.test_reset) {
                def contrast_rows: List<String> = java.nio.file.Files.readAllLines(audits[-1].contrast_stability)
                assert contrast_rows[1].split('\t', -1)[20] == '1'
                return 'DE stability resets and rebuilds a contrast streak after instability'
            }
            if (params.test_behavior == 'terminate') {
                if (params.test_stop_mode == 'missing') {
                    assert samples.every { sample -> !java.nio.file.Files.exists(sample.bam_dir.resolve('STOP')) }
                    assert rows.drop(1).every { row: String -> row.split('\t', -1)[12] == 'termination_disabled_no_run_id' }
                    assert !java.nio.file.Files.exists(fake_log)
                }
                else {
                    assert samples.every { sample -> java.nio.file.Files.isRegularFile(sample.bam_dir.resolve('STOP')) }
                    def expected_final_action: String = params.test_stop_mode == 'retry' ? 'stop_exists' : 'stop_created'
                    assert rows.drop(1).every { row: String -> row.split('\t', -1)[12] == expected_final_action }
                    def calls: List<String> = java.nio.file.Files.readAllLines(fake_log)
                    def expected_call_count: Integer = params.test_stop_mode == 'retry' ? 8 : 4
                    assert calls.size() == expected_call_count
                    assert calls.every { call: String ->
                        call.contains('stop --host minknow.local --port 9501') &&
                            call.contains('--client-certificate-path') &&
                            call.contains('--client-private-key-path') &&
                            call.contains('--ca-certificate-path') &&
                            call.contains('--run-id')
                    }
                    if (params.test_stop_mode == 'retry') {
                        def failed_rows: List<String> = java.nio.file.Files.readAllLines(audits[2].sample_stability)
                        def successful_rows: List<String> = java.nio.file.Files.readAllLines(audits[3].sample_stability)
                        assert failed_rows.drop(1).every { row: String -> row.split('\t', -1)[12] == 'stop_failed' }
                        assert successful_rows.drop(1).every { row: String -> row.split('\t', -1)[12] == 'stop_created' }
                    }
                }
            }
            else {
                assert samples.every { sample -> !java.nio.file.Files.exists(sample.bam_dir.resolve('STOP')) }
                assert rows.drop(1).every { row: String -> row.split('\t', -1)[12] == 'logged' }
            }
            java.nio.file.Files.createDirectories(
                output_root.resolve("differential_expression/batch_${final_index}")
            )
            def restored: Map = discover_stability_state(
                output_root,
                final_index + 1,
                settings.config,
            )
            def expected_restored_streak: Integer = params.test_stop_mode == 'retry' ? 4 : 2
            assert restored.streaks == ['group_treated_vs_control': expected_restored_streak]
            assert restored.eligible.values().every { value -> value }
            'DE stability reaches sample eligibility after two stable comparisons'
        }
        .view()
}
