#!/usr/bin/env nextflow

nextflow.enable.types = true

include { DifferentialExpressionResult ; Sample } from '../lib/sample.nf'
include { differential_stability ; discover_stability_state ; stability_config } from '../lib/stability.nf'

params {
    test_behavior: String = 'log'
    test_reset: Boolean = false
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
    def samples: List<Sample> = [
        record(name: 'c1', group: 'control', order: null, bam_dir: root.resolve('c1'), is_live: true),
        record(name: 'c2', group: 'control', order: null, bam_dir: root.resolve('c2'), is_live: true),
        record(name: 't1', group: 'treated', order: null, bam_dir: root.resolve('t1'), is_live: true),
        record(name: 't2', group: 'treated', order: null, bam_dir: root.resolve('t2'), is_live: true),
    ]
    samples.each { sample -> java.nio.file.Files.createDirectories(sample.bam_dir) }
    def final_index: Integer = params.test_reset ? 3 : 2
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
        max_call_churn_fraction: 0.1,
        max_lost_call_fraction: 0.1,
        max_fdr: 0.05,
        min_abs_lfc: 1.0,
        min_de_calls_for_fraction_metrics: 20,
        max_small_set_call_changes: 2,
    ]
    settings.config = stability_config(settings, samples)
    differential_stability(
        channel.fromList(results),
        samples,
        samples,
        params.test_behavior,
        settings,
        [previous_results: file("${projectDir}/../data/OPTIONAL_FILE"), streaks: [:], eligible: [:]],
    )
    differential_stability.out
        .collect()
        .map { audits ->
            assert audits.size() == final_index + 1
            def rows: List<String> = java.nio.file.Files.readAllLines(audits[-1].sample_stability)
            def expected_eligible: String = params.test_reset ? 'false' : 'true'
            assert rows.drop(1).every { row: String -> row.split('\t', -1)[7] == expected_eligible }
            if (params.test_reset) {
                def contrast_rows: List<String> = java.nio.file.Files.readAllLines(audits[-1].contrast_stability)
                assert contrast_rows[1].split('\t', -1)[20] == '1'
                return 'DE stability resets and rebuilds a contrast streak after instability'
            }
            if (params.test_behavior == 'terminate') {
                assert samples.every { sample -> java.nio.file.Files.isRegularFile(sample.bam_dir.resolve('STOP')) }
                assert rows.drop(1).every { row: String -> row.split('\t', -1)[10] == 'stop_created' }
            }
            else {
                assert samples.every { sample -> !java.nio.file.Files.exists(sample.bam_dir.resolve('STOP')) }
                assert rows.drop(1).every { row: String -> row.split('\t', -1)[10] == 'logged' }
            }
            java.nio.file.Files.createDirectories(
                output_root.resolve("differential_expression/batch_${final_index}")
            )
            def restored: Map = discover_stability_state(
                output_root,
                final_index + 1,
                settings.config,
            )
            assert restored.streaks == ['group_treated_vs_control': 2]
            assert restored.eligible.values().every { value -> value }
            'DE stability reaches sample eligibility after two stable comparisons'
        }
        .view()
}
