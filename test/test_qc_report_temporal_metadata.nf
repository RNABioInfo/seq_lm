#!/usr/bin/env nextflow

nextflow.enable.types = true

include { ChunkQCResult ; Sample } from '../lib/sample.nf'
include { parse_sample_order ; validate_timeline_samples } from '../lib/bam_ingress.nf'
include { qc_report_inputs_from_state } from '../modules/qc_report_helpers.nf'

def assert_temporal_error(action: Closure, expected_message: String) -> Void {
    try {
        action.call()
    }
    catch (exception: Exception) {
        def observed_message: String = exception.message ?: exception.cause?.message ?: "${exception}"
        assert observed_message.contains(expected_message):
            "Expected temporal validation error containing '${expected_message}', observed '${observed_message}'."
        return null
    }
    assert false: "Expected temporal validation to fail with: ${expected_message}"
}

workflow {
    def root: Path = java.nio.file.Files.createTempDirectory('seq-lm-temporal-report-')
    def bam_dir: Path = root.resolve('bams')
    java.nio.file.Files.createDirectories(bam_dir)

    def early_sample: Sample = record(
        name: 'early_rep',
        group: 'early',
        order: -5,
        bam_dir: bam_dir,
        is_live: false,
        protocol_run_id: null,
    )
    def late_sample: Sample = record(
        name: 'late_rep',
        group: 'late',
        order: 45,
        bam_dir: bam_dir,
        is_live: false,
        protocol_run_id: null,
    )
    def early_result: ChunkQCResult = record(
        batch_index: 2,
        sample: early_sample,
        bam: root.resolve('early.bam'),
        nanoplot_data: root.resolve('early.tsv.gz'),
        flagstat: root.resolve('early.tsv'),
    )
    def late_result: ChunkQCResult = record(
        batch_index: 2,
        sample: late_sample,
        bam: root.resolve('late.bam'),
        nanoplot_data: root.resolve('late.tsv.gz'),
        flagstat: root.resolve('late.tsv'),
    )
    def state: Map<String,List<ChunkQCResult>> = [
        early: [early_result],
        late: [late_result],
    ]

    def report_inputs: Map = qc_report_inputs_from_state(2, state)
    assert report_inputs.report_inputs_list*.sample*.order.flatten() == [-5, 45]
    assert report_inputs.rows == [
        'early_rep\tearly\t-5\t1\t2\tqc_results',
        'late_rep\tlate\t45\t1\t2\tqc_results',
    ].join('\n')

    assert parse_sample_order('-5', 'early_rep') == -5
    assert parse_sample_order('+45', 'late_rep') == 45
    assert parse_sample_order('', 'missing') == null
    validate_timeline_samples([
        [name: 'early_a', group: 'early', order: -5],
        [name: 'early_b', group: 'early', order: -5],
        [name: 'late_a', group: 'late', order: 45],
    ])
    assert_temporal_error(
        { -> parse_sample_order('1.5', 'fractional') },
        'Expected elapsed minutes as a signed integer',
    )
    assert_temporal_error(
        { -> validate_timeline_samples([
            [name: 'early', group: 'same', order: 0],
            [name: 'late', group: 'same', order: 15],
        ]) },
        'every group to use one elapsed minute',
    )
    assert_temporal_error(
        { -> validate_timeline_samples([
            [name: 'first', group: 'first', order: 0],
            [name: 'second', group: 'second', order: 0],
        ]) },
        'every elapsed minute to identify one group',
    )
    assert_temporal_error(
        { -> validate_timeline_samples([
            [name: 'only', group: 'only', order: 0],
        ]) },
        'at least two distinct elapsed minutes',
    )
    channel.of('QC report temporal metadata passed').view()
}
