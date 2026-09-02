#!/usr/bin/env nextflow

nextflow.enable.types = true

include {
    bam_ingress ;
    get_samples ;
    inspect_minknow_sample_sheet ;
    protocol_run_id_for_sample
} from '../lib/bam_ingress.nf'

workflow {
    def root: Path = java.nio.file.Files.createTempDirectory('seq-lm-minknow-metadata-')
    def run_dir: Path = root.resolve('run')
    def first_bam_dir: Path = run_dir.resolve('first-bams')
    def second_bam_dir: Path = run_dir.resolve('second-bams')
    java.nio.file.Files.createDirectories(first_bam_dir)
    java.nio.file.Files.createDirectories(second_bam_dir)
    java.nio.file.Files.writeString(first_bam_dir.resolve('first.bam'), 'bam')
    java.nio.file.Files.writeString(second_bam_dir.resolve('second.bam'), 'bam')
    java.nio.file.Files.writeString(
        run_dir.resolve('20260901_sample_sheet_run.csv'),
        '\uFEFFprotocol_run_id,position_id,flow_cell_id,sample_id\n run-first ,P1,F1, first \nrun-second,P2,F2,second\n',
    )
    def workflow_sheet: Path = root.resolve('workflow_samples.csv')
    java.nio.file.Files.writeString(
        workflow_sheet,
        "alias,group,bam_dir,is_live\nfirst,control,${first_bam_dir},false\nsecond,control,${second_bam_dir},false\n",
    )
    def ingress_args = record(
        live_analysis: false,
        timeline_analysis: false,
        sample_sheet_path: workflow_sheet,
        bam_poll_interval_ms: 100,
        bam_stability_polls: 3,
        termination_requested: true,
    )

    def missing_dir: Path = root.resolve('missing-sheet')
    java.nio.file.Files.createDirectories(missing_dir)
    assert inspect_minknow_sample_sheet(missing_dir).reason.contains('no file matching')
    def missing_first_bams: Path = missing_dir.resolve('first-bams')
    def missing_second_bams: Path = missing_dir.resolve('second-bams')
    java.nio.file.Files.createDirectories(missing_first_bams)
    java.nio.file.Files.createDirectories(missing_second_bams)
    def missing_workflow_sheet: Path = root.resolve('missing_workflow_samples.csv')
    java.nio.file.Files.writeString(
        missing_workflow_sheet,
        "alias,group,bam_dir\nmissing_first,control,${missing_first_bams}\nmissing_second,control,${missing_second_bams}\n",
    )
    def missing_samples: List = get_samples(record(
        live_analysis: true,
        timeline_analysis: false,
        sample_sheet_path: missing_workflow_sheet,
        bam_poll_interval_ms: 100,
        bam_stability_polls: 3,
        termination_requested: true,
    ))
    assert missing_samples*.protocol_run_id == [null, null]

    def nested_dir: Path = root.resolve('nested-only')
    java.nio.file.Files.createDirectories(nested_dir.resolve('child'))
    java.nio.file.Files.writeString(
        nested_dir.resolve('child/sample_sheet_nested.csv'),
        'sample_id,protocol_run_id\nfirst,nested-run\n',
    )
    assert inspect_minknow_sample_sheet(nested_dir).reason.contains('no file matching')

    def multiple_dir: Path = root.resolve('multiple')
    java.nio.file.Files.createDirectories(multiple_dir)
    java.nio.file.Files.writeString(multiple_dir.resolve('sample_sheet_a.csv'), 'sample_id,protocol_run_id\nfirst,a\n')
    java.nio.file.Files.writeString(multiple_dir.resolve('x_sample_sheet_b.csv'), 'sample_id,protocol_run_id\nfirst,b\n')
    assert inspect_minknow_sample_sheet(multiple_dir).reason.contains('multiple files')

    def missing_fields_dir: Path = root.resolve('missing-fields')
    java.nio.file.Files.createDirectories(missing_fields_dir)
    java.nio.file.Files.writeString(missing_fields_dir.resolve('sample_sheet.csv'), 'sample_id\nfirst\n')
    assert inspect_minknow_sample_sheet(missing_fields_dir).reason.contains('missing required fields')

    def malformed_dir: Path = root.resolve('malformed')
    java.nio.file.Files.createDirectories(malformed_dir)
    java.nio.file.Files.writeString(
        malformed_dir.resolve('sample_sheet.csv'),
        'sample_id,protocol_run_id\n"first,run-first\n',
    )
    assert inspect_minknow_sample_sheet(malformed_dir).reason != null

    def row_cases_dir: Path = root.resolve('row-cases')
    java.nio.file.Files.createDirectories(row_cases_dir)
    java.nio.file.Files.writeString(
        row_cases_dir.resolve('sample_sheet.csv'),
        'sample_id,protocol_run_id\nfirst,\nduplicate,one\nduplicate,two\n',
    )
    def row_cases: Map = inspect_minknow_sample_sheet(row_cases_dir)
    assert protocol_run_id_for_sample([name: 'missing'], row_cases).reason.contains('no sample_id')
    assert protocol_run_id_for_sample([name: 'duplicate'], row_cases).reason.contains('multiple sample_id')
    assert protocol_run_id_for_sample([name: 'first'], row_cases).reason.contains('blank protocol_run_id')

    bam_ingress(ingress_args)
        .map { batch ->
            assert batch.batch_index == 0
            assert batch.chunks*.sample*.protocol_run_id == ['run-first', 'run-second']
            assert batch.chunks*.sample*.name == ['first', 'second']
            'MinKNOW protocol run IDs are discovered and propagated with the sample batch'
        }
        .view()
}
