#!/usr/bin/env nextflow

nextflow.enable.types = true

include { finalize_qc_report } from '../modules/qc_report_helpers.nf'

workflow {
    def root: Path = java.nio.file.Files.createTempDirectory('seq-lm-final-report-')
    def report: Path = root.resolve('qc_report.html')
    def snapshot: Path = root.resolve('qc_report_snapshot_4.html')
    def state: Path = root.resolve('qc_report_state.json')
    java.nio.file.Files.writeString(
        report,
        '<html><body>Waiting for update.</body></html>',
    )
    java.nio.file.Files.writeString(
        snapshot,
        '<html><body>Final report</body></html>',
    )
    java.nio.file.Files.writeString(
        state,
        '{"latest_batch":"4","snapshot":"qc_report_snapshot_4.html"}\n',
    )

    finalize_qc_report(root)

    def finalized: String = java.nio.file.Files.readString(report)
    assert finalized.contains('Final report')
    assert !finalized.contains('Waiting for update.')
    assert java.nio.file.Files.isRegularFile(snapshot)
    assert java.nio.file.Files.isRegularFile(state)
    channel.of('final QC report replaces the live shell').view()
}
