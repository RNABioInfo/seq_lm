#!/usr/bin/env nextflow

nextflow.enable.types = true

include { publish_qc_report_snapshot ; next_qc_report_revision ; finalize_qc_report } from '../modules/qc_report_helpers.nf'

params {
    stale: Boolean = false
}

workflow {
    def root: Path = java.nio.file.Files.createTempDirectory('seq-lm-qc-report-publication-')
    (0..1).each { index: Integer ->
        def source: Path = java.nio.file.Files.createDirectories(root.resolve("source_${index}"))
        def snapshot_name: String = "qc_report_snapshot_revision_${index}.html"
        java.nio.file.Files.writeString(source.resolve('qc_report.html'), "shell ${index}")
        java.nio.file.Files.writeString(source.resolve(snapshot_name), "snapshot ${index}")
        java.nio.file.Files.writeString(
            source.resolve('qc_report_state.json'),
            new groovy.json.JsonBuilder([report_revision: index, latest_batch: "0", snapshot: snapshot_name, snapshot_bytes: "snapshot ${index}".bytes.length]).toString(),
        )
        publish_qc_report_snapshot(
            [source.resolve('qc_report_state.json'), source.resolve(snapshot_name), source.resolve('qc_report.html')],
            root.resolve('output'),
        )
    }
    def output: Path = root.resolve('output')
    def state: Map = new groovy.json.JsonSlurper().parse(output.resolve('qc_report_state.json').toFile()) as Map
    assert state.latest_batch == '0'
    assert state.report_revision == 1
    assert next_qc_report_revision(output) == 2
    assert java.nio.file.Files.readString(output.resolve('qc_report.html')) == 'shell 1'
    assert java.nio.file.Files.readString(output.resolve('qc_report_snapshot_revision_0.html')) == 'snapshot 0'
    assert java.nio.file.Files.readString(output.resolve('qc_report_snapshot_revision_1.html')) == 'snapshot 1'
    assert java.nio.file.Files.list(output).withCloseable { paths ->
        paths.noneMatch { path: Path -> path.fileName.toString().endsWith('.pending') }
    }
    if (params.stale) {
        def stale_source: Path = root.resolve('source_0')
        publish_qc_report_snapshot(
            [stale_source.resolve('qc_report_state.json'), stale_source.resolve('qc_report_snapshot_revision_0.html'), stale_source.resolve('qc_report.html')],
            output,
        )
    }
    finalize_qc_report(output)
    assert java.nio.file.Files.readString(output.resolve('qc_report.html')) == 'snapshot 1'
    // A snapshot left by an interrupted publication must never be reused.
    java.nio.file.Files.writeString(output.resolve('qc_report_snapshot_revision_5.html'), 'orphan')
    assert next_qc_report_revision(output) == 6
    channel.of('QC report publication passed').view()
}
