#!/usr/bin/env nextflow
nextflow.enable.types = true
params {
    terminate: Boolean = false
}

include { bam_ingress } from '../lib/bam_ingress.nf'

workflow {
    def root: Path = java.nio.file.Files.createTempDirectory('seq-lm-completion-mode-')
    def run: Path = root.resolve('20261012_1200_device_flowcell_12345678')
    def second: Path = root.resolve('second')
    java.nio.file.Files.createDirectories(run)
    java.nio.file.Files.createDirectories(second)
    java.nio.file.Files.writeString(run.resolve('startup.bam'), 'startup')
    java.nio.file.Files.writeString(second.resolve('startup.bam'), 'startup')
    java.nio.file.Files.writeString(run.resolve('final_summary.txt'), 'protocol_run_id=12345678-abcd-1234-abcd-123456789abc\nacquisition_run_id=abcdef123456\nacquisition_stopped=2026-10-12T12:00:00Z\nprocessing_stopped=2026-10-12T12:01:00Z\n')
    def sheet: Path = root.resolve('samples.csv')
    java.nio.file.Files.writeString(sheet, "alias,group,bam_dir,is_live\nfirst,control,${run},true\nsecond,control,${second},false\n")
    nextflow.Global.session.addIgniter {
        Thread.startDaemon('completion-mode-writer') {
            Thread.sleep(1600)
            java.nio.file.Files.writeString(run.resolve('late.bam'), 'complete')
            Thread.sleep(400)
            java.nio.file.Files.writeString(run.resolve('STOP.manual'), '')
        }
    }
    def args = record(live_analysis: true, timeline_analysis: false, sample_sheet_path: sheet,
        bam_poll_interval_ms: 100, bam_stability_polls: 3, termination_requested: params.terminate == true,
        detect_minknow_completion: params.terminate == true, minknow_completion_settle_seconds: 1)
    bam_ingress(args).collect(sort: false).map { batches ->
        assert batches*.batch_index == [0, 1]
        assert batches[1].chunks[0].bam_paths[0].fileName.toString() == 'late.bam'
        assert java.nio.file.Files.exists(run.resolve('STOP.manual'))
        'Opt-out/automated termination preserves manual STOP and accepts BAM after final summary settling time'
    }.view()
}
