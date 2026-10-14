#!/usr/bin/env nextflow
nextflow.enable.types = true
include { bam_ingress } from '../lib/bam_ingress.nf'

workflow {
    def root: Path = java.nio.file.Files.createTempDirectory('seq-lm-completion-ingress-')
    def run: Path = root.resolve('20261012_1200_device_flowcell_12345678')
    def first: Path = run.resolve('bam_pass/barcode01')
    def second: Path = run.resolve('bam_pass/barcode02')
    def other: Path = root.resolve('custom-run/bams')
    java.nio.file.Files.createDirectories(first)
    java.nio.file.Files.createDirectories(second)
    java.nio.file.Files.createDirectories(other)
    java.nio.file.Files.writeString(first.resolve('startup.bam'), 'startup')
    java.nio.file.Files.writeString(second.resolve('startup.bam'), 'startup')
    java.nio.file.Files.writeString(other.resolve('startup.bam'), 'startup')
    def valid: String = 'protocol_run_id=12345678-abcd-1234-abcd-123456789abc\nacquisition_run_id=abcdef123456\nacquisition_stopped=2026-10-12T12:00:00Z\nprocessing_stopped=2026-10-12T12:01:00Z\n'
    java.nio.file.Files.writeString(run.resolve('final_summary.txt'), valid)
    def sheet: Path = root.resolve('samples.csv')
    java.nio.file.Files.writeString(sheet, "alias,group,bam_dir,is_live,minknow_run_dir\nfirst,control,${first},true,\nsecond,control,${second},true,\nother,treated,${other},true,${other.parent}\n")
    // Keep writes tied to source activation rather than workflow compilation time.
    nextflow.Global.session.addIgniter {
        Thread.startDaemon('completion-test-writer') {
            Thread.sleep(400)
            java.nio.file.Files.writeString(first.resolve('late.bam'), 'partial')
            java.nio.file.Files.writeString(second.resolve('late.bam'), 'complete')
            java.nio.file.Files.writeString(other.resolve('late.bam'), 'complete')
            Thread.sleep(400)
            java.nio.file.Files.writeString(first.resolve('late.bam'), 'complete')
            Thread.sleep(1600)
            java.nio.file.Files.writeString(other.resolve('last.bam'), 'complete')
            java.nio.file.Files.writeString(other.parent.resolve('final_summary_12345678.txt'), valid)
        }
        Thread.startDaemon('completion-test-watchdog') {
            Thread.sleep(8000)
            [first, second, other].each { dir: Path -> java.nio.file.Files.writeString(dir.resolve('STOP.watchdog'), '') }
        }
    }
    def args = record(live_analysis: true, timeline_analysis: false, sample_sheet_path: sheet,
        bam_poll_interval_ms: 100, bam_stability_polls: 3, termination_requested: false,
        detect_minknow_completion: true, minknow_completion_settle_seconds: 1)
    bam_ingress(args).collect(sort: false).map { batches ->
        assert batches*.batch_index == [0, 1, 2]
        assert batches[1].chunks*.sample*.name == ['first', 'second', 'other']
        assert batches[1].chunks.every { chunk -> java.nio.file.Files.readString(chunk.bam_paths[0]) == 'complete' }
        assert batches[2].chunks*.sample*.name == ['other']
        assert batches[2].chunks[0].bam_paths[0].fileName.toString() == 'last.bam'
        assert [first, second, other].every { dir: Path -> !java.nio.file.Files.exists(dir.resolve('STOP.watchdog')) }
        'Shared barcode runs and independently completing override run drained all BAMs without STOP writes'
    }.view()
}
