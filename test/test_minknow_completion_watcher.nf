#!/usr/bin/env nextflow
nextflow.enable.types = true

workflow {
    def root: Path = java.nio.file.Files.createTempDirectory('seq-lm-watcher-')
    def run: Path = root.resolve('20261012_1200_device_flowcell_12345678')
    def bams: Path = run.resolve('bam_pass')
    java.nio.file.Files.createDirectories(bams)
    def summary: Path = run.resolve('final_summary.txt')
    def valid: String = 'protocol_run_id=12345678-abcd-1234-abcd-123456789abc\nacquisition_run_id=abcdef123456\nacquisition_stopped=2026-10-12T12:00:00Z\nprocessing_stopped=2026-10-12T12:01:00Z\n'
    java.nio.file.Files.writeString(summary, valid)
    def source = LiveBamWatcher.poll(bams, [], 100, 3,
        [enabled: true, run_dir: null, protocol_run_id: null, settle_seconds: 1],
        {
            Thread.startDaemon('watcher-summary-writer') {
                Thread.sleep(400)
                java.nio.file.Files.delete(summary)
                Thread.sleep(800)
                java.nio.file.Files.writeString(summary, 'partial')
                java.nio.file.Files.writeString(bams.resolve('first.bam'), 'complete')
                Thread.sleep(800)
                java.nio.file.Files.writeString(summary, valid)
                Thread.sleep(600)
                java.nio.file.Files.writeString(bams.resolve('last.bam'), 'partial')
                Thread.sleep(150)
                java.nio.file.Files.writeString(bams.resolve('last.bam'), 'complete')
            }
            Thread.startDaemon('watcher-test-watchdog') {
                Thread.sleep(7000)
                java.nio.file.Files.writeString(bams.resolve('STOP.watchdog'), '')
            }
        },
        { exception -> throw exception })
    source.collect(sort: false).map { events ->
        assert events*.source == ['poll', 'poll', 'minknow-final-summary']
        assert events.take(2)*.path*.fileName*.toString().toSorted() == ['first.bam', 'last.bam']
        assert events.last().path == summary
        assert events.findAll { event -> event.source == 'minknow-final-summary' }.size() == 1
        assert java.nio.file.Files.readString(events[1].path) == 'complete'
        assert !java.nio.file.Files.exists(bams.resolve('STOP.watchdog'))
        'Changing/invalid/disappearing evidence resets settling and emits exactly one stop after final BAM'
    }.view()
}
