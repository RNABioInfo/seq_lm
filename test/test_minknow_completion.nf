#!/usr/bin/env nextflow
nextflow.enable.types = true

workflow {
    def root: Path = java.nio.file.Files.createTempDirectory('seq-lm-summary-')
    def run: Path = root.resolve('20261012_1200_device_flowcell_12345678')
    java.nio.file.Files.createDirectories(run.resolve('bam_pass/barcode01'))
    java.nio.file.Files.createDirectories(run.resolve('bam_fail'))
    assert MinknowCompletion.resolveRunDirectory(run) == run
    assert MinknowCompletion.resolveRunDirectory(run.resolve('bam_pass')) == run
    assert MinknowCompletion.resolveRunDirectory(run.resolve('bam_fail')) == run
    assert MinknowCompletion.resolveRunDirectory(run.resolve('bam_pass/barcode01')) == run
    assert MinknowCompletion.resolveRunDirectory(root) == null
    assert MinknowCompletion.resolveRunDirectory(run.resolve('barcode01')) == null
    assert MinknowCompletion.resolveRunDirectory(run.resolve('arbitrary/subdir')) == null
    assert MinknowCompletion.resolveRunDirectory(root, run) == run
    assert MinknowCompletion.inspect(run).reason != null
    java.nio.file.Files.writeString(run.resolve('sequencing_summary.txt'), 'not a completion signal')
    assert MinknowCompletion.inspect(run).reason != null
    def summary: Path = run.resolve('final_summary.txt')
    def valid: String = 'protocol_run_id=12345678-abcd-1234-abcd-123456789abc\nacquisition_run_id=abcdef123456\nacquisition_stopped=2026-10-12T12:00:00+00:00\nprocessing_stopped=2026-10-12T12:01:00Z\n'
    def invalid: List<String> = ['', 'partial', 'protocol_run_id=12345678\n', valid.replace('12:01:00Z', 'invalid'), valid.replace('12:01:00Z', '11:00:00Z'), valid.replace('12345678-abcd', '87654321-abcd'), valid + 'protocol_run_id=duplicate\n']
    invalid.each { contents: String ->
        java.nio.file.Files.writeString(summary, contents)
        assert MinknowCompletion.inspect(run).reason != null
    }
    java.nio.file.Files.writeString(summary, valid)
    assert MinknowCompletion.inspect(run, 'wrong').reason != null
    def evidence: Map = MinknowCompletion.inspect(run, '12345678-abcd-1234-abcd-123456789abc')
    assert evidence.reason == null
    java.nio.file.Files.writeString(run.resolve('final_summary_other.txt'), valid)
    assert MinknowCompletion.inspect(run).reason != null
    java.nio.file.Files.delete(run.resolve('final_summary_other.txt'))
    def state: Map = [:]
    assert !MinknowCompletion.settled(state, evidence, [:], false, 2, 1000, 0)
    assert !MinknowCompletion.settled(state, evidence, [:], true, 2, 1000, 2_000_000_000)
    assert !MinknowCompletion.settled(state, evidence, ['late': [2, 3]], false, 2, 1000, 2_000_000_000)
    assert !MinknowCompletion.settled(state, evidence, ['late': [2, 4]], false, 2, 1000, 3_000_000_000)
    assert MinknowCompletion.settled(state, evidence, ['late': [2, 4]], false, 2, 1000, 4_000_000_000)
    assert !MinknowCompletion.settled(state, [reason: 'disappeared'], [:], false, 2, 1000, 5_000_000_000)
    assert state.isEmpty()
    java.nio.file.Files.writeString(summary, valid + 'extra=changed\n')
    assert MinknowCompletion.inspect(run).signature != evidence.signature
    java.nio.file.Files.delete(summary)
    assert MinknowCompletion.inspect(run).reason != null
    channel.of('MinKNOW layout, summary validation and deterministic settling checks passed').view()
}
