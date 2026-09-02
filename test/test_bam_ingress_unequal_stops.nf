#!/usr/bin/env nextflow

nextflow.enable.types = true

include { bam_ingress } from '../lib/bam_ingress.nf'

workflow {
    def root: Path = java.nio.file.Files.createTempDirectory('seq-lm-unequal-stops-')
    def control_dir: Path = root.resolve('control_1')
    def control_2_dir: Path = root.resolve('control_2')
    def treated_dir: Path = root.resolve('treated')
    def sample_sheet: Path = root.resolve('samples.csv')
    java.nio.file.Files.createDirectories(control_dir)
    java.nio.file.Files.createDirectories(control_2_dir)
    java.nio.file.Files.createDirectories(treated_dir)
    java.nio.file.Files.writeString(control_dir.resolve('startup.bam'), 'bam')
    java.nio.file.Files.writeString(control_2_dir.resolve('startup.bam'), 'bam')
    java.nio.file.Files.writeString(treated_dir.resolve('startup.bam'), 'bam')
    java.nio.file.Files.writeString(
        sample_sheet,
        "alias,group,bam_dir,is_live\ncontrol_1,control,${control_dir},true\ncontrol_2,control,${control_2_dir},true\ntreated_1,treated,${treated_dir},true\n",
    )

    Thread.startDaemon('unequal-stop-test-writer') {
        Thread.sleep(1200)
        java.nio.file.Files.writeString(control_dir.resolve('batch_1.bam'), 'bam')
        java.nio.file.Files.writeString(control_2_dir.resolve('batch_1.bam'), 'bam')
        java.nio.file.Files.writeString(treated_dir.resolve('batch_1.bam'), 'bam')
        Thread.sleep(300)
        java.nio.file.Files.writeString(treated_dir.resolve('STOP'), '')
        Thread.sleep(300)
        java.nio.file.Files.writeString(control_dir.resolve('batch_2.bam'), 'bam')
        java.nio.file.Files.writeString(control_2_dir.resolve('batch_2.bam'), 'bam')
        Thread.sleep(300)
        java.nio.file.Files.writeString(control_dir.resolve('STOP'), '')
        java.nio.file.Files.writeString(control_2_dir.resolve('STOP'), '')
    }

    def ingress_args = record(
        live_analysis: true,
        timeline_analysis: false,
        sample_sheet_path: sample_sheet,
        bam_poll_interval_ms: 100,
        bam_stability_polls: 3,
        termination_requested: false,
    )
    bam_ingress(ingress_args)
        .collect()
        .map { batches ->
            def ordered_batches = batches.toSorted { left, right -> left.batch_index <=> right.batch_index }
            assert ordered_batches*.batch_index == [0, 1, 2]
            assert ordered_batches[1].chunks*.sample*.name == ['control_1', 'control_2', 'treated_1']
            assert ordered_batches[2].chunks*.sample*.name == ['control_1', 'control_2']
            'Live ingress permits unequal sample STOP batch counts'
        }
        .view()
}
