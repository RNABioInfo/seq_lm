#!/usr/bin/env nextflow

nextflow.enable.types = true

include { bam_ingress } from '../lib/bam_ingress.nf'

workflow {
    def root: Path = java.nio.file.Files.createTempDirectory('seq-lm-bam-ingress-record-')
    def first_bam_dir: Path = root.resolve('first')
    def second_bam_dir: Path = root.resolve('second')
    def sample_sheet: Path = root.resolve('samples.csv')

    java.nio.file.Files.createDirectories(first_bam_dir)
    java.nio.file.Files.createDirectories(second_bam_dir)
    java.nio.file.Files.writeString(first_bam_dir.resolve('first.bam'), 'bam')
    java.nio.file.Files.writeString(second_bam_dir.resolve('second.bam'), 'bam')
    java.nio.file.Files.writeString(
        sample_sheet,
        "alias,group,bam_dir,is_live\nfirst,control,${first_bam_dir},false\nsecond,control,${second_bam_dir},false\n",
    )

    def ingress_args = record(
        live_analysis: false,
        timeline_analysis: false,
        sample_sheet_path: sample_sheet,
        bam_poll_interval_ms: 100,
        bam_stability_polls: 3,
        termination_requested: false,
    )
    bam_ingress(ingress_args)
        .map { batch ->
            assert batch.batch_index == 0
            assert batch.chunks*.sample*.name == ['first', 'second']
            'BAM ingress accepts record arguments'
        }
        .view()
}
