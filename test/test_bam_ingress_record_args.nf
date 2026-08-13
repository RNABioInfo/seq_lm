#!/usr/bin/env nextflow

nextflow.enable.types = true

include { get_samples } from '../lib/bam_ingress.nf'

workflow {
    def root: Path = java.nio.file.Files.createTempDirectory('seq-lm-bam-ingress-record-')
    def first_bam_dir: Path = root.resolve('first')
    def second_bam_dir: Path = root.resolve('second')
    def sample_sheet: Path = root.resolve('samples.csv')

    java.nio.file.Files.createDirectories(first_bam_dir)
    java.nio.file.Files.createDirectories(second_bam_dir)
    java.nio.file.Files.writeString(
        sample_sheet,
        "alias,group,bam_dir,is_live\nfirst,control,${first_bam_dir},false\nsecond,control,${second_bam_dir},false\n",
    )

    def ingress_args = record(
        live_analysis: false,
        timeline_analysis: false,
        sample_sheet_path: sample_sheet,
    )
    def samples = get_samples(ingress_args)

    assert samples*.name == ['first', 'second']
    channel.of('BAM ingress accepts record arguments').view()
}
