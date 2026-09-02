#!/usr/bin/env nextflow

nextflow.enable.types = true

include { bam_ingress } from '../lib/bam_ingress.nf'

workflow {
    def root: Path = java.nio.file.Files.createTempDirectory('seq-lm-stable-polling-')
    def live_dir: Path = root.resolve('live')
    def final_dir: Path = root.resolve('final')
    def sample_sheet: Path = root.resolve('samples.csv')
    java.nio.file.Files.createDirectories(live_dir)
    java.nio.file.Files.createDirectories(final_dir)
    java.nio.file.Files.writeString(live_dir.resolve('startup.bam'), 'startup')
    java.nio.file.Files.writeString(final_dir.resolve('startup.bam'), 'startup')
    java.nio.file.Files.writeString(
        sample_sheet,
        "alias,group,bam_dir,is_live\nlive,control,${live_dir},true\nfinal,control,${final_dir},false\n",
    )

    Thread.startDaemon('stable-polling-test-writer') {
        Thread.sleep(500)
        def next_bam: Path = live_dir.resolve('next.bam')
        java.nio.file.Files.writeString(next_bam, 'partial')
        Thread.sleep(150)
        java.nio.file.Files.writeString(next_bam, 'complete')
        Thread.sleep(500)
        java.nio.file.Files.writeString(live_dir.resolve('STOP'), '')
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
        .map { batch ->
            if (batch.batch_index == 1) {
                assert batch.chunks*.sample*.name == ['live']
                assert java.nio.file.Files.readString(batch.chunks[0].bam_paths[0]) == 'complete'
            }
            return batch.batch_index
        }
        .collect()
        .map { batch_indices ->
            assert batch_indices.toSorted() == [0, 1]
            'BAM ingress waits for three stable polls in a mixed live/final run'
        }
        .view()
}
