#!/usr/bin/env nextflow

nextflow.enable.types = true

include { ChunkQCResult ; QuantifiedSample ; Sample } from '../lib/sample.nf'
include { file_identity ; sha256_file ; write_sample_checkpoints } from '../lib/sample_checkpoints.nf'

workflow {
    def root: Path = java.nio.file.Files.createTempDirectory('seq-lm-checkpoint-writer-')
    def bam_dir: Path = root.resolve('bams')
    def genome: Path = root.resolve('genome.fa')
    def annotation: Path = root.resolve('annotation.gtf')
    def old_quant: Path = root.resolve('sample-old.quant')
    def quant: Path = root.resolve('sample.quant')
    def old_nanoplot: Path = root.resolve('nanoplot-old.tsv.gz')
    def nanoplot: Path = root.resolve('nanoplot.tsv.gz')
    def old_flagstat: Path = root.resolve('flagstat-old.tsv')
    def flagstat: Path = root.resolve('flagstat.tsv')
    java.nio.file.Files.createDirectories(bam_dir)
    java.nio.file.Files.writeString(bam_dir.resolve('sample.bam'), 'bam')
    java.nio.file.Files.writeString(genome, '>ref\nACGT\n')
    java.nio.file.Files.writeString(annotation, 'annotation')
    java.nio.file.Files.writeString(old_quant, 'tname\tlen\tnum_reads\ntx1\t4\t5\n')
    java.nio.file.Files.writeString(quant, 'tname\tlen\tnum_reads\ntx1\t4\t10\n')
    java.nio.file.Files.writeString(old_nanoplot, 'nanoplot-old')
    java.nio.file.Files.writeString(nanoplot, 'nanoplot')
    java.nio.file.Files.writeString(old_flagstat, 'flagstat-old')
    java.nio.file.Files.writeString(flagstat, 'flagstat')

    def sample: Sample = record(
        name: 'sample',
        group: 'control',
        order: null,
        bam_dir: bam_dir,
        is_live: true,
        protocol_run_id: null,
    )
    def old_quantified: QuantifiedSample = record(batch_index: 2, sample: sample, counts: old_quant)
    def quantified: QuantifiedSample = record(batch_index: 3, sample: sample, counts: quant)
    def old_qc: ChunkQCResult = record(
        batch_index: 2,
        sample: sample,
        bam: bam_dir.resolve('sample.bam'),
        nanoplot_data: old_nanoplot,
        flagstat: old_flagstat,
    )
    def qc: ChunkQCResult = record(
        batch_index: 3,
        sample: sample,
        bam: bam_dir.resolve('sample.bam'),
        nanoplot_data: nanoplot,
        flagstat: flagstat,
    )

    written = write_sample_checkpoints(
        channel.of(quantified, old_quantified).collect().map { collected ->
            collected.toList().toSorted { left, right -> left.batch_index <=> right.batch_index }
        },
        channel.of(qc, old_qc).collect().map { collected ->
            collected.toList().toSorted { left, right -> left.batch_index <=> right.batch_index }
        },
        channel.value([genome: file_identity(genome), annotation: file_identity(annotation)]),
    )
    written
        .map { outputs ->
            def sample_dir: Path = outputs instanceof Collection ? outputs.iterator().next() : outputs
            def final_marker: Path = sample_dir.resolve('FINAL')
            assert java.nio.file.Files.isRegularFile(final_marker)
            assert !java.nio.file.Files.exists(sample_dir.resolve('manifest.in.json'))
            def manifest: Map = new groovy.json.JsonSlurper().parse(final_marker) as Map
            assert java.nio.file.Files.readString(sample_dir.resolve('quantification/final.quant')).contains('\t10\n')
            assert java.nio.file.Files.readString(sample_dir.resolve('qc/nanoplot/chunk_2.tsv.gz')) == 'nanoplot-old'
            assert java.nio.file.Files.isRegularFile(sample_dir.resolve('qc/nanoplot/chunk_3.tsv.gz'))
            assert java.nio.file.Files.readString(sample_dir.resolve('qc/flagstat/chunk_2.tsv')) == 'flagstat-old'
            assert java.nio.file.Files.isRegularFile(sample_dir.resolve('qc/flagstat/chunk_3.tsv'))
            assert manifest.quantification.sha256 == sha256_file(
                sample_dir.resolve(manifest.quantification.path as String)
            )
            manifest.qc.each { manifest_qc ->
                assert manifest_qc.nanoplot_sha256 == sha256_file(
                    sample_dir.resolve(manifest_qc.nanoplot as String)
                )
                assert manifest_qc.flagstat_sha256 == sha256_file(
                    sample_dir.resolve(manifest_qc.flagstat as String)
                )
            }
            'sample checkpoint writer passed'
        }
        .view()
}
