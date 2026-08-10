#!/usr/bin/env nextflow

nextflow.enable.types = true

include { ChunkQCResult; QuantifiedSample; Sample } from '../lib/sample.nf'
include { file_identity; write_sample_checkpoints } from '../lib/sample_checkpoints.nf'

workflow {
    Path root = java.nio.file.Files.createTempDirectory('seq-lm-checkpoint-writer-')
    Path bam_dir = root.resolve('bams')
    Path genome = root.resolve('genome.fa')
    Path annotation = root.resolve('annotation.gtf')
    Path old_quant = root.resolve('sample-old.quant')
    Path quant = root.resolve('sample.quant')
    Path old_nanoplot = root.resolve('nanoplot-old.tsv.gz')
    Path nanoplot = root.resolve('nanoplot.tsv.gz')
    Path old_flagstat = root.resolve('flagstat-old.tsv')
    Path flagstat = root.resolve('flagstat.tsv')
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

    Sample sample = record(
        name: 'sample',
        group: 'control',
        order: null,
        bam_dir: bam_dir,
        is_live: true
    )
    QuantifiedSample old_quantified = record(batch_index: 2, sample: sample, counts: old_quant)
    QuantifiedSample quantified = record(batch_index: 3, sample: sample, counts: quant)
    ChunkQCResult old_qc = record(
        batch_index: 2,
        sample: sample,
        bam: bam_dir.resolve('sample.bam'),
        nanoplot_data: old_nanoplot,
        flagstat: old_flagstat
    )
    ChunkQCResult qc = record(
        batch_index: 3,
        sample: sample,
        bam: bam_dir.resolve('sample.bam'),
        nanoplot_data: nanoplot,
        flagstat: flagstat
    )

    written = write_sample_checkpoints(
        channel.of(quantified, old_quantified).collect().map { collected ->
            collected.toList().toSorted { left, right -> left.batch_index <=> right.batch_index }
        },
        channel.of(qc, old_qc).collect().map { collected ->
            collected.toList().toSorted { left, right -> left.batch_index <=> right.batch_index }
        },
        channel.value([genome: file_identity(genome), annotation: file_identity(annotation)])
    )
    written.map { outputs ->
        Path sample_dir = outputs instanceof Collection ? outputs.iterator().next() : outputs
        Path final_marker = sample_dir.resolve('FINAL')
        assert java.nio.file.Files.isRegularFile(final_marker)
        assert java.nio.file.Files.readString(sample_dir.resolve('quantification/final.quant')).contains('\t10\n')
        assert java.nio.file.Files.readString(sample_dir.resolve('qc/nanoplot/chunk_2.tsv.gz')) == 'nanoplot-old'
        assert java.nio.file.Files.isRegularFile(sample_dir.resolve('qc/nanoplot/chunk_3.tsv.gz'))
        assert java.nio.file.Files.readString(sample_dir.resolve('qc/flagstat/chunk_2.tsv')) == 'flagstat-old'
        assert java.nio.file.Files.isRegularFile(sample_dir.resolve('qc/flagstat/chunk_3.tsv'))
        'sample checkpoint writer passed'
    }.view()
}
