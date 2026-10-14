#!/usr/bin/env nextflow

nextflow.enable.types = true

include { oarfish_quant } from '../subworkflows/quantification.nf'
include { CumulativeCollatedBAM ; Sample } from '../lib/sample.nf'
include { oarfish_counts_file_name ; oarfish_out_name } from '../modules/generic_helpers.nf'

workflow {
    def root: Path = java.nio.file.Files.createTempDirectory('seq-lm-strand-process-')
    def genome: Path = root.resolve('genome.fa')
    def annotation: Path = root.resolve('annotation.gtf')
    def bam: Path = root.resolve('sample.bam')
    java.nio.file.Files.writeString(genome, '>ref\nACGT\n')
    java.nio.file.Files.writeString(annotation, 'annotation')
    java.nio.file.Files.writeString(bam, 'bam')
    def strand: String = (params as Map).get('strand_specificity', 'forward') as String
    def sample: Sample = record(
        name: 'sample', group: 'control', order: null, bam_dir: root,
        is_live: false, protocol_run_id: null,
    )
    def snapshot: CumulativeCollatedBAM = record(batch_index: 3, sample: sample, bam: bam)
    oarfish_quant(channel.of(snapshot), genome, annotation, strand)
        .map { result ->
            assert result.batch_index == 3
            assert result.sample.name == sample.name
            assert result.sample.group == sample.group
            assert result.sample.is_live == sample.is_live
            assert result.sample.bam_dir.toRealPath() == sample.bam_dir.toRealPath()
            assert result.counts.name == oarfish_counts_file_name(3, 'sample')
            assert java.nio.file.Files.readString(result.counts) == 'tname\tlen\tnum_reads\ntx1\t4\t10\n'
            def args: List<String> = java.nio.file.Files.readAllLines(
                result.counts.parent.resolve("${oarfish_out_name(3, 'sample')}.args")
            )
            def expected: Map<String,String> = [forward: 'fw', reverse: 'rc', unstranded: 'both']
            assert args == ['-j', '4', '--genome-alignments', 'sample.bam',
                '--annotation', 'annotation.gtf', '--genome-fasta', 'genome.fa',
                '--filter-group', 'no-filters', '-d', expected[strand],
                '--output', oarfish_out_name(3, 'sample')]: args
            "Oarfish ${strand} arguments and output contract passed"
        }
        .view()
}
