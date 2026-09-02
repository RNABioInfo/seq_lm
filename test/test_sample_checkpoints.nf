#!/usr/bin/env nextflow

nextflow.enable.types = true

include { Sample } from '../lib/sample.nf'
include {
    discover_sample_checkpoints ;
    file_identity ;
    next_analysis_snapshot_index ;
    sample_bam_inventory ;
    sample_checkpoint_schema ;
    sample_checkpoint_tools ;
    sha256_file
} from '../lib/sample_checkpoints.nf'

workflow {
    def root: Path = java.nio.file.Files.createTempDirectory('seq-lm-checkpoint-test-')
    def output_root: Path = root.resolve('chosen-output')
    def genome: Path = root.resolve('genome.fa')
    def annotation: Path = root.resolve('annotation.gtf')
    java.nio.file.Files.writeString(genome, '>ref\nACGT\n')
    java.nio.file.Files.writeString(annotation, 'ref\ttest\ttranscript\t1\t4\t.\t+\t.\ttranscript_id "tx1";\n')

    def cached_bam_dir: Path = root.resolve('cached-bams')
    def active_bam_dir: Path = root.resolve('active-bams')
    java.nio.file.Files.createDirectories(cached_bam_dir)
    java.nio.file.Files.createDirectories(active_bam_dir)
    java.nio.file.Files.writeString(cached_bam_dir.resolve('cached.bam'), 'bam-a')
    java.nio.file.Files.writeString(active_bam_dir.resolve('active.bam'), 'bam-b')

    def cached: Sample = record(
        name: 'cached',
        group: 'control',
        order: null,
        bam_dir: cached_bam_dir,
        is_live: true,
        protocol_run_id: null,
    )
    def active: Sample = record(
        name: 'active',
        group: 'treated',
        order: null,
        bam_dir: active_bam_dir,
        is_live: true,
        protocol_run_id: null,
    )

    def cached_output: Path = output_root.resolve('control/cached')
    def quant: Path = cached_output.resolve('quantification/final.quant')
    def nanoplot: Path = cached_output.resolve('qc/nanoplot/chunk_0.tsv.gz')
    def flagstat: Path = cached_output.resolve('qc/flagstat/chunk_0.tsv')
    java.nio.file.Files.createDirectories(quant.parent)
    java.nio.file.Files.createDirectories(nanoplot.parent)
    java.nio.file.Files.createDirectories(flagstat.parent)
    java.nio.file.Files.writeString(quant, 'tname\tlen\tnum_reads\ntx1\t4\t10\n')
    java.nio.file.Files.writeString(nanoplot, 'qc-data')
    java.nio.file.Files.writeString(flagstat, 'flagstat-data')

    def manifest: Map = [schema_version: sample_checkpoint_schema(), sample: [name: cached.name, group: cached.group], bams: sample_bam_inventory(cached), reference: [genome: file_identity(genome), annotation: file_identity(annotation)], tools: sample_checkpoint_tools(), quantification: [path: 'quantification/final.quant', sha256: sha256_file(quant)], qc: [[batch_index: 0, nanoplot: 'qc/nanoplot/chunk_0.tsv.gz', nanoplot_sha256: sha256_file(nanoplot), flagstat: 'qc/flagstat/chunk_0.tsv', flagstat_sha256: sha256_file(flagstat)]]]
    java.nio.file.Files.writeString(
        cached_output.resolve('FINAL'),
        groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(manifest)) + '\n',
    )

    if ("${(params as Map).get('stale') ?: false}".toBoolean()) {
        java.nio.file.Files.writeString(quant, 'corrupt')
    }

    def discovered: Map = discover_sample_checkpoints(
        [cached, active],
        output_root,
        genome,
        annotation,
    )
    assert discovered.restored*.sample*.name == ['cached']
    assert discovered.active*.name == ['active']
    assert discovered.restored[0].quantification.counts == quant
    assert discovered.restored[0].qc_results.size() == 1

    java.nio.file.Files.createDirectories(output_root.resolve('differential_expression/batch_2'))
    java.nio.file.Files.createDirectories(output_root.resolve('differential_expression/batch_5'))
    assert next_analysis_snapshot_index(output_root) == 6
    channel.of('sample checkpoint helpers passed').view()
}
