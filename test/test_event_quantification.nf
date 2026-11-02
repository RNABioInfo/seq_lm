#!/usr/bin/env nextflow
nextflow.enable.types = true
include { quantification } from '../subworkflows/quantification.nf'
include { Sample ; ChunkBAM ; SampleBatchSize } from '../lib/sample.nf'
workflow {
    def root: Path = java.nio.file.Files.createTempDirectory('seqlm-event-quantification-')
    def genome: Path = root.resolve('reference.fa')
    def annotation: Path = root.resolve('annotation.gtf')
    java.nio.file.Files.writeString(genome, '>chr1\n'+'A'*200+'\n')
    java.nio.file.Files.writeString(annotation,'chr1\ttest\ttranscript\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "tx1"; gene_biotype "protein_coding";\n')
    def samples: List<Sample> = ['fast','slow'].collect { name: String ->
        record(name:name,group:'control',order:null,bam_dir:root,is_live:true,protocol_run_id:null,minknow_run_dir:null)
    }
    def chunks: List<ChunkBAM> = (0..3).collect { index: Integer ->
        def sample = samples[index%2]
        def version: Integer = index.intdiv(2)
        def bam: Path = root.resolve("${sample.name}_${version}.bam")
        def sam: Path = root.resolve("${sample.name}_${version}.sam")
        java.nio.file.Files.writeString(sam,"@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:chr1\tLN:200\nread${index}\t4\t*\t0\t0\t*\t*\t0\t0\tAAAA\tIIII\n")
        def process: Process = new ProcessBuilder('samtools','view','-b','-o',bam.toString(),sam.toString()).start()
        assert process.waitFor()==0
        def key: String = "control\t${sample.name}"
        record(batch_index:index,sample:sample,bam:bam,
            observation:[started_at:100,observed_at:100+index,samples:[(key):[available_at:100+index,version:version,first_version:0]]])
    }
    def sizes: List<SampleBatchSize> = (0..3).collect { index: Integer -> record(batch_index:index,active_sample_count:1,experiment_sample_count:2) }
    outputs = quantification(channel.fromList(chunks),channel.fromList(sizes),[],genome,annotation,'unstranded',samples)
    def seen_fast: Integer = -1
    outputs.batches.map { batch ->
        def fast = batch.samples.find { quant -> quant.sample.name=='fast' }
        if (fast.batch_index >= 0) {
            def version: Integer = fast.observation.samples['control\tfast'].version
            assert version >= seen_fast
            seen_fast = version
        }
        return batch
    }.collect().map { batches ->
        assert batches.size()==4
        assert seen_fast==1
        def newest = batches.max { batch -> batch.report_sequence }
        assert newest.samples.every { quant -> quant.counts.text.contains('tx1\t100\t2') }
        assert newest.observation.samples.values().every { arrival -> arrival.version==1 }
        'Independent quantifications retain timestamps and reject late older versions'
    }.view()
}
