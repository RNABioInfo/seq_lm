#!/usr/bin/env nextflow
nextflow.enable.types = true
include { Sample ; QuantifiedSampleBatch ; GeneContrast } from '../lib/sample.nf'
include { differential_expression } from '../subworkflows/differential_expression.nf'
params {
    test_zero: Boolean = false
}
workflow {
    def root: Path = java.nio.file.Files.createTempDirectory('seqlm-r-deferral-')
    def annotation: Path = root.resolve('annotation.gtf')
    java.nio.file.Files.writeString(annotation,'chr\ttest\ttranscript\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "g1";\n')
    def optional: Path = root.resolve('OPTIONAL_FILE')
    java.nio.file.Files.writeString(optional,'')
    def samples: List = ['T0','T1'].collectMany { group: String ->
        (1..2).collect { index: Integer ->
            def sample: Sample = record(name:"${group}_${index}",group:group,order:null,bam_dir:root,is_live:false,protocol_run_id:null,minknow_run_dir:null)
            def path: Path = root.resolve("${group}_${index}.quant")
            java.nio.file.Files.writeString(path,'tname\tnum_reads\n'+(1..30).collect { gene: Integer -> "g${gene}\t${params.test_zero ? 0 : 1}\n" }.join(''))
            return record(batch_index:0,sample:sample,counts:path)
        }
    }
    def contrasts: List<GeneContrast> = [record(contrast_id:'pair',target_group:'T1',reference_group:'T0')]
    def batch: QuantifiedSampleBatch = record(batch_index:0,report_sequence:0,samples:samples)
    // Enrichment is requested, but must never receive an unavailable fit.
    outputs = differential_expression(channel.of(batch),contrasts,0,annotation,optional,true,1.0,1,2)
    outputs.report_batches.collect().map { reports ->
        assert reports.size()==1 && !reports[0].has_differential_results
        assert reports[0].differential_analysis_note
        return 'Unavailable fit is deferred before enrichment'
    }.view()
    outputs.precision_results.collect().map { snapshots ->
        assert snapshots.size()==1
        if (params.test_zero) assert snapshots[0].results.name=='OPTIONAL_FILE'
        else assert snapshots[0].results.resolve('analysis_deferred.txt').exists()
        return 'Exactly one deferred monitoring snapshot is emitted'
    }.view()
}
