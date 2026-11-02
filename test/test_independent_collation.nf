#!/usr/bin/env nextflow
nextflow.enable.types = true
include { independent_collated_batches } from '../subworkflows/quantification.nf'
workflow {
    def events: List = [[name:'slow',v:1,b:4],[name:'fast',v:0,b:1],[name:'fast',v:1,b:3],[name:'slow',v:0,b:0]].collect { row ->
        def key: String = "g\t${row.name}"
        record(batch_index:row.b,sample:record(name:row.name,group:'g'),bam:file("${row.name}_${row.v}.bam"),
            observation:[samples:[(key):[version:row.v,first_version:0,available_at:row.b+100]]])
    }
    def emitted: List = []
    independent_collated_batches(channel.fromList(events),2).map { row -> emitted.add(row.batch_index); return row }.collect().map { rows ->
        assert emitted==[1,3,0,4]
        assert rows.size()==4
        'Per-sample collation order does not wait on unrelated samples'
    }.view()
}
