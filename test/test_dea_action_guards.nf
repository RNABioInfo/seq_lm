#!/usr/bin/env nextflow
nextflow.enable.types = true
include { parse_guard_blocked } from '../lib/dea_monitoring.nf'
include { capture_dea_input } from '../lib/bam_ingress.nf'
workflow {
    def frontier: Map = [:]
    capture_dea_input(frontier,2,[started_at:100,samples:[k:[version:2,available_at:200]],closed_samples:['ended']])
    capture_dea_input(frontier,1,[started_at:100,samples:[k:[version:1,available_at:150]],closed_samples:[]])
    assert frontier._latest.k.version==2
    assert frontier._closed==['ended']
    def row: Map = [effectively_live:'true',processing_lag_seconds:'30',result_checked_at:'1000']
    def represented: Map = [available_at:990,version:2]
    assert !parse_guard_blocked(row,represented,[version:2],1000,1800,600)
    assert parse_guard_blocked(row,represented,[version:3],1000,1800,600)
    assert parse_guard_blocked(row,represented,[version:2],1600,1800,600)
    assert parse_guard_blocked(row,[available_at:-900,version:2],[version:2],1000,1800,600)
    assert parse_guard_blocked(row,null,null,1000,1800,600)
    assert !parse_guard_blocked(row+[effectively_live:'false'],null,null,1000,1800,600)
    channel.of('R actions recheck versions, ages, lag and missing input provenance').view()
}
