#!/usr/bin/env nextflow
nextflow.enable.types = true
include { represented_input_stale ; guard_precision_actions } from '../lib/precision_monitoring.nf'
include { capture_precision_input } from '../lib/bam_ingress.nf'
workflow {
    def frontier: Map = [:]
    capture_precision_input(frontier,2,[started_at:100,samples:[k:[version:2,available_at:200]],closed_samples:['ended']])
    capture_precision_input(frontier,1,[started_at:100,samples:[k:[version:1,available_at:150]],closed_samples:[]])
    assert frontier._latest.k.version==2
    assert frontier._closed==['ended']
    def row: Map = [effectively_live:'true']
    def represented: Map = [available_at:990,version:2]
    assert !represented_input_stale(row,represented,[version:2],1000,1800)
    assert represented_input_stale(row,represented,[version:3],1000,1800)
    assert represented_input_stale(row,[available_at:-900,version:2],[version:2],1000,1800)
    assert represented_input_stale(row,null,null,1000,1800)
    assert !represented_input_stale(row+[effectively_live:'false'],null,null,1000,1800)
    def first: Map = row + [group: 'target', sample: 'a', eligible: true, required_contrasts: 'pair']
    def fixed: Map = row + [group: 'control', sample: 'b', effectively_live: 'false', required_contrasts: 'pair']
    def other: Map = row + [group: 'other', sample: 'c', eligible: true, required_contrasts: 'unrelated']
    def arrivals: Map = ['target\ta': represented, 'control\tb': represented, 'other\tc': represented]
    def timing: Map = [now: 1000, checked_at: 990, horizon: 1800, maximum_lag: 600]
    guard_precision_actions([first, fixed, other], [], arrivals, ['control\tb': [version: 3]], timing, false)
    assert !first.eligible && first.reason == 'newer_input_or_stale_result'
    assert other.eligible

    first.eligible = true
    def lagged: List<Map> = [[contrast_id: 'unrelated', processing_lag_seconds: '595']]
    guard_precision_actions([first, fixed, other], lagged, arrivals, [:], timing, false)
    assert first.eligible && !other.eligible

    guard_precision_actions([first, fixed, other], [], arrivals, [:], timing, true)
    assert !first.eligible && first.reason == 'terminal_snapshot'
    channel.of('R actions recheck versions, ages, lag and missing input provenance').view()
}
