#!/usr/bin/env nextflow
nextflow.enable.types = true

include { join_ica_report_batches } from '../modules/ica_report_batches.nf'
include { join_report_batches } from '../modules/report_batches.nf'

params {
    kind: String = 'ica'
    scenario: String = 'valid'
}

workflow {
    def reports: List<Map> = [0, 1].collect { index: Integer ->
        [batch_index: index + 2, report_sequence: index, biotypes: file('biotypes.tsv'),
         differential_analysis_note: '', has_differential_results: false,
         results: file('OPTIONAL_FILE'), stability_results: file('OPTIONAL_FILE'),
         has_stability_results: false, ica_results: file('ica'), has_ica_results: true,
         ica_analysis_index: index + 10]
    }
    def snapshots: List<Map> = [0, 1].collect { index: Integer ->
        [batch_index: index + 2, report_sequence: index, analysis_index: index + 10, results: file('ica')]
    }
    if (params.scenario == 'missing') {
        snapshots.remove(1)
    }
    if (params.scenario == 'missing_report') {
        reports.remove(1)
    }
    if (params.scenario == 'duplicate_report') {
        reports.add(reports[0])
    }
    if (params.scenario == 'duplicate_snapshot') {
        snapshots.add(snapshots[0])
    }
    if (params.scenario == 'sequence_mismatch') {
        snapshots[0].report_sequence = 9
    }
    if (params.scenario == 'sequence_gap') {
        reports[1].report_sequence = 2
    }
    if (params.scenario == 'duplicate_sequence') {
        reports[1].report_sequence = 0
    }
    if (params.kind == 'ica') {
        joined = join_ica_report_batches(channel.fromList(reports), channel.fromList(snapshots))
    }
    else {
        trees = channel.fromList(snapshots).map { snapshot ->
            [qc_report_inputs: [latest_batch_index: snapshot.batch_index], qc_results: file('qc')]
        }
        joined = join_report_batches(channel.fromList(reports), trees)
    }
    joined.toList().map { values ->
        assert values.size() == 2
        'Report integrity join passed'
    }.view()
}
