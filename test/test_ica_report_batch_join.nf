#!/usr/bin/env nextflow

nextflow.enable.types = true

include { join_ica_report_batches } from '../modules/ica_report_batches.nf'

workflow {
    reports = channel.of(
        [batch_index: 4, report_sequence: 1, biotypes: file('b4.tsv'), differential_analysis_note: '', has_differential_results: false, results: file('OPTIONAL_FILE'), stability_results: file('OPTIONAL_FILE'), has_stability_results: false],
        [batch_index: 2, report_sequence: 0, biotypes: file('b2.tsv'), differential_analysis_note: '', has_differential_results: false, results: file('OPTIONAL_FILE'), stability_results: file('OPTIONAL_FILE'), has_stability_results: false],
    )
    snapshots = channel.of(
        [batch_index: 2, report_sequence: 0, analysis_index: 10, results: file('ica_10')],
        [batch_index: 4, report_sequence: 1, analysis_index: 11, results: file('ica_11')],
    )
    join_ica_report_batches(reports, snapshots)
    join_ica_report_batches.out
        .map { result ->
            assert result.has_ica_results
            assert result.ica_results.name == "ica_${result.ica_analysis_index}"
            result.batch_index
        }
        .collect()
        .map { indices ->
            assert indices.toSet() == [2, 4].toSet()
            'ICA report join passed'
        }
        .view()
}
