#!/usr/bin/env nextflow

nextflow.enable.types = true

include { join_report_batches } from '../modules/report_batches.nf'

workflow {
    /*
     * A failed later precondition may arrive before an earlier edgeR result.
     * QC inputs also arrive in the opposite order to prove the keyed join is
     * restored by the sequence assigned before readiness checking.
     */
    differential_results_ch = channel.of(
        [batch_index: 4, report_sequence: 1, read_depth_satisfied: false, results: file('OPTIONAL_FILE')],
        [batch_index: 2, report_sequence: 0, read_depth_satisfied: true, results: file('edgeR_batch_2')],
    )
    qc_report_trees_ch = channel.of(
        [qc_report_inputs: [latest_batch_index: 4, rows: 'batch 4'], qc_results: file('qc_batch_4')],
        [qc_report_inputs: [latest_batch_index: 2, rows: 'batch 2'], qc_results: file('qc_batch_2')],
    )

    joined_ch = join_report_batches(
        differential_results_ch,
        qc_report_trees_ch,
    )
    joined_ch
        .map { report -> report.batch_index }
        .collect()
        .map { batch_indices ->
            assert batch_indices.toList() == [2, 4]
            "report batch join passed: ${batch_indices.join(',')}"
        }
        .view()
}
