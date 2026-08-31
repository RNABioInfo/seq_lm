nextflow.enable.types = true

/**
 * Join ordered DEA readiness/results records to matching QC report inputs.
 *
 * A readiness record exists for every batch considered for DEA, including
 * batches that do not run edgeR. Its preassigned sequence preserves batch order
 * even when a fast failed precondition overtakes an earlier edgeR execution.
 */
workflow join_report_batches {
    take:
    differential_reports: Channel
    qc_report_trees: Channel

    main:
    sequenced_differential_results_ch = differential_reports.map { result ->
        tuple(
            result.batch_index,
            result.report_sequence,
            result.differential_analysis_note,
            result.has_differential_results,
            result.results,
            result.stability_results,
            result.has_stability_results,
        )
    }

    joined_report_batches_ch = sequenced_differential_results_ch
        .join(
            qc_report_trees.map { result ->
                tuple(
                    result.qc_report_inputs.latest_batch_index,
                    result.qc_report_inputs,
                    result.qc_results,
                )
            },
            by: 0
        )
        .map { batch_index: Integer, report_sequence: Integer, differential_analysis_note: String, has_differential_results: Boolean, differential_results_path, stability_results_path, has_stability_results: Boolean, report_inputs: Map, qc_results_path ->
            record(
                batch_index: batch_index,
                report_sequence: report_sequence,
                qc_report_inputs: report_inputs,
                qc_results: qc_results_path,
                differential_analysis_note: differential_analysis_note,
                has_differential_results: has_differential_results,
                differential_results: differential_results_path,
                stability_results: stability_results_path,
                has_stability_results: has_stability_results,
            )
        }

    def pending_reports: Map<Integer, ?> = [:]
    def next_sequence_to_emit: Integer = 0
    ordered_report_batches_ch = joined_report_batches_ch.flatMap { report_batch ->
        def sequence: Integer = report_batch.report_sequence
        if (pending_reports.containsKey(sequence)) {
            error("Received duplicate report sequence ${sequence}.")
        }
        pending_reports[sequence] = report_batch

        def ready_sequences: List<Integer> = pending_reports
            .keySet()
            .toSorted()
            .withIndex()
            .findAll { entry ->
                entry[0] == next_sequence_to_emit + entry[1]
            }
            .collect { entry -> entry[0] as Integer }
        def ready_reports: List = ready_sequences.collect { ready_sequence: Integer ->
            pending_reports.remove(ready_sequence)
        }
        next_sequence_to_emit += ready_sequences.size()
        return ready_reports
    }

    emit:
    ordered_report_batches_ch
}
