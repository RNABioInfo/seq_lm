nextflow.enable.types = true

/**
 * Join ordered quantification/analysis records to matching QC report inputs.
 *
 * Every reference-enabled report batch contains a transcript-biotype summary;
 * it may also carry matching DEA readiness/results and stability output. The
 * preassigned sequence preserves batch order when process runtimes differ.
 */
workflow join_report_batches {
    take:
    analysis_reports: Channel
    qc_report_trees: Channel

    main:
    sequenced_analysis_results_ch = analysis_reports.map { result ->
        tuple(
            result.batch_index,
            result.report_sequence,
            result.biotypes,
            result.differential_analysis_note,
            result.has_differential_results,
            result.results,
            result.stability_results,
            result.has_stability_results,
            result.ica_results,
            result.has_ica_results,
            result.ica_analysis_index,
        )
    }

    joined_report_batches_ch = sequenced_analysis_results_ch
        .join(
            qc_report_trees.map { result ->
                tuple(
                    result.qc_report_inputs.latest_batch_index,
                    result.qc_report_inputs,
                    result.qc_results,
                )
            },
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .map { batch_index: Integer, report_sequence: Integer, transcript_biotypes_path, differential_analysis_note: String, has_differential_results: Boolean, differential_results_path, stability_results_path, has_stability_results: Boolean, ica_results_path, has_ica_results: Boolean, ica_analysis_index: Integer, report_inputs: Map, qc_results_path ->
            record(
                batch_index: batch_index,
                report_sequence: report_sequence,
                qc_report_inputs: report_inputs + [report_sequence: report_sequence],
                qc_results: qc_results_path,
                transcript_biotypes: transcript_biotypes_path,
                differential_analysis_note: differential_analysis_note,
                has_differential_results: has_differential_results,
                differential_results: differential_results_path,
                stability_results: stability_results_path,
                has_stability_results: has_stability_results,
                ica_results: ica_results_path,
                has_ica_results: has_ica_results,
                ica_analysis_index: ica_analysis_index,
            )
        }

    def pending_reports: Map<Integer, ?> = [:]
    def next_sequence_to_emit: Integer = 0
    // The sentinel runs after channel completion, so gaps cannot silently leave
    // later reports buffered. It does not delay complete live reports.
    ordered_report_batches_ch = joined_report_batches_ch
        .map { report -> [report: report] }
        .concat(channel.of([complete: true]))
        .map { event ->
            if (event.complete) {
                if (!pending_reports.isEmpty()) {
                    error("Missing report sequence ${next_sequence_to_emit}; unpublished sequences: ${pending_reports.keySet().toSorted()}.")
                }
                return []
            }
            def report_batch = event.report
            def sequence: Integer = report_batch.report_sequence
            if (sequence < next_sequence_to_emit || pending_reports.containsKey(sequence)) {
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
        }.flatMap { reports -> reports }

    emit:
    ordered_report_batches_ch
}
