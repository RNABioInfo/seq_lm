nextflow.enable.types = true

/**
 * Join sparse, ordered edgeR results to matching QC report inputs.
 *
 * edgeR may begin after batch 0 when live samples initially have no data, so
 * batch_index itself cannot drive the existing contiguous order_batches helper.
 * Assign a contiguous sequence to edgeR emissions, then restore that sequence
 * after the keyed join in case QC inputs arrive out of order.
 */
workflow join_report_batches {
    take:
    differential_results: Channel
    qc_report_trees: Channel

    main:
    def next_report_sequence: Integer = 0
    sequenced_differential_results_ch = differential_results.map { result ->
        def report_sequence: Integer = next_report_sequence
        next_report_sequence += 1
        tuple(
            result.batch_index,
            report_sequence,
            result.results,
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
        .map { batch_index: Integer, report_sequence: Integer, differential_results_path, report_inputs: Map, qc_results_path ->
            record(
                batch_index: batch_index,
                report_sequence: report_sequence,
                qc_report_inputs: report_inputs,
                qc_results: qc_results_path,
                differential_results: differential_results_path,
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
