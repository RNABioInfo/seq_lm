nextflow.enable.types = true

/** Join every report record to the immutable ICA snapshot for the same batch. */
workflow join_ica_report_batches {
    take:
    analysis_reports: Channel
    ica_snapshots: Channel

    main:
    joined_ch = analysis_reports
        .map { report -> tuple(report.batch_index, report) }
        .join(
            ica_snapshots.map { snapshot -> tuple(snapshot.batch_index, snapshot) },
            by: 0,
            failOnDuplicate: true,
            failOnMismatch: true,
        )
        .map { batch_index: Integer, report, snapshot ->
            if (report.report_sequence != snapshot.report_sequence) {
                error("Inconsistent ICA report sequence for batch ${batch_index}.")
            }
            record(
                batch_index: batch_index,
                report_sequence: report.report_sequence,
                biotypes: report.biotypes,
                differential_analysis_note: report.differential_analysis_note,
                has_differential_results: report.has_differential_results,
                results: report.results,
                stability_results: report.stability_results,
                has_stability_results: report.has_stability_results,
                ica_results: snapshot.results,
                has_ica_results: true,
                ica_analysis_index: snapshot.analysis_index,
            )
        }

    emit:
    joined_ch
}
