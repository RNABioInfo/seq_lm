include { shell_quote; safe_name; } from './generic_helpers.nf'

/**
 * Sanitize table fields for the live QC report.
 */
String qc_report_field(Object value) {
    return "${value}".replaceAll(/[\t\r\n]+/, ' ')
}

/**
 * Root input directory staged next to the live QC report.
 */
String qc_report_root_dir() {
    return "qc_results"
}

String qc_report_nanoplot_sample_dir(Map report_inputs) {
    return "${qc_report_root_dir()}/nanoplot/${safe_name(report_inputs.sample.group)}/${safe_name(report_inputs.sample.name)}"
}

String qc_report_flagstat_sample_dir(Map report_inputs) {
    return "${qc_report_root_dir()}/flagstat/${safe_name(report_inputs.sample.group)}/${safe_name(report_inputs.sample.name)}"
}

/**
 * Rows for the live QC report sample table.
 */
String qc_report_rows(List report_inputs_list) {
    return report_inputs_list
        .collect { report_inputs ->
            [
                qc_report_field(report_inputs.sample.name),
                qc_report_field(report_inputs.sample.group),
                "${report_inputs.chunks.size()}",
                "${report_inputs.latest_batch_index}",
                qc_report_field(qc_report_root_dir())
            ].join('\t')
        }
        .join('\n')
}

/**
 * Flatten NanoPlot inputs into the same order used by qc_report_copy_commands.
 */
String qc_report_copy_commands(List report_inputs_list, List nanoplot_inputs, List flagstat_inputs) {
    List<String> commands = []
    Integer staged_index = 0

    report_inputs_list.each { report_inputs ->
        String nanoplot_dir = qc_report_nanoplot_sample_dir(report_inputs)
        String flagstat_dir = qc_report_flagstat_sample_dir(report_inputs)
        report_inputs.chunks.each { chunk ->
            String nanoplot_source = "${nanoplot_inputs[staged_index]}"
            String flagstat_source = "${flagstat_inputs[staged_index]}"
            String nanoplot_name = "${chunk.nanoplot_name}"
            String flagstat_name = "${chunk.flagstat_name}"

            commands << "mkdir -p ${shell_quote(nanoplot_dir)} ${shell_quote(flagstat_dir)}"
            commands << "cp ${shell_quote(nanoplot_source)} ${shell_quote("${nanoplot_dir}/${nanoplot_name}")}"
            commands << "cp ${shell_quote(flagstat_source)} ${shell_quote("${flagstat_dir}/${flagstat_name}")}"
            staged_index += 1
        }
    }

    return commands.join('\n')
}

/**
 * Convert cumulative QC state into the report input shape consumed by the
 * Python report process. The metadata intentionally contains no Path objects;
 * files are passed separately through typed process inputs.
 */
Map qc_report_inputs_from_state(Integer latest_batch_index, Map<String, List<ChunkQCResult>> state) {
    List<List<ChunkQCResult>> sorted_sample_chunks = state
        .values()
        .toList()
        .toSorted { List<ChunkQCResult> left_chunks, List<ChunkQCResult> right_chunks ->
            ChunkQCResult left = left_chunks[0]
            ChunkQCResult right = right_chunks[0]
            "${left.sample.group}/${left.sample.name}" <=>
                "${right.sample.group}/${right.sample.name}"
        }

    List report_inputs_list = sorted_sample_chunks.collect { List<ChunkQCResult> sample_chunks ->
        ChunkQCResult first = sample_chunks[0]
        List sorted_chunks = sample_chunks.toSorted { ChunkQCResult left, ChunkQCResult right ->
            left.batch_index <=> right.batch_index
        }

        [
            latest_batch_index: latest_batch_index,
            sample: [
                name: first.sample.name,
                group: first.sample.group
            ],
            chunks: sorted_chunks.collect { ChunkQCResult result ->
                [
                    batch_index: result.batch_index,
                    nanoplot_name: "nanoplot_data_chunk_${result.batch_index}.tsv.gz",
                    flagstat_name: "flagstat_data_chunk_${result.batch_index}.tsv"
                ]
            }
        ]
    }

    List<ChunkQCResult> sorted_all_chunks = sorted_sample_chunks.collectMany { List<ChunkQCResult> sample_chunks ->
        sample_chunks.toSorted { ChunkQCResult left, ChunkQCResult right ->
            left.batch_index <=> right.batch_index
        }
    }

    return [
        latest_batch_index: latest_batch_index,
        report_inputs_list: report_inputs_list,
        rows: qc_report_rows(report_inputs_list),
        nanoplot_inputs: sorted_all_chunks*.nanoplot_data,
        flagstat_inputs: sorted_all_chunks*.flagstat
    ]
}

Map accumulate_qc_report_chunk_state(
    Map<String, List<ChunkQCResult>> state,
    Integer batch_index,
    List<ChunkQCResult> chunk_results
) {
    chunk_results.each { ChunkQCResult result ->
        String key = "${result.sample.group}\t${result.sample.name}"
        List<ChunkQCResult> previous_chunks = state.containsKey(key) ? state[key] : []
        state[key] = (previous_chunks.findAll { ChunkQCResult chunk ->
            chunk.batch_index != result.batch_index
        } + [result]).toSorted { ChunkQCResult left, ChunkQCResult right ->
            left.batch_index <=> right.batch_index
        }
    }

    return qc_report_inputs_from_state(batch_index, state)
}
