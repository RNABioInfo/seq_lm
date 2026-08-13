nextflow.enable.types = true

include { ChunkQCResult } from '../lib/sample.nf'
include { shell_quote ; safe_name } from './generic_helpers.nf'

/**
 * Sanitize table fields for the live QC report.
 */
def qc_report_field(value: Object) -> String {
    return "${value}".replaceAll(/[\t\r\n]+/, ' ')
}

/**
 * Root input directory staged next to the live QC report.
 */
def qc_report_root_dir() -> String {
    return "qc_results"
}

def qc_report_nanoplot_sample_dir(report_inputs: Map) -> String {
    return "${qc_report_root_dir()}/nanoplot/${safe_name(report_inputs.sample.group)}/${safe_name(report_inputs.sample.name)}"
}

def qc_report_flagstat_sample_dir(report_inputs: Map) -> String {
    return "${qc_report_root_dir()}/flagstat/${safe_name(report_inputs.sample.group)}/${safe_name(report_inputs.sample.name)}"
}

/**
 * Rows for the live QC report sample table.
 */
def qc_report_rows(report_inputs_list: List) -> String {
    return report_inputs_list
        .collect { report_inputs ->
            [
            qc_report_field(report_inputs.sample.name), 
            qc_report_field(report_inputs.sample.group), 
            "${report_inputs.chunks.size()}", 
            "${report_inputs.latest_batch_index}", 
            qc_report_field(qc_report_root_dir())
            ]
            .join('\t')
        }
        .join('\n')
}

/**
 * Flatten NanoPlot inputs into the same order used by qc_report_copy_commands.
 */
def qc_report_copy_commands(report_inputs_list: List, nanoplot_inputs: List, flagstat_inputs: List) -> String {
    def commands: List<String> = []
    def staged_index: Integer = 0

    report_inputs_list.each { report_inputs ->
        def nanoplot_dir: String = qc_report_nanoplot_sample_dir(report_inputs)
        def flagstat_dir: String = qc_report_flagstat_sample_dir(report_inputs)
        report_inputs.chunks.each { chunk ->
            def nanoplot_source: String = "${nanoplot_inputs[staged_index]}"
            def flagstat_source: String = "${flagstat_inputs[staged_index]}"
            def nanoplot_name: String = "${chunk.nanoplot_name}"
            def flagstat_name: String = "${chunk.flagstat_name}"

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
def qc_report_inputs_from_state(latest_batch_index: Integer, state: Map<String,List<ChunkQCResult>>) -> Map {
    def sorted_sample_chunks: List<List<ChunkQCResult>> = state
        .values()
        .toList()
        .toSorted { left_chunks: List<ChunkQCResult>, right_chunks: List<ChunkQCResult> ->
            def left: ChunkQCResult = left_chunks[0]
            def right: ChunkQCResult = right_chunks[0]
            "${left.sample.group}/${left.sample.name}" <=> "${right.sample.group}/${right.sample.name}"
        }

    def report_inputs_list: List = sorted_sample_chunks.collect { sample_chunks: List<ChunkQCResult> ->
        def first: ChunkQCResult = sample_chunks[0]
        def sorted_chunks: List = sample_chunks.toSorted { left, right ->
            left.batch_index <=> right.batch_index
        }

        [latest_batch_index: latest_batch_index, sample: [name: first.sample.name, group: first.sample.group], chunks: sorted_chunks.collect { result ->
            [batch_index: result.batch_index, nanoplot_name: "nanoplot_data_chunk_${result.batch_index}.tsv.gz", flagstat_name: "flagstat_data_chunk_${result.batch_index}.tsv"]
        }]
    }

    def sorted_all_chunks: List<ChunkQCResult> = sorted_sample_chunks.collectMany { sample_chunks: List<ChunkQCResult> ->
        sample_chunks.toSorted { left, right ->
            left.batch_index <=> right.batch_index
        }
    }

    return [latest_batch_index: latest_batch_index, report_inputs_list: report_inputs_list, rows: qc_report_rows(report_inputs_list), nanoplot_inputs: sorted_all_chunks*.nanoplot_data, flagstat_inputs: sorted_all_chunks*.flagstat]
}

def accumulate_qc_report_chunk_state(state: Map<String,List<ChunkQCResult>>, batch_index: Integer, chunk_results: List<ChunkQCResult>) -> Map {
    chunk_results.each { result ->
        def key: String = "${result.sample.group}\t${result.sample.name}"
        def previous_chunks: List<ChunkQCResult> = state.containsKey(key) ? state[key] : []
        state[key] = (previous_chunks.findAll { chunk ->
            chunk.batch_index != result.batch_index
        } + [result]).toSorted { left, right ->
            left.batch_index <=> right.batch_index
        }
    }

    return qc_report_inputs_from_state(batch_index, state)
}
