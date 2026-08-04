nextflow.enable.types = true

include {
    ChunkQCResult;
    FlagstatQCResult;
    MergedIndexedChunkBAM;
    NanoPlotQCResult
} from '../lib/sample.nf'
include { flagstat_file_name; nanoplot_output_dir } from '../modules/generic_helpers.nf'

/**
 * Run chunk-level QC for each merged sample chunk.
 *
 * Each input is already sorted, merged, and indexed for exactly one sample and
 * one batch index. QC stays per chunk so live reports can be refreshed as new
 * sequencing chunks arrive instead of waiting for a final sample BAM.
 */
workflow quality_control {
    take:
        merged_bams: Channel<MergedIndexedChunkBAM>

    main:
        nanoplot_qc_ch = nanoplot_qc(merged_bams)
        flagstat_qc_ch = samtools_flagstat_qc(merged_bams)
        qc_results_ch = nanoplot_qc_ch
            .map { NanoPlotQCResult result ->
                tuple(qc_result_key(result), result)
            }
            .join(
                flagstat_qc_ch.map { FlagstatQCResult result ->
                    tuple(qc_result_key(result), result)
                },
                by: 0
            )
            .map { _key, NanoPlotQCResult nanoplot_result, FlagstatQCResult flagstat_result ->
                record(
                    batch_index: nanoplot_result.batch_index,
                    sample: nanoplot_result.sample,
                    bam: nanoplot_result.bam,
                    bam_index: nanoplot_result.bam_index,
                    nanoplot_data: nanoplot_result.nanoplot_data,
                    flagstat: flagstat_result.flagstat
                )
            }

    emit:
        qc_results_ch
}

String qc_result_key(result) {
    return "${result.batch_index}\t${result.sample.group}\t${result.sample.name}"
}

/**
 * Extract the NanoPlot-compatible per-read table for one merged sample chunk.
 * The report only consumes this table, so skip NanoPlot's unused plot and HTML
 * rendering.
 */
process nanoplot_qc {
    label 'seq_lm_qc'
    container 'rnabioinfo/seq_lm_quality_control:v1.0.0'
    cpus 1

    input:
        merged_bam: MergedIndexedChunkBAM

    output:
        record(
            batch_index: merged_bam.batch_index,
            sample: merged_bam.sample,
            bam: merged_bam.bam,
            bam_index: merged_bam.bam_index,
            nanoplot_data: file("${nanoplot_output_dir(merged_bam)}/NanoPlot-data.tsv.gz")
        )

    script:
        String output_dir = nanoplot_output_dir(merged_bam)
        """
        bam-qc-table \
            --threads ${task.cpus} \
            --bam ${merged_bam.bam} \
            --output ${output_dir}/NanoPlot-data.tsv.gz
        """
}

/**
 * Run samtools flagstat independently for one merged sample chunk. NanoPlot
 * and flagstat are joined after both complete so neither blocks the other.
 */
process samtools_flagstat_qc {
    label 'seq_lm_qc'
    container 'rnabioinfo/seq_lm_samtools:v1.0.0'
    cpus 2

    input:
        merged_bam: MergedIndexedChunkBAM

    output:
        record(
            batch_index: merged_bam.batch_index,
            sample: merged_bam.sample,
            bam: merged_bam.bam,
            bam_index: merged_bam.bam_index,
            flagstat: file(flagstat_file_name(merged_bam.batch_index, merged_bam.sample.name))
        )

    script:
        String flagstat = flagstat_file_name(merged_bam.batch_index, merged_bam.sample.name)
        """
        samtools flagstat -@ ${task.cpus - 1} -O tsv ${merged_bam.bam} > ${flagstat}
        """
}
