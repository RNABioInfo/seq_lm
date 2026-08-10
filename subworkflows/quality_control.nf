nextflow.enable.types = true

include {
    ChunkBAM;
    ChunkQCResult;
    FlagstatQCResult;
    NanoPlotQCResult
} from '../lib/sample.nf'
include { flagstat_file_name; nanoplot_output_dir } from '../modules/generic_helpers.nf'

/**
 * Run chunk-level QC for each merged sample chunk.
 *
 * Each input contains exactly one sample and one batch index. QC scans the BAM
 * sequentially, so coordinate sorting and indexing are unnecessary. QC stays
 * per chunk so live reports can refresh as sequencing chunks arrive.
 */
workflow quality_control {
    take:
        merged_bams: Channel<ChunkBAM>

    main:
        nanoplot_qc_ch = nanoplot_qc(merged_bams)
        flagstat_qc_ch = samtools_flagstat_qc(merged_bams)
        qc_results_ch = nanoplot_qc_ch
            .map { result ->
                record(
                    qc_key: qc_result_key(result),
                    nanoplot_result: result
                )
            }
            .join(flagstat_qc_ch.map { result ->
                    record(
                        qc_key: qc_result_key(result),
                        flagstat_result: result
                    )
                },
                by: 'qc_key'
            )
            .map { joined ->
                record(
                    batch_index: joined.nanoplot_result.batch_index,
                    sample: joined.nanoplot_result.sample,
                    bam: joined.nanoplot_result.bam,
                    nanoplot_data: joined.nanoplot_result.nanoplot_data,
                    flagstat: joined.flagstat_result.flagstat
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
        merged_bam: ChunkBAM

    output:
        record(
            batch_index: merged_bam.batch_index,
            sample: merged_bam.sample,
            bam: merged_bam.bam,
            nanoplot_data: file("${nanoplot_output_dir(merged_bam)}/NanoPlot-data.tsv.gz")
        )

    script:
        def output_dir: String = nanoplot_output_dir(merged_bam)
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
        merged_bam: ChunkBAM

    output:
        record(
            batch_index: merged_bam.batch_index,
            sample: merged_bam.sample,
            bam: merged_bam.bam,
            flagstat: file(flagstat_file_name(merged_bam.batch_index, merged_bam.sample.name))
        )

    script:
        def flagstat: String = flagstat_file_name(merged_bam.batch_index, merged_bam.sample.name)
        """
        samtools flagstat -@ ${task.cpus - 1} -O tsv ${merged_bam.bam} > ${flagstat}
        """
}
