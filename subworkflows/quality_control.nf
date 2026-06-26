nextflow.enable.types = true

include {
    ChunkQCResult;
    MergedChunkBAM;
    NanoPlotQCResult
} from '../lib/sample.nf'

/**
 * Return a filesystem-safe token for process-local QC output names.
 */
String safeName(String value) {
    return value.replaceAll(/[^A-Za-z0-9._-]/, '_')
}

/**
 * Build the stable NanoPlot output directory for one sample chunk.
 */
String nanoplotOutputDir(MergedChunkBAM merged_bam) {
    return "${safeName(merged_bam.sample.id)}_${merged_bam.batch_index}.nanoplot"
}

/**
 * Build the stable flagstat TSV filename for one sample chunk.
 */
String flagstatFileName(Integer batch_index, String sample_id) {
    return "${safeName(sample_id)}_${batch_index}.flagstat.tsv"
}

/**
 * Run chunk-level QC for each merged sample chunk.
 *
 * Each input is already sorted, merged, and indexed for exactly one sample and
 * one batch index. QC stays per chunk so live reports can be refreshed as new
 * sequencing chunks arrive instead of waiting for a final sample BAM.
 */
workflow quality_control {
    take:
        merged_bams: Channel<MergedChunkBAM>

    main:
        nanoplot_qc_ch = nanoplot_qc(merged_bams)

    emit:
        samtools_flagstat_qc(nanoplot_qc_ch)
}

/**
 * Run NanoPlot for one merged sample chunk and emit the compressed raw
 * `NanoPlot-data.tsv.gz` table used by downstream QC reporting.
 */
process nanoplot_qc {
    label 'seqLM_qc'
    container 'seqlm/quality_control'
    cpus 4

    input:
        merged_bam: MergedChunkBAM

    output:
        record(
            batch_index: merged_bam.batch_index,
            sample: merged_bam.sample,
            bam: merged_bam.bam,
            bam_index: merged_bam.bam_index,
            nanoplot_data: file("${nanoplotOutputDir(merged_bam)}/NanoPlot-data.tsv.gz")
        )

    script:
        String output_dir = nanoplotOutputDir(merged_bam)
        """
        NanoPlot -t ${task.cpus} \
            --raw \
            -o ${output_dir} \
            --no_static \
            --tsv_stats \
            --drop_outliers \
            --loglength \
            --title ${merged_bam.bam.name} \
            --bam ${merged_bam.bam}
        """
}

/**
 * Run samtools flagstat after NanoPlot for the same merged sample chunk and
 * emit the combined chunk QC record. This is intentionally chunk-level so
 * report inputs can update after every live batch.
 */
process samtools_flagstat_qc {
    label 'seqLM_qc'
    container 'seqlm/samtools'
    cpus 4

    input:
        nanoplot_result: NanoPlotQCResult

    output:
        record(
            batch_index: nanoplot_result.batch_index,
            sample: nanoplot_result.sample,
            bam: nanoplot_result.bam,
            bam_index: nanoplot_result.bam_index,
            nanoplot_data: nanoplot_result.nanoplot_data,
            flagstat: file(flagstatFileName(nanoplot_result.batch_index, nanoplot_result.sample.id))
        )

    script:
        String flagstat = flagstatFileName(nanoplot_result.batch_index, nanoplot_result.sample.id)
        """
        samtools flagstat -@ ${task.cpus} -O tsv ${nanoplot_result.bam} > ${flagstat}
        """
}
