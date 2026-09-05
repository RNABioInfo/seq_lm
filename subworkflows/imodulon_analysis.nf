#!/usr/bin/env nextflow
nextflow.enable.types = true

include { QuantifiedSampleBatch ; ICABatchResult } from '../lib/sample.nf'
include { shell_quote } from '../modules/generic_helpers.nf'
include { publish_ica_snapshot } from '../lib/ica_publication.nf'

workflow imodulon_analysis {
    take:
    batches: Channel<QuantifiedSampleBatch>
    matrix: Path
    annotation: Path
    gene_map: Path
    settings: Map
    first_index: Integer
    output_root: Path

    main:
    prepared = prepare_ica_model(matrix, annotation, gene_map, settings)
    def results: Channel<ICABatchResult> = infer_ica_activities(batches, prepared, settings, first_index)
    // A synchronous channel operator serializes filesystem publication. Native
    // process submission can reorder side effects even with fair/maxForks.
    published = results.map { result ->
        def destination: Path = publish_ica_snapshot(result, output_root)
        record(
            batch_index: result.batch_index,
            report_sequence: result.report_sequence,
            analysis_index: result.analysis_index,
            results: destination,
        )
    }

    emit:
    snapshots = published
}

process prepare_ica_model {
    label 'seq_lm_ica'
    container 'rnabioinfo/seq_lm_ica:v1.0.0'
    cpus 1

    input:
    matrix: Path
    annotation: Path
    gene_map: Path
    settings: Map

    stage:
    stageAs matrix, 'model_weights/input.csv'
    stageAs annotation, 'annotation/input.gtf'
    stageAs gene_map, 'gene_map/input.tsv'

    output:
    file('prepared_ica')

    script:
    def map_args: String = settings.has_gene_map ? '--gene-map gene_map/input.tsv' : ''
    """
    imodulon-analysis prepare \\
        --matrix model_weights/input.csv \\
        --annotation annotation/input.gtf \\
        ${map_args} \\
        --min-gene-coverage ${settings.min_gene_coverage} \\
        --output prepared_ica
    """
}

process infer_ica_activities {
    label 'seq_lm_ica'
    container 'rnabioinfo/seq_lm_ica:v1.0.0'
    cpus 1
    // Preserve input report-sequence order even when computations finish out of order.
    fair true

    input:
    batch: QuantifiedSampleBatch
    prepared: Path
    settings: Map
    first_index: Integer

    stage:
    stageAs batch.samples*.counts, 'quant/input?.quant'
    stageAs prepared, 'prepared_ica'

    output:
    record(
        batch_index: batch.batch_index,
        report_sequence: batch.report_sequence,
        analysis_index: first_index + batch.report_sequence,
        results: file('ica_results'),
    )

    script:
    def rows: String = batch.samples.withIndex().collect { sample, index: Integer ->
        [sample.sample.name, sample.sample.group, sample.sample.order == null ? '' : sample.sample.order,
         "input${index + 1}.quant", sample.batch_index].join('\t')
    }.join('\n')
    """
    printf 'name\\tgroup\\torder\\tcount_file\\tsource_batch_index\\n' > quant_manifest.tsv
    printf '%s\\n' ${shell_quote(rows)} >> quant_manifest.tsv
    imodulon-analysis analyze \\
        --prepared prepared_ica \\
        --manifest quant_manifest.tsv \\
        --counts-dir quant \\
        --log-base ${settings.log_base} \\
        --pseudocount ${settings.pseudocount} \\
        --min-reads ${settings.min_read_count} \\
        --cutoff ${settings.padj_cutoff} \\
        --batch-index ${batch.batch_index} \\
        --report-sequence ${batch.report_sequence} \\
        --analysis-index ${first_index + batch.report_sequence} \\
        --output ica_results
    """
}

