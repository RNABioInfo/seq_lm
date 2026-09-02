#!/usr/bin/env nextflow

nextflow.enable.types = true

include { DifferentialExpressionBatch ; QuantifiedSampleBatch } from '../lib/sample.nf'

include {
    differential_expression_results_dir ;
    optional_file ;
    shell_quote
} from '../modules/generic_helpers.nf'

/** Run differential expression for every complete quantification batch. */
workflow differential_expression {
    take:
    quantified_batches: Channel<QuantifiedSampleBatch>
    first_analysis_index: Integer
    annotation: Path
    gene_sets: Path
    gene_set_enrichment: Boolean
    de_lfc_cutoff: Number
    min_read_count: Integer
    min_replicate_sample_count: Integer

    main:
    def next_analysis_index: Integer = first_analysis_index
    differential_batches_ch = quantified_batches.map { quant_batch ->
        def differential_batch: DifferentialExpressionBatch = record(
            batch_index: quant_batch.batch_index,
            analysis_index: next_analysis_index,
            report_sequence: quant_batch.report_sequence,
            samples: quant_batch.samples,
        )
        next_analysis_index += 1
        return differential_batch
    }

    checked_quantified_sample_batches_ch = check_differential_analysis_readiness(
        differential_batches_ch,
        min_read_count,
        min_replicate_sample_count,
    )
    readiness_status_ch = checked_quantified_sample_batches_ch.map { checked ->
        def readiness_text: String = checked.readiness.text.trim().toLowerCase()
        def readiness_notices: Map<String, String> = [
            insufficient_read_depth: 'For DEA, the required read depth is not yet satisfied.',
            no_matching_feature_ids: 'For DEA, the sample quantifications contain no matching feature IDs.',
        ]
        if (readiness_text != 'ready' && !readiness_notices.containsKey(readiness_text)) {
            error(
                "Unexpected differential-analysis readiness result for batch " +
                "${checked.batch_index}: '${readiness_text}'."
            )
        }
        tuple(checked.batch_index, readiness_text == 'ready', readiness_notices[readiness_text] ?: '')
    }
    readiness_decisions_ch = differential_batches_ch
        .map { quant_batch -> tuple(quant_batch.batch_index, quant_batch) }
        .join(readiness_status_ch, by: 0)
        .map { _batch_index: Integer, quant_batch, differential_analysis_ready: Boolean, differential_analysis_note: String ->
            record(
                quant_batch: quant_batch,
                differential_analysis_ready: differential_analysis_ready,
                differential_analysis_note: differential_analysis_note,
            )
        }
    ready_quantified_sample_batches_ch = readiness_decisions_ch
        .filter { decision -> decision.differential_analysis_ready }
        .map { decision -> decision.quant_batch }
    deferred_report_batches_ch = readiness_decisions_ch
        .filter { decision -> !decision.differential_analysis_ready }
        .map { decision ->
            record(
                batch_index: decision.quant_batch.batch_index,
                report_sequence: decision.quant_batch.report_sequence,
                differential_analysis_note: decision.differential_analysis_note,
                has_differential_results: false,
                results: optional_file(),
            )
        }

    edgeR_results_ch = run_differential_expression_edgeR(
        ready_quantified_sample_batches_ch,
        gene_sets,
        annotation,
        de_lfc_cutoff,
    )
    if (gene_set_enrichment) {
        analysis_results_ch = run_gene_set_variation_analysis(edgeR_results_ch)
    }
    else {
        analysis_results_ch = edgeR_results_ch
    }
    completed_report_batches_ch = analysis_results_ch.map { result ->
        record(
            batch_index: result.batch_index,
            report_sequence: result.report_sequence,
            differential_analysis_note: '',
            has_differential_results: true,
            results: result.results,
        )
    }

    emit:
    results = analysis_results_ch
    report_batches = completed_report_batches_ch.mix(deferred_report_batches_ch)
}

/** Check read depth and shared feature IDs before launching edgeR. */
process check_differential_analysis_readiness {
    label 'seq_lm_qc'
    container 'rnabioinfo/seq_lm_report:v1.0.0'
    cpus 1

    input:
    quant_batch: DifferentialExpressionBatch
    min_read_count: Integer
    min_replicate_sample_count: Integer

    stage:
    stageAs quant_batch.samples*.counts, 'quant/input?.quant'

    output:
    record(
        batch_index: quant_batch.batch_index,
        readiness: file('differential_analysis_ready.txt'),
    )

    script:
    def manifest_rows: String = quant_batch.samples
        .withIndex()
        .collect { sample, index: Integer ->
            [de_manifest_field(sample.sample.name), de_manifest_field(sample.sample.group), "input${index + 1}.quant"].join('\t')
        }
        .join('\n')
    def quoted_manifest_rows: String = shell_quote(manifest_rows)
    """
        printf 'name\\tgroup\\tcount_file\\n' > quant_manifest.tsv
        printf '%s\\n' ${quoted_manifest_rows} >> quant_manifest.tsv

        stability_read_count \
            --min_read_count ${min_read_count} \
            --min_replicate_sample_count ${min_replicate_sample_count} \
            --metadata quant_manifest.tsv \
            --counts_dir quant \
            > differential_analysis_ready.txt
        """
}

/** Rebuild the full count matrix and rerun edgeR for one live batch. */
process run_differential_expression_edgeR {
    label 'seq_lm_dea'
    container 'rnabioinfo/seq_lm_dea_r:v1.0.0'
    cpus 1
    maxForks 1

    input:
    quant_batch: DifferentialExpressionBatch
    gene_sets: Path
    annotation: Path
    de_lfc_cutoff: Number

    stage:
    stageAs quant_batch.samples*.counts, 'quant/input?.quant'

    output:
    record(
        batch_index: quant_batch.batch_index,
        analysis_index: quant_batch.analysis_index,
        report_sequence: quant_batch.report_sequence,
        results: file(differential_expression_results_dir(quant_batch.analysis_index)),
    )

    script:
    def results_dir: String = differential_expression_results_dir(quant_batch.analysis_index)
    def gene_set_args: String = gene_sets.name == optional_file().name ? '' : "--gene_sets ${gene_sets}"
    def manifest_rows: String = quant_batch.samples
        .withIndex()
        .collect { sample, index: Integer ->
            [de_manifest_field(sample.sample.name), de_manifest_field(sample.sample.group), "quant/input${index + 1}.quant"].join('\t')
        }
        .join('\n')
    def quoted_manifest_rows: String = shell_quote(manifest_rows)
    """
        printf 'name\\tgroup\\tcount_file\\n' > quant_manifest.tsv
        printf '%s\\n' ${quoted_manifest_rows} >> quant_manifest.tsv

        mkdir ${results_dir}

        edgeR-analysis \
            --quant_manifest quant_manifest.tsv \
            --output_dir ${results_dir} \
            ${gene_set_args} \
            --annotation ${annotation} \
            --lfc ${de_lfc_cutoff}
        """
}

/** Add GSVA outputs to one complete edgeR result tree. */
process run_gene_set_variation_analysis {
    label 'seq_lm_gsva'
    container 'rnabioinfo/seq_lm_gsva:v1.1.0'
    cpus 1
    maxForks 1

    input:
    differential_result: Map

    stage:
    stageAs differential_result.results, 'edgeR_results'

    output:
    record(
        batch_index: differential_result.batch_index,
        analysis_index: differential_result.analysis_index,
        report_sequence: differential_result.report_sequence,
        results: file(differential_expression_results_dir(differential_result.analysis_index)),
    )

    script:
    def results_dir: String = differential_expression_results_dir(differential_result.analysis_index)
    """
        mkdir ${results_dir}
        cp -R edgeR_results/. ${results_dir}/

        gsva-analysis \
            --feature_counts ${results_dir}/feature_counts.tsv \
            --sample_metadata ${results_dir}/sample_metadata.tsv \
            --gene_set_resolution ${results_dir}/gene_set_resolution.tsv \
            --output_dir ${results_dir}
        """
}

def de_manifest_field(value: Object) -> String {
    return "${value}".replaceAll(/[\t\r\n]+/, ' ')
}
