#!/usr/bin/env nextflow

nextflow.enable.types = true

include { bam_ingress ; get_samples ; validate_samples } from './lib/bam_ingress.nf'
include {
    ChunkQCResult ;
    Sample ;
    SampleBatchSize ;
    SampleChunkBAMGroup
} from './lib/sample.nf'
include {
    chunk_bam_name ;
    optional_file ;
    shell_quote ;
    safe_name ;
    output ;
    publish_differential_results
} from './modules/generic_helpers.nf'
include {
    qc_report_copy_commands ;
    qc_report_root_dir ;
    accumulate_qc_report_chunk_state ;
    qc_report_inputs_from_state ;
    finalize_qc_report
} from './modules/qc_report_helpers.nf'
include { join_report_batches } from './modules/report_batches.nf'
include { quality_control } from './subworkflows/quality_control.nf'
include { quantification } from './subworkflows/quantification.nf'
include { differential_expression } from './subworkflows/differential_expression.nf'
include { order_batches } from './lib/util.nf'
include {
    differential_stability ;
    discover_stability_state ;
    stability_config
} from './lib/stability.nf'
include {
    discover_sample_checkpoints ;
    file_identity ;
    next_analysis_snapshot_index ;
    sample_checkpoint_key ;
    write_sample_checkpoints
} from './lib/sample_checkpoints.nf'
include { validate_parameters } from './lib/validation.nf'

process make_report {
    label 'seq_lm'

    input:
    metadata: Map
    per_read_stats: Path
    versions: Path
    params_json: Path

    stage:
    stageAs versions, 'versions/*'
    stageAs params_json, 'params.json'

    output:
    file('wf-template-*.html')

    script:
    def report_name: String = 'wf-template-report.html'
    def metadataJSON: String = new groovy.json.JsonBuilder(metadata).toPrettyString()
    def stats_args: String = (per_read_stats.name == optional_file().name) ? '' : "--stats ${per_read_stats}"
    """
        echo '${metadataJSON}' > metadata.json
        workflow-glue report ${report_name} \
            --versions versions \
            ${stats_args} \
            --params params.json \
            --metadata metadata.json
        """
}

process prepare_chunk_bam {
    label 'seq_lm'
    container 'rnabioinfo/seq_lm_samtools:v1.0.0'
    cpus 1

    input:
    input_group: SampleChunkBAMGroup

    stage:
    stageAs input_group.bams, 'bam?'

    output:
    record(
        batch_index: input_group.batch_index,
        sample: input_group.sample,
        bam: file(chunk_bam_name(input_group.batch_index, input_group.sample)),
    )

    script:
    def chunk_bam: String = chunk_bam_name(input_group.batch_index, input_group.sample)
    def chunk_bam_arg: String = shell_quote(chunk_bam)
    def bam_args: String = input_group.bams.collect { bam: Path -> shell_quote(bam.toString()) }.join(' ')
    def prepare_chunk: String = input_group.bams.size() == 1
        ? """
            ln -s -- ${bam_args} ${chunk_bam_arg}
        """
        : """
            : > bams.txt
            printf '%s\\n' ${bam_args} > bams.txt

            first_bam=1
            while IFS= read -r input_bam; do
                samtools view -H "\$input_bam" | awk '\$1 == "@SQ"' > current.sq
                if [ "\$first_bam" -eq 1 ]; then
                    cp current.sq expected.sq
                    first_bam=0
                elif ! cmp -s expected.sq current.sq; then
                    printf 'Incompatible BAM sequence dictionaries in sample %s batch %s: %s\\n' \
                        ${shell_quote(input_group.sample.name)} \
                        ${input_group.batch_index} \
                        "\$input_bam" >&2
                    exit 1
                fi
            done < bams.txt

            samtools cat -o ${chunk_bam_arg} -b bams.txt
        """
    """
        ${prepare_chunk}
        """
}

process qc_report_input_tree {
    debug true
    label 'seq_lm_qc'
    cpus 1
    maxForks 1

    input:
    qc_report_inputs: Map
    nanoplot_inputs: List<Path>
    flagstat_inputs: List<Path>

    stage:
    stageAs nanoplot_inputs, 'qc_report_sources/nanoplot/input?/NanoPlot-data.tsv.gz'
    stageAs flagstat_inputs, 'qc_report_sources/flagstat/input?/*'

    output:
    record(
        qc_report_inputs: qc_report_inputs,
        qc_results: file(qc_report_root_dir()),
    )

    script:
    def copy_commands: String = qc_report_copy_commands(
        qc_report_inputs.report_inputs_list,
        nanoplot_inputs,
        flagstat_inputs,
    )
    """
        ${copy_commands}
        """
}

/**
 * EPI2ME-displayable live QC and differential-analysis report.
 */
process qc_report {
    debug true
    label 'seq_lm_qc'
    container 'rnabioinfo/seq_lm_report:v1.0.0'
    cpus 1
    maxForks 1

    publishDir params.out_dir, mode: 'copy', pattern: 'qc_report*', overwrite: true

    input:
    qc_report_inputs: Map
    qc_results: Path
    transcript_biotypes: Path
    has_transcript_biotypes: Boolean
    differential_results: Path
    differential_analysis_note: String
    stability_results: Path
    has_differential_results: Boolean
    has_stability_results: Boolean

    stage:
    stageAs qc_results, 'qc_results'
    stageAs transcript_biotypes, 'transcript_biotypes.tsv'
    stageAs differential_results, 'differential_results'
    stageAs stability_results, 'stability_results.tsv'

    output:
    report_files = files('qc_report*', arity: '3')

    script:
    def rows: String = qc_report_inputs.rows
    def quoted_rows: String = shell_quote(rows)
    def params_json: String = new groovy.json.JsonBuilder(params).toPrettyString()
    def quoted_params_json: String = shell_quote(params_json)
    def differential_args: String = has_differential_results
        ? '--differential-results differential_results'
        : ''
    def transcript_biotype_args: String = has_transcript_biotypes
        ? '--transcript-biotypes transcript_biotypes.tsv'
        : ''
    def gene_set_args: String = params.gene_set_enrichment && has_differential_results ? '--gene-set-enrichment' : ''
    def readiness_args: String = differential_analysis_note
        ? "--dea-readiness-notice ${shell_quote(differential_analysis_note)}"
        : ''
    def stability_args: String = has_stability_results
        ? '--stability-results stability_results.tsv'
        : ''
    """
        printf 'name\\tgroup\\tchunks_seen\\tlatest_batch_index\\tqc_dir\\n' > report_samples.tsv
        printf '%s\\n' ${quoted_rows} >> report_samples.tsv

        mkdir versions
        printf 'qc_report,workflow\\n' > versions/versions.txt
        printf '%s\\n' ${quoted_params_json} > params.json

        export MPLCONFIGDIR="\$PWD/.matplotlib"
        mkdir -p "\$MPLCONFIGDIR"

        workflow-glue qc_report qc_report.html \
            --samples report_samples.tsv \
            --versions versions \
            --params params.json \
            --latest-batch ${qc_report_inputs.latest_batch_index} \
            ${transcript_biotype_args} \
            ${differential_args} \
            ${gene_set_args} \
            ${readiness_args} \
            ${stability_args} \
            --stability-behavior ${params.stability_analysis_behavior} \
            --lfc-cutoff ${params.de_lfc_cutoff} \
            --padj-cutoff ${params.de_padj_cutoff}
        """
}

// workflow module
workflow sample_pipeline {
    take:
    sample_batches: Channel
    restored_quantifications: List
    restored_qc_results: List
    fresh_sample_count: Integer
    first_analysis_index: Integer
    reference_genome: Path
    reference_annotation: Path
    gene_sets: Path
    quantification_enabled: Boolean
    differential_expression_enabled: Boolean
    gene_set_enrichment_enabled: Boolean
    de_lfc_cutoff: Number
    min_read_count: Integer
    min_replicate_sample_count: Integer
    all_samples: List<Sample>
    active_samples: List<Sample>
    stability_behavior: String
    stability_settings: Map
    initial_stability_state: Map
    minknow_connection: Map

    main:
    /*
         * `bam_ingress` emits synchronized sample batches. The pipeline spreads
         * each non-empty batch into one record per sample chunk, prepares a
         * sequential BAM for QC, and then refreshes the live QC report after
         * each complete batch. Keep this flow in channel operators; only use
         * Groovy helpers for naming and report-state formatting.
         */
    sample_batch_size_ch = sample_batches.map { batch ->
        record(
            batch_index: batch.batch_index,
            active_sample_count: batch.chunks.count { chunk -> !chunk.bam_paths.isEmpty() },
            experiment_sample_count: batch.experiment_sample_count,
        )
    }

    sample_chunk_bam_group_ch = sample_batches.flatMap { batch ->
        batch.chunks
            .findAll { chunk -> !chunk.bam_paths.isEmpty() }
            .collect { chunk ->
                record(
                    batch_index: batch.batch_index,
                    sample: chunk.sample,
                    bams: chunk.bam_paths,
                )
            }
    }

    merged_chunk_bam_ch = prepare_chunk_bam(sample_chunk_bam_group_ch)
    qc_result_ch = quality_control(merged_chunk_bam_ch)
    if (quantification_enabled) {
        quantification(
            merged_chunk_bam_ch,
            sample_batch_size_ch,
            restored_quantifications,
            reference_genome,
            reference_annotation,
        )
        quantified_samples_ch = quantification.out.quantifications
        quantified_sample_batches_ch = quantification.out.batches
        quantification_report_batches_ch = quantification.out.report_batches

        quantification_output_ch = quantified_samples_ch.map { quantified_sample ->
            tuple(
                quantified_sample.counts,
                "${safe_name(quantified_sample.sample.group)}/${safe_name(quantified_sample.sample.name)}/quantification",
            )
        }
        output(quantification_output_ch)

        if (differential_expression_enabled) {
            differential_expression(
                quantified_sample_batches_ch,
                first_analysis_index,
                reference_annotation,
                gene_sets,
                gene_set_enrichment_enabled,
                de_lfc_cutoff,
                min_read_count,
                min_replicate_sample_count,
            )
            differential_results_ch = differential_expression.out.results
            differential_report_batches_ch = differential_expression.out.report_batches

            publish_differential_results(
                differential_results_ch.map { result ->
                    tuple(result.batch_index, result.analysis_index, result.results)
                }
            )
            if (stability_behavior != 'disabled') {
                differential_stability(
                    differential_results_ch,
                    all_samples,
                    active_samples,
                    stability_behavior,
                    stability_settings,
                    initial_stability_state,
                    minknow_connection,
                )
                stability_audits_ch = differential_stability.out
                differential_reports_with_stability_ch = differential_report_batches_ch
                    .filter { report_batch -> report_batch.has_differential_results }
                    .map { report_batch -> tuple(report_batch.batch_index, report_batch) }
                    .join(
                        stability_audits_ch.map { audit ->
                            tuple(audit.batch_index, audit.sample_stability)
                        },
                        by: 0,
                    )
                    .map { _batch_index: Integer, report_batch, sample_stability_path ->
                        record(
                            batch_index: report_batch.batch_index,
                            report_sequence: report_batch.report_sequence,
                            differential_analysis_note: report_batch.differential_analysis_note,
                            has_differential_results: true,
                            results: report_batch.results,
                            stability_results: sample_stability_path,
                            has_stability_results: true,
                        )
                    }
                    .mix(
                        differential_report_batches_ch
                            .filter { report_batch -> !report_batch.has_differential_results }
                            .map { report_batch ->
                                record(
                                    batch_index: report_batch.batch_index,
                                    report_sequence: report_batch.report_sequence,
                                    differential_analysis_note: report_batch.differential_analysis_note,
                                    has_differential_results: false,
                                    results: report_batch.results,
                                    stability_results: optional_file(),
                                    has_stability_results: false,
                                )
                            }
                    )
            }
            else {
                differential_reports_with_stability_ch = differential_report_batches_ch.map { report_batch ->
                    record(
                        batch_index: report_batch.batch_index,
                        report_sequence: report_batch.report_sequence,
                        differential_analysis_note: report_batch.differential_analysis_note,
                        has_differential_results: report_batch.has_differential_results,
                        results: report_batch.results,
                        stability_results: optional_file(),
                        has_stability_results: false,
                    )
                }
            }

            analysis_report_batches_ch = quantification_report_batches_ch
                .map { report_batch -> tuple(report_batch.batch_index, report_batch) }
                .join(
                    differential_reports_with_stability_ch.map { report_batch ->
                        tuple(report_batch.batch_index, report_batch)
                    },
                    by: 0,
                )
                .map { batch_index: Integer, quantification_report, differential_report ->
                    if (quantification_report.report_sequence != differential_report.report_sequence) {
                        error("Inconsistent report sequence for batch ${batch_index}.")
                    }
                    record(
                        batch_index: batch_index,
                        report_sequence: quantification_report.report_sequence,
                        biotypes: quantification_report.biotypes,
                        differential_analysis_note: differential_report.differential_analysis_note,
                        has_differential_results: differential_report.has_differential_results,
                        results: differential_report.results,
                        stability_results: differential_report.stability_results,
                        has_stability_results: differential_report.has_stability_results,
                    )
                }
        }
        else {
            analysis_report_batches_ch = quantification_report_batches_ch.map { report_batch ->
                record(
                    batch_index: report_batch.batch_index,
                    report_sequence: report_batch.report_sequence,
                    biotypes: report_batch.biotypes,
                    differential_analysis_note: '',
                    has_differential_results: false,
                    results: optional_file(),
                    stability_results: optional_file(),
                    has_stability_results: false,
                )
            }
        }
    }

    if (quantification_enabled && fresh_sample_count > 0) {
        checkpoint_quantifications_ch = quantified_samples_ch
            .collect()
            .map { collected_quantifications ->
                collected_quantifications
                    .toList()
                    .toSorted { left, right ->
                        sample_checkpoint_key(left.sample) <=> sample_checkpoint_key(right.sample) ?: left.batch_index <=> right.batch_index
                    }
            }
        checkpoint_qc_results_ch = qc_result_ch
            .collect()
            .map { collected_qc_results ->
                collected_qc_results
                    .toList()
                    .toSorted { left, right ->
                        sample_checkpoint_key(left.sample) <=> sample_checkpoint_key(right.sample) ?: left.batch_index <=> right.batch_index
                    }
            }
        write_sample_checkpoints(
            checkpoint_quantifications_ch,
            checkpoint_qc_results_ch,
            [genome: file_identity(reference_genome), annotation: file_identity(reference_annotation)],
        )
    }

    def pending_qc_batches: Map<Integer, List<ChunkQCResult>> = [:]
    nonempty_qc_report_chunk_result_batches_ch = qc_result_ch
        .join(sample_batch_size_ch, by: 'batch_index')
        .map { joined ->
            def result: ChunkQCResult = record(
                batch_index: joined.batch_index,
                sample: joined.sample,
                bam: joined.bam,
                nanoplot_data: joined.nanoplot_data,
                flagstat: joined.flagstat,
            )
            def chunk_results: List<ChunkQCResult> = pending_qc_batches.computeIfAbsent(
                joined.batch_index
            ) { [] }
            chunk_results.add(result)
            if (chunk_results.size() < joined.active_sample_count) {
                return null
            }
            if (chunk_results.size() > joined.active_sample_count) {
                error(
                    "QC batch ${joined.batch_index} received more than " + "${joined.active_sample_count} result(s)."
                )
            }
            pending_qc_batches.remove(joined.batch_index)
            return record(
                batch_index: joined.batch_index,
                chunk_results: chunk_results,
            )
        }
        .filter { batch -> batch != null }

    empty_qc_report_chunk_result_batches_ch = sample_batch_size_ch
        .filter { batch_size -> batch_size.active_sample_count == 0 }
        .map { batch_size ->
            record(batch_index: batch_size.batch_index, chunk_results: [])
        }

    qc_report_chunk_result_batches_ch = order_batches(
        nonempty_qc_report_chunk_result_batches_ch.mix(empty_qc_report_chunk_result_batches_ch)
    )
    def qc_report_state: Map<String, List<ChunkQCResult>> = [:]
    restored_qc_results.each { qc_result ->
        def key: String = sample_checkpoint_key(qc_result.sample)
        def sample_results: List<ChunkQCResult> = qc_report_state.containsKey(key) ? qc_report_state[key] : []
        qc_report_state[key] = (sample_results + [qc_result]).toSorted { left, right ->
            left.batch_index <=> right.batch_index
        }
    }
    qc_report_inputs_ch = qc_report_chunk_result_batches_ch.flatMap { batch ->
        if (batch.chunk_results.empty) {
            return qc_report_state.empty
                ? []
                : [qc_report_inputs_from_state(
                    batch.batch_index,
                    qc_report_state,
                )]
        }
        return [accumulate_qc_report_chunk_state(
            qc_report_state,
            batch.batch_index,
            batch.chunk_results,
        )]
    }
    qc_report_metadata_ch = qc_report_inputs_ch.map { report_inputs: Map ->
        def report_metadata: Map = [latest_batch_index: report_inputs.latest_batch_index, report_inputs_list: report_inputs.report_inputs_list, rows: report_inputs.rows]
        report_metadata
    }
    qc_report_nanoplot_inputs_ch = qc_report_inputs_ch.map { report_inputs: Map ->
        report_inputs.nanoplot_inputs
    }
    qc_report_flagstat_inputs_ch = qc_report_inputs_ch.map { report_inputs: Map ->
        report_inputs.flagstat_inputs
    }
    qc_report_input_tree_ch = qc_report_input_tree(
        qc_report_metadata_ch,
        qc_report_nanoplot_inputs_ch,
        qc_report_flagstat_inputs_ch,
    )

    if (quantification_enabled) {
        qc_report_ready_ch = join_report_batches(
            analysis_report_batches_ch,
            qc_report_input_tree_ch,
        )
        qc_report(
            qc_report_ready_ch.map { result -> result.qc_report_inputs },
            qc_report_ready_ch.map { result -> result.qc_results },
            qc_report_ready_ch.map { result -> result.transcript_biotypes },
            true,
            qc_report_ready_ch.map { result -> result.differential_results },
            qc_report_ready_ch.map { result -> result.differential_analysis_note },
            qc_report_ready_ch.map { result -> result.stability_results },
            qc_report_ready_ch.map { result -> result.has_differential_results },
            qc_report_ready_ch.map { result -> result.has_stability_results },
        )
    }
    else {
        qc_report(
            qc_report_input_tree_ch.map { result -> result.qc_report_inputs },
            qc_report_input_tree_ch.map { result -> result.qc_results },
            optional_file(),
            false,
            optional_file(),
            '',
            optional_file(),
            false,
            false,
        )
    }
}

// Entrypoint workflow
workflow {
    main:
    WorkflowMain.initialise(workflow, params, log)
    validate_parameters(params)

    if (params.disable_ping == false) {
        Pinguscript.ping_post(workflow, 'start', 'none', params.out_dir, params)
    }

    reference_genome = params.reference_genome
        ? file(params.reference_genome, checkIfExists: true)
        : optional_file()
    reference_annotation = params.reference_annotation
        ? file(params.reference_annotation, checkIfExists: true)
        : optional_file()
    def quantification_enabled: Boolean = params.reference_genome != null &&
        "${params.reference_genome}".trim() != '' &&
        params.reference_annotation != null &&
        "${params.reference_annotation}".trim() != ''
    gene_sets = params.gene_set_enrichment
        ? file(params.gene_sets, checkIfExists: true)
        : optional_file()
    output_root = file(params.out_dir).toAbsolutePath().normalize()

    ingress_args = record(
        live_analysis: params.live_analysis,
        timeline_analysis: params.timeline_analysis,
        sample_sheet_path: params.sample_sheet ? file(params.sample_sheet) : null,
        bam_poll_interval_ms: (params.bam_poll_interval_seconds as Integer) * 1000,
        bam_stability_polls: params.bam_stability_polls as Integer,
        termination_requested: params.stability_analysis_behavior == 'terminate',
    )
    all_samples = get_samples(ingress_args)
    validate_samples(all_samples)
    checkpoint_state = quantification_enabled
        ? discover_sample_checkpoints(
            all_samples,
            output_root,
            reference_genome,
            reference_annotation,
        )
        : [restored: [], active: all_samples]
    restored_quantifications = checkpoint_state.restored*.quantification
    restored_qc_results = checkpoint_state.restored.collectMany { restored ->
        restored.qc_results
    }
    first_analysis_index = params.differential_expression
        ? next_analysis_snapshot_index(output_root)
        : 0
    stability_parameter_values = [
        behavior: params.stability_analysis_behavior as String,
        num_stable_batches: params.num_stable_batches as Integer,
        max_feature_diff_fraction: params.stability_max_feature_diff_fraction,
        max_median_abs_lfc_delta: params.stability_max_median_abs_lfc_delta,
        min_jaccard_similarity: params.stability_min_jaccard_similarity,
        max_call_churn_fraction: params.stability_max_call_churn_fraction,
        max_lost_call_fraction: params.stability_max_lost_call_fraction,
        max_fdr: params.de_padj_cutoff,
        min_abs_lfc: params.de_lfc_cutoff,
        min_de_calls_for_fraction_metrics: params.stability_min_de_calls_for_fraction_metrics as Integer,
        max_small_set_call_changes: params.stability_max_small_set_call_changes as Integer,
    ]
    stability_parameter_values.config = stability_config(stability_parameter_values, all_samples)
    initial_stability_state = params.stability_analysis_behavior == 'disabled'
        ? [previous_results: optional_file(), streaks: [:], eligible: [:]]
        : discover_stability_state(
            output_root,
            first_analysis_index,
            stability_parameter_values.config,
        )
    minknow_connection = [
        command: 'seq-run-manager',
        host: params.minknow_host as String,
        port: params.minknow_port as Integer,
        client_certificate: params.minknow_client_certificate
            ? file(params.minknow_client_certificate).toAbsolutePath().normalize()
            : null,
        client_private_key: params.minknow_client_private_key
            ? file(params.minknow_client_private_key).toAbsolutePath().normalize()
            : null,
        ca_certificate: params.minknow_ca_certificate
            ? file(params.minknow_ca_certificate).toAbsolutePath().normalize()
            : null,
    ]

    // Finalized samples are restored from the CLI output directory and never
    // enter BAM preparation, QC, collation, or Oarfish. Ingress polls only
    // samples without a valid FINAL checkpoint while retaining the full
    // experiment size for downstream analysis readiness.
    sample_batch_ch = bam_ingress(
        checkpoint_state.active,
        ingress_args,
        all_samples.size(),
    )
    sample_pipeline(
        sample_batch_ch,
        restored_quantifications,
        restored_qc_results,
        checkpoint_state.active.size(),
        first_analysis_index,
        reference_genome,
        reference_annotation,
        gene_sets,
        quantification_enabled,
        params.differential_expression,
        params.gene_set_enrichment,
        params.de_lfc_cutoff,
        params.min_read_count,
        params.min_replicate_sample_count,
        all_samples,
        checkpoint_state.active,
        params.stability_analysis_behavior,
        stability_parameter_values,
        initial_stability_state,
        minknow_connection,
    )

    onComplete:
    if (workflow.success) {
        finalize_qc_report(file(params.out_dir).toAbsolutePath().normalize())
    }
    if (params.disable_ping == false) {
        Pinguscript.ping_post(workflow, 'end', 'none', params.out_dir, params)
    }
}
