#!/usr/bin/env nextflow

nextflow.enable.types = true

include {
    bam_ingress ; capture_precision_input ;
    get_samples ;
    sample_sheet_has_order_column ;
    validate_samples ;
    validate_timeline_samples
} from './lib/bam_ingress.nf'
include {
    ChunkQCResult ;
    Sample ;
    GeneContrast ;
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
    finalize_qc_report ;
    publish_qc_report_snapshot ;
    next_qc_report_revision
} from './modules/qc_report_helpers.nf'
include { join_report_batches } from './modules/report_batches.nf'
include { join_ica_report_batches } from './modules/ica_report_batches.nf'
include { quality_control } from './subworkflows/quality_control.nf'
include { quantification } from './subworkflows/quantification.nf'
include { differential_expression } from './subworkflows/differential_expression.nf'
include { resolve_gene_contrasts ; gene_contrast_groups ; resolve_ica_groups ; select_analysis_batch } from './lib/gene_contrasts.nf'
include { order_batches } from './lib/util.nf'
include { precision_monitoring as dea_monitoring ; precision_monitoring as ica_monitoring } from './lib/precision_monitoring.nf'
include {
    discover_sample_checkpoints ;
    file_identity ;
    sha256_file ;
    next_analysis_snapshot_index ;
    sample_checkpoint_key ;
    write_sample_checkpoints
} from './lib/sample_checkpoints.nf'
include { imodulon_analysis } from './subworkflows/imodulon_analysis.nf'
include { next_ica_snapshot_index } from './lib/ica_publication.nf'
include { validate_parameters; ica_analysis_enabled } from './lib/validation.nf'

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
        observation: input_group.observation,
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
    fair true

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
    fair true

    input:
    qc_report_inputs: Map
    qc_results: Path
    transcript_biotypes: Path
    has_transcript_biotypes: Boolean
    differential_results: Path
    differential_analysis_note: String
    monitoring_results: Path
    ica_results: Path
    has_differential_results: Boolean
    gene_set_timecourse_enabled: Boolean
    has_monitoring_results: Boolean
    has_ica_results: Boolean
    ica_analysis_index: Integer
    first_report_revision: Integer

    stage:
    stageAs qc_results, 'qc_results'
    stageAs transcript_biotypes, 'transcript_biotypes.tsv'
    stageAs differential_results, 'differential_results'
    stageAs monitoring_results, 'monitoring_results'
    stageAs ica_results, 'ica_results'

    output:
    files('report{.html,_state.json,_snapshot_revision_*.html}', arity: '3')

    script:
    def rows: String = qc_report_inputs.rows
    def quoted_rows: String = shell_quote(rows)
    def report_revision: Integer = first_report_revision + (qc_report_inputs.report_sequence != null ? qc_report_inputs.report_sequence : qc_report_inputs.latest_batch_index) as Integer
    def params_json: String = new groovy.json.JsonBuilder(params).toPrettyString()
    def quoted_params_json: String = shell_quote(params_json)
    def differential_args: String = has_differential_results
        ? '--differential-results differential_results'
        : ''
    def transcript_biotype_args: String = has_transcript_biotypes
        ? '--transcript-biotypes transcript_biotypes.tsv'
        : ''
    def gene_set_args: String = params.gene_set_enrichment && has_differential_results ? '--gene-set-enrichment' : ''
    def temporal_args: String = gene_set_timecourse_enabled && has_differential_results
        ? '--temporal-analysis'
        : ''
    def readiness_args: String = differential_analysis_note
        ? "--dea-readiness-notice ${shell_quote(differential_analysis_note)}"
        : ''
    def stability_args: String = has_monitoring_results
        ? '--monitoring-results monitoring_results'
        : ''
    def ica_args: String = has_ica_results
        ? "--imodulon-results ica_results --imodulon-batch ${qc_report_inputs.latest_batch_index} --imodulon-sequence ${qc_report_inputs.report_sequence} --imodulon-analysis-index ${ica_analysis_index}"
        : ''
    """
        printf 'name\\tgroup\\torder\\tchunks_seen\\tlatest_batch_index\\tqc_dir\\n' > report_samples.tsv
        printf '%s\\n' ${quoted_rows} >> report_samples.tsv

        mkdir versions
        printf 'qc_report,workflow\\n' > versions/versions.txt
        printf '%s\\n' ${quoted_params_json} > params.json

        export MPLCONFIGDIR="\$PWD/.matplotlib"
        mkdir -p "\$MPLCONFIGDIR"

        workflow-glue qc_report report.html \
            --samples report_samples.tsv \
            --versions versions \
            --params params.json \
            --latest-batch ${qc_report_inputs.latest_batch_index} \
            --report-revision ${report_revision} \
            ${transcript_biotype_args} \
            ${differential_args} \
            ${gene_set_args} \
            ${temporal_args} \
            ${readiness_args} \
            ${stability_args} \
            ${ica_args} \
            --monitoring-behavior ${params.monitoring_behavior} \
            --stopping-method ${params.stopping_method} \
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
    strand_specificity: String
    gene_sets: Path
    dea_stopping_genes: Path
    gene_contrasts: List<GeneContrast>
    ica_groups: Set<String>
    quantification_enabled: Boolean
    differential_expression_enabled: Boolean
    gene_set_enrichment_enabled: Boolean
    gene_set_timecourse_enabled: Boolean
    de_lfc_cutoff: Number
    min_read_count: Integer
    min_replicate_sample_count: Integer
    all_samples: List<Sample>
    active_samples: List<Sample>
    stability_behavior: String
    stopping_settings: Map
    precision_observations: Map
    minknow_connection: Map
    ica_matrix: Path
    ica_gene_map: Path
    ica_imodulon_table: Path
    ica_settings: Map
    first_ica_index: Integer
    ica_output_root: Path

    main:
    def first_report_revision: Integer = next_qc_report_revision(ica_output_root)
    /*
         * `bam_ingress` emits synchronized sample batches. The pipeline spreads
         * each non-empty batch into one record per sample chunk, prepares a
         * sequential BAM for QC, and then refreshes the live QC report after
         * each complete batch. Keep this flow in channel operators; only use
         * Groovy helpers for naming and report-state formatting.
         */
    observed_sample_batches_ch = sample_batches.map { batch ->
        if (batch.get('observation') != null) capture_precision_input(precision_observations,batch.batch_index,batch.observation)
        return batch
    }
    sample_batch_size_ch = observed_sample_batches_ch.map { batch ->
        record(
            batch_index: batch.batch_index,
            active_sample_count: batch.chunks.count { chunk -> !chunk.bam_paths.isEmpty() },
            experiment_sample_count: batch.experiment_sample_count,
        )
    }

    sample_chunk_bam_group_ch = observed_sample_batches_ch.flatMap { batch ->
        batch.chunks
            .findAll { chunk -> !chunk.bam_paths.isEmpty() }
            .collect { chunk ->
                record(
                    batch_index: batch.batch_index,
                    sample: chunk.sample,
                    bams: chunk.bam_paths,
                    observation: batch.observation,
                )
            }
    }

    merged_chunk_bam_ch = prepare_chunk_bam(sample_chunk_bam_group_ch)
    qc_result_ch = quality_control(merged_chunk_bam_ch)
    
    if (quantification_enabled) {
        quantification_outputs = quantification(
            merged_chunk_bam_ch,
            sample_batch_size_ch,
            restored_quantifications,
            reference_genome,
            reference_annotation,
            strand_specificity,
            all_samples,
        )
        quantified_samples_ch = quantification_outputs.quantifications
        quantified_sample_batches_ch = quantification_outputs.batches
        quantification_report_batches_ch = quantification_outputs.report_batches

        quantification_output_ch = quantified_samples_ch.map { quantified_sample ->
            tuple(
                quantified_sample.counts,
                "${safe_name(quantified_sample.sample.group)}/${safe_name(quantified_sample.sample.name)}/quantification",
            )
        }
        output(quantification_output_ch)

        if (ica_settings.enabled) {
            ica_batches_ch = quantified_sample_batches_ch.map { batch -> select_analysis_batch(batch, ica_groups) }
            ica_snapshots_ch = imodulon_analysis(
                ica_batches_ch,
                ica_matrix,
                reference_annotation,
                ica_gene_map,
                ica_imodulon_table,
                ica_settings,
                first_ica_index,
                ica_output_root,
            )
        }

        if (differential_expression_enabled) {
            differential_outputs = differential_expression(
                quantified_sample_batches_ch.map { batch -> select_analysis_batch(batch, gene_contrast_groups(gene_contrasts)) },
                gene_contrasts,
                first_analysis_index,
                reference_annotation,
                gene_sets,
                dea_stopping_genes,
                gene_set_enrichment_enabled,
                de_lfc_cutoff,
                min_read_count,
                min_replicate_sample_count,
            )
            differential_results_ch = differential_outputs.results
            differential_report_batches_ch = differential_outputs.report_batches

            publish_differential_results(
                differential_results_ch.map { result ->
                    tuple(result.batch_index, result.analysis_index, result.results)
                }
            )
            if (stability_behavior != 'disabled' && stopping_settings.method == 'dea') {
                dea_snapshots_ch = differential_outputs.precision_results
                dea_audits_ch = dea_monitoring(
                    dea_snapshots_ch, all_samples, active_samples, stability_behavior,
                    stopping_settings, precision_observations, minknow_connection,
                )
                differential_reports_with_stability_ch = differential_report_batches_ch
                    .map { report_batch -> tuple(report_batch.batch_index, report_batch) }
                    .join(dea_audits_ch.map { audit -> tuple(audit.batch_index, audit.monitoring_results) }, by: 0)
                    .map { _batch_index: Integer, report_batch, monitoring_path ->
                        record(
                            batch_index: report_batch.batch_index,
                            report_sequence: report_batch.report_sequence,
                            differential_analysis_note: report_batch.differential_analysis_note,
                            has_differential_results: report_batch.has_differential_results,
                            results: report_batch.results,
                            monitoring_results: monitoring_path,
                            has_monitoring_results: true,
                        )
                    }
            }
            else {
                differential_reports_with_stability_ch = differential_report_batches_ch.map { report_batch ->
                    record(
                        batch_index: report_batch.batch_index,
                        report_sequence: report_batch.report_sequence,
                        differential_analysis_note: report_batch.differential_analysis_note,
                        has_differential_results: report_batch.has_differential_results,
                        results: report_batch.results,
                        monitoring_results: optional_file(),
                        has_monitoring_results: false,
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
                        monitoring_results: differential_report.monitoring_results,
                        has_monitoring_results: differential_report.has_monitoring_results,
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
                    monitoring_results: optional_file(),
                    has_monitoring_results: false,
                )
            }
        }

        if (ica_settings.enabled) {
            ica_report_batches_ch = join_ica_report_batches(
                analysis_report_batches_ch,
                ica_snapshots_ch,
            )
            analysis_report_batches_ch = ica_report_batches_ch
            if (stability_behavior != 'disabled' && stopping_settings.method == 'ica') {
                ica_audits_ch = ica_monitoring(
                    ica_snapshots_ch, all_samples, active_samples, stability_behavior,
                    stopping_settings, precision_observations, minknow_connection,
                )
                analysis_report_batches_ch = analysis_report_batches_ch
                    .map { report -> tuple(report.batch_index, report) }
                    .join(
                        ica_audits_ch.map { audit -> tuple(audit.batch_index, audit) },
                        by: 0, failOnDuplicate: true, failOnMismatch: true,
                    )
                    .map { batch_index: Integer, report, audit ->
                        if (report.ica_analysis_index != audit.analysis_index) {
                            error("Inconsistent ICA monitoring analysis index for batch ${batch_index}.")
                        }
                        record(
                            batch_index: batch_index, report_sequence: report.report_sequence,
                            biotypes: report.biotypes, differential_analysis_note: report.differential_analysis_note,
                            has_differential_results: report.has_differential_results, results: report.results,
                            monitoring_results: audit.monitoring_results, has_monitoring_results: true,
                            ica_results: report.ica_results, has_ica_results: true,
                            ica_analysis_index: report.ica_analysis_index,
                        )
                    }
            }
        }
        else {
            analysis_report_batches_ch = analysis_report_batches_ch.map { report ->
                record(
                    batch_index: report.batch_index,
                    report_sequence: report.report_sequence,
                    biotypes: report.biotypes,
                    differential_analysis_note: report.differential_analysis_note,
                    has_differential_results: report.has_differential_results,
                    results: report.results,
                    monitoring_results: report.monitoring_results,
                    has_monitoring_results: report.has_monitoring_results,
                    ica_results: optional_file(),
                    has_ica_results: false,
                    ica_analysis_index: 0,
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
            strand_specificity,
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
        qc_report_files_ch = qc_report(
            qc_report_ready_ch.map { result -> result.qc_report_inputs },
            qc_report_ready_ch.map { result -> result.qc_results },
            qc_report_ready_ch.map { result -> result.transcript_biotypes },
            true,
            qc_report_ready_ch.map { result -> result.differential_results },
            qc_report_ready_ch.map { result -> result.differential_analysis_note },
            qc_report_ready_ch.map { result -> result.monitoring_results },
            qc_report_ready_ch.map { result -> result.ica_results },
            qc_report_ready_ch.map { result -> result.has_differential_results },
            gene_set_timecourse_enabled,
            qc_report_ready_ch.map { result -> result.has_monitoring_results },
            qc_report_ready_ch.map { result -> result.has_ica_results },
            qc_report_ready_ch.map { result -> result.ica_analysis_index },
            first_report_revision,
        )
    }
    else {
        qc_report_files_ch = qc_report(
            qc_report_input_tree_ch.map { result -> result.qc_report_inputs },
            qc_report_input_tree_ch.map { result -> result.qc_results },
            optional_file(),
            false,
            optional_file(),
            '',
            optional_file(),
            optional_file(),
            false,
            false,
            false,
            false,
            0,
            first_report_revision,
        )
    }
    qc_report_files_ch
        .map { report_files: Collection<Path> ->
            def ordered_files: List<Path> = report_files.toList()
            publish_qc_report_snapshot(ordered_files, ica_output_root)
        }
        .collect()
}

// Entrypoint workflow
workflow {
    main:
    WorkflowMain.initialise(workflow, params, log)
    validate_parameters(params)
    def ica_enabled: Boolean = ica_analysis_enabled(params)
    def monitoring_enabled: Boolean = params.monitoring_behavior != 'disabled'
    def dea_monitoring_enabled: Boolean = monitoring_enabled && params.stopping_method == 'dea'
    def ica_monitoring_enabled: Boolean = monitoring_enabled && params.stopping_method == 'ica'
    def termination_dir: String = "termination/${params.stopping_method}"

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
    def dea_stopping_genes: Path = dea_monitoring_enabled && params.dea_stopping_genes
        ? file(params.dea_stopping_genes, checkIfExists: true)
        : optional_file()
    output_root = file(params.out_dir).toAbsolutePath().normalize()
    def sample_sheet_path: Path? = params.sample_sheet ? file(params.sample_sheet) : null
    def has_order_column: Boolean = (params.gene_set_enrichment || ica_enabled) && sample_sheet_path != null &&
        sample_sheet_has_order_column(sample_sheet_path)
    def gene_set_timecourse_enabled: Boolean = params.gene_set_enrichment && has_order_column
    def ica_timecourse_enabled: Boolean = ica_enabled && has_order_column

    ingress_args = record(
        live_analysis: params.live_analysis,
        sample_sheet_path: sample_sheet_path,
        bam_poll_interval_ms: (params.bam_poll_interval_seconds as Integer) * 1000,
        bam_stability_polls: params.bam_stability_polls as Integer,
        termination_requested: params.monitoring_behavior == 'terminate',
        detect_minknow_completion: params.detect_minknow_completion,
        minknow_completion_settle_seconds: params.minknow_completion_settle_seconds,
        precision_event_mode: params.live_analysis && monitoring_enabled,
        precision_state_directory: params.live_analysis && monitoring_enabled
            ? output_root.resolve(termination_dir).toString() : null,
    )
    all_samples = get_samples(ingress_args)
    validate_samples(all_samples, params.differential_expression && !params.gene_contrast_file && !ica_enabled)
    def gene_contrasts: List<GeneContrast> = params.differential_expression || params.gene_contrast_file
        ? resolve_gene_contrasts(all_samples, params.gene_contrast_file ? file(params.gene_contrast_file) : null)
        : []
    def ica_groups: Set<String> = ica_enabled
        ? resolve_ica_groups(all_samples, params.ica_control_group as String, params.ica_condition_groups as String)
        : [] as Set<String>
    if (gene_set_timecourse_enabled) {
        def contrast_groups: Set<String> = gene_contrast_groups(gene_contrasts)
        validate_timeline_samples(all_samples.findAll { sample -> contrast_groups.contains(sample.group) })
    }
    if (ica_timecourse_enabled) {
        validate_timeline_samples(all_samples.findAll { sample -> ica_groups.contains(sample.group) })
    }
    checkpoint_state = quantification_enabled
        ? discover_sample_checkpoints(
            all_samples,
            output_root,
            reference_genome,
            reference_annotation,
            params.strand_specificity,
        )
        : [restored: [], active: all_samples]
    restored_quantifications = checkpoint_state.restored*.quantification
    restored_qc_results = checkpoint_state.restored.collectMany { restored ->
        restored.qc_results
    }
    first_analysis_index = params.differential_expression
        ? next_analysis_snapshot_index(output_root)
        : 0
    def ica_control: String = (params.ica_control_group as String).trim()
    def stopping_contrasts: List<GeneContrast> = gene_contrasts
    if (ica_monitoring_enabled) {
        def targets: List<String> = ica_groups.findAll { group: String -> group != ica_control }.toList().toSorted()
        if (targets.empty) {
            error('ICA stopping requires at least one condition group besides the ICA control.')
        }
        stopping_contrasts = targets.withIndex().collect { target: String, index: Integer ->
            record(contrast_id: "contrast_${index + 1}" as String, target_group: target, reference_group: ica_control)
        }
    }

    def precision_sources: List<String> = [
        'bin/precision_termination/engine.py',
        'bin/precision_termination/forecast.py',
        'bin/precision_termination/cli.py',
        'bin/workflow_glue/readiness.py',
        'lib/LiveBamWatcher.groovy',
        'lib/bam_ingress.nf',
        'lib/precision_monitoring.nf',
        'lib/acquisition_actions.nf',
        'subworkflows/quantification.nf',
    ]
    def method_sources: List<String> = ica_monitoring_enabled
        ? ['bin/precision_termination/ica.py', 'bin/imodulon_analysis/analysis.py', 'subworkflows/imodulon_analysis.nf']
        : ['bin/de_analysis/profile_intervals.R', 'bin/de_analysis/edgeR_analysis.R',
           'bin/de_analysis/read_gmt.R', 'bin/check-de-readiness', 'subworkflows/differential_expression.nf']
    def algorithm_identity: Map = (precision_sources + method_sources).collectEntries { source: String ->
        [(source): sha256_file(projectDir.resolve(source))]
    }

    stopping_settings = [
        contrasts: stopping_contrasts,
        method: params.stopping_method as String,
        live_analysis: params.live_analysis,
        output_root: output_root.toString(),
        command: projectDir.resolve('bin/precision-termination').toString(),
        container_engine: workflow.containerEngine,
        identities: [
            model_settings: [
                min_read_count: params.min_read_count,
                min_replicate_sample_count: params.min_replicate_sample_count,
                de_lfc_cutoff: params.de_lfc_cutoff,
                de_padj_cutoff: params.de_padj_cutoff,
                strand_specificity: params.strand_specificity,
            ],
            genome: reference_genome ? sha256_file(reference_genome) : null,
            annotation: reference_annotation ? sha256_file(reference_annotation) : null,
            dea_stopping_genes: dea_stopping_genes.name != optional_file().name ? sha256_file(dea_stopping_genes) : null,
            controls: restored_quantifications.collect { quant ->
                [key: "${quant.sample.group}/${quant.sample.name}", sha256: sha256_file(quant.counts)]
            },
            algorithm: algorithm_identity,
        ],
        precision_config: monitoring_enabled ? [
            gain_cutoff: params.stopping_gain_cutoff,
            information_radius_scale: params.stopping_information_radius_scale,
            horizon_seconds: (params.stopping_forecast_minutes as Number) * 60,
            max_processing_lag_seconds: params.stopping_max_processing_lag_seconds,
            min_points: params.stopping_min_points,
            earliest_seconds: params.stopping_min_seconds,
            confirmations: params.stopping_confirmations,
            interruption_fraction: params.stopping_interruption_fraction,
            recovery_fraction: params.stopping_recovery_fraction,
            recovery_intervals: params.stopping_recovery_intervals,
        ] : [:],
    ]
    precision_observations = [:]
    minknow_connection = [
        command: 'seq-run-manager',
        container_engine: workflow.containerEngine,
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
        params.strand_specificity,
        gene_sets,
        dea_stopping_genes,
        gene_contrasts,
        ica_groups,
        quantification_enabled,
        params.differential_expression,
        params.gene_set_enrichment,
        gene_set_timecourse_enabled,
        params.de_lfc_cutoff,
        params.min_read_count,
        params.min_replicate_sample_count,
        all_samples,
        checkpoint_state.active,
        params.monitoring_behavior,
        stopping_settings,
        precision_observations,
        minknow_connection,
        ica_enabled ? file(params.ica_matrix, checkIfExists: true) : optional_file(),
        ica_enabled && params.ica_gene_map ? file(params.ica_gene_map, checkIfExists: true) : optional_file(),
        ica_enabled && params.ica_imodulon_table ? file(params.ica_imodulon_table, checkIfExists: true) : optional_file(),
        [enabled: ica_enabled,
         monitoring: ica_monitoring_enabled,
         stopping_components: ica_monitoring_enabled && params.ica_stopping_components
             ? file(params.ica_stopping_components, checkIfExists: true) : null,
         has_imodulon_table: params.ica_imodulon_table != null && "${params.ica_imodulon_table}".trim() != '',
         has_gene_map: params.ica_gene_map != null && "${params.ica_gene_map}".trim() != '',
         control_group: ica_control,
         condition_groups: ica_groups.toList().toSorted(),
         log_base: params.ica_log_base, pseudocount: params.ica_pseudocount,
         min_gene_coverage: params.ica_min_gene_coverage,
         min_read_count: params.min_read_count, min_replicate_sample_count: params.min_replicate_sample_count,
         padj_cutoff: params.ica_padj_cutoff],
        ica_enabled ? next_ica_snapshot_index(output_root) : 0,
        output_root,
    )

    onComplete:
    if (workflow.success) {
        finalize_qc_report(file(params.out_dir).toAbsolutePath().normalize())
    }
    if (params.disable_ping == false) {
        Pinguscript.ping_post(workflow, 'end', 'none', params.out_dir, params)
    }
}
