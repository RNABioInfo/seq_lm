#!/usr/bin/env nextflow

nextflow.enable.types = true

include { bam_ingress } from './lib/bam_ingress.nf'
include {
    ChunkQCResult;
    Sample;
    SampleBatchSize;
    SampleChunkBAMGroup
} from './lib/sample.nf'
include {
    chunk_bam_name;
    optional_file;
    shell_quote
} from './modules/generic_helpers.nf'
include { output } from './modules/generic_helpers.nf'
include { qc_report_copy_commands; qc_report_root_dir; accumulate_qc_report_chunk_state } from './modules/qc_report_helpers.nf'
include { join_report_batches } from './modules/report_batches.nf'
include { quality_control } from './subworkflows/quality_control.nf'
include { differential_expression } from './subworkflows/differential_expression.nf'
include { order_batches } from './lib/util.nf'


process write_config {
    label 'seq_lm'
    cpus 1
    exec:
        log.info 'Writing config file...'

        // Writing experiment config file should only happen at the first run of the experiment
        def configOut = new File("${params.out_dir}/experiment.config")

        configOut.withWriter { w ->
            w << 'params {\n'
            params.each { k, v ->
                if (k.startsWith('ex')) {
                    if (k == 'ex_run_number') {
                        v = v + 1
                    }
                    def line = ''
                    if (v instanceof String) {
                        line = "\t${k} = \"${v}\"\n"
                    } else {
                        line = "\t${k} = ${v}\n"
                    }
                    w << line
                }
            }
            w << '}\n'
        }
}

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
        String report_name = 'wf-template-report.html'
        String metadataJSON = new groovy.json.JsonBuilder(metadata).toPrettyString()
        String stats_args = \
            (per_read_stats.name == optional_file().name) ? '' : "--stats $per_read_stats"
        """
        echo '${metadataJSON}' > metadata.json
        workflow-glue report $report_name \
            --versions versions \
            $stats_args \
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
            bam: file(chunk_bam_name(input_group.batch_index, input_group.sample))
        )

    script:
        String chunk_bam = chunk_bam_name(input_group.batch_index, input_group.sample)
        String chunk_bam_arg = shell_quote(chunk_bam)
        String bam_args = input_group.bams.collect { Path bam -> shell_quote(bam.toString()) }.join(' ')
        String prepare_chunk = input_group.bams.size() == 1 ? """
            ln -s -- ${bam_args} ${chunk_bam_arg}
        """ : """
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
            qc_results: file(qc_report_root_dir())
        )

    script:
        String copy_commands = qc_report_copy_commands(
            qc_report_inputs.report_inputs_list,
            nanoplot_inputs,
            flagstat_inputs
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

    publishDir(
        params.out_dir,
        mode: 'copy',
        pattern: 'qc_report*',
        overwrite: true
    )

    input:
        qc_report_inputs: Map
        qc_results: Path
        differential_results: Path

    stage:
        stageAs qc_results, 'qc_results'
        stageAs differential_results, 'differential_results'

    output:
        report_files = files('qc_report*', arity: '3')

    script:
        String rows = qc_report_inputs.rows
        String quoted_rows = shell_quote(rows)
        String params_json = new groovy.json.JsonBuilder(params).toPrettyString()
        String quoted_params_json = shell_quote(params_json)
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
            --differential-results differential_results \
            --lfc-cutoff ${params.de_lfc_cutoff} \
            --padj-cutoff ${params.de_padj_cutoff}
        """
}

// workflow module
workflow sample_pipeline {
    take:
        sample_batches: Channel
        reference_genome: Path
        reference_annotation: Path
        gene_sets: Path
    main:
        /*
         * `bam_ingress` emits synchronized sample batches. The pipeline spreads
         * each non-empty batch into one record per sample chunk, prepares a
         * sequential BAM for QC, and then refreshes the live QC report after
         * each complete batch. Keep this flow in channel operators; only use
         * Groovy helpers for naming and report-state formatting.
         */
        sample_batch_size_ch = sample_batches
            .map { batch ->
                record(
                    batch_index: batch.batch_index,
                    active_sample_count: batch.chunks.count { chunk -> !chunk.bam_paths.isEmpty() },
                    experiment_sample_count: batch.experiment_sample_count
                )
            }

        sample_chunk_bam_group_ch = sample_batches
            .flatMap { batch ->
                batch.chunks
                    .findAll { chunk -> !chunk.bam_paths.isEmpty() }
                    .collect { chunk ->
                        record(
                            batch_index: batch.batch_index,
                            sample: chunk.sample,
                            bams: chunk.bam_paths
                        )
                    }
            }

        merged_chunk_bam_ch = prepare_chunk_bam(sample_chunk_bam_group_ch)
        qc_result_ch = quality_control(merged_chunk_bam_ch)
        differential_expression(
            merged_chunk_bam_ch,
            sample_batch_size_ch,
            reference_genome,
            reference_annotation,
            gene_sets
        )

        differential_output_ch = differential_expression.out.quantifications
            .map { quantified_sample ->
                tuple(
                    quantified_sample.counts,
                    "${quantified_sample.sample.group}/${quantified_sample.sample.name}/quantification"
                )
            }
            .mix(
                differential_expression.out.results.map { differential_result ->
                    tuple(differential_result.results, 'differential_expression')
                }
            )
        output(differential_output_ch)

        Map<Integer, List<ChunkQCResult>> pending_qc_batches = [:]
        nonempty_qc_report_chunk_result_batches_ch = qc_result_ch
            .join(sample_batch_size_ch, by: 'batch_index')
            .map { joined ->
                ChunkQCResult result = record(
                    batch_index: joined.batch_index,
                    sample: joined.sample,
                    bam: joined.bam,
                    nanoplot_data: joined.nanoplot_data,
                    flagstat: joined.flagstat
                )
                List<ChunkQCResult> chunk_results = pending_qc_batches.computeIfAbsent(
                    joined.batch_index
                ) { [] }
                chunk_results.add(result)
                if (chunk_results.size() < joined.active_sample_count) {
                    return null
                }
                if (chunk_results.size() > joined.active_sample_count) {
                    error(
                        "QC batch ${joined.batch_index} received more than " +
                        "${joined.active_sample_count} result(s)."
                    )
                }
                pending_qc_batches.remove(joined.batch_index)
                return record(
                    batch_index: joined.batch_index,
                    chunk_results: chunk_results
                )
            }
            .filter { batch -> batch != null }

        empty_qc_report_chunk_result_batches_ch = sample_batch_size_ch
            .filter { SampleBatchSize batch_size -> batch_size.active_sample_count == 0 }
            .map { SampleBatchSize batch_size ->
                record(batch_index: batch_size.batch_index, chunk_results: [])
            }

        qc_report_chunk_result_batches_ch = order_batches(
            nonempty_qc_report_chunk_result_batches_ch
                .mix(empty_qc_report_chunk_result_batches_ch)
        )
        Map<String, List<ChunkQCResult>> qc_report_state = [:]
        qc_report_inputs_ch = qc_report_chunk_result_batches_ch.flatMap { batch ->
                if (batch.chunk_results.empty) {
                    return []
                }
                return [accumulate_qc_report_chunk_state(
                    qc_report_state,
                    batch.batch_index,
                    batch.chunk_results
                )]
            }
        qc_report_metadata_ch = qc_report_inputs_ch.map { Map report_inputs ->
                Map report_metadata = [
                    latest_batch_index: report_inputs.latest_batch_index,
                    report_inputs_list: report_inputs.report_inputs_list,
                    rows: report_inputs.rows
                ]
                report_metadata
            }
        qc_report_nanoplot_inputs_ch = qc_report_inputs_ch.map { Map report_inputs ->
            report_inputs.nanoplot_inputs
        }
        qc_report_flagstat_inputs_ch = qc_report_inputs_ch.map { Map report_inputs ->
            report_inputs.flagstat_inputs
        }
        qc_report_input_tree_ch = qc_report_input_tree(
            qc_report_metadata_ch,
            qc_report_nanoplot_inputs_ch,
            qc_report_flagstat_inputs_ch
        )

        qc_report_ready_ch = join_report_batches(
            differential_expression.out.results,
            qc_report_input_tree_ch
        )
        qc_report(
            qc_report_ready_ch.map { result -> result.qc_report_inputs },
            qc_report_ready_ch.map { result -> result.qc_results },
            qc_report_ready_ch.map { result -> result.differential_results }
        )
}

Map prepare_run(String experiment_dir, Integer _run_number, Integer replicate_count) {
    def runName = "run_${params.ex_run_number}"
    def runDir = file("${params.out_dir}/${runName}")

    def metadataFile = new File("${experiment_dir}/metadata.tsv")

    // Add header to metadata file if file is empty
    if (metadataFile.length() == 0) {
        metadataFile.withWriter { w ->
            w << 'run_number\treplicate_number\treplicate_dir\n'
        }
    }

    // Create run directories for each replicate
    (1..replicate_count).each { count ->
        def replicateName = "replicate_${count}"
        def replicateDir = file("${runDir}/${replicateName}")
        replicateDir.mkdirs()

        metadataFile << "${params.ex_run_number}\t${count}\t${replicateDir}\n"
    }

    return [runName: runName, runDir: runDir]
}

// Entrypoint workflow
workflow {
    main:
    WorkflowMain.initialise(workflow, params, log)

    if (params.de_lfc_cutoff < 0) {
        error('--de_lfc_cutoff must be nonnegative.')
    }
    if (params.de_padj_cutoff <= 0 || params.de_padj_cutoff > 1) {
        error('--de_padj_cutoff must be greater than 0 and at most 1.')
    }

    if (params.disable_ping == false) {
        Pinguscript.ping_post(workflow, 'start', 'none', params.out_dir, params)
    }

    // TODO: Implement parameter validation
    // validate_experiment_dir(params.out_dir, params.ex_run_number)

    // TODO: Implement sequencing setup checks

    // Config is stored in order to fetch parameters in subsequent runs
    // write_config()

    // // Setup the run
    // runInfo = prepare_run(params.out_dir, params.ex_run_number, params.ex_replicate_count)
    // runDir = runInfo.runDir

    // // Start the sequencing run
    // metadataFile = channel.fromPath("${params.out_dir}/metadata.tsv")
    // keyFile = channel.fromPath(params.ex_mk_key)
    // certificateFile = channel.fromPath(params.ex_mk_cert)
    // sequencingArgs = get_sequencing_arguments(runDir)
    // startSequencing(sequencingArgs, keyFile, certificateFile, metadataFile)

    // The ingress emits synchronized sample batches. QC remains chunk-local;
    // differential-expression quantification accumulates every sample chunk.
    sample_batch_ch = bam_ingress(
        record(
            live_analysis: params.live_analysis,
            timeline_analysis: params.timeline_analysis,
            sample_sheet_path: params.sample_sheet ? file(params.sample_sheet) : null
        )
    )

    if (!params.reference_genome) {
        error('Differential expression requires --ex_reference_genome.')
    }
    if (!params.reference_annotation) {
        error('Differential expression requires --reference_annotation.')
    }

    reference_genome = file(params.reference_genome, checkIfExists: true)
    reference_annotation = file(params.reference_annotation, checkIfExists: true)
    gene_sets = file(params.gene_sets, checkIfExists: true)
    sample_pipeline(sample_batch_ch, reference_genome, reference_annotation, gene_sets)

    onComplete:
    if (params.disable_ping == false) {
        Pinguscript.ping_post(workflow, 'end', 'none', params.out_dir, params)
    }
}
