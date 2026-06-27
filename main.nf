#!/usr/bin/env nextflow

nextflow.enable.types = true

include { bam_ingress } from './lib/bam_ingress.nf'
include {
    ChunkQCResult;
    IndexedChunkBAMGroup;
    Sample;
    SampleChunkBAM;
    SampleQCReportInputs
} from './lib/sample.nf'
include { getSamplePath; getSeqSummaryFile; getSequencingArguments } from './lib/util.nf'
include { validateExperimentDir } from './lib/validation.nf'
include { getVersions; getParams; output } from './modules/generic_helpers.nf'
include { quality_control } from './subworkflows/quality_control.nf'

Path optionalFile() {
    return file("$projectDir/data/OPTIONAL_FILE")
}

String safeName(String value) {
    return value.replaceAll(/[^A-Za-z0-9._-]/, '_')
}

/**
 * Quote arbitrary text as one POSIX shell argument for process scripts.
 */
String shellQuote(String value) {
    return "'" + value.replace("'", "'\"'\"'") + "'"
}

String sortedChunkBamName(SampleChunkBAM input_bam) {
    return "${safeName(input_bam.sample.id)}_${input_bam.batch_index}_${input_bam.bam_index_in_chunk}.sorted.bam"
}

String mergedChunkBamName(Integer batch_index, Sample sample) {
    return "${safeName(sample.id)}_${batch_index}.merged.bam"
}

/**
 * Stable key for accumulating live QC state by biological sample.
 */
String sampleKey(Sample sample) {
    return sample.id
}

/**
 * Published directory for one chunk-level QC metric group.
 */
String sampleQCChunkDir(ChunkQCResult result, String metric_name) {
    return "${result.sample.group}/${result.sample.alias}/qc/chunk_${result.batch_index}/${metric_name}"
}

/**
 * Published directory for placeholder report-input manifests for one sample.
 */
String sampleQCReportInputDir(Sample sample) {
    return "${sample.group}/${sample.alias}/qc/report_inputs"
}

/**
 * Manifest filename for one placeholder report-input update.
 */
String qcReportInputManifestName(SampleQCReportInputs report_inputs) {
    return "${safeName(report_inputs.sample.id)}_chunk_${report_inputs.latest_batch_index}_qc_report_inputs.tsv"
}

/**
 * Append one chunk QC result to the mutable per-sample live state and emit the
 * future report-generator input shape for that sample.
 */
SampleQCReportInputs accumulateSampleQCResult(Map<String, List<ChunkQCResult>> state, ChunkQCResult result) {
    String key = sampleKey(result.sample)
    List<ChunkQCResult> previous_chunks = state.containsKey(key) ? state[key] : []
    List<ChunkQCResult> sorted_chunks = (previous_chunks + [result])
        .toSorted { ChunkQCResult left, ChunkQCResult right ->
            left.batch_index <=> right.batch_index
        }

    state[key] = sorted_chunks

    return record(
        latest_batch_index: result.batch_index,
        sample: result.sample,
        chunks: sorted_chunks
    )
}

/**
 * Tabular content logged by the placeholder report-input process.
 */
String qcReportManifestRows(SampleQCReportInputs report_inputs) {
    return report_inputs.chunks
        .collectMany { ChunkQCResult chunk ->
            [
                "${report_inputs.sample.id}\t${chunk.batch_index}\tnanoplot\t${chunk.nanoplot_data}",
                "${report_inputs.sample.id}\t${chunk.batch_index}\tsamtools_flagstat\t${chunk.flagstat}"
            ]
        }
        .join('\n')
}

/**
 * TEMPORARY: sanitize table fields for the throwaway live QC report.
 */
String temporary_qc_report_field(Object value) {
    return "${value}".replaceAll(/[\t\r\n]+/, ' ')
}

/**
 * TEMPORARY: rows for the throwaway live QC report sample table.
 */
String temporary_qc_report_rows(List<SampleQCReportInputs> report_inputs_list) {
    return report_inputs_list
        .collect { SampleQCReportInputs report_inputs ->
            [
                temporary_qc_report_field(report_inputs.sample.id),
                temporary_qc_report_field(report_inputs.sample.alias),
                temporary_qc_report_field(report_inputs.sample.group),
                temporary_qc_report_field(report_inputs.sample.type),
                "${report_inputs.chunks.size()}",
                "${report_inputs.latest_batch_index}"
            ].join('\t')
        }
        .join('\n')
}

/**
 * TEMPORARY: accumulate all sample rows for the throwaway live QC report.
 */
Map accumulate_temporary_qc_report_state(
    Map<String, SampleQCReportInputs> state,
    SampleQCReportInputs report_inputs
) {
    state[sampleKey(report_inputs.sample)] = report_inputs

    List<SampleQCReportInputs> report_inputs_list = state
        .values()
        .toList()
        .toSorted { SampleQCReportInputs left, SampleQCReportInputs right ->
            "${left.sample.group}/${left.sample.alias}/${left.sample.id}" <=>
                "${right.sample.group}/${right.sample.alias}/${right.sample.id}"
        }
    Integer latest_batch_index = report_inputs_list
        .collect { SampleQCReportInputs input -> input.latest_batch_index }
        .max()

    return [
        latest_batch_index: latest_batch_index,
        rows: temporary_qc_report_rows(report_inputs_list)
    ]
}

process writeConfig {
    label 'seq_lm'
    cpus 1
    exec:
        log.info 'Writing config file...'

        // Writing experiment config file should only happen at the first run of the experiment
        def configOut = new File("${params.ex_dir}/experiment.config")

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

process makeReport {
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
            (per_read_stats.name == optionalFile().name) ? '' : "--stats $per_read_stats"
        """
        echo '${metadataJSON}' > metadata.json
        workflow-glue report $report_name \
            --versions versions \
            $stats_args \
            --params params.json \
            --metadata metadata.json
        """
}

process featureCounts {
    label 'seq_lm'
    container 'pegi3s/feature-counts'
    cpus params.threads
    input:
        tuple(meta: Map, bam: Path)
        ref_annotation: Path
    output:
        tuple(meta, countsFile: Path)
    script:
        def runName = meta['runName']
        def replicateName = meta['replicateName']
        def countsFile = "${runName}_${replicateName}_counts.txt"
        """
        featureCounts -a $ref_annotation --fracOverlap 0.9 -M -s 1 -T $task.cpus -L --largestOverlap -o $countsFile $bam
        echo "\$(tail -n +2 $countsFile)" > $countsFile
        """
}

process mergeFeatureCounts {
    debug true
    label 'seq_lm'
    container 'seq_lm_dea'
    cpus 1
    input:
        tuple(meta: Map, newCounts: Path, allCounts: Path)
    stage:
        stageAs newCounts, 'new_counts.txt'
        stageAs allCounts, 'all_counts.txt'
    output:
        tuple(meta, file("${meta['runName']}_${meta['replicateName']}_counts.txt"))
    script:
        def runName = meta['runName']
        def replicateName = meta['replicateName']
        def mergedCountsFile = "${runName}_${replicateName}_counts.txt"
        """
        workflow-glue merge_feature_counts -n $newCounts -a $allCounts -o $mergedCountsFile
        """
}

process bamQC {
    debug true
    label 'seq_lm'
    container 'seq_lm_qualitycontrol'
    cpus 1
    input:
        tuple(meta: Map, bam: Path, bam_index: Path)
        seq_summary: Path

    output:
        tuple(meta, file("run_${params.ex_run_number}_qc.html", optional: true))

    script:
        if (seq_summary.name != optionalFile().name) {
            def qcHTML = "run_${params.ex_run_number}_qc.html"
            def qcTitle = "Run ${params.ex_run_number} QC Report"
            """
            pycoQC --summary_file ${seq_summary} --bam_file ${bam} --report_title "${qcTitle}" --html_outfile ${qcHTML}
            """
        } else {
            log.info 'No seq_summary.txt file found, skipping QC report.'
            '''
            '''
        }
}

process bamIndex {
    label 'seq_lm'
    container 'seq_lm/samtools'
    errorStrategy 'ignore'
    cpus 1
    input:
        tuple(meta: Map, bam: Path)
    output:
        tuple(meta: Map, file(bam), file("${bam}.bai"))
    script:
        """
        samtools sort -@ $task.cpus -o ${bam} ${bam}
        samtools index ${bam}
        """
}

process bam_sort_index {
    label 'seq_lm'
    container 'seq_lm/samtools'
    cpus 4

    input:
        input_bam: SampleChunkBAM

    output:
        record(
            batch_index: input_bam.batch_index,
            sample: input_bam.sample,
            bam_index_in_chunk: input_bam.bam_index_in_chunk,
            bam_count: input_bam.bam_count,
            bam: file(sortedChunkBamName(input_bam)),
            bam_index: file("${sortedChunkBamName(input_bam)}.bai")
        )

    script:
        String sorted_bam = sortedChunkBamName(input_bam)
        """
        samtools sort -o ${sorted_bam} -@ ${task.cpus} ${input_bam.bam}
        samtools index -@ ${task.cpus} ${sorted_bam}
        """
}

process bam_merge_index {
    label 'seq_lm'
    container 'seq_lm/samtools'
    cpus 4

    input:
        input_group: IndexedChunkBAMGroup

    stage:
        stageAs input_group.bams*.bam, 'bam?'
        stageAs input_group.bams*.bam_index, 'bam?.bai'

    output:
        record(
            batch_index: input_group.batch_index,
            sample: input_group.sample,
            bam: file(mergedChunkBamName(input_group.batch_index, input_group.sample)),
            bam_index: file("${mergedChunkBamName(input_group.batch_index, input_group.sample)}.bai")
        )

    script:
        String merged_bam = mergedChunkBamName(input_group.batch_index, input_group.sample)
        String bam_args = input_group.bams*.bam.join(' ')
        """
        printf '%s\\n' ${bam_args} > bams.txt
        samtools merge -o ${merged_bam} -@ ${task.cpus} -b bams.txt
        samtools index -@ ${task.cpus} ${merged_bam}
        """
}

/**
 * Placeholder for the future QC report generator.
 *
 * For now this process only logs and writes the accumulated list of chunk QC
 * TSVs for one sample whenever a new chunk completes.
 */
process qc_report_input_log {
    debug true
    label 'seq_lm_qc'
    cpus 1

    input:
        report_inputs: SampleQCReportInputs

    stage:
        stageAs report_inputs.chunks*.nanoplot_data, 'qc_inputs/chunk?/nanoplot/NanoPlot-data.tsv.gz'
        stageAs report_inputs.chunks*.flagstat, 'qc_inputs/chunk?/samtools_flagstat/*'

    output:
        record(
            latest_batch_index: report_inputs.latest_batch_index,
            sample: report_inputs.sample,
            manifest: file(qcReportInputManifestName(report_inputs))
        )

    script:
        String manifest = qcReportInputManifestName(report_inputs)
        String manifest_rows = qcReportManifestRows(report_inputs)
        String quoted_manifest_rows = shellQuote(manifest_rows)
        """
        printf 'sample_id\\tbatch_index\\tmetric\\tqc_tsv\\n' > ${manifest}
        printf '%s\\n' ${quoted_manifest_rows} >> ${manifest}

        echo "QC report inputs for sample ${report_inputs.sample.id} through chunk ${report_inputs.latest_batch_index}:"
        cat ${manifest}
        """
}

/**
 * TEMPORARY: EPI2ME-displayable live QC placeholder report.
 *
 * This writes only the current list of samples with QC results and must be
 * removed when the permanent live QC report is implemented.
 */
process temporary_qc_report {
    debug true
    label 'seq_lm_qc'
    container 'seq_lm/seq_lm_dea'
    cpus 1
    maxForks 1

    publishDir(
        params.ex_dir,
        mode: 'copy',
        pattern: 'temporary_qc_report.html',
        overwrite: true
    )

    input:
        temporary_qc_report_inputs: Map

    output:
        file('temporary_qc_report.html')

    script:
        String rows = temporary_qc_report_inputs.rows
        String quoted_rows = shellQuote(rows)
        String params_json = new groovy.json.JsonBuilder(params).toPrettyString()
        String quoted_params_json = shellQuote(params_json)
        """
        printf 'sample_id\\talias\\tgroup\\ttype\\tchunks_seen\\tlatest_batch_index\\n' > temporary_qc_report_samples.tsv
        printf '%s\\n' ${quoted_rows} >> temporary_qc_report_samples.tsv

        mkdir versions
        printf 'temporary_qc_report,temporary\\n' > versions/versions.txt
        printf '%s\\n' ${quoted_params_json} > params.json

        export MPLCONFIGDIR="\$PWD/.matplotlib"
        mkdir -p "\$MPLCONFIGDIR"

        workflow-glue temporary_qc_report temporary_qc_report.html \
            --samples temporary_qc_report_samples.tsv \
            --versions versions \
            --params params.json \
            --latest-batch ${temporary_qc_report_inputs.latest_batch_index}
        """
}

process differentiaExpression {
    container 'seq_lm/dea'
    cpus params.threads

    input:
        quantSF: Path

    script:
        """
        deseq-analysis -q ${quantSF} -t "${task.cpus}"
        """
}

// workflow module
workflow sample_pipeline {
    take:
        sampleChunk: Channel
    main:
        getVersions()
        workflow_params = getParams()

        newBamIn = sampleChunk.map { meta, newBam, _allBam -> tuple(meta, newBam) }

        bamIndexResult = bamIndex(newBamIn)

        //Quality control of aligned reads
        bamQCRes = bamQC(bamIndexResult, optionalFile()).map { meta, res ->
            def samplePath = getSamplePath(meta)
            return tuple(res, "${samplePath}/qc")
        }

        featureCountsRes = featureCounts(newBamIn, params.ex_reference_annotation).map { meta, newCounts ->
            def samplePath = getSamplePath(meta)
            def allCounts = file("${params.ex_dir}/${samplePath}/*counts.txt")
            return tuple(meta, newCounts, allCounts)
        }.branch { meta, newCounts, allCounts ->
            initial: allCounts.empty
                def samplePath = getSamplePath(meta)
                return tuple(newCounts, samplePath)
            merge: true
        }

        mergedCountsRes = mergeFeatureCounts(featureCountsRes.merge).map { meta, mergedCounts ->
            def samplePath = getSamplePath(meta)
            return tuple(mergedCounts, samplePath)
        }

        output(featureCountsRes.initial.mix(mergedCountsRes).mix(bamQCRes))

    emit:
        workflow_params as Value
}

Map prepareRun(String experiment_dir, Integer _run_number, Integer replicate_count) {
    def runName = "run_${params.ex_run_number}"
    def runDir = file("${params.ex_dir}/${runName}")

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

    if (params.disable_ping == false) {
        Pinguscript.ping_post(workflow, 'start', 'none', params.ex_dir, params)
    }

    // TODO: Implement parameter validation
    // validateExperimentDir(params.ex_dir, params.ex_run_number)

    // TODO: Implement sequencing setup checks

    // Config is stored in order to fetch parameters in subsequent runs
    // writeConfig()

    // // Setup the run
    // runInfo = prepareRun(params.ex_dir, params.ex_run_number, params.ex_replicate_count)
    // runDir = runInfo.runDir

    // // Start the sequencing run
    // metadataFile = channel.fromPath("${params.ex_dir}/metadata.tsv")
    // keyFile = channel.fromPath(params.ex_mk_key)
    // certificateFile = channel.fromPath(params.ex_mk_cert)
    // sequencingArgs = getSequencingArguments(runDir)
    // startSequencing(sequencingArgs, keyFile, certificateFile, metadataFile)

    // The ingress emits synchronized sample batches. Each non-empty sample chunk
    // is handled independently so live QC does not accumulate previous chunks.
    sample_batch_ch = bam_ingress(
        record(
        live_analysis: params.live_analysis,
        timeline_analysis: params.timeline_analysis,
        sample_sheet_path: params.sample_sheet_path ? file(params.sample_sheet_path) : null
        )
    )

    sample_chunk_bam_ch = sample_batch_ch
        .flatMap { batch ->
            batch.chunks
                .findAll { chunk -> !chunk.bam_paths.isEmpty() }
                .collectMany { chunk ->
                    chunk.bam_paths.withIndex().collect { Path bam_path, Integer bam_index ->
                        record(
                            batch_index: batch.batch_index,
                            sample: chunk.sample,
                            bam_index_in_chunk: bam_index,
                            bam_count: chunk.bam_paths.size(),
                            bam: bam_path
                        )
                    }
                }
        }

    indexed_chunk_bam_ch = bam_sort_index(sample_chunk_bam_ch)

    indexed_chunk_bam_group_ch = indexed_chunk_bam_ch
        .map { indexed_bam ->
            tuple("${indexed_bam.batch_index}:${indexed_bam.sample.id}", indexed_bam.bam_count, indexed_bam)
        }
        .groupBy()
        .map { String _group_key, Bag grouped_bams ->
            List sorted_bams = grouped_bams.toSorted { left, right ->
                left.bam_index_in_chunk <=> right.bam_index_in_chunk
            }
            def first_bam = sorted_bams[0]
            return record(
                batch_index: first_bam.batch_index,
                sample: first_bam.sample,
                bams: sorted_bams
            )
        }

    merged_chunk_bam_ch = bam_merge_index(indexed_chunk_bam_group_ch)
    qc_result_ch = quality_control(merged_chunk_bam_ch)

    Map<String, List<ChunkQCResult>> sample_qc_report_state = [:]
    sample_qc_report_inputs_ch = qc_result_ch.map { result ->
        accumulateSampleQCResult(sample_qc_report_state, result)
    }
    qc_report_input_log_ch = qc_report_input_log(sample_qc_report_inputs_ch)

    Map<String, SampleQCReportInputs> temporary_qc_report_state = [:]
    temporary_qc_report_inputs_ch = sample_qc_report_inputs_ch.map { report_inputs ->
        accumulate_temporary_qc_report_state(temporary_qc_report_state, report_inputs)
    }
    temporary_qc_report(temporary_qc_report_inputs_ch)

    qc_publish_ch = qc_result_ch.flatMap { result ->
        [
            tuple(result.nanoplot_data, sampleQCChunkDir(result, 'nanoplot')),
            tuple(result.flagstat, sampleQCChunkDir(result, 'samtools_flagstat'))
        ]
    }

    qc_report_publish_ch = qc_report_input_log_ch.map { result ->
        tuple(result.manifest, sampleQCReportInputDir(result.sample))
    }

    output(
        qc_publish_ch.mix(qc_report_publish_ch)
    )

    // sample_pipeline(sample_pipeline_input)

    // Start differential expression analysis if there is more than one run
    // if (params.ex_run_number > 1) {
    //     quantResults = channel.watchPath("$runDir/**counts.txt")
    //     .until { result -> result.name.startsWith('STOP') }
    //     .filter { result -> result.name.endsWith('counts.txt') }
    //     .map { _result -> file("$params.ex_dir/**counts.txt") }
    //     // Waits until there are count files for all replicates
    //     .filter { result -> result.size() == (params.ex_run_number * params.ex_replicate_count) }

    //     differentiaExpression(quantResults)
    // }

    onComplete:
    if (params.disable_ping == false) {
        Pinguscript.ping_post(workflow, 'end', 'none', params.ex_dir, params)
    }
}
