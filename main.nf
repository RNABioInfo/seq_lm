#!/usr/bin/env nextflow

nextflow.enable.types = true

include { bamIngress } from './lib/bamIngress.nf'
include { getSamplePath; getSeqSummaryFile; getSequencingArguments } from './lib/util.nf'
include { validateExperimentDir } from './lib/validation.nf'
include { getVersions; getParams; output } from './modules/generic_helpers.nf'

Path optionalFile() {
    return file("$projectDir/data/OPTIONAL_FILE")
}

process writeConfig {
    label 'seqLM'
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
    label 'seqLM'
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
    label 'seqLM'
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
    label 'seqLM'
    container 'seqlm_dea'
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
    label 'seqLM'
    container 'seqlm_quality'
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
    label 'seqLM'
    container 'seqlm/samtools'
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

process differentiaExpression {
    container 'seqlm/dea'
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

        newBamIn = sampleChunk.map { meta, newBam, _allBam -> [meta, newBam] }

        bamIndexResult = bamIndex(newBamIn)

        //Quality control of aligned reads
        bamQCRes = bamQC(bamIndexResult, optionalFile()).map { meta, res ->
            def samplePath = getSamplePath(meta)
            return [res, "${samplePath}/qc"]
        }

        featureCountsRes = featureCounts(newBamIn, params.ex_reference_annotation).map { meta, newCounts ->
            def samplePath = getSamplePath(meta)
            def allCounts = file("${params.ex_dir}/${samplePath}/*counts.txt")
            return [meta, newCounts, allCounts]
        }.branch { meta, newCounts, allCounts ->
            initial: allCounts.empty
                def samplePath = getSamplePath(meta)
                return [newCounts, samplePath]
            merge: true
        }

        mergedCountsRes = mergeFeatureCounts(featureCountsRes.merge).map { meta, mergedCounts ->
            def samplePath = getSamplePath(meta)
            return [mergedCounts, samplePath]
        }

        output(featureCountsRes.initial.mix(mergedCountsRes).mix(bamQCRes))

    emit:
        workflow_params as Value
}

process startSequencing {
    container 'seqlm_seq'
    debug true
    label 'seqLM'
    cpus 1
    input:
        argumentMap: Map
        keyFile: Path
        certificateFile: Path
        metadataFile: Path
    script:
        def argumentString = argumentMap.collect { k, v -> "--${k} '${v}'" }.join(' ')
        println argumentString
        """
        seq-run-manager ${argumentString} --key_path ${keyFile} --certificate_path ${certificateFile} --metadata ${metadataFile}
        """
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
    validateExperimentDir(params.ex_dir, params.ex_run_number)

    // TODO: Implement sequencing setup checks

    // Config is stored in order to fetch parameters in subsequent runs
    writeConfig()

    // Setup the run
    runInfo = prepareRun(params.ex_dir, params.ex_run_number, params.ex_replicate_count)
    runName = runInfo.runName
    runDir = runInfo.runDir

    // Start the sequencing run
    metadataFile = channel.fromPath("${params.ex_dir}/metadata.tsv")
    keyFile = channel.fromPath(params.ex_mk_key)
    certificateFile = channel.fromPath(params.ex_mk_cert)
    sequencingArgs = getSequencingArguments(runDir)
    startSequencing(sequencingArgs, keyFile, certificateFile, metadataFile)

    // Sample chunk is [map[runName, replicateName], newBam, [allBam]]
    sampleChunk = bamIngress([
    'input':runDir,
    'runName':runName,
    'bam_stats': params.wf.bam_stats,
    'watch_path': params.watch_path])

    sample_pipeline(sampleChunk)

    // Start differential expression analysis if there is more than one run
    if (params.ex_run_number > 1) {
        quantResults = channel.watchPath("$runDir/**counts.txt")
        .until { result -> result.name.startsWith('STOP') }
        .filter { result -> result.name.endsWith('counts.txt') }
        .map { _result -> file("$params.ex_dir/**counts.txt") }
        // Waits until there are count files for all replicates
        .filter { result -> result.size() == (params.ex_run_number * params.ex_replicate_count) }

        differentiaExpression(quantResults)
    }

    onComplete:
    if (params.disable_ping == false) {
        Pinguscript.ping_post(workflow, 'end', 'none', params.ex_dir, params)
    }

    onError:
    if (params.disable_ping == false) {
        Pinguscript.ping_post(workflow, 'error', "$workflow.errorMessage", params.ex_dir, params)
    }
}
