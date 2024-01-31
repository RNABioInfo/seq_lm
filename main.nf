#!/usr/bin/env nextflow

// Developer notes
//
// This template workflow provides a basic structure to copy in order
// to create a new workflow. Current recommended practices are:
//     i) create a simple command-line interface.
//    ii) include an abstract workflow scope named "real_time_pipeline" to be used
//        in a module fashion
//   iii) a second concrete, but anonymous, workflow scope to be used
//        as an entry point when using this workflow in isolation.

import groovy.json.JsonBuilder
import java.text.SimpleDateFormat

nextflow.enable.dsl = 2

include { bam_ingress } from './lib/bamIngress.nf'
include { getSamplePath } from './lib/util.nf'
include { getSeqSummaryFile } from './lib/util.nf'

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")

process getVersions {
    label 'preproc'
    cpus 1
    output:
        path 'versions.txt'
    script:
    """
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    bamstats --version | sed 's/^/bamstats,/' >> versions.txt
    """
}

process getParams {
    label 'seqLM'
    cpus 1
    output:
        path 'params.json'
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}

// See https://github.com/nextflow-io/nextflow/issues/1636. This is the only way to
// publish files from a workflow whilst decoupling the publish from the process steps.
// The process takes a tuple containing the filename and the name of a sub-directory to
// put the file into. If the latter is `null`, puts it into the top-level directory.
process output {
    // publish inputs to output directory
    debug true
    label 'seqLM'
    publishDir(
        params.ex_dir,
        mode: 'copy',
        saveAs: { dirname ? "$dirname/$fname" : fname }
    )
    input:
        tuple path(fname), val(dirname)
    output:
        path fname
    """
    """
}

process writeConfig {
    label 'seqLM'
    cpus 1
    exec:
        log.info 'Writing config file...'

        // Writing experiment config file should only happen at the first run of the experiment
        configOut = new File("${params.ex_dir}/experiment.config")

        configOut.withWriter { w ->
            w << 'params {\n'
            params.each { k, v ->
                if (k.startsWith('ex')) {
                    if ( k == 'ex_run_number' ) {
                        v = v + 1
                    }
                    line = ""
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
        val metadata
        path per_read_stats
        path 'versions/*'
        path 'params.json'
    output:
        path 'wf-template-*.html'
    script:
        String report_name = 'wf-template-report.html'
        String metadata = new JsonBuilder(metadata).toPrettyString()
        String stats_args = \
            (per_read_stats.name == OPTIONAL_FILE.name) ? '' : "--stats $per_read_stats"
    """
    echo '${metadata}' > metadata.json
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
        tuple val(meta), path(bam)
        path ref_annotation
    output:
        tuple val(meta), path(countsFile)
    script:
        runName = meta['runName']
        replicateName = meta['replicateName']
        countsFile = "${runName}_${replicateName}_counts.txt"
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
        tuple val(meta), path(newCounts, name: "new_counts.txt"), path(allCounts, name: "all_counts.txt")
    output:
        tuple val(meta), path(mergedCountsFile)
    script:
        runName = meta['runName']
        replicateName = meta['replicateName']
        mergedCountsFile = "${runName}_${replicateName}_counts.txt"
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
        tuple val(meta), path(bam), path(bam_index)
        path seq_summary

    output:
        tuple val(meta), path(qcHTML), optional: true

    script:
        if (seq_summary.name != OPTIONAL_FILE.name) {
            qcHTML = "run_${params.ex_run_number}_qc.html"
            qcTitle = "Run ${params.ex_run_number} QC Report"
            """
            pycoQC --summary_file ${seq_summary} --bam_file ${bam} --report_title "${qcTitle}" --html_outfile ${qcHTML}
            """
        } else {
            log.info 'No seq_summary.txt file found, skipping QC report.'
            """
            """
        }
}

process bamIndex {
    label 'seqLM'
    container 'staphb/samtools'
    errorStrategy 'ignore'
    cpus 1
    input:
        tuple val(meta), path(bam)
    output:
        tuple val(meta), path(bam), path(bam_index)
    script:
        bam_index = "${bam}.bai"
        """
        samtools sort -@ $task.cpus -o ${bam} ${bam}
        samtools index ${bam}
        """
}

process differentiaExpression {
    container "seqlm_dea"
    cpus params.threads

    input:
        path quantSF

    script:
        """
        workflow-glue deseq -q ${quantSF} -t "${task.cpus}"
        """
}

// workflow module
workflow sample_pipeline {
    take:
        sampleChunk
    main:
        software_versions = getVersions()
        workflow_params = getParams()
        
        newBamIn = sampleChunk.map { meta, newBam, allBam -> [meta, newBam] }

        bamIndexResult = bamIndex(newBamIn)

        //Quality control of aligned reads
        bamQCRes = bamQC(bamIndexResult, OPTIONAL_FILE) |
        map { meta, res -> 
            samplePath = getSamplePath(meta)
            return [res, "${samplePath}/qc"] 
        }

        featureCountsRes = featureCounts(newBamIn, params.ref_annotation) |
        map { meta, newCounts -> 
            samplePath = getSamplePath(meta)
            allCounts = file("${params.ex_dir}/${samplePath}/*counts.txt")
            return [meta, newCounts, allCounts] } | 
        branch { meta, newCounts, allCounts -> 
            initial: allCounts.empty
                samplePath = getSamplePath(meta)
                return [newCounts, samplePath]
            merge: true }

        mergedCountsRes = mergeFeatureCounts(featureCountsRes.merge) |
        map { meta, mergedCounts -> 
            samplePath = getSamplePath(meta)
            return [mergedCounts, samplePath] }

        featureCountsRes.initial |
        mix(mergedCountsRes) |
        mix(bamQCRes) | 
        output

    emit:
        workflow_params
}

// Entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {
    if (params.disable_ping == false) {
        Pinguscript.ping_post(workflow, 'start', 'none', params.ex_dir, params)
    }

    writeConfig()

    runName = "run_${params.ex_run_number}"
    runDir = file("${params.ex_dir}/${runName}")

    (1..params.ex_replicate_number).each { replicate_number ->
        replicateDir = file("${runDir}/replicate_${replicate_number}")
        replicateDir.mkdirs()
    }

    // Sample chunk is [map[runName, replicateName], newBam, [allBam]]
    sampleChunk = bam_ingress([
    'input':runDir,
    'runName':runName,
    'bam_stats': params.wf.bam_stats,
    'watch_path': params.watch_path])

    sample_pipeline(sampleChunk)

    if (params.ex_run_number > 1) {
        quantResults = Channel.watchPath("$runDir/**counts.txt")
        .until { it.name.startsWith("STOP") }
        .filter { it.name.endsWith("counts.txt") }
        .map { file("$params.ex_dir/**counts.txt") }
        // Waits until there are count files for all replicates
        .filter { it.size() == (params.ex_run_number * params.ex_replicate_number) }

        differentiaExpression(quantResults)
    }
}

if (params.disable_ping == false) {
    workflow.onComplete {
        Pinguscript.ping_post(workflow, 'end', 'none', params.ex_dir, params)
    }

    workflow.onError {
        Pinguscript.ping_post(workflow, 'error', "$workflow.errorMessage", params.ex_dir, params)
    }
}
