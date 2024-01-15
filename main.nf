#!/usr/bin/env nextflow

// Developer notes
//
// This template workflow provides a basic structure to copy in order
// to create a new workflow. Current recommended practices are:
//     i) create a simple command-line interface.
//    ii) include an abstract workflow scope named "pipeline" to be used
//        in a module fashion
//   iii) a second concrete, but anonymous, workflow scope to be used
//        as an entry point when using this workflow in isolation.

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { bam_ingress } from './lib/bamIngress.nf'

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")

process count_transcripts {
    container = "combinelab/salmon:latest"
    // Count transcripts using Salmon.
    // library type is specified as forward stranded (-l SF) as it should have either been through pychopper or come from direct RNA reads.
    label "isoforms"
    cpus params.threads
    input:
        tuple val(meta), path(bam), path(ref_transcriptome)
    output:
        path "*transcript_counts.tsv", emit: counts
        path "*seqkit.stats", emit: seqkit_stats
    """
    salmon quant --noErrorModel -p "${task.cpus}" -t "${ref_transcriptome}" -l SF -a "${bam}" -o counts
    mv counts/quant.sf "${meta.alias}.transcript_counts.tsv"
    seqkit bam  "${bam}" 2>  "${meta.alias}.seqkit.stats"
    """
}


process getVersions {
    label "seqLM"
    cpus 1
    output:
        path "versions.txt"
    script:
    """
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    bamstats --version | sed 's/^/bamstats,/' >> versions.txt
    """
}


process getParams {
    label "seqLM"
    cpus 1
    output:
        path "params.json"
    script:
        String paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}


process makeReport {
    label "seqLM"
    input:
        val metadata
        path per_read_stats
        path "versions/*"
        path "params.json"
    output:
        path "wf-template-*.html"
    script:
        String report_name = "wf-template-report.html"
        String metadata = new JsonBuilder(metadata).toPrettyString()
        String stats_args = \
            (per_read_stats.name == OPTIONAL_FILE.name) ? "" : "--stats $per_read_stats"
    """
    echo '${metadata}' > metadata.json
    workflow-glue report $report_name \
        --versions versions \
        $stats_args \
        --params params.json \
        --metadata metadata.json
    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636. This is the only way to
// publish files from a workflow whilst decoupling the publish from the process steps.
// The process takes a tuple containing the filename and the name of a sub-directory to
// put the file into. If the latter is `null`, puts it into the top-level directory.
process output {
    // publish inputs to output directory
    label "seqLM"
    publishDir (
        params.out_dir,
        mode: "copy",
        saveAs: { dirname ? "$dirname/$fname" : fname }
    )
    input:
        tuple path(fname), val(dirname)
    output:
        path fname
    """
    """
}

// Creates a new directory named after the sample alias and moves the fastcat results
// into it.
process collectFastqIngressResultsInDir {
    label "seqLM"
    input:
        // both the fastcat seqs as well as stats might be `OPTIONAL_FILE` --> stage in
        // different sub-directories to avoid name collisions
        tuple val(meta), path(concat_seqs, stageAs: "seqs/*"), path(fastcat_stats,
            stageAs: "stats/*")
    output:
        // use sub-dir to avoid name clashes (in the unlikely event of a sample alias
        // being `seq` or `stats`)
        path "out/*"
    script:
    String outdir = "out/${meta["alias"]}"
    String metaJson = new JsonBuilder(meta).toPrettyString()
    String concat_seqs = \
        (concat_seqs.fileName.name == OPTIONAL_FILE.name) ? "" : concat_seqs
    String fastcat_stats = \
        (fastcat_stats.fileName.name == OPTIONAL_FILE.name) ? "" : fastcat_stats
    """
    mkdir -p $outdir
    echo '$metaJson' > metamap.json
    mv metamap.json $concat_seqs $fastcat_stats $outdir
    """
}

// workflow module
workflow pipeline {
    take:
        bam_ingress_results
    main:
        software_versions = getVersions()
        workflow_params = getParams()
        metadata = bam_ingress_results.map { it[0] }.toList()
        metadata.view()
    emit:
        workflow_params
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    if (params.disable_ping == false) {
        Pinguscript.ping_post(workflow, "start", "none", params.out_dir, params)
    }

    samples = bam_ingress([
        "input":params.bam,
        "sample":params.sample,
        "bam_stats": params.wf.bam_stats,
        "watch_path": params.watch_path])

    pipeline(samples)
}

if (params.disable_ping == false) {
    workflow.onComplete {
        Pinguscript.ping_post(workflow, "end", "none", params.out_dir, params)
    }

    workflow.onError {
        Pinguscript.ping_post(workflow, "error", "$workflow.errorMessage", params.out_dir, params)
    }
}
