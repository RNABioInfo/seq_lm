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

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")

process count_transcripts {
    container = 'combinelab/salmon:latest'
    // Count transcripts using Salmon.
    // library type is specified as forward stranded (-l SF) as it should have either been through pychopper or come from direct RNA reads.
    label 'isoforms'
    cpus params.threads
    input:
        tuple val(meta), path(bam), path(ref_transcriptome)
    output:
        path '*transcript_counts.tsv', emit: counts
        path '*seqkit.stats', emit: seqkit_stats
    """
    salmon quant --noErrorModel -p "${task.cpus}" -t "${ref_transcriptome}" -l SF -a "${bam}" -o counts
    mv counts/quant.sf "${meta.alias}.transcript_counts.tsv"
    seqkit bam  "${bam}" 2>  "${meta.alias}.seqkit.stats"
    """
}

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

process writeConfig {
    label 'seqLM'
    cpus 1
    exec:
        log.info 'Writing config file'

        // Writing experiment config file should only happen at the first run of the experiment
        configOut = new File("${params.ex_dir}/experiment.config")

        if (configOut.isEmpty()) {
            configOut << 'params {\n'
            params.each { k, v ->
                if (k.startsWith('ex')) {
                    String line = "\t${k} = ${v}\n"
                    configOut << line
                }
            }
            configOut << '}\n'
        } else {
            println 'Config file already exists, skipping.'
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

// See https://github.com/nextflow-io/nextflow/issues/1636. This is the only way to
// publish files from a workflow whilst decoupling the publish from the process steps.
// The process takes a tuple containing the filename and the name of a sub-directory to
// put the file into. If the latter is `null`, puts it into the top-level directory.
process output {
    // publish inputs to output directory
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
    '''
    '''
}

// Creates a new directory named after the sample alias and moves the fastcat results
// into it.
process collectFastqIngressResultsInDir {
    label 'seqLM'
    input:
        // both the fastcat seqs as well as stats might be `OPTIONAL_FILE` --> stage in
        // different sub-directories to avoid name collisions
        tuple val(meta), path(concat_seqs, stageAs: 'seqs/*'), path(fastcat_stats,
            stageAs: 'stats/*')
    output:
        // use sub-dir to avoid name clashes (in the unlikely event of a sample alias
        // being `seq` or `stats`)
        path 'out/*'
    script:
    String outdir = "out/${meta['alias']}"
    String metaJson = new JsonBuilder(meta).toPrettyString()
    String concat_seqs = \
        (concat_seqs.fileName.name == OPTIONAL_FILE.name) ? '' : concat_seqs
    String fastcat_stats = \
        (fastcat_stats.fileName.name == OPTIONAL_FILE.name) ? '' : fastcat_stats
    """
    mkdir -p $outdir
    echo '$metaJson' > metamap.json
    mv metamap.json $concat_seqs $fastcat_stats $outdir
    """
}

process salmonQuant {
    label 'seqLM'
    errorStrategy 'ignore'
    cpus params.threads
    container 'combinelab/salmon:latest'

    input:
        path bam
        path refTranscriptome
    output:
        path quantSF

    script:
    quantDir = 'quantification'
    quantSF = "${quantDir}/quant*"
    """
    salmon quant --ont -p "$task.cpus" -t "$refTranscriptome" -l SF -a $bam -o quantification
    """
}

process bamQC {
    label 'seqLM'
    container 'seqlm_quality'
    cpus 1
    input:
        tuple path(bam), path(bam_index)
        path seq_summary

    output:
        path qcHTML

    script:
        qcHTML = "run_${params.ex_run_number}_qc.html"
        qcTitle = "Run ${params.ex_run_number} QC Report"
        """
        pycoQC --summary_file ${seq_summary} --bam_file ${bam} --report_title "${qcTitle}" --html_outfile ${qcHTML}
        """
}

process bamIndex {
    label 'seqLM'
    container 'staphb/samtools'
    cpus params.threads
    input:
        path bam
    output:
        tuple path(bam), path(bam_index)
    script:
        bam_index = "${bam}.bai"
        """
        samtools sort -@ $task.cpus -o ${bam} ${bam}
        samtools index ${bam}
        """
}

// workflow module
workflow real_time_pipeline {
    take:
        newBam
        seqSummary
        allBam
        runName
    main:
        software_versions = getVersions()
        workflow_params = getParams()

        bamIndexResult = bamIndex(newBam)

        // Quality control of aligned reads
        bamQCRes = bamQC(bamIndexResult, seqSummary) |
        map { qc_result -> [qc_result, "${runName}/qc"] }

        // Count transcripts using Salmon
        salmonQuantRes = salmonQuant(allBam, params.ref_transcriptome) |
        map { quant_result -> [quant_result, 'quant'] }

        bamQCRes | 
        concat(salmonQuantRes) | 
        output

    emit:
        workflow_params
}

// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {
    if (params.disable_ping == false) {
        Pinguscript.ping_post(workflow, 'start', 'none', params.ex_dir, params)
    }

    if (params.watch_path) {
        writeConfig()

        run_name = "run_${params.ex_run_number}"
        runDir = file("${params.ex_dir}/${run_name}")
        runDir.mkdirs()

        // Emit a new map of bam files whenever a new bam file is added to the run directory
        input = Channel.watchPath("$runDir/*.bam")
        .until { it.name.startsWith('STOP') }
        .multiMap { newBam -> 
            newBam: newBam
            seqSummary: file(newBam.parent / "seq_summary.txt", checkIfExists: true)
            allBam: file("$runDir/*.bam")
        }

        real_time_pipeline(
            input.newBam, 
            input.seqSummary, 
            input.allBam, 
            run_name)


        // sample_chunk = bam_ingress([
        //     'input':runDir,
        //     'sample':params.ex_name,
        //     'bam_stats': params.wf.bam_stats,
        //     'watch_path': params.watch_path])

    } else {
        error 'Retrospective analysis not yet implemented.'
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
