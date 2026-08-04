nextflow.enable.types = true

include {
    CumulativeSampleBAMGroup;
    MergedIndexedChunkBAM;
    Sample
} from '../lib/sample.nf'

Path optional_file() {
    return file("$projectDir/data/OPTIONAL_FILE")
}

String safe_name(String value) {
    return value.replaceAll(/[^A-Za-z0-9._-]/, '_')
}

/**
 * Quote arbitrary text as one POSIX shell argument for process scripts.
 */
String shell_quote(String value) {
    return "'" + value.replace("'", "'\"'\"'") + "'"
}

String merged_chunk_bam_name(Integer batch_index, Sample sample) {
    return "${safe_name(sample.name)}_${batch_index}.merged.bam"
}

String cumulative_collated_bam_name(CumulativeSampleBAMGroup input_group) {
    return "${safe_name(input_group.sample.name)}_${input_group.batch_index}.cumulative.collated.bam"
}

String nanoplot_output_dir(MergedIndexedChunkBAM merged_bam) {
    return "${safe_name(merged_bam.sample.name)}_${merged_bam.batch_index}.nanoplot"
}

String flagstat_file_name(Integer batch_index, String sample_name) {
    return "${safe_name(sample_name)}_${batch_index}.flagstat.tsv"
}

String oarfish_out_name(Integer batch_index, String sample_name) {
    return "${safe_name(sample_name)}_${batch_index}"
}

String oarfish_counts_file_name(Integer batch_index, String sample_name) {
    return "${oarfish_out_name(batch_index, sample_name)}.quant"
}

String differential_expression_results_dir(Integer batch_index) {
    return "batch_${batch_index}"
}

String get_sample_path(Map meta) {
    return "${meta['runName']}/${meta['replicateName']}"
}

Path get_seq_summary_file(Path bam_file) {
    Path summary_file = file("${bam_file.parent}/seq_summary.txt")
    if (summary_file.exists()) {
        return summary_file
    }
    return optional_file()
}

Map get_sequencing_arguments(Path _run_dir) {
    Map args = [:]
    args['experiment_id'] = params.ex_name
    args['run_id'] = params.ex_run_number
    args['kit'] = params.ex_kit
    if (!params.ex_special_alignment) {
        args['reference_genome'] = params.reference_genome
    }
    args['basecall_config'] = params.ex_basecall_config
    return args
}

def validate_experiment_dir(String experiment_path, int run_number) {
    def experiment_dir = new File(experiment_path)
    if (!experiment_dir.exists()) {
        experiment_dir.mkdirs()
    } else if (!experiment_dir.isDirectory()) {
        throw new RuntimeException("Experiment directory ${experiment_path} is not a directory.")
    } else if (!experiment_dir.canWrite()) {
        throw new RuntimeException("Experiment directory ${experiment_path} is not writable.")
    } else if (!experiment_dir.canRead()) {
        throw new RuntimeException("Experiment directory ${experiment_path} is not readable.")
    }

    String config_file_name = "experiment.config"
    File config_file = new File(experiment_path, config_file_name)
    boolean config_file_exists = config_file.exists()
    if (!config_file_exists && run_number > 1 ) {
        throw new RuntimeException("Experiment directory ${experiment_path} does not contain a configuration file.")
    } else if ( config_file_exists && run_number == 1 ) {
        throw new RuntimeException("Experiment directory ${experiment_path} already contains a configuration file.")
    }
}

Boolean is_empty(path) {
    java.nio.file.Files.newDirectoryStream(path).withCloseable { directory ->
        !directory.iterator().hasNext()
    }
}

process get_versions {
    label 'preproc'
    cpus 1
    output:
        file('versions.txt')
    script:
        """
        python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
        """
}

process get_params {
    label 'seq_lm'
    cpus 1
    output:
        file('params.json')
    script:
        def paramsJSON = new groovy.json.JsonBuilder(params).toPrettyString()
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
    label 'seq_lm'
    publishDir(
        params.out_dir,
        mode: 'copy',
        saveAs: { dirname ? "$dirname/$fname" : fname }
    )
    input:
        tuple(fname: Path, dirname: String?)
    output:
        file(fname.name)
    script:
        '''
        '''
}
