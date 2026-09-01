nextflow.enable.types = true

include {
    ChunkBAM ;
    CumulativeCollatedBAMGroup ;
    Sample
} from '../lib/sample.nf'

def optional_file() -> Path {
    return file("${projectDir}/data/OPTIONAL_FILE")
}

def safe_name(value: String) -> String {
    return value.replaceAll(/[^A-Za-z0-9._-]/, '_')
}

/**
 * Quote arbitrary text as one POSIX shell argument for process scripts.
 */
def shell_quote(value: String) -> String {
    return "'" + value.replace("'", "'\"'\"'") + "'"
}

def chunk_bam_name(batch_index: Integer, sample) -> String {
    return "${safe_name(sample.name)}_${batch_index}.chunk.bam"
}

def collated_chunk_bam_name(batch_index: Integer, sample) -> String {
    return "${safe_name(sample.name)}_${batch_index}.chunk.collated.bam"
}

def cumulative_collated_bam_name(input_group) -> String {
    return "${safe_name(input_group.sample.name)}_${input_group.batch_index}.cumulative.collated.bam"
}

def nanoplot_output_dir(merged_bam) -> String {
    return "${safe_name(merged_bam.sample.name)}_${merged_bam.batch_index}.nanoplot"
}

def flagstat_file_name(batch_index: Integer, sample_name: String) -> String {
    return "${safe_name(sample_name)}_${batch_index}.flagstat.tsv"
}

def oarfish_out_name(batch_index: Integer, sample_name: String) -> String {
    return "${safe_name(sample_name)}_${batch_index}"
}

def oarfish_counts_file_name(batch_index: Integer, sample_name: String) -> String {
    return "${oarfish_out_name(batch_index, sample_name)}.quant"
}

def differential_expression_results_dir(batch_index: Integer) -> String {
    return "batch_${batch_index}"
}

def get_sample_path(meta: Map) -> String {
    return "${meta['runName']}/${meta['replicateName']}"
}

def get_seq_summary_file(bam_file: Path) -> Path {
    def summary_file: Path = file("${bam_file.parent}/seq_summary.txt")
    if (summary_file.exists()) {
        return summary_file
    }
    return optional_file()
}

def get_sequencing_arguments(_run_dir: Path) -> Map {
    def args: Map = [:]
    args['experiment_id'] = params.ex_name
    args['run_id'] = params.ex_run_number
    args['kit'] = params.ex_kit
    if (!params.ex_special_alignment) {
        args['reference_genome'] = params.reference_genome
    }
    args['basecall_model'] = params.ex_basecall_model
    return args
}

def validate_experiment_dir(experiment_path: Path, run_number: Integer) {
    def experiment_dir = new File(experiment_path)
    if (!experiment_dir.exists()) {
        experiment_dir.mkdirs()
    }
    else if (!experiment_dir.isDirectory()) {
        throw new RuntimeException("Experiment directory ${experiment_path} is not a directory.")
    }
    else if (!experiment_dir.canWrite()) {
        throw new RuntimeException("Experiment directory ${experiment_path} is not writable.")
    }
    else if (!experiment_dir.canRead()) {
        throw new RuntimeException("Experiment directory ${experiment_path} is not readable.")
    }

    def config_file_name: String = "experiment.config"
    def config_file = new File(experiment_path, config_file_name)
    def config_file_exists = config_file.exists()
    if (!config_file_exists && run_number > 1) {
        throw new RuntimeException("Experiment directory ${experiment_path} does not contain a configuration file.")
    }
    else if (config_file_exists && run_number == 1) {
        throw new RuntimeException("Experiment directory ${experiment_path} already contains a configuration file.")
    }
}

def is_empty(path) -> Boolean {
    java.nio.file.Files
        .newDirectoryStream(path)
        .withCloseable { directory ->
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
    echo '${paramsJSON}' > params.json
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
    publishDir params.out_dir, mode: 'copy', saveAs: { dirname ? "${dirname}/${fname}" : fname }

    input:
    tuple(fname: Path, dirname: String?)

    output:
    file(fname.name)

    script:
    '''
        '''
}

/**
 * Publish complete differential-analysis trees under the user-selected output
 * root. The immutable global batch is retained and `latest` is refreshed in
 * serialized analysis order.
 */
process publish_differential_results {
    label 'seq_lm'
    cpus 1
    maxForks 1
    fair true

    publishDir "${params.out_dir}/differential_expression", mode: 'copy', overwrite: false, saveAs: { _fname -> "batch_${analysis_index}" }
    publishDir "${params.out_dir}/differential_expression", mode: 'copy', overwrite: true, saveAs: { _fname -> 'latest' }

    input:
    tuple(local_batch_index: Integer, analysis_index: Integer, results: Path)

    output:
    record(
        batch_index: local_batch_index,
        analysis_index: analysis_index,
        results: file(results.name),
    )

    script:
    '''
        '''
}
