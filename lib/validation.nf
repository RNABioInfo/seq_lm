nextflow.enable.types = true

def validate_experiment_dir(experiment_path: String, run_number: Integer) {
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
    def config_file: File = new File(experiment_path, config_file_name)
    def config_file_exists: Boolean = config_file.exists()
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
