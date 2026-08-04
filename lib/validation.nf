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
