def validateExperimentDir(String experimentPath, int runNumber) {
    def experiment_dir = new File(experimentPath)
    if (!experiment_dir.exists()) {
        experiment_dir.mkdirs()
    } else if (!experiment_dir.isDirectory()) {
        throw new RuntimeException("Experiment directory ${experimentPath} is not a directory.")
    } else if (!experiment_dir.canWrite()) {
        throw new RuntimeException("Experiment directory ${experimentPath} is not writable.")
    } else if (!experiment_dir.canRead()) {
        throw new RuntimeException("Experiment directory ${experimentPath} is not readable.")
    }

    String configFileName = "experiment.config"
    File configFile = new File(experimentPath, configFileName)
    boolean configFileExists = configFile.exists()
    if (!configFileExists && runNumber > 1 ) {
        throw new RuntimeException("Experiment directory ${experimentPath} does not contain a configuration file.")
    } else if ( configFileExists && runNumber == 1 ) {
        throw new RuntimeException("Experiment directory ${experimentPath} already contains a configuration file.")
    }
}

Boolean isEmpty(path) {
    java.nio.file.Files.newDirectoryStream(path).withCloseable { directory ->
        !directory.iterator().hasNext()
    }
}
