OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")

def getSamplePath(Map meta) {
    return "${meta['runName']}/${meta['replicateName']}"
}

def getSeqSummaryFile(Path bamFile) {
    summaryFile = file("${bamFile.parent}/seq_summary.txt")
    if (summaryFile.exists()) {
        return summaryFile
    } else {
        return OPTIONAL_FILE
    }
}

Map getSequencingArguments(runDir) {
    Map args = [:]
    args['experiment_id'] = params.ex_name
    args['run_id'] = params.ex_run_number
    args['kit'] = params.ex_kit
    if (!params.ex_special_alignment) {
        args['reference_genome'] = params.ex_reference_genome
    }
    args['basecall_config'] = params.ex_basecall_config
    return args
}
