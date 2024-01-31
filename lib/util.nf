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
