EXPECTED_COLUMNS = [
    "name",
    "group",
    "chunks_seen",
    "latest_batch_index",
    "qc_dir",
]
OPTIONAL_COLUMNS = ["order"]
DISPLAY_COLUMNS = {
    "name": "Name",
    "group": "Group",
    "order": "Time (min)",
    "chunks_seen": "QC chunks seen",
    "latest_batch_index": "Latest batch index",
    "qc_dir": "QC input directory",
}
NANOPLOT_COLUMNS = [
    "readIDs",
    "quals",
    "aligned_quals",
    "lengths",
    "aligned_lengths",
    "mapQ",
    "percentIdentity",
]
NANOPLOT_METRICS = [
    "Number of mapped reads",
    "Number of bases",
    "Number of aligned bases",
    "Median read length",
    "Mean read length",
    "Median read quals",
    "Mean read quals",
    "Median MapQ",
    "Mean MapQ",
]
