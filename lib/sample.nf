nextflow.enable.types = true

record Sample {
    name: String
    group: String
    order: Integer?
    bam_dir: Path
    is_live: Boolean
}

record SampleChunk {
    sample: Sample
    bam_paths: List<Path>
}

record SampleBatch {
    batch_index: Integer
    chunks: List<SampleChunk>
    experiment_sample_count: Integer
}

record SampleBatchSize {
    batch_index: Integer
    active_sample_count: Integer
    experiment_sample_count: Integer
}

record SampleChunkBAMGroup {
    batch_index: Integer
    sample: Sample
    bams: List<Path>
}

record MergedIndexedChunkBAM {
    batch_index: Integer
    sample: Sample
    bam: Path
    bam_index: Path
}

record MergedChunkBAM {
    batch_index: Integer
    sample: Sample
    bam: Path
}

record MergedIndexedChunkBAMBatch {
    batch_index: Integer
    bams: List<MergedIndexedChunkBAM>
    experiment_sample_count: Integer
}

record CumulativeSampleBAMGroup {
    batch_index: Integer
    sample: Sample
    bams: List<MergedIndexedChunkBAM>
}

record CumulativeSampleBAMSnapshot {
    batch_index: Integer
    sample_bams: List<CumulativeSampleBAMGroup>
    experiment_sample_count: Integer
}

record QuantifiedSample {
    batch_index: Integer
    sample: Sample
    counts: Path
}

record QuantifiedSampleBatch {
    batch_index: Integer
    samples: List<QuantifiedSample>
}

record QuantifiedSampleUpdateBatch {
    batch_index: Integer
    samples: List<QuantifiedSample>
    experiment_sample_count: Integer
}

record DifferentialExpressionResult {
    batch_index: Integer
    results: Path
}

/**
 * NanoPlot output for one merged sample chunk.
 *
 * This is chunk-level QC state: each record describes the NanoPlot tabular
 * read-level data for one sample at one batch index.
 */
record NanoPlotQCResult {
    batch_index: Integer
    sample: Sample
    bam: Path
    bam_index: Path
    nanoplot_data: Path
}

/**
 * Samtools flagstat output for one merged sample chunk.
 *
 * This is chunk-level QC state: each record describes the alignment summary
 * table for one sample at one batch index.
 */
record FlagstatQCResult {
    batch_index: Integer
    sample: Sample
    bam: Path
    bam_index: Path
    flagstat: Path
}

/**
 * Combined QC outputs for one merged sample chunk.
 *
 * The main workflow publishes these files and also folds them into the
 * accumulated per-sample report input state.
 */
record ChunkQCResult {
    batch_index: Integer
    sample: Sample
    bam: Path
    bam_index: Path
    nanoplot_data: Path
    flagstat: Path
}

/**
 * Accumulated QC inputs for one sample after a chunk finishes.
 *
 * This is sample-level live state: every emission contains all chunk QC
 * results seen so far for the sample and is the future report-generator input.
 */
record SampleQCReportInputs {
    latest_batch_index: Integer
    sample: Sample
    chunks: List<ChunkQCResult>
}
