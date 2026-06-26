nextflow.enable.types = true

record Sample {
    id: String
    alias: String
    type: SampleType
    group: String
    order: Integer?
    bam_dir: Path
}

record SampleChunk {
    sample: Sample
    bam_paths: List<Path>
}

record SampleBatch {
    batch_index: Integer
    chunks: List<SampleChunk>
}

record SampleChunkBAM {
    batch_index: Integer
    sample: Sample
    bam_index_in_chunk: Integer
    bam_count: Integer
    bam: Path
}

record IndexedChunkBAM {
    batch_index: Integer
    sample: Sample
    bam_index_in_chunk: Integer
    bam_count: Integer
    bam: Path
    bam_index: Path
}

record IndexedChunkBAMGroup {
    batch_index: Integer
    sample: Sample
    bams: List<IndexedChunkBAM>
}

record MergedChunkBAM {
    batch_index: Integer
    sample: Sample
    bam: Path
    bam_index: Path
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

enum SampleType {
    CONTROL,
    CONDITION
}
