include { Sample; SampleType; SampleChunk; SampleBatch } from './sample.nf'

nextflow.enable.types = true

record BamIngressArgs {
    live_analysis: Boolean
    timeline_analysis: Boolean?
    sample_sheet_path: Path?
}

record SampleChunkEvent {
    batch_index: Integer
    sample_index: Integer
    chunk: SampleChunk
}

/**
 * Watch per-sample BAM directories and return synchronized sample chunks.
 *
 * The first emitted item contains one `SampleChunk` for every sample, with all
 * BAMs that already exist under each sample directory. Later items are emitted
 * only after every sample has produced one new BAM.
 */
def bam_ingress(BamIngressArgs ingress_args) {
    def samples: List<Sample> = get_samples(ingress_args)
    validate_samples(samples)

    if (ingress_args.live_analysis) {
        return live_analysis_ingress(samples)
    }

    error('Retrospective analysis not yet implemented.')
}

def get_samples(BamIngressArgs ingress_args) -> List<Sample> {
    if (ingress_args.sample_sheet_path == null) {
        error("BAM ingress requires sample_sheet_path.")
    }
    if (!ingress_args.sample_sheet_path.exists()) {
        error("BAM ingress sample sheet does not exist: ${ingress_args.sample_sheet_path}")
    }
    if (!ingress_args.sample_sheet_path.isFile()) {
        error("BAM ingress sample sheet is not a file: ${ingress_args.sample_sheet_path}")
    }

    return parse_sample_sheet(ingress_args)
}

def parse_sample_type(String value) -> SampleType {
    def normalized = value.trim().toUpperCase()

    if (normalized == 'CONTROL') {
        return SampleType.CONTROL
    }
    if (normalized == 'CONDITION') {
        return SampleType.CONDITION
    }

    error("Invalid sample type '${value}'. Expected one of: CONTROL, CONDITION")
}

def parse_sample_sheet(BamIngressArgs ingress_args) -> List<Sample> {
    def required_fields: Set<String> = ["id", "alias", "type", "group", "bam_dir"].toSet()

    return ingress_args.sample_sheet_path
        .splitCsv(header: true, strip: true)
        .collect { row ->
            def keys: Set<String> = row.keySet() as Set<String>
            def missing = required_fields - keys

            if (missing) {
                error("Sample sheet is missing required fields: ${missing}")
            }

            def sample: Sample = record(
                id: row.id,
                alias: row.alias,
                type: parse_sample_type(row.type),
                group: row.group,
                order: row.order ? row.order.toInteger() : null,
                bam_dir: file(row.bam_dir)
            )

            if (ingress_args.timeline_analysis && sample.order == null) {
                error("You need to provide a sample order when running in timeline mode.")
            }

            return sample
        }
}

def validate_samples(List<Sample> samples) {
    if (samples.empty) {
        error("BAM ingress requires at least one sample.")
    }

    def ids = samples.collect { sample -> sample.id }
    def duplicate_ids = ids.findAll { id -> ids.count(id) > 1 }.unique()
    if (duplicate_ids) {
        error("Sample ids must be unique. Duplicate ids: ${duplicate_ids}")
    }

    def group_aliases = samples.collect { sample -> "${sample.group}\t${sample.alias}" }
    def duplicate_group_aliases = group_aliases
        .findAll { group_alias -> group_aliases.count(group_alias) > 1 }
        .unique()
        .collect { group_alias ->
            def fields = group_alias.split('\t', 2)
            "${fields[0]}/${fields[1]}"
        }
    if (duplicate_group_aliases) {
        error("Sample aliases must be unique within each group. Duplicate group/alias pairs: ${duplicate_group_aliases}")
    }

    samples.each { sample ->
        if (!sample.bam_dir.exists()) {
            error("BAM directory does not exist for sample '${sample.id}': ${sample.bam_dir}")
        }
        if (!sample.bam_dir.isDirectory()) {
            error("BAM path is not a directory for sample '${sample.id}': ${sample.bam_dir}")
        }
    }
}

def live_analysis_ingress(List<Sample> samples) {
    log.info("Running workflow in live analysis mode.")

    List<SampleChunkEvent> existing_events = samples.withIndex().collect { Sample sample, Integer sample_index ->
        make_sample_chunk_event(
            0,
            sample_index,
            make_sample_chunk(sample, get_bam_files_in_dir(sample.bam_dir))
        )
    }

    def ch_events = channel.fromList(existing_events)
    samples.withIndex().each { Sample sample, Integer sample_index ->
        ch_events = ch_events.mix(make_live_analysis_ingress_channel(sample, sample_index))
    }

    return ch_events
        .map { SampleChunkEvent event -> tuple(event.batch_index, samples.size(), event) }
        .groupBy()
        .map { Integer batch_index, Bag<SampleChunkEvent> events ->
            def chunks: List<SampleChunk> = events
                .toSorted { SampleChunkEvent left, SampleChunkEvent right -> left.sample_index <=> right.sample_index }
                .collect { SampleChunkEvent event -> event.chunk }
            def batch: SampleBatch = record(batch_index: batch_index, chunks: chunks)
            return batch
        }
}

def make_live_analysis_ingress_channel(Sample sample, Integer sample_index) {
    Integer batch_index = 0

    return channel.watchPath("${sample.bam_dir}/**", "create")
        .until { Path watched_path -> watched_path.name.startsWith('STOP') }
        .filter { Path watched_path -> watched_path.name.endsWith('.bam') }
        .map { Path watched_path ->
            batch_index += 1
            return make_sample_chunk_event(
                batch_index,
                sample_index,
                make_sample_chunk(sample, [watched_path])
            )
        }
}

def make_sample_chunk_event(Integer batch_index, Integer sample_index, SampleChunk chunk) -> SampleChunkEvent {
    return record(batch_index: batch_index, sample_index: sample_index, chunk: chunk)
}

def make_sample_chunk(Sample sample, List<Path> bam_paths) -> SampleChunk {
    return record(sample: sample, bam_paths: sort_bam_paths(bam_paths))
}

/**
 * Convert a synchronized sample batch to the legacy per-BAM tuple shape consumed
 * by `sample_pipeline`.
 */
def sampleChunksToBamInputs(SampleBatch batch) -> List {
    return batch.chunks.collectMany { chunk ->
        def meta = sample_to_meta(chunk.sample)

        chunk.bam_paths.collect { bam_path ->
            tuple(meta, bam_path, chunk.bam_paths)
        }
    }
}

def sample_to_meta(Sample sample) -> Map {
    return [
        runName: sample.group,
        replicateName: sample.alias,
        sampleId: sample.id,
        alias: sample.alias,
        sampleType: "${sample.type}",
        group: sample.group,
        order: sample.order
    ]
}

def get_bam_files_in_dir(Path dir) -> List<Path> {
    return sort_bam_paths(files("${dir}/**.bam", type: 'file'))
}

def sort_bam_paths(Iterable<Path> bam_paths) -> List<Path> {
    return bam_paths.toSorted { Path left, Path right -> "${left}" <=> "${right}" }
}
