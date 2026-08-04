include { Sample; SampleChunk; SampleBatch } from './sample.nf'

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

record LivePathEvent {
    source: String
    path: Path
}

/**
 * Watch per-sample BAM directories and return synchronized sample chunks.
 *
 * The first emitted item contains one `SampleChunk` for every sample, with all
 * BAMs that already exist under each sample directory. Later items contain the
 * live subset and are emitted only after every live sample produces one BAM.
 */
def bam_ingress(BamIngressArgs ingress_args) {
    def samples: List<Sample> = get_samples(ingress_args)
    validate_samples(samples)
    return mixed_analysis_ingress(samples, ingress_args.live_analysis)
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

def parse_sample_sheet(BamIngressArgs ingress_args) -> List<Sample> {
    def required_fields: Set<String> = ["alias", "group", "bam_dir"].toSet()

    return ingress_args.sample_sheet_path
        .splitCsv(header: true, strip: true)
        .collect { row ->
            def keys: Set<String> = row.keySet() as Set<String>
            def missing = required_fields - keys

            if (missing) {
                error("Sample sheet is missing required fields: ${missing}")
            }
            def removed_fields = keys.intersect(["id", "type"].toSet())
            if (removed_fields) {
                error("Sample sheet contains removed fields: ${removed_fields}. Use alias, group, bam_dir, and optional order.")
            }

            def sample: Sample = record(
                name: row.alias,
                group: row.group,
                order: row.order ? row.order.toInteger() : null,
                bam_dir: file(row.bam_dir),
                is_live: parse_is_live(row.is_live, row.alias)
            )

            if (!sample.name?.trim()) {
                error("Sample sheet contains a row with an empty alias.")
            }
            if (!sample.group?.trim()) {
                error("Sample sheet contains a row with an empty group for alias '${sample.name}'.")
            }
            if (!row.bam_dir?.trim()) {
                error("Sample sheet contains a row with an empty bam_dir for sample '${sample_label(sample)}'.")
            }

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

    def control_count = samples.count { sample -> sample.group.trim().equalsIgnoreCase('control') }
    if (control_count < 2) {
        error("Sample sheet requires at least two control samples in the group column.")
    }

    def group_names = samples.collect { sample -> sample_key(sample) }
    def duplicate_group_names = group_names
        .findAll { group_name -> group_names.count(group_name) > 1 }
        .unique()
        .collect { group_name ->
            def fields = group_name.split('\t', 2)
            "${fields[0]}/${fields[1]}"
        }
    if (duplicate_group_names) {
        error("Sample names must be unique within each group. Duplicate group/name pairs: ${duplicate_group_names}")
    }

    samples.each { sample ->
        if (!sample.bam_dir.exists()) {
            error("BAM directory does not exist for sample '${sample_label(sample)}': ${sample.bam_dir}")
        }
        if (!sample.bam_dir.isDirectory()) {
            error("BAM path is not a directory for sample '${sample_label(sample)}': ${sample.bam_dir}")
        }
    }
}

Boolean parse_is_live(Object value, Object alias) {
    String normalized = value == null ? '' : "${value}".trim().toLowerCase()
    if (!normalized || normalized == 'true') {
        return true
    }
    if (normalized == 'false') {
        return false
    }
    error(
        "Invalid is_live value '${value}' for sample '${alias}'. " +
        "Expected true, false, or a blank value."
    )
}

String sample_key(Sample sample) {
    return "${sample.group}\t${sample.name}"
}

String sample_label(Sample sample) {
    return "${sample.group}/${sample.name}"
}

def mixed_analysis_ingress(List<Sample> samples, Boolean live_analysis) {
    List<Sample> live_samples = samples.findAll { Sample sample ->
        live_analysis && sample.is_live
    }
    log.info(
        "Running workflow with ${live_samples.size()} live sample(s) and " +
        "${samples.size() - live_samples.size()} final sample(s)."
    )

    Map<String, List<Path>> existing_bams_by_sample = [:]
    List<SampleChunk> existing_chunks = samples.collect { Sample sample ->
        List<Path> existing_bams = get_bam_files_in_dir(sample.bam_dir)
        existing_bams_by_sample[sample_key(sample)] = existing_bams
        Boolean effectively_live = live_analysis && sample.is_live
        log.info(
            "BAM ingress startup scan for sample '${sample_label(sample)}' " +
            "found ${existing_bams.size()} existing BAM file(s) in ${sample.bam_dir}; " +
            "effective mode=${effectively_live ? 'live' : 'final'}."
        )
        if (!effectively_live && existing_bams.empty) {
            error(
                "Final sample '${sample_label(sample)}' must contain at least one BAM file at startup: " +
                "${sample.bam_dir}"
            )
        }
        make_sample_chunk(sample, existing_bams)
    }

    SampleBatch startup_batch = record(
        batch_index: 0,
        chunks: existing_chunks,
        experiment_sample_count: samples.size()
    )
    if (live_samples.empty) {
        log.info('No effective live samples remain; BAM ingress will emit startup batch 0 and close.')
        return channel.of(startup_batch)
    }

    Set<String> stopped_samples = java.util.Collections.synchronizedSet(new LinkedHashSet<String>())
    def stopped_sample_count = new java.util.concurrent.atomic.AtomicInteger(0)
    Map<String, Integer> stopped_batch_counts = java.util.Collections.synchronizedMap([:])
    def ch_events = channel.empty()
    live_samples.each { Sample sample ->
        Integer sample_index = samples.indexOf(sample)
        ch_events = ch_events.mix(make_live_analysis_ingress_channel(
            sample,
            sample_index,
            existing_bams_by_sample[sample_key(sample)],
            stopped_samples,
            stopped_sample_count,
            stopped_batch_counts,
            live_samples.size()
        ))
    }

    Map<Integer, List<SampleChunkEvent>> pending_batch_events = [:]
    def live_batches_ch = ch_events
        .map { SampleChunkEvent event ->
            List<SampleChunkEvent> events = pending_batch_events[event.batch_index]
            if (events == null) {
                events = []
                pending_batch_events[event.batch_index] = events
            }
            events.add(event)
            if (events.size() < live_samples.size()) {
                return null
            }
            if (events.size() > live_samples.size()) {
                error(
                    "Live BAM batch ${event.batch_index} received more than ${live_samples.size()} sample event(s)."
                )
            }
            pending_batch_events.remove(event.batch_index)
            def chunks: List<SampleChunk> = events
                .toSorted { SampleChunkEvent left, SampleChunkEvent right -> left.sample_index <=> right.sample_index }
                .collect { SampleChunkEvent grouped_event -> grouped_event.chunk }
            def batch: SampleBatch = record(
                batch_index: event.batch_index,
                chunks: chunks,
                experiment_sample_count: samples.size()
            )
            log.info(
                "Emitting synchronized live BAM batch ${batch.batch_index}: " +
                batch.chunks.collect { chunk -> "${sample_label(chunk.sample)}=${chunk.bam_paths.size()} BAM(s)" }.join(", ")
            )
            return batch
        }
        .filter { SampleBatch batch -> batch != null }

    return channel.of(startup_batch).concat(live_batches_ch)
}

def make_live_analysis_ingress_channel(
    Sample sample,
    Integer sample_index,
    List<Path> initial_bams,
    Set<String> stopped_samples,
    java.util.concurrent.atomic.AtomicInteger stopped_sample_count,
    Map<String, Integer> stopped_batch_counts,
    Integer live_sample_count
) {
    List<Path> existing_stop_files = get_stop_files_in_dir(sample.bam_dir)
    Integer batch_index = 0
    Set<String> accepted_bam_paths = new LinkedHashSet<String>()
    initial_bams.each { Path bam_path -> accepted_bam_paths.add(path_key(bam_path)) }

    if (!existing_stop_files.empty) {
        log.warn(
            "Sample '${sample_label(sample)}' has existing STOP file(s) " +
            "at startup: ${existing_stop_files}. No future BAM watcher will be opened for this sample."
        )
        log_sample_stop(
            sample,
            existing_stop_files[0],
            "existing",
            batch_index,
            stopped_samples,
            stopped_sample_count,
            stopped_batch_counts,
            live_sample_count
        )
        return channel.empty()
    }

    log.info(
        "Starting live BAM watcher for sample '${sample_label(sample)}' " +
        "in ${sample.bam_dir}; waiting for BAM create events and STOP create/modify events."
    )

    return LiveBamWatcher.watch(
        sample.bam_dir,
        { Integer watched_dir_count ->
            log.info(
                "Live BAM watcher is active for sample '${sample_label(sample)}'; " +
                "watching ${watched_dir_count} director${watched_dir_count == 1 ? 'y' : 'ies'} under ${sample.bam_dir}."
            )
        },
        { Throwable error ->
            log.error(
                "Live BAM watcher failed for sample '${sample_label(sample)}'; " +
                "closing this sample stream.",
                error
            )
        }
    )
        .map { event -> make_live_path_event(event.source as String, event.path as Path) }
        .flatMap { LivePathEvent event ->
            if (is_stop_path(event.path)) {
                List<Path> pending_bams = get_bam_files_in_dir(sample.bam_dir)
                    .findAll { Path bam_path -> !accepted_bam_paths.contains(path_key(bam_path)) }
                if (!pending_bams.empty) {
                    log.info(
                        "Draining ${pending_bams.size()} pending BAM file(s) for sample '${sample_label(sample)}' " +
                        "before honoring STOP at ${event.path}."
                    )
                }

                List<SampleChunkEvent> drained_events = pending_bams.collect { Path pending_bam ->
                    accepted_bam_paths.add(path_key(pending_bam))
                    batch_index += 1
                    log.info(
                        "Accepted live BAM event for sample '${sample_label(sample)}': " +
                        "batch_index=${batch_index}, source=drain, path=${pending_bam}."
                    )
                    return make_sample_chunk_event(
                        batch_index,
                        sample_index,
                        make_sample_chunk(sample, [pending_bam])
                    )
                }
                log_sample_stop(
                    sample,
                    event.path,
                    event.source,
                    batch_index,
                    stopped_samples,
                    stopped_sample_count,
                    stopped_batch_counts,
                    live_sample_count
                )
                return drained_events
            }
            if (event.source != "create") {
                log.info(
                    "Ignoring non-create live BAM path event for sample '${sample_label(sample)}': " +
                    "source=${event.source}, path=${event.path}."
                )
                return []
            }
            if (!is_bam_path(event.path)) {
                log.info(
                    "Ignoring non-BAM live path event for sample '${sample_label(sample)}': " +
                    "source=${event.source}, path=${event.path}."
                )
                return []
            }
            if (!accepted_bam_paths.add(path_key(event.path))) {
                log.info(
                    "Ignoring duplicate live BAM path event for sample '${sample_label(sample)}': " +
                    "source=${event.source}, path=${event.path}."
                )
                return []
            }

            batch_index += 1
            log.info(
                "Accepted live BAM event for sample '${sample_label(sample)}': " +
                "batch_index=${batch_index}, source=${event.source}, path=${event.path}."
            )
            return [make_sample_chunk_event(
                batch_index,
                sample_index,
                make_sample_chunk(sample, [event.path])
            )]
        }
}

def make_live_path_event(String source, Path path) -> LivePathEvent {
    return record(source: source, path: path)
}

def log_sample_stop(
    Sample sample,
    Path stop_path,
    String source,
    Integer batch_count,
    Set<String> stopped_samples,
    java.util.concurrent.atomic.AtomicInteger stopped_sample_count,
    Map<String, Integer> stopped_batch_counts,
    Integer live_sample_count
) {
    if (!stopped_samples.add(sample_key(sample))) {
        log.warn(
            "Duplicate STOP event for already stopped sample '${sample_label(sample)}': " +
            "source=${source}, path=${stop_path}."
        )
        return
    }

    stopped_batch_counts[sample_key(sample)] = batch_count
    Integer stopped_count = stopped_sample_count.incrementAndGet()
    log.info(
        "Observed STOP for sample '${sample_label(sample)}': " +
        "source=${source}, path=${stop_path}. Stopped ${stopped_count}/${live_sample_count} sample watcher(s)."
    )
    if (stopped_count == live_sample_count) {
        List<Integer> batch_counts = stopped_batch_counts.values().toList()
        if (batch_counts.unique().size() != 1) {
            String counts = stopped_batch_counts.collect { key, count ->
                "${key.replace('\t', '/')}=${count}"
            }.sort().join(', ')
            throw new IllegalStateException(
                "Live samples stopped with unequal synchronized BAM batch counts: ${counts}. " +
                "Every live sample must contribute one BAM to every post-startup batch."
            )
        }
        log.info(
            "All ${live_sample_count} live sample watcher(s) have observed STOP after " +
            "${batch_counts[0]} synchronized batch(es); " +
            "live ingress is closed and downstream tasks are draining."
        )
    }
}

Boolean is_stop_path(Path path) {
    return path.name.startsWith('STOP')
}

Boolean is_bam_path(Path path) {
    return path.name.endsWith('.bam')
}

def make_sample_chunk_event(Integer batch_index, Integer sample_index, SampleChunk chunk) -> SampleChunkEvent {
    return record(batch_index: batch_index, sample_index: sample_index, chunk: chunk)
}

def make_sample_chunk(Sample sample, List<Path> bam_paths) -> SampleChunk {
    return record(sample: sample, bam_paths: sort_bam_paths(bam_paths))
}

def get_bam_files_in_dir(Path dir) -> List<Path> {
    return sort_bam_paths(files("${dir}/**.bam", type: 'file'))
}

def get_stop_files_in_dir(Path dir) -> List<Path> {
    return files("${dir}/**", type: 'file')
        .findAll { Path path -> is_stop_path(path) }
        .toSorted { Path left, Path right -> "${left}" <=> "${right}" }
}

def sort_bam_paths(Iterable<Path> bam_paths) -> List<Path> {
    return bam_paths.toSorted { Path left, Path right -> "${left}" <=> "${right}" }
}

String path_key(Path path) {
    return "${path.toAbsolutePath().normalize()}"
}
