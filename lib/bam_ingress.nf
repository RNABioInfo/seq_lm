include { Sample ; SampleChunk ; SampleBatch } from './sample.nf'

nextflow.enable.types = true

record BamIngressArgs {
    live_analysis: Boolean
    timeline_analysis: Boolean?
    sample_sheet_path: Path?
    bam_poll_interval_ms: Integer
    bam_stability_polls: Integer
}

record LiveSampleEvent {
    batch_index: Integer
    sample_index: Integer
    sample: Sample
    chunk: SampleChunk?
    is_stop: Boolean
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
def bam_ingress(ingress_args) {
    def samples: List<Sample> = get_samples(ingress_args)
    validate_samples(samples)
    return mixed_analysis_ingress(samples, ingress_args, samples.size())
}

/**
 * Start ingress for only the samples that still require upstream processing.
 * experiment_sample_count retains the complete sample-sheet size so restored
 * quantifications can satisfy the remaining differential-analysis inputs.
 */
def bam_ingress(samples: List<Sample>, ingress_args, experiment_sample_count: Integer) {
    return mixed_analysis_ingress(samples, ingress_args, experiment_sample_count)
}

def get_samples(ingress_args) -> List<Sample> {
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

def parse_sample_sheet(ingress_args) -> List<Sample> {
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
                is_live: parse_is_live(row.is_live, row.alias),
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

def validate_samples(samples: List<Sample>) {
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

def parse_is_live(value: Object, alias: Object) -> Boolean {
    def normalized: String = value == null ? '' : "${value}".trim().toLowerCase()
    if (!normalized || normalized == 'true') {
        return true
    }
    if (normalized == 'false') {
        return false
    }
    error(
        "Invalid is_live value '${value}' for sample '${alias}'. " + "Expected true, false, or a blank value."
    )
}

def sample_key(sample) -> String {
    return "${sample.group}\t${sample.name}"
}

def sample_label(sample) -> String {
    return "${sample.group}/${sample.name}"
}

def mixed_analysis_ingress(samples: List<Sample>, ingress_args, experiment_sample_count: Integer) {
    def live_samples: List<Sample> = samples.findAll { sample ->
        ingress_args.live_analysis && sample.is_live
    }
    log.info(
        "Running workflow with ${live_samples.size()} live sample(s) and " + "${samples.size() - live_samples.size()} final sample(s)."
    )

    def existing_bams_by_sample: Map<String,List<Path>> = [:]
    def existing_chunks: List<SampleChunk> = samples.collect { sample ->
        def existing_bams: List<Path> = get_bam_files_in_dir(sample.bam_dir)
        existing_bams_by_sample[sample_key(sample)] = existing_bams
        def effectively_live: Boolean = ingress_args.live_analysis && sample.is_live
        log.info(
            "BAM ingress startup scan for sample '${sample_label(sample)}' " + "found ${existing_bams.size()} existing BAM file(s) in ${sample.bam_dir}; " + "effective mode=${effectively_live ? 'live' : 'final'}."
        )
        if (!effectively_live && existing_bams.empty) {
            error(
                "Final sample '${sample_label(sample)}' must contain at least one BAM file at startup: " + "${sample.bam_dir}"
            )
        }
        make_sample_chunk(sample, existing_bams)
    }

    def startup_batch: SampleBatch = record(
        batch_index: 0,
        chunks: existing_chunks,
        experiment_sample_count: experiment_sample_count,
    )
    if (live_samples.empty) {
        log.info('No effective live samples remain; BAM ingress will emit startup batch 0 and close.')
        return channel.of(startup_batch)
    }

    def ch_events = channel.empty()
    live_samples.each { sample ->
        def sample_index: Integer = samples.indexOf(sample)
        ch_events = ch_events.mix(
            make_live_analysis_ingress_channel(
                sample,
                sample_index,
                existing_bams_by_sample[sample_key(sample)],
                live_samples.size(),
                ingress_args.bam_poll_interval_ms,
                ingress_args.bam_stability_polls,
            )
        )
    }

    def live_sample_indices: List<Integer> = live_samples.collect { sample -> samples.indexOf(sample) }
    def pending_batch_events: Map<Integer,Map<Integer,LiveSampleEvent>> = [:]
    def stopped_final_batches: Map<Integer,Integer> = [:]
    def next_batch_index: Integer = 1
    def live_batches_ch = ch_events
        .flatMap { event ->
            if (event.is_stop) {
                if (stopped_final_batches.containsKey(event.sample_index)) {
                    error("Received duplicate terminal event for sample '${sample_label(event.sample)}'.")
                }
                stopped_final_batches[event.sample_index] = event.batch_index
                log.info(
                    "Sample '${sample_label(event.sample)}' left live batch synchronization after local batch ${event.batch_index}."
                )
            }
            else {
                def events: Map<Integer,LiveSampleEvent> = pending_batch_events.computeIfAbsent(event.batch_index) { [:] }
                if (events.put(event.sample_index, event) != null) {
                    error("Live BAM batch ${event.batch_index} received a duplicate event for sample '${sample_label(event.sample)}'.")
                }
            }

            def highest_pending_index: Integer = pending_batch_events.empty ? next_batch_index - 1 : pending_batch_events.keySet().max()
            def candidate_indices: List<Integer> = highest_pending_index < next_batch_index
                ? []
                : (next_batch_index..highest_pending_index).toList()
            def readiness: Map = candidate_indices.inject([open: true, indices: []]) { state: Map, candidate_index: Integer ->
                if (!state.open) {
                    return state
                }
                def required_indices: List<Integer> = live_sample_indices.findAll { sample_index: Integer ->
                    !stopped_final_batches.containsKey(sample_index) || stopped_final_batches[sample_index] >= candidate_index
                }
                def events: Map<Integer,LiveSampleEvent> = pending_batch_events[candidate_index] ?: [:]
                if (required_indices.empty || !required_indices.every { sample_index: Integer -> events.containsKey(sample_index) }) {
                    state.open = false
                    return state
                }
                state.indices.add(candidate_index)
                return state
            }
            def ready_batches: List<SampleBatch> = (readiness.indices as List<Integer>).collect { ready_index: Integer ->
                def required_indices: List<Integer> = live_sample_indices.findAll { sample_index: Integer ->
                    !stopped_final_batches.containsKey(sample_index) || stopped_final_batches[sample_index] >= ready_index
                }
                def events: Map<Integer,LiveSampleEvent> = pending_batch_events[ready_index]
                def chunks: List<SampleChunk> = required_indices
                    .toSorted()
                    .collect { sample_index: Integer -> events[sample_index].chunk }
                pending_batch_events.remove(ready_index)
                def batch: SampleBatch = record(
                    batch_index: ready_index,
                    chunks: chunks,
                    experiment_sample_count: experiment_sample_count,
                )
                log.info(
                    "Emitting live BAM batch ${batch.batch_index} from ${chunks.size()} remaining sample(s): " + chunks.collect { chunk -> "${sample_label(chunk.sample)}=${chunk.bam_paths.size()} BAM(s)" }.join(', ')
                )
                return batch
            }
            next_batch_index += ready_batches.size()
            return ready_batches
        }

    return channel.of(startup_batch).concat(live_batches_ch)
}

def make_live_analysis_ingress_channel(
    sample,
    sample_index: Integer,
    initial_bams: List<Path>,
    live_sample_count: Integer,
    poll_interval_ms: Integer,
    stability_polls: Integer
) {
    def existing_stop_files: List<Path> = get_stop_files_in_dir(sample.bam_dir)
    def batch_index: Integer = 0
    def accepted_bam_paths: Set<String> = new LinkedHashSet<String>()
    initial_bams.each { bam_path: Path -> accepted_bam_paths.add(path_key(bam_path)) }

    if (!existing_stop_files.empty) {
        log.warn(
            "Sample '${sample_label(sample)}' has existing STOP file(s) " + "at startup: ${existing_stop_files}. No future BAM poller will be opened for this sample."
        )
        log_sample_stop(sample, existing_stop_files[0], 'existing', batch_index, live_sample_count)
        return channel.of(make_live_stop_event(batch_index, sample_index, sample))
    }

    log.info(
        "Starting live BAM poller for sample '${sample_label(sample)}' in ${sample.bam_dir}; " +
            "interval=${poll_interval_ms} ms, required stable polls=${stability_polls}."
    )

    return LiveBamWatcher
        .poll(
            sample.bam_dir,
            initial_bams,
            poll_interval_ms as Long,
            stability_polls,
            {
                log.info(
                    "Live BAM poller is active for sample '${sample_label(sample)}' under ${sample.bam_dir}."
                )
            },
        ) { error: Throwable ->
            log.error(
                "Live BAM poller failed for sample '${sample_label(sample)}'; " + "closing this sample stream.",
                error,
            )
        }
        .map { event -> make_live_path_event(event.source as String, event.path as Path) }
        .flatMap { event ->
            if (is_stop_path(event.path)) {
                log_sample_stop(sample, event.path, event.source, batch_index, live_sample_count)
                return [make_live_stop_event(batch_index, sample_index, sample)]
            }
            if (event.source != "poll") {
                log.info(
                    "Ignoring non-poll live BAM path event for sample '${sample_label(sample)}': " + "source=${event.source}, path=${event.path}."
                )
                return []
            }
            if (!is_bam_path(event.path)) {
                log.info(
                    "Ignoring non-BAM live path event for sample '${sample_label(sample)}': " + "source=${event.source}, path=${event.path}."
                )
                return []
            }
            if (!accepted_bam_paths.add(path_key(event.path))) {
                log.info(
                    "Ignoring duplicate live BAM path event for sample '${sample_label(sample)}': " + "source=${event.source}, path=${event.path}."
                )
                return []
            }

            batch_index += 1
            log.info(
                "Accepted live BAM event for sample '${sample_label(sample)}': " + "batch_index=${batch_index}, source=${event.source}, path=${event.path}."
            )
            return [make_sample_chunk_event(
                batch_index,
                sample_index,
                make_sample_chunk(sample, [event.path]),
            )]
        }
}

def make_live_path_event(source: String, path: Path) {
    return record(source: source, path: path)
}

def log_sample_stop(sample, stop_path: Path, source: String, batch_count: Integer, live_sample_count: Integer) {
    log.info(
        "Observed STOP for sample '${sample_label(sample)}': " + "source=${source}, path=${stop_path}, final_local_batch=${batch_count}; this poller is closing (${live_sample_count} initially live sample(s))."
    )
}

def is_stop_path(path: Path) -> Boolean {
    return path.name.startsWith('STOP')
}

def is_bam_path(path: Path) -> Boolean {
    return path.name.endsWith('.bam')
}

def make_sample_chunk_event(batch_index: Integer, sample_index: Integer, chunk) -> LiveSampleEvent {
    return record(batch_index: batch_index, sample_index: sample_index, sample: chunk.sample, chunk: chunk, is_stop: false)
}

def make_live_stop_event(batch_index: Integer, sample_index: Integer, sample) -> LiveSampleEvent {
    return record(batch_index: batch_index, sample_index: sample_index, sample: sample, chunk: null, is_stop: true)
}

def make_sample_chunk(sample, bam_paths: List<Path>) {
    return record(sample: sample, bam_paths: sort_bam_paths(bam_paths))
}

def get_bam_files_in_dir(dir: Path) -> List<Path> {
    return sort_bam_paths(files("${dir}/**.bam", type: 'file'))
}

def get_stop_files_in_dir(dir: Path) -> List<Path> {
    return files("${dir}/**", type: 'file')
        .findAll { path: Path -> is_stop_path(path) }
        .toSorted { left: Path, right: Path -> "${left}" <=> "${right}" }
}

def sort_bam_paths(bam_paths: Iterable<Path>) -> List<Path> {
    return bam_paths.toSorted { left: Path, right: Path -> "${left}" <=> "${right}" }
}

def path_key(path: Path) -> String {
    return "${path.toAbsolutePath().normalize()}"
}
