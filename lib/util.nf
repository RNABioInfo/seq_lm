nextflow.enable.types = true

/**
 * Restore ingress order after asynchronous processes complete whole batches.
 *
 * Every batch index, including batches with no process inputs, must be present
 * on the input channel. Empty placeholders allow the buffer to advance without
 * letting a later batch update cumulative state too early.
 */
workflow order_batches {
    take:
    batches: Channel

    main:
    def pending_batches: Map<Integer,Object> = [:]
    def next_batch_index: Integer = 0
    ordered_batches_ch = batches.flatMap { batch ->
        if (pending_batches.containsKey(batch.batch_index)) {
            error("Received duplicate completed batch index ${batch.batch_index}.")
        }
        pending_batches[batch.batch_index] = batch

        def ready_indices: List<Integer> = pending_batches
            .keySet()
            .toSorted()
            .withIndex()
            .findAll { entry -> entry[0] == next_batch_index + entry[1] }
            .collect { entry -> entry[0] as Integer }
        def ready_batches: List = ready_indices.collect { batch_index: Integer ->
            pending_batches.remove(batch_index)
        }
        next_batch_index += ready_indices.size()
        return ready_batches
    }

    emit:
    ordered_batches_ch
}

def get_sample_path(meta: Map) -> String {
    return "${meta['runName']}/${meta['replicateName']}"
}

def get_seq_summary_file(bam_file: Path) -> Path {
    def summary_file: Path = file("${bam_file.parent}/seq_summary.txt")
    if (summary_file.exists()) {
        return summary_file
    }
    return file("${projectDir}/data/OPTIONAL_FILE")
}

def get_sequencing_arguments(_run_dir: Path) -> Map {
    def args: Map = [:]
    args['experiment_id'] = params.ex_name
    args['run_id'] = params.ex_run_number
    args['kit'] = params.ex_kit
    if (!params.ex_special_alignment) {
        args['reference_genome'] = params.reference_genome
    }
    args['basecall_config'] = params.ex_basecall_config
    return args
}
