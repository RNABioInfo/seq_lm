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
        Map<Integer, Object> pending_batches = [:]
        Integer next_batch_index = 0
        ordered_batches_ch = batches.flatMap { batch ->
            if (pending_batches.containsKey(batch.batch_index)) {
                error("Received duplicate completed batch index ${batch.batch_index}.")
            }
            pending_batches[batch.batch_index] = batch

            List<Integer> ready_indices = pending_batches
                .keySet()
                .toSorted()
                .withIndex()
                .findAll { entry -> entry[0] == next_batch_index + entry[1] }
                .collect { entry -> entry[0] as Integer }
            List ready_batches = ready_indices.collect { Integer batch_index ->
                pending_batches.remove(batch_index)
            }
            next_batch_index += ready_indices.size()
            return ready_batches
        }

    emit:
        ordered_batches_ch
}

String get_sample_path(Map meta) {
    return "${meta['runName']}/${meta['replicateName']}"
}

Path get_seq_summary_file(Path bam_file) {
    Path summary_file = file("${bam_file.parent}/seq_summary.txt")
    if (summary_file.exists()) {
        return summary_file
    }
    return file("$projectDir/data/OPTIONAL_FILE")
}

Map get_sequencing_arguments(Path _run_dir) {
    Map args = [:]
    args['experiment_id'] = params.ex_name
    args['run_id'] = params.ex_run_number
    args['kit'] = params.ex_kit
    if (!params.ex_special_alignment) {
        args['reference_genome'] = params.reference_genome
    }
    args['basecall_config'] = params.ex_basecall_config
    return args
}
