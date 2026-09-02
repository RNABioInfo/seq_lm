#!/usr/bin/env nextflow

nextflow.enable.types = true

include {
    ChunkBAM ;
    CollatedChunkBAM ;
    CumulativeCollatedBAM ;
    CumulativeCollatedBAMGroup ;
    QuantifiedSample ;
    QuantifiedSampleBatch ;
    SampleBatchSize
} from '../lib/sample.nf'

include { order_batches } from '../lib/util.nf'

include {
    collated_chunk_bam_name ;
    cumulative_collated_bam_name ;
    oarfish_counts_file_name ;
    oarfish_out_name ;
    shell_quote
} from '../modules/generic_helpers.nf'

/**
 * Build cumulative per-sample Oarfish quantifications independently of DEA.
 * Complete batches contain the newest quantification for every experiment
 * sample, including restored or already stopped samples.
 */
workflow quantification {
    take:
    merged_bams: Channel<ChunkBAM>
    batch_sizes: Channel<SampleBatchSize>
    restored_quantifications: List<QuantifiedSample>
    genome: Path
    annotation: Path

    main:
    collated_bams = collate_chunk_bam(merged_bams)
    def pending_merged_bam_batches: Map<Integer, Map> = [:]
    nonempty_merged_bam_batches_ch = collated_bams
        .map { merged_bam ->
            record(
                batch_index: merged_bam.batch_index,
                merged_bam: merged_bam,
            )
        }
        .join(batch_sizes, by: 'batch_index')
        .map { joined ->
            def pending_batch: Map = pending_merged_bam_batches[joined.batch_index]
            if (pending_batch == null) {
                pending_batch = [expected_count: joined.active_sample_count, experiment_sample_count: joined.experiment_sample_count, bams: []]
                pending_merged_bam_batches[joined.batch_index] = pending_batch
            }
            if (pending_batch.expected_count != joined.active_sample_count || pending_batch.experiment_sample_count != joined.experiment_sample_count) {
                error("Inconsistent size metadata for merged BAM batch ${joined.batch_index}.")
            }
            pending_batch.bams.add(joined.merged_bam)
            if (pending_batch.bams.size() < pending_batch.expected_count) {
                return null
            }
            if (pending_batch.bams.size() > pending_batch.expected_count) {
                error(
                    "Merged BAM batch ${joined.batch_index} received more than " + "${pending_batch.expected_count} BAM(s)."
                )
            }
            pending_merged_bam_batches.remove(joined.batch_index)
            return record(
                batch_index: joined.batch_index,
                bams: pending_batch.bams.toSorted { left, right ->
                    "${left.sample.group}/${left.sample.name}" <=> "${right.sample.group}/${right.sample.name}"
                },
                experiment_sample_count: pending_batch.experiment_sample_count,
            )
        }
        .filter { batch -> batch != null }

    empty_merged_bam_batches_ch = batch_sizes
        .filter { batch_size -> batch_size.active_sample_count == 0 }
        .map { batch_size ->
            record(
                batch_index: batch_size.batch_index,
                bams: [],
                experiment_sample_count: batch_size.experiment_sample_count,
            )
        }

    merged_bam_batches_ch = order_batches(
        nonempty_merged_bam_batches_ch.mix(empty_merged_bam_batches_ch)
    )

    def cumulative_bam_state: Map<String, List<CollatedChunkBAM>> = [:]
    cumulative_snapshots_ch = merged_bam_batches_ch.map { batch ->
        accumulate_cumulative_bam_state(cumulative_bam_state, batch)
    }

    snapshot_sizes_ch = cumulative_snapshots_ch.map { snapshot ->
        record(
            batch_index: snapshot.batch_index,
            active_sample_count: snapshot.sample_bams.size(),
            experiment_sample_count: snapshot.experiment_sample_count,
        )
    }

    def cumulative_sample_bams_ch: Channel<CumulativeCollatedBAMGroup> = cumulative_snapshots_ch.flatMap { snapshot ->
        snapshot.sample_bams
    }
    quantified_samples_ch = oarfish_quant(
        assemble_cumulative_collated_bam(cumulative_sample_bams_ch),
        genome,
        annotation,
    )

    def pending_quantified_sample_batches: Map<Integer, Map> = [:]
    nonempty_quantified_sample_updates_ch = quantified_samples_ch
        .join(snapshot_sizes_ch, by: 'batch_index')
        .map { joined ->
            def quantified_sample = record(
                batch_index: joined.batch_index,
                sample: joined.sample,
                counts: joined.counts,
            )
            def pending_batch: Map = pending_quantified_sample_batches[joined.batch_index]
            if (pending_batch == null) {
                pending_batch = [expected_count: joined.active_sample_count, experiment_sample_count: joined.experiment_sample_count, samples: []]
                pending_quantified_sample_batches[joined.batch_index] = pending_batch
            }
            if (pending_batch.expected_count != joined.active_sample_count || pending_batch.experiment_sample_count != joined.experiment_sample_count) {
                error(
                    "Inconsistent size metadata for quantified sample batch " + "${joined.batch_index}."
                )
            }
            pending_batch.samples.add(quantified_sample)
            if (pending_batch.samples.size() < pending_batch.expected_count) {
                return null
            }
            if (pending_batch.samples.size() > pending_batch.expected_count) {
                error(
                    "Quantified sample batch ${joined.batch_index} received more than " + "${pending_batch.expected_count} sample(s)."
                )
            }
            pending_quantified_sample_batches.remove(joined.batch_index)
            return record(
                batch_index: joined.batch_index,
                samples: pending_batch.samples.toSorted { left, right ->
                    "${left.sample.group}/${left.sample.name}" <=> "${right.sample.group}/${right.sample.name}"
                },
                experiment_sample_count: pending_batch.experiment_sample_count,
            )
        }
        .filter { batch -> batch != null }

    empty_quantified_sample_updates_ch = snapshot_sizes_ch
        .filter { batch_size -> batch_size.active_sample_count == 0 }
        .map { batch_size ->
            record(
                batch_index: batch_size.batch_index,
                samples: [],
                experiment_sample_count: batch_size.experiment_sample_count,
            )
        }

    ordered_quantified_sample_updates_ch = order_batches(
        nonempty_quantified_sample_updates_ch.mix(empty_quantified_sample_updates_ch)
    )

    def latest_quantifications: Map<String, QuantifiedSample> = restored_quantifications.collectEntries { quantified_sample ->
        [(sample_key(quantified_sample.sample)): quantified_sample]
    }
    def next_report_sequence: Integer = 0
    quantified_sample_batches_ch = ordered_quantified_sample_updates_ch.flatMap { update_batch ->
        update_batch.samples.each { quantified_sample ->
            latest_quantifications[sample_key(quantified_sample.sample)] = quantified_sample
        }
        if (latest_quantifications.size() < update_batch.experiment_sample_count) {
            log.info(
                "Deferring quantification summary for batch ${update_batch.batch_index}: " + "${latest_quantifications.size()} of ${update_batch.experiment_sample_count} sample(s) " + "have quantifications."
            )
            return []
        }

        def quant_batch = record(
            batch_index: update_batch.batch_index,
            report_sequence: next_report_sequence,
            samples: latest_quantifications.values().toList().toSorted { left, right ->
                "${left.sample.group}/${left.sample.name}" <=> "${right.sample.group}/${right.sample.name}"
            },
        )
        next_report_sequence += 1
        return [quant_batch]
    }

    biotype_map_ch = extract_transcript_biotype_map(annotation)
    transcript_biotype_report_batches_ch = summarize_transcript_biotypes(
        quantified_sample_batches_ch,
        biotype_map_ch,
    )

    emit:
    quantifications = quantified_samples_ch
    batches = quantified_sample_batches_ch
    report_batches = transcript_biotype_report_batches_ch
}

/** Fold one synchronized batch into the cumulative per-sample BAM state. */
def accumulate_cumulative_bam_state(state: Map<String, List<CollatedChunkBAM>>, batch) {
    batch.bams.each { merged_bam ->
        def key: String = sample_key(merged_bam.sample)
        def previous_bams: List<CollatedChunkBAM> = state.containsKey(key) ? state[key] : []
        state[key] = (previous_bams.findAll { previous_bam ->
            previous_bam.batch_index != merged_bam.batch_index
        } + [merged_bam]).toSorted { left, right ->
            left.batch_index <=> right.batch_index
        }
    }

    def sample_bams: List<CumulativeCollatedBAMGroup> = batch.bams
        .collect { updated_bam ->
            def bams: List<CollatedChunkBAM> = state[sample_key(updated_bam.sample)]
            def latest_bam = bams[-1]
            record(
                batch_index: batch.batch_index,
                sample: latest_bam.sample,
                bams: bams,
            )
        }
        .toSorted { left, right ->
            "${left.sample.group}/${left.sample.name}" <=> "${right.sample.group}/${right.sample.name}"
        }

    return record(
        batch_index: batch.batch_index,
        sample_bams: sample_bams,
        experiment_sample_count: batch.experiment_sample_count,
    )
}

def sample_key(sample) -> String {
    return "${sample.group}\t${sample.name}"
}

/** Strip unsupported tags and collate every newly arrived chunk once. */
process collate_chunk_bam {
    label 'seq_lm_dea'
    container 'rnabioinfo/seq_lm_samtools:v1.0.0'
    cpus 2

    input:
    chunk_bam: ChunkBAM

    output:
    record(
        batch_index: chunk_bam.batch_index,
        sample: chunk_bam.sample,
        bam: file(collated_chunk_bam_name(chunk_bam.batch_index, chunk_bam.sample)),
    )

    script:
    def collated_bam: String = collated_chunk_bam_name(chunk_bam.batch_index, chunk_bam.sample)
    """
        samtools view -u -x ts ${chunk_bam.bam} |
            samtools collate \
                -o ${collated_bam} \
                -@ ${task.cpus - 1} \
                -
        """
}

/** Assemble already-collated chunks into one cumulative sample BAM. */
process assemble_cumulative_collated_bam {
    label 'seq_lm_dea'
    container 'rnabioinfo/seq_lm_samtools:v1.0.0'
    cpus 1

    input:
    input_group: CumulativeCollatedBAMGroup

    stage:
    stageAs input_group.bams*.bam, 'collated_chunk?.bam'

    output:
    record(
        batch_index: input_group.batch_index,
        sample: input_group.sample,
        bam: file(cumulative_collated_bam_name(input_group)),
    )

    script:
    def cumulative_bam: String = cumulative_collated_bam_name(input_group)
    def cumulative_bam_arg: String = shell_quote(cumulative_bam)
    def bam_args: String = input_group.bams*.bam.collect { bam -> shell_quote(bam.toString()) }.join(' ')

    if (input_group.bams.size() == 1) {
        return """
            ln -s -- ${bam_args} ${cumulative_bam_arg}
            """
    }

    return """
        printf '%s\\n' ${bam_args} > bams.txt
        samtools cat -o ${cumulative_bam_arg} -b bams.txt
        """
}

/** Quantify one cumulative sample snapshot from a name-collated genome BAM. */
process oarfish_quant {
    label 'seq_lm_dea'
    container 'rnabioinfo/seq_lm_oarfish:v1.0.0'
    cpus 4

    input:
    merged_bam: CumulativeCollatedBAM
    genome: Path
    annotation: Path

    output:
    record(
        batch_index: merged_bam.batch_index,
        sample: merged_bam.sample,
        counts: file(oarfish_counts_file_name(merged_bam.batch_index, merged_bam.sample.name)),
    )

    script:
    def output_prefix: String = oarfish_out_name(merged_bam.batch_index, merged_bam.sample.name)
    """
        oarfish \
            -j ${task.cpus} \
            --genome-alignments ${merged_bam.bam} \
            --annotation ${annotation} \
            --genome-fasta ${genome} \
            --filter-group no-filters \
            -d fw \
            --output ${output_prefix}
        """
}

/** Extract the transcript/alias-to-biotype map once per annotation. */
process extract_transcript_biotype_map {
    label 'seq_lm_qc'
    container 'rnabioinfo/seq_lm_report:v1.0.0'
    cpus 1

    input:
    annotation: Path

    output:
    biotype_map = file('transcript_biotype_map.tsv')

    script:
    """
        workflow-glue transcript_biotypes \
            --mode map \
            --annotation ${annotation} \
            --output transcript_biotype_map.tsv
        """
}

/** Summarize the latest complete Oarfish cohort into fixed biotype classes. */
process summarize_transcript_biotypes {
    label 'seq_lm_qc'
    container 'rnabioinfo/seq_lm_report:v1.0.0'
    cpus 1

    input:
    quant_batch: QuantifiedSampleBatch
    biotype_map: Path

    stage:
    stageAs quant_batch.samples*.counts, 'quant/input?.quant'

    output:
    record(
        batch_index: quant_batch.batch_index,
        report_sequence: quant_batch.report_sequence,
        biotypes: file('transcript_biotype_fractions.tsv'),
    )

    script:
    def manifest_rows: String = quant_batch.samples
        .withIndex()
        .collect { sample, index: Integer ->
            [quant_manifest_field(sample.sample.name), quant_manifest_field(sample.sample.group), "input${index + 1}.quant"].join('\t')
        }
        .join('\n')
    def quoted_manifest_rows: String = shell_quote(manifest_rows)
    """
        printf 'name\\tgroup\\tcount_file\\n' > quant_manifest.tsv
        printf '%s\\n' ${quoted_manifest_rows} >> quant_manifest.tsv

        workflow-glue transcript_biotypes \
            --mode summarize \
            --manifest quant_manifest.tsv \
            --counts-dir quant \
            --mapping ${biotype_map} \
            --output transcript_biotype_fractions.tsv
        """
}

def quant_manifest_field(value: Object) -> String {
    return "${value}".replaceAll(/[\t\r\n]+/, ' ')
}
