#!/usr/bin/env nextflow

nextflow.enable.types = true

include {
    ChunkBAM;
    CollatedChunkBAM;
    CollatedChunkBAMBatch;
    CumulativeCollatedBAM;
    CumulativeCollatedBAMGroup;
    CumulativeCollatedBAMSnapshot;
    QuantifiedSample;
    QuantifiedSampleBatch;
    QuantifiedSampleUpdateBatch;
    SampleBatchSize
} from '../lib/sample.nf'

include { order_batches } from '../lib/util.nf'

include {
    collated_chunk_bam_name;
    cumulative_collated_bam_name;
    differential_expression_results_dir;
    oarfish_counts_file_name;
    oarfish_out_name;
    shell_quote
} from '../modules/generic_helpers.nf'

/**
 * Quantify the cumulative BAM snapshot and rerun differential expression for
 * every complete live batch.
 *
 * A complete batch contains one newly collated chunk for each active sample.
 * Active samples are accumulated and requantified; unchanged final samples
 * reuse their startup quantification. edgeR starts only after every sample in
 * the experiment has a current quantification.
 */
workflow differential_expression {
    take:
        merged_bams: Channel<ChunkBAM>
        batch_sizes: Channel<SampleBatchSize>
        genome: Path
        annotation: Path
        gene_sets: Path

    main:
        collated_bams = collate_chunk_bam(merged_bams)
        Map<Integer, Map> pending_merged_bam_batches = [:]
        nonempty_merged_bam_batches_ch = collated_bams
            .map { CollatedChunkBAM merged_bam ->
                record(
                    batch_index: merged_bam.batch_index,
                    merged_bam: merged_bam
                )
            }
            .join(batch_sizes, by: 'batch_index')
            .map { joined ->
                Map pending_batch = pending_merged_bam_batches[joined.batch_index]
                if (pending_batch == null) {
                    pending_batch = [
                        expected_count: joined.active_sample_count,
                        experiment_sample_count: joined.experiment_sample_count,
                        bams: []
                    ]
                    pending_merged_bam_batches[joined.batch_index] = pending_batch
                }
                if (
                    pending_batch.expected_count != joined.active_sample_count ||
                    pending_batch.experiment_sample_count != joined.experiment_sample_count
                ) {
                    error("Inconsistent size metadata for merged BAM batch ${joined.batch_index}.")
                }
                pending_batch.bams.add(joined.merged_bam)
                if (pending_batch.bams.size() < pending_batch.expected_count) {
                    return null
                }
                if (pending_batch.bams.size() > pending_batch.expected_count) {
                    error(
                        "Merged BAM batch ${joined.batch_index} received more than " +
                        "${pending_batch.expected_count} BAM(s)."
                    )
                }
                pending_merged_bam_batches.remove(joined.batch_index)
                return record(
                    batch_index: joined.batch_index,
                    bams: pending_batch.bams.toSorted { left, right ->
                        "${left.sample.group}/${left.sample.name}" <=>
                            "${right.sample.group}/${right.sample.name}"
                    },
                    experiment_sample_count: pending_batch.experiment_sample_count
                )
            }
            .filter { batch -> batch != null }

        empty_merged_bam_batches_ch = batch_sizes
            .filter { SampleBatchSize batch_size -> batch_size.active_sample_count == 0 }
            .map { SampleBatchSize batch_size ->
                record(
                    batch_index: batch_size.batch_index,
                    bams: [],
                    experiment_sample_count: batch_size.experiment_sample_count
                )
            }

        merged_bam_batches_ch = order_batches(
            nonempty_merged_bam_batches_ch.mix(empty_merged_bam_batches_ch)
        )

        Map<String, List<CollatedChunkBAM>> cumulative_bam_state = [:]
        cumulative_snapshots_ch = merged_bam_batches_ch.map { CollatedChunkBAMBatch batch ->
            accumulate_cumulative_bam_state(cumulative_bam_state, batch)
        }

        snapshot_sizes_ch = cumulative_snapshots_ch.map { CumulativeCollatedBAMSnapshot snapshot ->
                record(
                    batch_index: snapshot.batch_index,
                    active_sample_count: snapshot.sample_bams.size(),
                    experiment_sample_count: snapshot.experiment_sample_count
                )
            }

        Channel<CumulativeCollatedBAMGroup> cumulative_sample_bams_ch = cumulative_snapshots_ch.flatMap { CumulativeCollatedBAMSnapshot snapshot ->
            snapshot.sample_bams
        }
        quantified_samples_ch = oarfish_quant(
            assemble_cumulative_collated_bam(cumulative_sample_bams_ch),
            genome,
            annotation
        )

        Map<Integer, Map> pending_quantified_sample_batches = [:]
        nonempty_quantified_sample_updates_ch = quantified_samples_ch
            .join(snapshot_sizes_ch, by: 'batch_index')
            .map { joined ->
                QuantifiedSample quantified_sample = record(
                    batch_index: joined.batch_index,
                    sample: joined.sample,
                    counts: joined.counts
                )
                Map pending_batch = pending_quantified_sample_batches[joined.batch_index]
                if (pending_batch == null) {
                    pending_batch = [
                        expected_count: joined.active_sample_count,
                        experiment_sample_count: joined.experiment_sample_count,
                        samples: []
                    ]
                    pending_quantified_sample_batches[joined.batch_index] = pending_batch
                }
                if (
                    pending_batch.expected_count != joined.active_sample_count ||
                    pending_batch.experiment_sample_count != joined.experiment_sample_count
                ) {
                    error(
                        "Inconsistent size metadata for quantified sample batch " +
                        "${joined.batch_index}."
                    )
                }
                pending_batch.samples.add(quantified_sample)
                if (pending_batch.samples.size() < pending_batch.expected_count) {
                    return null
                }
                if (pending_batch.samples.size() > pending_batch.expected_count) {
                    error(
                        "Quantified sample batch ${joined.batch_index} received more than " +
                        "${pending_batch.expected_count} sample(s)."
                    )
                }
                pending_quantified_sample_batches.remove(joined.batch_index)
                return record(
                    batch_index: joined.batch_index,
                    samples: pending_batch.samples.toSorted { left, right ->
                        "${left.sample.group}/${left.sample.name}" <=>
                            "${right.sample.group}/${right.sample.name}"
                    },
                    experiment_sample_count: pending_batch.experiment_sample_count
                )
            }
            .filter { batch -> batch != null }

        empty_quantified_sample_updates_ch = snapshot_sizes_ch
            .filter { SampleBatchSize batch_size -> batch_size.active_sample_count == 0 }
            .map { SampleBatchSize batch_size ->
                record(
                    batch_index: batch_size.batch_index,
                    samples: [],
                    experiment_sample_count: batch_size.experiment_sample_count
                )
            }

        ordered_quantified_sample_updates_ch = order_batches(
            nonempty_quantified_sample_updates_ch.mix(empty_quantified_sample_updates_ch)
        )

        Map<String, QuantifiedSample> latest_quantifications = [:]
        quantified_sample_batches_ch = ordered_quantified_sample_updates_ch
            .flatMap { QuantifiedSampleUpdateBatch update_batch ->
                update_batch.samples.each { QuantifiedSample quantified_sample ->
                    latest_quantifications[sample_key(quantified_sample.sample)] = quantified_sample
                }
                if (latest_quantifications.size() < update_batch.experiment_sample_count) {
                    log.info(
                        "Deferring differential expression for batch ${update_batch.batch_index}: " +
                        "${latest_quantifications.size()} of ${update_batch.experiment_sample_count} sample(s) " +
                        "have quantifications."
                    )
                    return []
                }

                QuantifiedSampleBatch quant_batch = record(
                    batch_index: update_batch.batch_index,
                    samples: latest_quantifications.values().toList().toSorted { left, right ->
                        "${left.sample.group}/${left.sample.name}" <=>
                            "${right.sample.group}/${right.sample.name}"
                    }
                )
                return [quant_batch]
            }

        differential_results_ch = run_differential_expression_edgeR(
            quantified_sample_batches_ch,
            gene_sets,
            annotation
        )

    emit:
        quantifications = quantified_samples_ch
        results = differential_results_ch
}

/**
 * Fold one synchronized batch into the cumulative per-sample BAM state.
 */
CumulativeCollatedBAMSnapshot accumulate_cumulative_bam_state(
    Map<String, List<CollatedChunkBAM>> state,
    CollatedChunkBAMBatch batch
) {
    batch.bams.each { CollatedChunkBAM merged_bam ->
        String sample_key = "${merged_bam.sample.group}\t${merged_bam.sample.name}"
        List<CollatedChunkBAM> previous_bams = state.containsKey(sample_key) ? state[sample_key] : []
        state[sample_key] = (previous_bams.findAll { previous_bam ->
            previous_bam.batch_index != merged_bam.batch_index
        } + [merged_bam]).toSorted { left, right ->
            left.batch_index <=> right.batch_index
        }
    }

    List<CumulativeCollatedBAMGroup> sample_bams = batch.bams
        .collect { CollatedChunkBAM updated_bam ->
            List<CollatedChunkBAM> bams = state[sample_key(updated_bam.sample)]
            CollatedChunkBAM latest_bam = bams[-1]
            record(
                batch_index: batch.batch_index,
                sample: latest_bam.sample,
                bams: bams
            )
        }
        .toSorted { left, right ->
            "${left.sample.group}/${left.sample.name}" <=>
                "${right.sample.group}/${right.sample.name}"
        }

    return record(
        batch_index: batch.batch_index,
        sample_bams: sample_bams,
        experiment_sample_count: batch.experiment_sample_count
    )
}

String sample_key(sample) {
    return "${sample.group}\t${sample.name}"
}

/**
 * Strip Oarfish-unsupported tags and collate each newly arrived chunk exactly
 * once. Cumulative snapshots can then be assembled with block-level BAM
 * concatenation instead of repeatedly collating old alignments.
 */
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
            bam: file(collated_chunk_bam_name(chunk_bam.batch_index, chunk_bam.sample))
        )

    script:
        String collated_bam = collated_chunk_bam_name(chunk_bam.batch_index, chunk_bam.sample)
        """
        samtools view -u -x ts ${chunk_bam.bam} |
            samtools collate \
                -o ${collated_bam} \
                -@ ${task.cpus - 1} \
                -
        """
}

/**
 * Assemble all already-collated chunks for one cumulative sample snapshot.
 * samtools cat copies BAM blocks without decompressing and recollating the
 * historical chunks. Read names must be globally unique across chunks.
 */
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
            bam: file(cumulative_collated_bam_name(input_group))
        )

    script:
        String cumulative_bam = cumulative_collated_bam_name(input_group)
        String cumulative_bam_arg = shell_quote(cumulative_bam)
        String bam_args = input_group.bams*.bam.collect { bam -> shell_quote(bam.toString()) }.join(' ')

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

/**
 * Quantify one cumulative sample snapshot from a name-collated genome BAM.
 */
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
            counts: file(oarfish_counts_file_name(merged_bam.batch_index, merged_bam.sample.name))
        )

    script:
        String output_prefix = oarfish_out_name(merged_bam.batch_index, merged_bam.sample.name)
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

/**
 * Rebuild the full count matrix and rerun edgeR for one live batch.
 */
process run_differential_expression_edgeR {
    label 'seq_lm_dea'
    container 'rnabioinfo/seq_lm_dea_r:v1.0.0'
    cpus 1
    maxForks 1

    input:
        quant_batch: QuantifiedSampleBatch
        gene_sets: Path
        annotation: Path

    stage:
        stageAs quant_batch.samples*.counts, 'quant/input?.quant'

    output:
        record(
            batch_index: quant_batch.batch_index,
            results: file(differential_expression_results_dir(quant_batch.batch_index))
        )

    script:
        String results_dir = differential_expression_results_dir(quant_batch.batch_index)
        String manifest_rows = quant_batch.samples.withIndex().collect { QuantifiedSample sample, Integer index ->
            [
                de_manifest_field(sample.sample.name),
                de_manifest_field(sample.sample.group),
                "quant/input${index + 1}.quant"
            ].join('\t')
        }.join('\n')
        String quoted_manifest_rows = shell_quote(manifest_rows)
        """
        printf 'name\\tgroup\\tcount_file\\n' > quant_manifest.tsv
        printf '%s\\n' ${quoted_manifest_rows} >> quant_manifest.tsv

        mkdir ${results_dir}

        edgeR-analysis \
            --quant_manifest quant_manifest.tsv \
            --output_dir ${results_dir} \
            --gene_sets ${gene_sets} \
            --annotation ${annotation} \
            --lfc ${params.de_lfc_cutoff}
        """
}

String de_manifest_field(Object value) {
    return "${value}".replaceAll(/[\t\r\n]+/, ' ')
}
