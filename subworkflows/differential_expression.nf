#!/usr/bin/env nextflow

nextflow.enable.types = true

include {
    ChunkBAM ;
    CollatedChunkBAM ;
    CollatedChunkBAMBatch ;
    CumulativeCollatedBAM ;
    CumulativeCollatedBAMGroup ;
    CumulativeCollatedBAMSnapshot ;
    QuantifiedSample ;
    QuantifiedSampleBatch ;
    QuantifiedSampleUpdateBatch ;
    SampleBatchSize
} from '../lib/sample.nf'

include { order_batches } from '../lib/util.nf'

include {
    collated_chunk_bam_name ;
    cumulative_collated_bam_name ;
    differential_expression_results_dir ;
    oarfish_counts_file_name ;
    oarfish_out_name ;
    optional_file ;
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
    restored_quantifications: List<QuantifiedSample>
    first_analysis_index: Integer
    genome: Path
    annotation: Path
    gene_sets: Path
    gene_set_enrichment: Boolean

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
    def next_analysis_index: Integer = first_analysis_index
    quantified_sample_batches_ch = ordered_quantified_sample_updates_ch.flatMap { update_batch ->
        update_batch.samples.each { quantified_sample ->
            latest_quantifications[sample_key(quantified_sample.sample)] = quantified_sample
        }
        if (latest_quantifications.size() < update_batch.experiment_sample_count) {
            log.info(
                "Deferring differential expression for batch ${update_batch.batch_index}: " + "${latest_quantifications.size()} of ${update_batch.experiment_sample_count} sample(s) " + "have quantifications."
            )
            return []
        }

        def quant_batch = record(
            batch_index: update_batch.batch_index,
            analysis_index: next_analysis_index,
            samples: latest_quantifications.values().toList().toSorted { left, right ->
                "${left.sample.group}/${left.sample.name}" <=> "${right.sample.group}/${right.sample.name}"
            },
        )
        next_analysis_index += 1
        return [quant_batch]
    }

    edgeR_results_ch = run_differential_expression_edgeR(
        quantified_sample_batches_ch,
        gene_sets,
        annotation,
    )
    if (gene_set_enrichment) {
        analysis_results_ch = run_gene_set_variation_analysis(edgeR_results_ch)
    }
    else {
        analysis_results_ch = edgeR_results_ch
    }

    emit:
    quantifications = quantified_samples_ch
    results = analysis_results_ch
}

/**
 * Fold one synchronized batch into the cumulative per-sample BAM state.
 */
def accumulate_cumulative_bam_state(state: Map<String, List<CollatedChunkBAM>>, batch) {
    batch.bams.each { merged_bam ->
        def sample_key: String = "${merged_bam.sample.group}\t${merged_bam.sample.name}"
        def previous_bams: List<CollatedChunkBAM> = state.containsKey(sample_key) ? state[sample_key] : []
        state[sample_key] = (previous_bams.findAll { previous_bam ->
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
        analysis_index: quant_batch.analysis_index,
        results: file(differential_expression_results_dir(quant_batch.analysis_index)),
    )

    script:
    def results_dir: String = differential_expression_results_dir(quant_batch.analysis_index)
    def gene_set_args: String = gene_sets.name == optional_file().name ? '' : "--gene_sets ${gene_sets}"
    def manifest_rows: String = quant_batch.samples
        .withIndex()
        .collect { sample, index: Integer ->
            [de_manifest_field(sample.sample.name), de_manifest_field(sample.sample.group), "quant/input${index + 1}.quant"].join('\t')
        }
        .join('\n')
    def quoted_manifest_rows: String = shell_quote(manifest_rows)
    """
        printf 'name\\tgroup\\tcount_file\\n' > quant_manifest.tsv
        printf '%s\\n' ${quoted_manifest_rows} >> quant_manifest.tsv

        mkdir ${results_dir}

        edgeR-analysis \
            --quant_manifest quant_manifest.tsv \
            --output_dir ${results_dir} \
            ${gene_set_args} \
            --annotation ${annotation} \
            --lfc ${params.de_lfc_cutoff}
        """
}

/**
 * Score resolved gene sets per sample and test score differences after edgeR.
 * The copied edgeR tree and GSVA additions form one complete batch result.
 */
process run_gene_set_variation_analysis {
    label 'seq_lm_gsva'
    container 'rnabioinfo/seq_lm_gsva:v1.1.0'
    cpus 1
    maxForks 1

    input:
    differential_result: Map

    stage:
    stageAs differential_result.results, 'edgeR_results'

    output:
    record(
        batch_index: differential_result.batch_index,
        analysis_index: differential_result.analysis_index,
        results: file(differential_expression_results_dir(differential_result.analysis_index)),
    )

    script:
    def results_dir: String = differential_expression_results_dir(differential_result.analysis_index)
    """
        mkdir ${results_dir}
        cp -R edgeR_results/. ${results_dir}/

        gsva-analysis \
            --feature_counts ${results_dir}/feature_counts.tsv \
            --sample_metadata ${results_dir}/sample_metadata.tsv \
            --gene_set_resolution ${results_dir}/gene_set_resolution.tsv \
            --output_dir ${results_dir}
        """
}

def de_manifest_field(value: Object) -> String {
    return "${value}".replaceAll(/[\t\r\n]+/, ' ')
}
