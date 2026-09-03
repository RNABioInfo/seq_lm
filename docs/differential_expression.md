# Differential Expression Workflow

Quantification, differential expression, and gene-set enrichment have separate
activation rules:

* Supplying both `--reference_genome` and `--reference_annotation` enables
  cumulative Oarfish quantification and transcript-biotype QC. Supplying only
  one reference input is an error.
* `--differential_expression` enables edgeR analysis of completed Oarfish
  quantification batches and therefore requires both reference inputs.
* `--gene_set_enrichment` enables fry and GSVA after edgeR. It requires both
  `--differential_expression` and `--gene_sets`.

Both analysis switches default to `true`, preserving the complete analysis
workflow when references are provided. Set both to `false` and omit both
references for quality control only. Set `--differential_expression=false`
while retaining both references for QC plus quantification and transcript-
biotype composition. The GMT input is optional whenever gene-set enrichment is
disabled.

The standalone quantification workflow refreshes after every complete live BAM
batch. It uses the same chunk stream as quality control, but quantification is
cumulative rather than chunk-local:

1. BAMs belonging to the same sample chunk are concatenated without sorting or
   indexing. A single input BAM is reused directly.
2. Each new chunk is stripped of its `ts` tags and collated by read name exactly
   once with `samtools collate`.
3. For each live sample, the already-collated chunks seen through the current
   batch are assembled into a temporary cumulative BAM with block-level
   `samtools cat`. Final samples are processed only at startup.
4. Oarfish 0.10 quantifies the cumulative genome alignments against
   `--reference_annotation`, using `--reference_genome` for soft-clip
   rescue.
5. The latest quantification is retained for every sample. Each live batch
   replaces entries for samples contributing a new chunk, reuses entries for
   final or already stopped samples, and emits a complete experiment batch.
   If differential expression is enabled, edgeR rebuilds the full count matrix
   from that batch and tests each non-control group versus the control group.

The BAM sample sheet must contain `alias`, `group`, and `bam_dir`; `is_live` is
optional. It accepts case-insensitive `true` or `false`, and missing or blank
values default to `true`. A sample is live only when both `--live_analysis` and
its `is_live` value are true. Consequently, `--live_analysis=false` processes
every sample once and exits. Every final sample must have at least one BAM at
startup, while live samples may begin empty.

## Extending an experiment across workflow invocations

Completed sample-level work is persisted below the directory supplied through
`--out_dir`. After a sample's input stream closes normally (including a
sample-level `STOP`), the workflow publishes its final Oarfish counts and raw
QC inputs, then writes `FINAL` as the last file in the sample directory:

```text
<out_dir>/<group>/<alias>/
  FINAL
  quantification/final.quant
  qc/nanoplot/chunk_<batch_index>.tsv.gz
  qc/flagstat/chunk_<batch_index>.tsv
```

`FINAL` is a JSON manifest containing the sample identity, BAM paths, sizes and
modification times, reference identities, upstream container versions, and
checksums for the derived files. On a later invocation, a valid finalized
sample is restored from `--out_dir` and does not enter BAM concatenation, QC,
collation, or Oarfish. New sample-sheet rows are processed normally, after
which transcript-biotype QC and, when enabled, edgeR, fry, GSVA, and the
integrated report are rebuilt using both the restored and new quantifications.
This restoration is independent of Nextflow's `-resume` cache and remains
active when differential expression is disabled.

If any finalized BAM, reference identity, manifest field, or derived artifact
has changed, the workflow stops rather than silently reuse stale results. To
recompute a finalized sample, delete its complete
`<out_dir>/<group>/<alias>/` directory (or delete the complete output
directory) and launch the workflow again. A partial sample directory without
`FINAL` is not trusted and is recomputed.

For example:

```csv
alias,group,bam_dir,is_live
control_1,control,/data/control_1,false
control_2,control,/data/control_2,false
treated_1,treated,/data/treated_1,true
treated_2,treated,/data/treated_2,true
```

Post-startup batches initially wait for one new BAM from every live sample.
After a sample receives `STOP`, its pending BAMs are drained and it leaves all
later synchronization barriers; remaining samples may therefore continue for
different numbers of batches. At least two samples must use `control` as their
group. Each other group is tested separately against the control group.
Oarfish's EM-estimated `num_reads` values are consumed by the edgeR analysis,
which rebuilds the current count matrix for every complete live batch.

### Differential-expression stability

`--monitoring_behavior` controls automatic depth decisions and defaults
to `disabled`. `log` performs a dry run and records when a sample would stop.
`terminate` uses the sample's discovered MinKNOW `protocol_run_id` to call
`seq-run-manager stop`, then atomically creates `STOP` in that sample's
`bam_dir`. It requires `--minknow_client_certificate`,
`--minknow_client_private_key`, and `--minknow_ca_certificate`; the MinKNOW
manager defaults to `host.docker.internal:9501` and can be changed with
`--minknow_host` and `--minknow_port`.

For every sample, startup discovery inspects the direct parent of `bam_dir` for
exactly one `*sample_sheet*.csv` file and matches `alias` to its `sample_id`.
Missing, ambiguous, malformed, or unmatched metadata disables termination only
for that sample and produces a warning. Failed MinKNOW stop requests also warn,
do not create `STOP`, and retry after the next stable batch.

Every successful edgeR snapshot is compared with the previous successful
snapshot. The checks cover filtered-feature identity, median absolute logFC
change, and DE-call Jaccard similarity. Deferred DEA-readiness batches
neither increment nor reset stability. The first successful snapshot is only a
baseline. With the default `--num_stable_batches 3`, at least four successful
edgeR snapshots are therefore required.

Streaks are tracked independently for every non-control-versus-control
contrast. A treated sample becomes eligible when its own group contrast has
reached the threshold. A control sample becomes eligible only when every
contrast involving the control group has reached it. Only effectively live,
non-restored samples receive actions.

Immutable audit files are published for every analyzed snapshot:

```text
<out_dir>/stability/batch_<analysis_index>/
  config.json
  contrast_stability.tsv
  sample_stability.tsv
```

`sample_stability.tsv` includes the optional `protocol_run_id` and the result
of each attempted termination action.

They contain metric values, individual checks, consecutive streaks, required
sample contrasts, eligibility transitions, and action results. A compatible
audit and DE snapshot seed the next invocation; changing the behavior,
thresholds, or sample structure resets the streak with a warning.

Read names must be globally unique across the BAM chunks for a sample. This is
normally guaranteed by Nanopore UUID read identifiers. It ensures that all
alignment records for one read remain adjacent when independently collated
chunks are concatenated. Oarfish rereads the complete temporary cumulative BAM
on every update, but historical chunks are not decompressed and recollated.

### NCBI prokaryotic annotations

NCBI prokaryotic GTF files normally describe protein-coding genes with `gene`
and `CDS` records but do not emit the `transcript` and `exon` records required
as Oarfish projection targets. Convert such an annotation before passing it to
`--reference_annotation`:

```bash
bin/oarfish-gtf-convert genomic.gtf genomic.oarfish.gtf
```

The converter preserves existing RNA transcript/exon models and synthesizes one
single-exon transcript per protein-coding gene. It uses the complete `gene`
span, including stop-codon bases that NCBI may exclude from the corresponding
`CDS`, and retains the CDS `transcript_id`, protein ID, gene name, and locus tag.
It rejects missing or ambiguous CDS-to-transcript mappings and refuses to
overwrite an output unless `--force` is supplied.

This representation supports gene-sized bacterial quantification targets; it
does not infer operons, transcription start/termination sites, or untranslated
regions. Use an experimentally curated transcript annotation when those units
are required.

Quantification behavior and annotation-driven biotype classification are
documented separately in
[Quantification and transcript-biotype QC](quantification.md).

Differential-expression checkpoints are published exclusively under the
CLI-selected output directory:

```text
<ex_dir>/differential_expression/batch_<batch_index>/
<ex_dir>/differential_expression/latest/
```

Snapshot indices continue from the highest existing `batch_<N>` across
separate invocations. Snapshots are immutable; `latest` is refreshed with the
newest complete edgeR, fry, and GSVA result tree.

Each checkpoint contains edgeR differential-expression results, fry gene-set
tests, and GSVA sample-level scores and limma contrasts. The integrated report
consumes the edgeR, fry, and GSVA result tables. Results produced while
sequencing is active are provisional and are replaced statistically by later,
better-powered checkpoints.

## Integrated differential-analysis report

The stable `qc_report.html` artifact exposes **Quality Control**, **Differential
Analysis**, and **Gene Set Enrichment** as primary tabs once matching QC and
differential outputs exist for a batch. Within Differential Analysis, the
Overview tab contains a PCA using `log2(CPM + 1)` for all samples, an edgeR
leading-logFC MDS plot, and a biological coefficient of variation (BCV) plot.
The MDS plot uses coordinates and axis metadata exported in
`edgeR_mds_data.tsv`; the BCV plot uses `edgeR_bcv_data.tsv`, showing tagwise
BCVs as feature-level points with trended and common BCVs overlaid. Contrast
and plot-type subtabs separate:

* an edgeR `logFC` versus `logCPM` plot;
* a volcano plot of `logFC` versus `-log10(FDR)`;
* a heatmap containing up to 20 significant genes in each direction.

The adjacent **Result Stability** tab reports every sample represented in the
DE snapshot. Its second column reports the number of consecutive stable
comparisons for that sample. For a sample requiring multiple contrasts, this
is the minimum streak across those contrasts, because every required contrast
must reach the threshold before the sample becomes eligible. The remaining
columns report group, live state, required contrasts, eligibility, action
result, and analyzed batch. The report waits for the matching stability audit
before publishing a successful DE snapshot, so the table and plots always
describe the same batch.

The heatmap selects genes by effect size after significance filtering, shows
only the two compared conditions, and displays row z-scores of `log2(CPM + 1)`.
Rows are ordered by average-linkage hierarchical clustering of those displayed
profiles; the report does not draw a dendrogram. Blue represents low and red
represents high row-standardized expression.
Condition colors are shared by the PCA, significant-gene directions, and the
heatmap sample annotation. Heatmap cell colors encode expression z-scores
rather than condition.

`--de_lfc_cutoff` defaults to `1.0` and is passed to edgeR `glmTreat` as well as
the report. `--de_padj_cutoff` defaults to `0.05` and is applied to the edgeR
`FDR` column. A feature is colored and eligible for the heatmap only when it
passes both thresholds.

These plots report differential-expression associations. They do not establish
that a reported gene causally drives the experimental phenotype. Results from
live sequencing batches remain provisional until the final batch is analyzed.

The Gene Set Enrichment primary tab separates **GSVA scores**, **GSVA
differential**, and **fry enrichment** analyses. The fry contrast and plot-type
subtabs include a signed summary displaying up to 30 gene sets with the smallest
directional fry FDR values. Positive bars represent coordinated expression
toward the target condition and negative bars represent the reference
condition. Significant bars use the corresponding condition color; gene sets
that do not pass `--de_padj_cutoff` are gray. Mixed P-values and FDRs remain
available in the hover metadata but are not assigned an up/down sign.

Every tested gene set is available in the same Bootstrap dropdown-tab control
used for read-flow samples. The selector and signed FDR plot use the gene-set
identifier from the first GMT column; the description and fry statistics are
available in the collapsible **Gene-set details** block. Genes are ordered from
negative to positive edgeR `logFC`; colored ticks locate the selected set's
retained members, and the neutral worm shows their relative local density using
limma's tricube moving average with a span of 0.45.

## edgeR `fry` gene-set analysis

The `edgeR-analysis` command runs differential expression and rotation-based
gene-set testing from feature-count tables. It accepts a tab-delimited manifest
with `name`, `group`, and `count_file` columns, an optional GTF or GFF3
annotation, and any valid GMT collection:

```bash
edgeR-analysis \
    --quant_manifest quant_manifest.tsv \
    --annotation genomic.oarfish.gtf \
    --gene_sets pathways.gmt \
    --output_dir edgeR_results
```

Count files may be tab- or comma-delimited. Common feature columns (`tname`,
`target_id`, `transcript_id`, `gene_id`, `feature_id`, `name`, or `id`) and
count columns (`num_reads`, `expected_count`, `est_counts`, `count`, `counts`,
or `read_count`) are detected case-insensitively. For other layouts, provide
`--quant_id_column` and `--quant_count_column`. Relative `count_file` paths are
resolved relative to the manifest rather than the current working directory.

GMT members can match count-table feature IDs directly, in which case no
annotation is needed. Otherwise, the command builds an identifier map from
identifier-like GTF/GFF3 attributes, including gene, transcript, protein,
locus-tag, `ID`, `Parent`, `Name`, and `Alias` fields. Use
`--annotation_id_attributes key1,key2` for custom attributes and
`--strip_id_versions` when one input includes terminal numeric versions such as
`.1` and the other does not. Namespaced GFF3 identifiers such as
`gene:GENE001` and `transcript:TX001` are handled automatically.

The script uses `control` as the reference group unless
`--control_group <label>` is supplied. It filters low-expression transcripts,
applies TMM normalization, estimates robust dispersions, fits an edgeR
quasi-likelihood model, and tests each non-control group against the reference.
Gene-level differential-expression results use `glmTreat`; its minimum effect
size defaults to one log2 fold and can be changed with `--lfc`.

At the analysis level, the command writes `feature_counts.tsv`,
`sample_metadata.tsv`, `edgeR_mds_data.tsv`, and `edgeR_bcv_data.tsv`. The MDS
table contains the sample coordinates and edgeR axis metadata; the BCV table
contains average log CPM plus tagwise, trended, and common dispersion and BCV
values for every retained feature. For each contrast the command writes
`edgeR_results.tsv`, `fry_results.tsv`, and `fry_signed_significance.png`. The
plot preserves arbitrary GMT set names and shows the top 30 sets by directional
FDR by default; change this with `--plot_top_n`. The first column of
`edgeR_results.tsv` is the generic count-table `feature_id`; annotation
identifiers are added as columns when an annotation is supplied.

At analysis level, `gene_set_resolution.tsv` records every GMT member, its
resolved count-table feature ID, and the matching route. The
`gene_set_coverage.tsv` summary distinguishes matched GMT members from unique
count-matrix features and records how many survive expression filtering. A set
requires at least two retained features to be tested. These files should be
checked before interpreting a null enrichment result.

In `fry_results.tsv`, `Direction`, `PValue`, and `FDR` describe the directional
test for coordinated up- or downregulation. `PValue.Mixed` and `FDR.Mixed`
describe the non-directional test for any coordinated differential expression.
The signed plot displays `-log10(FDR)`, positive for upregulated sets and
negative for downregulated sets.

The static PNG remains an edgeR command output. The integrated HTML report
recreates the signed view from `fry_results.tsv` and uses
`gene_set_resolution.tsv` to construct interactive, single-set barcode plots.
The report validates that each set's unique resolved features retained in
`edgeR_results.tsv` agree with fry's `NGenes` value.

Gene-set enrichment is a self-contained expression association in the tested
organism, conditions, and sampling context. It does not by itself demonstrate
pathway activity, biological mechanism, or causality. Live-analysis enrichment
results are provisional and may change as additional reads alter expression
estimates and statistical power.

## GSVA sample-level gene-set analysis

After edgeR completes, `gsva-analysis` reuses `feature_counts.tsv`,
`sample_metadata.tsv`, and `gene_set_resolution.tsv`. It transforms the
edgeR-filtered TMM-normalized CPM matrix with `log2(CPM + 1)`, removes features
that do not vary across the current samples, and runs the standard GSVA method
with a Gaussian kernel. Gene sets require at least two variable retained
features. The standalone command is:

```bash
gsva-analysis \
    --feature_counts edgeR_results/feature_counts.tsv \
    --sample_metadata edgeR_results/sample_metadata.tsv \
    --gene_set_resolution edgeR_results/gene_set_resolution.tsv \
    --output_dir edgeR_results
```

The batch root receives:

* `gsva_scores.tsv`, a gene-set-by-sample score matrix;
* `gsva_scores_long.tsv`, a Python-friendly tidy score table with sample groups;
* `gsva_gene_set_coverage.tsv`, including variable-feature filtering status;
* `gsva_parameters.tsv`, recording the fixed scoring configuration;
* `gsva_score_heatmap.png`, with row-standardized colors for display only;
* `gsva_group_boxplots.pdf`, with raw sample scores by group.

Each existing `group_<target>_vs_<control>/` directory also receives
`gsva_limma_results.tsv`. Its `effect_size` is the target-minus-control
difference in GSVA scores, and its `adjusted_p_value` is Benjamini-Hochberg
adjusted within that contrast. The score tables retain raw GSVA scores; only
the heatmap is row-standardized.

GSVA scores summarize relative expression patterns within this dataset. They
are not direct biochemical measurements of pathway activity and do not
establish mechanism or causality. The HTML report displays a variance-ranked,
row-standardized score heatmap, raw score distributions selected through a
gene-set dropdown, and the coverage table. GSVA heatmap rows retain their
selection order, while sample columns retain the order defined by the metadata.
Blue represents low and red represents high scores.
Plots and the displayed coverage table use gene-set identifiers rather than
descriptions. Its limma views include a
multi-contrast effect/significance dot plot plus per-contrast volcano and
significant-score heatmap subtabs. Volcano x-axes are labeled as target-minus-
control **GSVA score differences**, not log fold changes; BH-adjusted p-values
control the report significance threshold. The across-contrast score-difference
scale follows the same convention: negative/low values are blue and
positive/high values are red.

## Temporal gene-set analysis

Enable `--timeline_analysis` to add the primary **Temporal Analysis** report
tab. This mode requires differential expression, gene-set enrichment, and a
signed integer `order` value for every sample-sheet row. `order` is elapsed
time in minutes. The first version represents one trajectory: samples in the
same group are independent biological replicates at one minute, every group
has one minute, and every minute has one group.

The gene-set dropdown controls two stacked figures. The upper figure connects
mean raw GSVA scores and shows one-sample-SD whiskers plus the individual
sample scores. The lower figure is a gene-by-time heatmap of gene-wise z-scores
calculated from mean `log2(TMM-normalized CPM + 1)` at each time point. Genes
with similar standardized trajectories are adjacent after average-linkage
clustering. Hover fields retain the absolute mean logCPM, SD, group, and
replicate count; SD is unavailable where a time point has one replicate.

These figures are descriptive and do not perform temporal regression,
smoothing, repeated-measures modeling, or a treatment-by-time test. Connecting
lines do not estimate unmeasured intermediate states. GSVA is a dataset-relative
expression summary rather than a biochemical pathway-activity assay, and bulk
temporal patterns can reflect composition, batch, or another variable
confounded with time.

The workflow uses `rnabioinfo/seq_lm_gsva:v1.1.0`. Build and publish that image
before deploying this workflow version:

```bash
docker_containers/helper_scripts/build_and_publish_docker_image.sh \
    seq_lm_gsva rnabioinfo
```

Organism- and strain-specific carbon-stress GMT files are documented in
[`data/gene_sets/README.md`](../data/gene_sets/README.md).
