# Fixed-matrix iModulon analysis

Supply `--ica_matrix` to infer component activities from each complete cumulative
Oarfish batch. Every non-control group is also compared with the control group
using independent two-sided Welch tests and Benjamini–Hochberg (BH) correction.
ICA runs alongside differential expression and does not require edgeR, GMT,
fry, or GSVA. Active ICA results are included in a dedicated report tab.

The `ica_analysis` toggle defaults to `true`, preserving matrix-based activation.
Set `ica_analysis` to `false` in your parameter file (or `params.ica_analysis = false`
in Nextflow config) to disable ICA even when a matrix is configured. Disabled ICA
skips matrix/map file loading and snapshot publication. Without `ica_matrix`, ICA
is skipped regardless of the toggle. Other analysis branches keep their own settings.

ICA has no separate time-course toggle. When active ICA sees an `order` column in
the sample sheet, it automatically adds elapsed-minute activity trajectories.
Every row must then provide a signed integer `order`; each group must identify one
time, each time must identify one group, and at least two times are required.
Validation happens before analysis tasks start. Samples at a time are independent
biological replicates. Unequal intervals and negative times are supported.
Inference remains the existing Welch comparison of each non-control time point
against the shared controls; no trend or repeated-measures model is fitted.
Without an `order` column, ICA uses the comparison layout.

## Run

Build the local analysis image once from the `seq_lm_wf` directory. The image is
not automatically downloaded or published by this change:

```bash
docker build -t rnabioinfo/seq_lm_ica:v1.0.0 ../docker_containers/seq_lm_ica
```

The image follows the other `docker_containers` examples: it provides pinned
runtime dependencies, while Nextflow supplies the CLI from the pipeline's `bin`
directory. Its `docker_version.json` matches the workflow's `v1.0.0` image tag.

An ICA-only analysis (plus existing QC/quantification) can use an `ica.json`
parameter file. JSON preserves boolean and numeric types for schema validation:

```json
{
  "sample_sheet": "samples.csv",
  "reference_genome": "genome.fasta",
  "reference_annotation": "transcripts.gtf",
  "ica_analysis": true,
  "ica_matrix": "weights.tsv",
  "differential_expression": false,
  "gene_set_enrichment": false,
  "live_analysis": false,
  "out_dir": "output"
}
```

```bash
nextflow run . -params-file ica.json
```

Both reference inputs are required. Existing sample-sheet requirements apply:
`alias`, `group`, and `bam_dir`, with at least two control samples. Aliases may
repeat across groups. ICA recognizes one case-insensitive spelling of `control`;
using both `control` and `Control` as distinct labels is an error. Use lowercase
`control` for compatibility with the existing differential-expression branch.
Adding the optional `order` column enables ICA time-course analysis automatically;
omit the column for comparison-only ICA. Sample rows represent independent
biological replicates, not paired observations.

## Matrix and mapping

A comma- or tab-delimited matrix must have genes in rows, a header row, and one
column per component after the first identifier column:

```text
model_gene_id	carbon_program	stress_program
gene_A	0.35	-0.08
gene_B	-0.20	0.42
```

The first header may be blank, as in common exported `M.csv` or `S.csv` files.
Gene IDs and component headers must be nonempty and unique. Weights must be
finite numbers. Orientation is never guessed. All components are fitted jointly;
their names, order, signs, and scales are preserved. The original training activity
matrix, PyModulon objects, module-membership thresholds, and metadata sidecars
are not required. Component column names serve as display identifiers.

The run annotation resolves Oarfish `tname` targets to genes. Exact gene IDs and
unambiguous locus-tag aliases are matched to model rows. GTF transcript IDs and
GFF3 transcript-parent relationships are supported. Multiple transcripts of one
gene are summed before logarithmic transformation.

For another identifier namespace or a validated cross-strain mapping, provide
`--ica_gene_map mapping.tsv`:

```text
gene_id	model_gene_id
run_gene_1	gene_A
run_gene_2	gene_B
```

`gene_id` is the canonical annotation gene ID (`gene_id` in GTF, or the gene's
`ID` in GFF3 when no `gene_id` is supplied). Mappings must be one-to-one and refer
to known annotation genes with transcript targets and known model genes. Explicit
maps are authoritative: automatic alias matching is not added. Ambiguous
parentage, aliases, duplicate destinations, and duplicate source genes are errors.
No orthologs or model downloads are generated automatically.

## Parameters and numerical contract

| Parameter | Default | Meaning |
| --- | --- | --- |
| `ica_analysis` | true | Allow ICA when a matrix is supplied; false disables the branch |
| `ica_matrix` | unset | Weights CSV/TSV; enables analysis when the toggle is true |
| `ica_gene_map` | unset | Explicit one-to-one gene map; requires the matrix |
| `ica_log_base` | 2.0 | Finite logarithm base greater than 1 |
| `ica_pseudocount` | 1.0 | Positive finite value added to million-scaled abundance |
| `ica_min_gene_coverage` | 1.0 | Required model-gene fraction, greater than 0 and at most 1 |
| `ica_min_read_count` | 10000 | Nonnegative assigned-abundance threshold for every sample |
| `ica_padj_cutoff` | 0.05 | BH significance threshold, greater than 0 and at most 1 |

All parameters are declared in the Nextflow config and schema. They are independent
of DEA read-count and significance parameters. Full gene coverage is required by
default. To permit a partial projection deliberately, lower
`--ica_min_gene_coverage`, for example to `0.9`, after assessing the mapping.
Every component remains in the fit and the retained shared matrix must have full
column rank. Missing genes are excluded, never represented as zero expression.
Expected mapped transcript rows must exist in every quantification, including
rows with a measured zero.

For each sample, divide aggregated gene abundance by **all Oarfish-assigned
abundance**, including targets outside the model, and multiply by one million.
There is no length correction, biotype exclusion, TMM normalization, or
expression-based gene filtering. Then:

```text
L = log_base(million_scaled_gene_abundance + pseudocount)
reference[g] = mean(L[g, control_samples])
X = L - reference
A = pseudoinverse(M_shared) @ X_shared
residual = X_shared - M_shared @ A
```

The annotation mapping, shared genes, and pseudoinverse are prepared once per
invocation. Numerical rank uses the SVD cutoff
`max(M_shared.shape) * machine_epsilon * largest_singular_value`. The same
pseudoinverse is applied to all samples. Live updates recompute the reference and
all activities, including activities for stopped or restored samples.

The supplied weights must be compatible with this abundance/transform convention.
Changing log base or pseudocount does not make a model trained on standardized
expression or another measurement protocol compatible. Compatibility cannot be
established from the matrix alone. Coverage, retained squared-weight fractions,
and residuals are diagnostic measurements, not probabilities of biological validity.
Component activities are model-scale estimates relative to this experiment's
controls; annotations do not establish causal regulation or pathway flux.

## Differential activity

Each component is tested for target-minus-control mean activity difference,
allowing unequal group variances. Tests require two or more biological samples
in both groups. Individual samples remain available when a group is too small
for a test. SD is unavailable for a single sample. Tests with zero variance in
both groups are unavailable; one zero-variance group is allowed if the other
supplies estimable variance.

Tests use uncentered projections of logged expression. Subtracting the common
control baseline leaves the group difference and within-group variances
unchanged. Published activities and descriptive group means are control-centered.
The effect is an **activity difference**, not a log fold change. Confidence
intervals are nominal two-sided 95% intervals, not simultaneous intervals.

BH correction is performed separately for each contrast and snapshot, with every
supplied component in the family. Untestable components enter the correction as
p=1 internally but retain unavailable p-values, adjusted p-values, and significance
in the output. A valid test is significant when adjusted p-value is at most
`ica_padj_cutoff`. No biological effect-size threshold is imposed across arbitrarily
scaled components.

Live batches are repeated estimates from accumulating reads, not independent
replicates. Snapshot p-values do not control error from repeated checking or
optional stopping, and ICA does not trigger sequencing termination.

## Outputs and live status

Snapshots are immutable directories under `output/ica/batch_<analysis_index>/`.
Indices continue from existing ICA snapshots independently of DEA indices.
`ica/latest.json` is updated atomically after the newest complete snapshot is
published, including deferred snapshots. Read its `path` relative to `ica/` and
check `status` before consuming activity tables. Do not run concurrent workflow
invocations against the same output directory.

Each ready snapshot contains:

| File | Schema/content |
| --- | --- |
| `sample_metadata.tsv` | `sample_id`, `alias`, `group`, `order`, `count_file`, `source_batch_index`, `assigned_abundance`, `ready` |
| `activities.tsv` | `component_id` followed by stable sample-ID columns |
| `activities_long.tsv` | `component_id`, `sample_id`, `activity`, `alias`, `group`, `order` |
| `activity_summary.tsv` | `component_id`, `group`, `mean`, `sd`, `n` |
| `differential_activity.tsv` | Component/contrast statistics described below |
| `gene_mapping.tsv` | `model_gene_id`, `gene_id`, `transcript_id`, `method`, `status` |
| `component_coverage.tsv` | `component_id`, `gene_coverage`, `retained_squared_weight_fraction` |
| `projection_qc.tsv` | `sample_id`, residual sum of squares, centered sum of squares, residual RMSE, normalized residual |
| `centered_expression.tsv` | `model_gene_id` followed by stable sample-ID columns |
| `reference_expression.tsv` | `model_gene_id`, `reference_expression`, JSON array of `control_sample_ids` |
| `provenance.json` | Schema version 1; hashes, model mapping/diagnostics, settings, control identities, software versions, local batch/report sequence/global analysis indices |
| `status.json` | `ready` or `deferred`, statistical availability and test-status counts or deferred sample IDs/reason |

`differential_activity.tsv` columns are `component_id`, `target_group`,
`control_group`, `activity_difference`, `target_mean`, `control_mean`, `target_sd`,
`control_sd`, `target_n`, `control_n`, `standard_error`, `degrees_of_freedom`,
`t_statistic`, `ci_lower`, `ci_upper`, `p_value`, `adjusted_p_value`, `significant`,
and `status`. Test statuses are `tested`, `insufficient_replicates`, or
`zero_variance`. An experiment containing controls only has an empty table with
these headers. Missing numeric/statistical values are empty TSV fields; JSON
uses null rather than NaN.

Sample IDs are SHA-256 identifiers of the JSON-encoded full group/alias pair.
Labels remain separate and no component-derived filesystem names are created.
Normalized residual means squared residual norm divided by squared centered
expression norm; it is unavailable if that denominator is zero.

A batch is deferred if any sample has zero assigned abundance or is below
`ica_min_read_count`. Its directory contains only sample metadata, provenance,
and status. It never reuses earlier activity tables. Malformed inputs, invalid
mappings, insufficient model coverage, and deficient rank stop analysis. Mapping
and rank diagnostics prepared before a failure are available in the Nextflow
model-preparation work directory. Final low-depth samples remain explicitly
deferred; lower the threshold and rerun if appropriate.

## Report visualization

When ICA is active, the matching immutable ICA snapshot is required before its
QC report batch is generated. The main **iModulon Analysis** tab includes an
overview, control-centered activity heatmaps, per-component sample and effect
views, differential-activity plots and tables, diagnostics, and interpretation
notes. One contrast is shown directly; multiple contrasts use contrast tabs and
an across-contrast effect matrix. Deferred snapshots show assigned abundance and
readiness without stale plots.

Report joins reject duplicate batches, unmatched inputs, inconsistent snapshot
identities, and missing publication sequences. A missing batch stops the run
instead of silently dropping reports or leaving later snapshots buffered.
The loader also rejects incomplete component/sample or component/contrast tables,
conflicting time metadata, malformed numbers, and inconsistent statistical
availability. Empty statistics remain valid only where inference is unavailable.

Live reports use a persistent `report_revision`, independently of local batch
numbers and ICA analysis indices. Each HTML snapshot is named
`qc_report_snapshot_revision_<revision>.html`; `qc_report_state.json` identifies
the revision and source batch. Revision allocation continues across restarts and
skips snapshots left by interrupted publication. Older reports cannot replace
newer reports, and restarting the batch counter does not prevent an open viewer
from refreshing. The complete HTML is published before the state pointer advances.
The final report remains self-contained. An already-open viewer created by an
older pipeline version needs one manual reload after upgrading.

Volcano plots retain the published adjusted p-values in tooltips. Only exact
zeros are replaced with the smallest positive representable value for plotting;
positive p-values are unchanged. The dashed horizontal line marks the configured
adjusted p-value cutoff, and the vertical dotted line marks zero activity difference.

When the sample sheet contains `order`, activities and target-minus-control
differences are also shown against the numeric elapsed-minute scale. Connected
lines join group means only and do not imply paired samples. Activities are not
z-scored, smoothed, clustered, or interpolated; signs and scales remain those of
the supplied model.

Changes to the model, mapping, or transform recompute downstream ICA from restored
quantifications without changing upstream checkpoints. New hashes and settings
are stored with each snapshot.

## Standalone CLI and verification

The numerical CLI can also consume a quantification manifest directly:

```bash
imodulon-analysis prepare --matrix weights.tsv --annotation transcripts.gtf \
    --output prepared_ica
imodulon-analysis analyze --prepared prepared_ica --manifest quant_manifest.tsv \
    --counts-dir quant --output ica_results
```

The manifest requires `name`, `group`, and `count_file`, and optionally `order`
and `source_batch_index`. Outputs must name new directories. `prepare --gene-map`
and `--min-gene-coverage`, and `analyze --log-base`, `--pseudocount`, `--min-reads`,
`--cutoff`, `--batch-index`, `--analysis-index`, and `--report-sequence` expose the
same analysis settings and provenance.

`containers/imodulon/smoke_test.sh` builds the image and exercises both CLI stages.
The Python tests are under `bin/imodulon_analysis/tests`. The Nextflow integration
fixture is `test/test_imodulon_analysis.nf`; its local config uses the Python
environment on PATH. To run the same pinned numerical environment locally:

```bash
python3 -m venv /tmp/ica-tests
/tmp/ica-tests/bin/pip install -r containers/imodulon/requirements.txt pytest==8.3.5
/tmp/ica-tests/bin/python -m pytest bin/imodulon_analysis/tests -q
PATH="/tmp/ica-tests/bin:$PATH" bash test/run_imodulon_test.sh
bash test/run_report_integrity_test.sh
nextflow lint .
```

Report regressions are in `bin/workflow_glue/tests/test_imodulon_plots.py` and
`test_qc_report.py`, using the report Python environment. They cover malformed
snapshots, unavailable tests, actual numerical CLI outputs, and revision-based
refresh. The JavaScript refresh regression uses Node with mocked DOM/network
events. The Nextflow report integrity script exercises rejected joins and stale
publication as expected failures, plus restart revision allocation and final HTML.

The integration runner creates synthetic fixtures in a temporary directory and
runs twice, deliberately delaying the first batch to test ordering and checking
restart index allocation. No biological reference download is required.
