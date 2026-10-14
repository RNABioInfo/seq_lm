# Fixed-matrix iModulon analysis

Enable `--ica_analysis=true` and supply `--ica_matrix` to infer component activities
from each complete cumulative Oarfish batch. seqLM performs fixed-model activity
inference: it does not train ICA, discover new components, or determine iModulon
gene-membership thresholds. Each selected non-control group is compared with
`ica_control_group` using a two-sided Welch test for each component, followed by
Benjamini–Hochberg (BH) correction within that contrast and snapshot. Biological
samples are assumed independent; tests can be dependent, including contrasts
that share controls. This branch can run with or without differential expression
and does not require edgeR, GMT, fry, or GSVA. Active results are included in a
dedicated report tab.

The repository CLI `imodulon-analysis` implements this analysis using NumPy,
SciPy, and pandas. PyModulon is not a runtime dependency and neither
`pymodulon.util.infer_activities` nor `pymodulon.util.dima` is called.

Enable ICA explicitly in your run parameters.
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

Nextflow automatically fetches the published `rnabioinfo/seq_lm_ica:v1.0.0`
image from Docker Hub through the selected container runtime. The image supplies
pinned runtime dependencies, while Nextflow supplies the CLI from the workflow's
`bin` directory.

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
`alias`, `group`, and `bam_dir`, with at least two samples in
`--ica_control_group` (default `control`). Group names match exactly after
whitespace is trimmed. `--ica_condition_groups early,late` selects those
conditions plus all samples in the configured ICA control group. Unset or
whitespace-only values select all samples. Empty entries, duplicate names, and
unknown groups are errors. The gene contrast CSV does not affect ICA.

The selected cohort determines ICA readiness, centering, activity comparisons,
time-course plots, and monitoring. Excluded groups remain in QC; they can also
appear in DEA if selected by gene contrasts. ICA provenance records the control
and selected groups, and changing either resets monitoring history. Unselected
live samples still block stopping a shared MinKNOW run.
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
are not required. Component column names remain the stable identifiers.

To show curated names and descriptions, supply `--ica_imodulon_table iM_table.csv`
from the **same model** as the weight matrix. PyModulon exports commonly call this
`imodulon_table.csv`; iModulonDB uses `iM_table.csv`. This is model-specific
annotation, not a universal mapping from component numbers to biological functions.
The CSV/TSV must contain one row for every matrix component, keyed by `component_id`,
`k`, or its first exported index column, and a `name` column. Optional columns are
`function` (or `description`), `regulator` (or `regulator_readable`), and `category`.
Unknown, duplicate, and missing component IDs are rejected. The normalized table
and its source SHA-256 are recorded in model provenance and carried into snapshots.
Names appear alongside IDs in plot axes, titles, and component selectors; functions
appear in heatmap/volcano hovers and component details, with annotations in report
tables. Unannotated models retain their original IDs without inferred names.

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
| `ica_analysis` | false | Enable ICA when a matrix is supplied; false disables the branch |
| `ica_control_group` | `control` | ICA centering and activity-comparison reference; at least two samples |
| `ica_condition_groups` | unset | Comma-separated conditions; includes ICA control automatically; unset uses all samples |
| `ica_matrix` | unset | Weights CSV/TSV; enables analysis when the toggle is true |
| `ica_imodulon_table` | unset | Matching component names/functions CSV/TSV; requires the matrix |
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
expression-based gene filtering. This is million-scaled assigned abundance;
it should not be described as length-normalized TPM. Then:

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

More explicitly, let `P = pseudoinverse(M_shared)` and `B = P @ L_shared`.
For component `j`, target samples `T`, and control samples `C`, the implementation
calculates:

```text
activity_difference[j] = mean(B[j, T]) - mean(B[j, C])
                      = mean(A[j, T]) - mean(A[j, C])
                      = (P @ (mean(L_shared[:, T], axis=1)
                            - mean(L_shared[:, C], axis=1)))[j]
```

Samples contribute equally to each mean, without read-depth weighting. The
published control mean is zero up to floating-point precision. Positive
differences mean higher activity in the supplied component's sign convention;
they do not necessarily mean activation of every member gene.

The test uses sample variances `vT` and `vC` with `ddof=1` (exactly constant
groups are explicitly assigned zero variance), and sample counts `nT` and `nC`:

```text
SE = sqrt(vT/nT + vC/nC)
df = (vT/nT + vC/nC)^2 / ((vT/nT)^2/(nT-1) + (vC/nC)^2/(nC-1))
t = activity_difference / SE
p = 2 * scipy.stats.t.sf(abs(t), df)
CI95 = activity_difference +/- scipy.stats.t.ppf(0.975, df) * SE
```

No paired, covariate-adjusted, or repeated-measures model is fitted. Two samples
per group are the implementation's minimum for an estimable test, not a validated
sample-size recommendation. Very small groups provide unstable variance
estimates; the test does not account for uncertainty in the supplied ICA matrix.

BH correction is performed separately for each contrast and snapshot, with every
supplied component in the family. Untestable components enter the correction as
p=1 internally but retain unavailable p-values, adjusted p-values, and significance
in the output. A valid test is significant when adjusted p-value is at most
`ica_padj_cutoff`. No biological effect-size threshold is imposed across arbitrarily
scaled components.

BH adjustment is the implemented procedure, not evidence of unconditional FDR
control: its guarantees depend on valid input p-values and suitable dependence
conditions. Independence of ICA components in a training compendium does not
establish independence of the contrast tests in this new experiment.

Live batches are repeated estimates from accumulating reads, not independent
replicates. Snapshot p-values do not control error from repeated checking or
optional stopping. The convergence metric below does not use these p-values to
decide termination.

## Published basis and distinction from PyModulon DiMA

The projection follows the least-squares calculation published by
[Sastry et al. (2019), Nature Communications](https://www.nature.com/articles/s41467-019-13483-w)
and the activity-difference projection described by
[Sastry et al. (2021), PLOS Computational Biology](https://pmc.ncbi.nlm.nih.gov/articles/PMC7888660/).
It also uses the same mathematical operation as
[PyModulon's `infer_activities`](https://pymodulon.readthedocs.io/en/stable/_modules/pymodulon/util.html#infer_activities).
Numerical agreement requires identical gene selection, expression preparation,
and compatible pseudoinverse tolerances; seqLM rejects deficient rank explicitly.
These references establish the projection method, not validation of seqLM's
specific Oarfish normalization or its complete analysis pipeline.

seqLM combines that projection with
[Welch's unequal-variance test (1947)](https://doi.org/10.1093/biomet/34.1-2.28)
and [Benjamini–Hochberg correction (1995)](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x).
This is a separate statistical choice from the DiMA method described in the
2019 paper. That method fits a component-specific log-normal distribution to
absolute within-condition biological-replicate activity differences, then
evaluates the absolute difference between condition means against that
distribution and applies FDR and an activity-difference threshold. The
[documented PyModulon 0.2.1 `dima` source](https://pymodulon.readthedocs.io/en/stable/_modules/pymodulon/util.html#dima)
uses the supplied activity matrix and project/condition metadata; it does not
automatically obtain training-compendium noise from a weight matrix.

Welch testing estimates noise from the two groups being compared and needs no
replicate compendium. This is a practical rationale for the current interface,
not evidence that Welch is superior to DiMA. The repository does not establish
comparative calibration or superiority on representative biological experiments.
Existing numerical tests compare the implementation with SciPy's Welch test;
this verifies the calculation, not biological calibration.

DiMA has no documented universal minimum sample count. A contrast may contain
one sample per condition if suitable replicate-noise data are supplied elsewhere,
but a few replicate differences cannot justify a reliable fitted noise
distribution. In the documented 0.2.1 source, pairs overwrite the same
project/condition column, retaining only the last pair per group; this is a
version-specific implementation detail to verify before adopting that package.
Neither a universal DiMA sample minimum nor improved significance calibration
should be inferred from package usage. Historical noise must be assessed for
compatibility with this experiment's protocol and sequencing depth.

## Implementation packages and tools

| Tool/package | Role |
| --- | --- |
| Nextflow (host environment: >=26.04.3), OpenJDK 17 | Workflow orchestration |
| Docker or Singularity/Apptainer, according to the run profile | Process containers |
| SAMtools | Remove `ts` tags, collate reads, assemble cumulative BAMs |
| Oarfish | Transcript `num_reads` estimates from genome-aligned reads and annotation |
| `imodulon-analysis` 1.0.0, Python 3.12.8 | Repository numerical CLI and ICA container runtime |
| NumPy 1.26.4 | SVD, pseudoinverse, projection, variances, local BH implementation |
| SciPy 1.13.1 (`scipy.stats`) | Student-t probabilities and confidence intervals |
| pandas 2.2.3 | Input/output tables, sample joins, group summaries |
| Bokeh, ezcharts, colorcet (report container) | Interactive plots, report layout/tables, colors |

Numerical dependencies are pinned in `../docker_containers/seq_lm_imodulon/requirements.txt`;
the standalone CLI runtime is defined in `../docker_containers/seq_lm_imodulon/Dockerfile`. Standard-library
helpers handle argument parsing, annotation/CSV/JSON parsing, paths, and hashes.
The ICA container is `rnabioinfo/seq_lm_ica:v1.0.0`. Container tags identify images;
they are not the SAMtools or Oarfish executable version numbers. Actual numerical
package versions are recorded in each snapshot's `provenance.json`.

## Sequencing termination

ICA can control sequencing termination independently of differential expression:

```groovy
ica_analysis = true
ica_matrix = '/path/to/model.csv'
stopping_method = 'ica'
monitoring_behavior = 'log'  // use 'terminate' for MinKNOW termination
ica_max_fold_change = 1.05
ica_min_read_growth_fraction = 0.1
num_stable_batches = 3
```

Supply the reference genome and annotation as usual. For ICA-only analysis,
disable both `differential_expression` and `gene_set_enrichment`. Termination
requires live analysis and the same MinKNOW credentials as DEA monitoring.
Both analyses may run, but only `stopping_method` controls termination; ICA
monitoring does not require significant differential activity.

For each sample and **every component in the supplied model**, the evaluator
retains `num_stable_batches + 1` informative cumulative checkpoints. A checkpoint
is informative when its assigned abundance has increased by at least
`ica_min_read_growth_fraction` since the last accepted checkpoint. The default
therefore requires three depth increases of at least 10%, spanning four snapshots.
The first ready snapshot establishes a baseline and cannot trigger stopping.

For component j, let R be its highest activity minus its lowest activity in the
window, and W the largest absolute gene weight in its full model column. The
reported multiplicative change is `ica_log_base ** (R * W)`. Every component must
be at or below `ica_max_fold_change`. Computation uses logarithms to assess the
threshold safely when exponentiation would overflow. The metric is unchanged by
rescaling or reversing the sign of an ICA column together with its activities.
It works for zero and negative activities without dividing by those activities.

For example, an activity range of 0.10 and a maximum absolute weight of 0.50
give `2 ** 0.05 = 1.035`, which passes a 1.05 tolerance. Decreases and increases
are treated symmetrically as multiplicative factors of at least one. This is
each component's contribution to reconstructed **logged, normalized expression
with the pseudocount**, not an observed gene fold change, total reconstruction
error, biological activity fold change, or confidence bound. Validate the
tolerance, depth floor, and window on completed runs before relying on automatic
termination. Minimum-read readiness and model mapping/rank checks still apply.

Unchanged or insufficiently increased depth never advances the window, but its
current activities are still checked: a moving control baseline can revoke
eligibility. Decreasing depth restarts that sample's window. Deferred or invalid
snapshots reset the windows; missing/nonfinite activities and missing model
scales block stopping. Changes to model, sample/control identities, transformation,
software, or monitoring settings also reset history. Each new workflow invocation,
including resume, starts a fresh monitoring window rather than restoring an old
stop decision.

Audits are published atomically under `stability/ica/batch_<ICA analysis index>/`:
`component_stability.tsv` records activity ranges, weight scales, multiplicative
changes, and pass/fail results; `sample_stability.tsv` records passing/blocking
components, window readiness, run blockers, and stop actions; `config.json`
records the monitoring settings. The matching HTML report includes an **ICA
Monitoring** tab even when DEA is disabled. Every snapshot must include valid model
weight scales.

A live sample is eligible only if all its components pass the full window.
All required active live samples sharing a MinKNOW acquisition must be eligible
before that acquisition stops. Stop acknowledgement and STOP-marker handling
use the same coordinator as DEA. Remaining BAMs are drained and final results
reported using the existing ingress/finalization path. In `log` mode the
coordinator records when it would stop, without issuing a stop request.

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
QC report batch is generated. The main **iModulon Analysis** tab includes
control-centered activity heatmaps, per-component sample and effect
views, differential-activity plots and tables, and diagnostics. One contrast is shown directly; multiple contrasts use contrast tabs and
an across-contrast effect matrix. The mean-activity heatmap uses ordered discrete time points with equally spaced,
adjacent columns; trajectory plots retain elapsed-minute spacing. Overview and
Methods sub-tabs are omitted. Complete `order` metadata selects **Time course details**;
without it, **Component details** shows sample distributions and group effects.
Only one of these detail tabs is displayed. **Diagnostics** contains only the
Component annotations table. Snapshot readiness, projection QC, mapping, and rank
statistics remain available in the analysis output files. Deferred snapshots show
a brief status message and component annotations without stale plots.

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
The final report remains self-contained.

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

The manifest requires `name`, `group`, `count_file`, `order`, and
`source_batch_index` columns. Leave `order` values blank when time-course analysis
is not needed. Outputs must name new directories. `prepare --gene-map`
and `--min-gene-coverage`, and `analyze --log-base`, `--pseudocount`, `--min-reads`,
`--cutoff`, `--batch-index`, `--analysis-index`, and `--report-sequence` expose the
same analysis settings and provenance.

The Python tests are under `bin/imodulon_analysis/tests`. The Nextflow integration
fixture is `test/test_imodulon_analysis.nf`; its local config uses the Python
environment on PATH. To run the same pinned numerical environment locally:

```bash
python3 -m venv /tmp/ica-tests
/tmp/ica-tests/bin/pip install -r ../docker_containers/seq_lm_imodulon/requirements.txt pytest==8.3.5
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
