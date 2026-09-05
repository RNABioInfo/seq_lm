# Existing ICA model integration assessment

Assessment date: 2026-09-04. Pipeline inspected at commit `66cd122`, including
the working tree. This is a design assessment; the interfaces below are
proposed and are not implemented.

The recommended first version is an optional, organism-independent branch that
projects complete quantification batches onto a supplied, fixed ICA basis and
adds an **ICA / iModulons** tab to the existing `qc_report.html`. Its outputs are
sample activities, descriptive group/time summaries, and projection diagnostics.
Training ICA models, discovering new modules, differential-activity significance
testing, and ICA-driven sequencing termination are outside this version.

## Verified fit with the existing pipeline

| Existing implementation | Integration consequence |
| --- | --- |
| `subworkflows/quantification.nf` emits `QuantifiedSampleBatch` with `batch_index`, `report_sequence`, and the latest quantification for every sample. Restored and stopped samples are retained. | Consume `quantification.out.batches`; do not rerun Oarfish or consume individual asynchronous sample updates. |
| `bin/de_analysis/edgeR_analysis.R` applies `filterByExpr`, `normLibSizes`, then exports CPM as `feature_counts.tsv`. | Prepare ICA expression directly from the complete Oarfish outputs. The DE matrix has the wrong transformation and a changing feature universe. |
| The inspected Oarfish `.quant` files contain `tname`, `len`, and `num_reads`. | Add explicit transcript-to-gene aggregation and a documented abundance transformation. There is no ready-made model-compatible TPM matrix in these files. |
| `modules/report_batches.nf` joins on `batch_index` and restores `report_sequence` order. `main.nf:qc_report` serializes report generation. | Carry ICA results/status through the same joins and ordering mechanism. |
| `bin/workflow_glue/qc_report.py` builds a live shell and self-contained final snapshots. Plot modules use Bokeh/ezcharts. | Extend this report with an independent ICA loader and tab, retaining offline final-report behavior. |
| `qc_report_types/temporal_plots.py` depends on `DifferentialResult`, GSVA membership, and log-CPM. | Reuse display patterns and extract neutral metadata helpers; do not pass ICA activities through the GSVA result type. |
| `lib/validation.nf` currently requires gene-set enrichment for `timeline_analysis`. | Allow timeline analysis when either gene-set enrichment or ICA is enabled; apply the corresponding backend checks independently. |

The source search found no existing ICA inference stage or model bundle in the
pipeline. Existing carbon-stress GMT collections are gene sets, not ICA models:
they do not contain the continuous signed weights required for projection.

## Pseudomonas model and example compatibility

The shared analysis identifies the published KT2440 `putidaPRECISE321` model.
Inspection of the [SBRG repository](https://github.com/SBRG/modulome_ppu/tree/f63a0dfab124ea9022123ce6227d132f39de0108)
at commit `f63a0dfab124ea9022123ce6227d132f39de0108` confirmed:

- `data/interim/ica_runs/220/S.csv` is the gene-weight matrix used by the
  characterization notebook: **5,564 genes × 84 components**, with `PP_####`
  locus IDs. Here `S` is the matrix called `M` in PyModulon.
- `data/interim/ica_runs/220/A.csv` contains the corresponding training
  activities; `data/raw_data/sample_table.csv` has 321 sample records.
- `data/processed_data/imodulon_table.csv` contains 84 curated module records.
- The characterization notebook renames components and applies manually
  adjusted membership thresholds. The importer must preserve the published
  component identity, sign, scale, names, and curated memberships.
- The notebook writes `putidaPRECISE321.json`, but that generated JSON is absent
  from the inspected Git tree. Do not assume that filename is downloadable from
  the repository. Prefer a verified curated export, or build and validate a
  deterministic adapter from the pinned matrices and curation records.

These file roles and curation steps are documented in the
[characterization notebook](https://github.com/SBRG/modulome_ppu/blob/f63a0dfab124ea9022123ce6227d132f39de0108/notebooks/5_iModulon_characterization.ipynb).
Membership thresholds are for annotation/display; projection uses the full
weight matrix, including weights below those thresholds.

The repository's supplied example reference is **NBRC 14164**, whereas this
model is **KT2440**. The existing
`data/gene_sets/p_putida_kt2440_to_nbrc14164_mapping.tsv` contains 334 source
proteins, of which 307 map to 305 target loci. Direct intersection with the
downloaded model covers **307/5,564 genes = 5.52%**. This is coverage of the
available mapping file, not a measure of genome-wide conservation between the
strains. It is insufficient to treat this example as a validated projection.

For KT2440 data, resolve annotation aliases to the model locus IDs. For NBRC
14164, supply a genome-wide orthology map and label results as transferred from
KT2440. Default to unambiguous one-to-one orthologs; exclude and report collisions
rather than copying one target measurement into several model genes. The
existing mapping includes target collisions and cannot be reused unchanged.
Sequence correspondence alone also does not establish conserved regulation.

## Expression preparation and inference

Use a model-declared transformation rather than hardcoding a normalization for
all organisms and models. The published
[normalization notebook](https://github.com/SBRG/modulome_ppu/blob/f63a0dfab124ea9022123ce6227d132f39de0108/notebooks/2_expression_visualization.ipynb)
subtracts each project's reference-condition mean **in log-expression space**.
It does not center over all experimental samples. For this pipeline, explicitly
select baseline samples, normally the pre-starvation/control replicates.

For matched genes `G`, the proposed calculation is:

```text
B[g,s] = model-compatible gene abundance
L[g,s] = model-declared log transform of B[g,s]
X[g,s] = L[g,s] - mean(L[g,r] for reference samples r)
A       = pinv(M[G,:]) @ X[G,:]
E       = X[G,:] - M[G,:] @ A
```

The pseudoinverse over shared genes matches
[PyModulon's inference implementation](https://pymodulon.readthedocs.io/en/latest/_modules/pymodulon/util.html#infer_activities).
It does not perform normalization, validate strain compatibility, or assess
projection quality for the caller.

The abundance adapter is a necessary validation step. Oarfish estimates
long-read transcript abundance; applying short-read count/length normalization
automatically would be inappropriate. The
[Oarfish paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12261437/) models read
counts as proportional to transcript copies without the additional transcript
length factor used for fragmented short reads. A candidate adapter is
million-scaled transcript abundance, aggregated to genes before the model's
log transform. Verify its suitability for this bacterial annotation and library
preparation, including the normalization denominator and RNA-biotype exclusions.
The exact PRECISE321 log base and pseudocount must be confirmed and recorded
before freezing that adapter; the inspected normalization notebook begins with
already logged values. Do not promise numerical equivalence merely because two
tables are both called TPM.

Required behavior:

- Resolve quantification features through the run annotation, then through an
  explicit model-ID map. For gene-level models, aggregate compatible transcript
  abundances before logging. Reject ambiguous parentage. Additional species may
  need another supported expression adapter; unsupported conventions should
  produce a clear error.
- Keep the model gene intersection fixed for the run, independent of edgeR
  filtering and observed batch expression. Distinguish quantified zero from
  absent or unmappable features. Never fill absent model genes with zero.
- Normalize over the declared expression universe before selecting the model
  intersection, avoiding a denominator based only on successfully mapped genes.
- Fit all model components jointly; selecting a module in the report only
  changes display. Do not refit a small subset or normalize component columns.
- Record matrix rank, conditioning, shared-gene coverage, per-component retained
  squared-weight fraction, curated-member coverage, and per-sample residual
  error. A rank-deficient basis must not produce apparently identifiable
  individual activities. Quality thresholds require calibration, not arbitrary
  universal percentages. Return explicit unavailable/low-confidence states.
- Keep baseline sample identities fixed. If their live quantifications change,
  recompute all sample activities together and record the baseline snapshot.
  Finalized baseline samples give a more stable live reference. Report scores
  as relative to this experiment's reference, not directly interchangeable
  with original compendium activities.

## Proposed model and result contracts

Accept one local, versioned model bundle per invocation. No runtime network
lookup or automatic model selection is needed. Pseudomonas is an example bundle,
not a species-specific code path.

```text
model/
  manifest.json       # schema/model version, source, organism/strain/assembly,
                      # ID namespace, feature level, transform, checksums
  weights.tsv         # model_gene_id × stable component_id; signed full weights
  components.tsv      # component_id, display name, annotation/category, source
  memberships.tsv     # published membership/thresholds, when available
```

The manifest must define supported abundance preparation, log base, pseudocount,
centering convention, and any model-specific scaling. Plain TSV/JSON keeps the
core independent of PyModulon serialization. A PyModulon/iModulonDB importer can
produce this contract. Training `A`/`X` are useful validation resources, but are
not runtime requirements for fixed-basis inference. Generic ICA bundles may lack
curated memberships; those should still support activities, with gene loadings
shown without invented membership calls.

Proposed switches:

| Parameter | Behavior |
| --- | --- |
| `--ica_model` | Optional bundle path; presence enables ICA and requires both reference inputs. |
| `--ica_gene_map` | Optional explicit run-gene → model-gene mapping; required for transferred namespaces/strains that annotation aliases cannot resolve. |
| `--ica_reference_group` | Baseline group, initially `control`; validate its existence and usable abundance independently of DEA readiness. |
| `--timeline_analysis` | Existing opt-in timeline switch, extended to permit ICA without GMT/GSVA. |

No separate enable flag is necessary initially. Setting the existing DE and
gene-set switches to false should permit quantification + ICA + report. A single
baseline sample permits descriptive projection, but baseline uncertainty and
the absence of a replicate SD must be visible. Current sample-sheet minimum
control requirements must also be checked for ICA-only operation.

Proposed per-snapshot outputs:

```text
ica/<invocation_id>/batch_<batch_index>/
  activities.tsv          # component_id, sample_id, group, order, activity
  activity_summary.tsv    # component/group mean, SD when n>1, n
  gene_mapping.tsv        # feature/run gene/model gene, mapping status/evidence
  component_coverage.tsv  # retained weights and member coverage
  projection_qc.tsv       # residual and numerical diagnostics
  expression.tsv          # matched, transformed, centered expression for drill-down
  provenance.json         # model/map/reference hashes, settings, batch identity
  status.json             # ready/deferred/unavailable with a specific reason
```

Use a stable sample key derived from group and alias, with alias as a display
label, because aliases need only be unique within a group in the base workflow.
Namespace outputs across invocations so restarting batch numbering does not
overwrite prior results. Rebuild ICA from restored quantifications when the
model, map, transform, or baseline selection changes; these downstream changes
should not invalidate completed alignment/quantification work.

## Workflow and report changes

Add `subworkflows/ica_projection.nf` and a small Python inference/validation CLI.
Validate the model and annotation map once, then consume each complete
quantification batch. Keep this branch independent of the edgeR readiness gate,
GMT input, fry, and GSVA. Numerical projection can use NumPy/SciPy; no training
framework or model-training process is required.

Extend records in `lib/sample.nf`, orchestration in `main.nf`, and
`modules/report_batches.nf` to carry ICA result paths and status. Emit exactly
one ICA result/status for every input report sequence, including deferred
batches, so report joins do not drop or indefinitely wait for a sequence.
Malformed models/mappings are configuration errors; insufficient live data
should yield a reportable deferred state. Preserve the existing serialized
publishing and final-snapshot mechanism.

Add `qc_report_types/ica_plots.py`, `--ica-results` to the report CLI, and an
independent primary tab with:

1. **Overview:** model/strain, baseline, mapping coverage, fit diagnostics, and
   an activity heatmap across samples. Show raw model-scale activities by
   default; any row-standardized display must be explicitly labeled.
2. **Component view:** searchable component selector; individual sample values,
   group means, SD and replicate counts; curated annotation and signed gene
   weights; optional links to existing gene-level DE results when available.
3. **Time course:** when enabled, activity versus sample-sheet `order`, showing
   individual biological replicates and group summaries, plus a member-gene
   expression view using ICA's own transformed data. The current metadata
   supports one group per time and independent replicates, not paired subjects
   or multiple treatment trajectories. Preserve that scope and state it.

Keep biological time separate from cumulative sequencing batch. Batches are
repeated estimates from accumulating reads, not additional biological
replicates. Render unavailable diagnostics as unavailable, including residual
ratios with a zero denominator. The final HTML should embed the plotted data
and work without iModulonDB or an external service.

Use the term **inferred component activity**. A regulatory annotation is an
interpretive association; projection does not identify a causal driver, measure
TF protein activity, or establish biochemical pathway flux. Do not import
PyModulon DIMA significance thresholds or add p-values in this first version.

Update `nextflow.config`, `nextflow_schema.json`, parameter/sample-sheet
validation, runtime dependencies/container packaging, and user documentation.
All changed `.nf` code must follow the repository's Nextflow 26.04 typing rules
and pass `nextflow lint .`.

## Acceptance criteria and remaining decisions

Before considering the implementation usable:

- Verify a pinned model import, component renaming, signs, and curated
  memberships. Check algebraic recovery on `X = M @ A_known`, then quantify
  projection agreement/residuals using published expression and activities;
  real-data agreement need not be exact because ICA is an approximation.
- Verify row-order invariance, annotation aggregation, ambiguous mappings,
  missing versus zero features, duplicate IDs, nonfinite values, deficient rank,
  absent/zero-depth references, and baseline recomputation.
- Exercise ICA-only and combined analyses, deferred live batches, out-of-order
  completion, stopped/restored samples, and restart output identities. Confirm
  each report contains matching quantification and ICA snapshots.
- Verify report selectors, negative activities, missing SD, optional memberships,
  optional DE, timeline validation, and final offline behavior. Reuse current
  report-batch/finalization test infrastructure for those guarantees.
- Demonstrate a second model with non-Pseudomonas identifiers to detect hidden
  species-specific assumptions.
- For the NBRC 14164 example, obtain a genome-wide one-to-one map and measure
  per-component coverage/conditioning. Benchmark the Oarfish expression adapter
  against compatible reference data and assess sensitivity to missing genes
  before interpreting biological activity patterns.

The main engineering work is the expression/model adapter and batch-safe report
integration. The matrix calculation itself is small. The key unresolved
scientific requirements are the exact published transform, long-read-to-model
compatibility, and genome-wide mapping for the NBRC example. These should be
resolved as explicit adapter validation, while keeping the pipeline interface
generic and the first release limited to inference and reporting.
