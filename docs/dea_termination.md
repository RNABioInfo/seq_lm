# Curve-based R stopping after BAM-triggered differential expression

The existing live BAM trigger discovers stable new files. Each affected sample is
quantified, DEA fits the configured target/reference pair, and the R check consumes
that fit immediately. It does not wait for GSVA, ICA or report rendering. There is
no separate R watcher, periodic analysis clock, or stream of empty time batches.
File-discovery polling remains for reliable operation on mounted filesystems.

```groovy
params {
    stopping_method = 'dea'              // alternative: 'ica'
    monitoring_behavior = 'log'         // disabled | log | terminate
    dea_stopping_genes = null           // optional GMT file; defaults to all features
    dea_gain_cutoff = 0.01
    dea_information_radius_scale = 1
    dea_forecast_minutes = 30
    dea_max_processing_lag_seconds = 600
    dea_min_points = 6
    dea_min_seconds = 3600
    dea_confirmations = 3
    dea_interruption_fraction = 0.25
    dea_recovery_fraction = 0.5
    dea_recovery_intervals = 2
}
```

The default means **less than 1% predicted reduction in aggregate R over the next
30 minutes**. For hourly BAM output use `--dea_forecast_minutes 60`; the same cutoff
then means 1% over an hour. Changing the horizon changes the scientific question,
not just the update frequency. At least six informative observations and one hour
are always required, so hourly files require several hours before eligibility.

## Input timing and independent contrasts

Every accepted sample update has its own timestamp and monotonic version. Several
BAMs accepted for one sample in a discovery pass are grouped. Collation preserves
that sample's input order without waiting on other samples. Completed
quantifications replace older versions only; a slow old task cannot roll back the
latest sample state. DEA snapshots use the latest completed inputs, including
reusable fixed full-depth controls.

A contrast adds an R observation only when every still-live member has more
assigned reads than at its previous evaluated snapshot. Partial updates accumulate
until this holds. Unrelated samples and repeated fits do not advance or reset its
confirmation streak. Natural completion keeps a library's collected depth in the
model and resets the affected contrast's confirmations.

File availability and analysis completion are distinct. Batch indices are IDs,
never clocks. Existing startup BAMs establish initial depth without historical
arrival observations. The observation origin persists across restarts; startup
versions after a restart establish fresh rate baselines rather than interpreting
the downtime as zero yield.

The checks deliberately use file-availability intervals. They estimate delivery
rates, which can differ from physical sequencing rates during basecalling or
transfer delays. No scan or pore-check metadata is required. No output means
waiting, not evidence that enough information has been collected.

## Model-based gene intervals and aggregate R

Each configured contrast uses only its target and reference groups, with the same
edgeR fit for DEA, fry and interval calculations. Readiness is evaluated per pair;
unavailable or insufficient-depth groups are deferred and do not enter GSVA as
biological zero-expression samples. The unfiltered quantification feature IDs
freeze the denominator before expression filtering.

The staged `profile_intervals.R` solver inverts the fitted gene's conditional
negative-binomial profile test at 95%, using fitted offsets, effective dispersion
and empirical-Bayes parameters. For unshrunk log2 fold change beta and interval
[L,U], the gene radius is:

```
R_g = max(beta - L, U - beta)
c_g = 1 / (1 + (R_g / s)^2)
I   = mean(c_g) over the fixed feature universe
R_information = s * sqrt(1/I - 1)
```

Filtered, missing and unbounded intervals contribute zero while remaining in the
denominator. Zero contribution does not mean zero uncertainty. If I is zero,
aggregate uncertainty is unavailable/unbounded and cannot justify a stop. The
configurable scale s defaults to 1 log2FC unit and is not estimated from the data.

### Optional stopping gene sets

Use `--dea_stopping_genes targets.gmt` to calculate the stopping uncertainty score
over the union of all gene sets in a GMT file. Overlapping members and aliases
that resolve to the same feature count once. This input is independent of
`--gene_sets` and works with `--gene_set_enrichment false`. It is consumed only
when `stopping_method=dea` and monitoring is enabled.

The existing GMT parser and annotation mapping resolve members against the
unfiltered count-matrix feature IDs. Gene names can resolve to multiple transcript
features; all resolved features enter the union, without gene-level aggregation.
Malformed GMT records, duplicate set names, empty selections and selections with
no matching features fail. Partially matched collections use the resolved union
and warn about unmatched members. Inspect `dea_stopping_resolution.tsv` and the
deduplicated `dea_stopping_features.tsv` in each fitted analysis output.

Only the monitoring universe changes. The full edgeR fit, interval calculations,
DEA/enrichment outputs and library read depths are retained. Selected features
filtered out by expression or lacking finite intervals remain in the denominator
with zero information contribution. Every contrast uses the same fixed union;
unexpected changes fail, and changing the GMT contents resets monitoring history
and eligibility. Fitted snapshots must contain the selection export; deferred
analyses remain deferred. Monitoring audits and reports identify the uncertainty
scope and feature count.

The stopping rule still measures diminishing expected improvement, without an
absolute precision requirement or a significance filter. A selected-feature stop
does not guarantee sufficient precision for every selected feature or for features
outside the selection. Omit `--dea_stopping_genes` to retain all-feature monitoring.

## Forecast over the configured horizon

Let H be `dea_forecast_minutes` in seconds. Each live sample's rate intervals are
consecutive observed cumulative-depth differences divided by actual elapsed time.
Rolling totals integrate these rates, clipping intervals at the window boundaries.
The forecast is:

```
additional_reads = H * max(
    average rate over the preceding H,
    duration-weighted 75th percentile of interval rates over the preceding 3H
)
```

The weighted percentile is the smallest rate reaching 75% of total interval
duration. At least three valid intervals are required. Intervals longer than H are
excluded from the productive-rate estimate; missing output is never zero-imputed.
Automatic stopping requires adequate coverage over the preceding H and latest data
no older than H for every live member. Hourly output therefore cannot qualify with
a 30-minute horizon. Available coarse forecasts remain reportable, but do not
bypass this guard.

Fit the nonnegative curve `R_information^2 = a + b/(D/1e6)` to all informative
observations, where D includes only monitored libraries' depth, not fixed controls.
Keep the collected depth of naturally completed libraries. The forecast uses only
this depth curve: `curve_gain = 1 - fitted_R(D + additional_reads) / fitted_R(D)`.
The comparison is between two points on the fitted curve; observed current R is
used separately for fit-discrepancy and observed-improvement checks. There is no
gene-model depth projection or second forecast branch.

Add a margin of at least 0.002 (0.2 percentage points) to the curve gain, increased
by the 90th percentile of positive underprediction errors from the last six
completed curve forecasts. Each error is the later observed relative R reduction
minus that forecast's unadjusted curve gain. This calibration is specific to the
curve-only predictor. A forecast matures only after observed R points
bracket its deadline. Horizon comparisons interpolate log R between observed
points at most H apart; they never extrapolate. Later observations update forecast
validation, never earlier decisions.

The nonnegative fit, elbow, fit-discrepancy, observed-improvement and raw-R stability
checks remain. Gain and observed improvement must be strictly below the cutoff.
Three consecutive qualifying fresh observations are required by default, spanning
at least 2H/3 (20 minutes at the default horizon). A burst of files cannot confirm
stopping in seconds.

## Interruptions, processing lag and actions

An observed horizon rate below 25% of the preserved healthy reference enters an
interruption state. A long low-yield interval can also expose an interruption when
files resume. Preserve that reference throughout the interruption. Recovery needs
two consecutive valid intervals at least 50% of the reference, together spanning
at least H/3. Failed checks and inadequate cadence reset confirmations.

Processing lag is independent of output cadence: if analysis completion is more
than `dea_max_processing_lag_seconds` after the snapshot's latest represented live
data, automatic stopping is withheld. Data age is checked for every live member.
There is no low-yield timeout or metadata-only stopping fallback.

Before issuing an acquisition stop, recheck current eligibility, input versions,
data age and processing lag. A known newer BAM not represented in the fit blocks
the action. Every required live member of a shared MinKNOW acquisition must qualify.
Log mode issues no stop requests. Stop failures are recorded and retried on later
eligible results. STOP markers and final-output draining retain their existing
behavior. ICA uses the same acquisition actions without depending on R.

## Outputs, resume and migration

Precision monitoring now uses `stopping_method=dea`, `dea_` parameters and
`dea-termination` / `dea_termination` entry points. Former method and parameter names are not aliases. State is written under `termination/dea`; earlier
`termination/r` state and reports remain untouched and are not imported.


`termination/dea/clock.json` records the persistent observation origin and version
inventory. `state.json` uses **schema version 3**, with `forecast_method=depth_curve`, per-sample arrival histories,
productive references, contrast observations, forecast deadlines/errors and
confirmation spans. State and per-analysis audits are committed atomically.
Configuration, cohort, reference or algorithm identity changes reset eligibility.
Earlier monitoring state is incompatible; use a fresh output directory while reusing
verified fixed-control inputs.

Each `termination/dea/batch_<analysis_index>/` contains `audit.json`,
`contrast_monitoring.tsv`, `sample_monitoring.tsv`, `read_rates.tsv`,
`forecast_history.tsv`, `read_rate_history.tsv`, `R_history.tsv` and its state
snapshot. Monitoring tables include the raw `curve_gain`, `forecast_margin` and
their sum, `predicted_conservative_gain`. DEA output retains `precision.tsv`, `precision_samples.tsv`,
`precision_methods.tsv` and `precision_model.rds`. Reports show the configured
horizon, R against time/depth, forecasts against later observed gain, fractions
below R=0.5 and R=1, throughput, gaps, coverage, lag and confirmation spans.
The precision table contains feature IDs, interval bounds, unshrunk log2 fold
changes and radii. Fitted means and projection-only dispersion columns are not
exported. The saved RDS fit remains available for statistical reproducibility.

Only current settings are validated. DEA stopping parameters are checked when
`stopping_method=dea` and monitoring is enabled; ICA stopping parameters, including
`num_stable_batches`, are checked only for enabled ICA monitoring. ICA analysis
parameters remain independently validated whenever ICA analysis is enabled.
Obsolete parameter names have no custom migration checks or aliases.

The aggregate, forecast, margin, interpolation and stopping rule are **heuristics**,
not accuracy guarantees or sequential confidence bounds. Per-gene intervals are
conditional model-based intervals, not a guarantee for each gene after optional
stopping. A 1% aggregate improvement cutoff does not mean 99% correct gene calls.

## Validation and containers

Run `bash test/run_monitoring_tests.sh`, the Python regression suite, the native R
profile tests and `nextflow lint .`. The saved-data evaluator is
`test/evaluate_saved_dea_checkpoints.py REFERENCE_DIRECTORY OUTPUT`; it does not run
the sequencing replay. For the native out-of-order quantification test, run
`bash test/run_event_quantification_test.sh` with samtools and the report Python
dependencies available. It uses a synthetic quantifier and does not benchmark
biological quantification. Synthetic timing tests separately cover delivery delays,
bursts and interruptions that saved checkpoint times cannot reconstruct.

The live Python engine uses only the standard library. Interval calculations use
workflow-staged scripts and existing edgeR functionality; report dependencies are
already in the report image. No Docker image rebuild or publication is required.
