# Quality control and reports

> **_Author Note:_**  This documentation was generated using OpenAI GPT-6.1 and manually checked for correctness.

The workflow writes `report.html` to the output directory. The **Quality
Control** tab presents read counts, alignment summaries, and distributions
of read length, read quality, and mapping quality for each sample. When
both reference files are available, **Transcript Biotypes** shows the
distribution of assigned reads across transcript biotypes.
[Quantification guide](quantification.md)

Additional tabs appear when results from the enabled analyses are available.

- **Differential Expression** presents an overview of the samples and
  results for each group comparison. [Differential expression guide](differential_expression.md)
- **Gene Set Enrichment** shows GSVA scores and fry enrichment results.
  When the sample sheet includes time points, **Time Course** appears under
  **GSVA scores**, after **Grouped Distributions**.
- **Independent Components** shows iModulon activities and comparisons with
  the control group. [iModulon guide](imodulon_analysis.md)
- **Automatic Termination** appears after **Quality Control** when stopping
  checks are enabled. It shows the selected method, how much collecting more
  reads is expected to improve the precision of the results, whether the
  stopping rule was met, and what action was taken. The tab shows a waiting
  message until the first check is available.
  [Automated stopping guide](dea_termination.md)

## Reasons in the termination tab

| Reason | What it means |
| --- | --- |
| Precision check passed | Recent and predicted improvements are below the cutoff. Check **Stopping rule met** to see whether enough checks passed over the required time. |
| More reads may improve precision | Recent or predicted improvement is still at or above the cutoff. |
| Waiting for a stable uncertainty trend | The forecast is unusable, the trend has yet to flatten, or recent uncertainty is unstable. |
| Waiting for more sequencing data | The run is warming up, read arrivals are sparse or slow, or a library is waiting for new reads. |
| Waiting for usable analysis results | Results are pending, unavailable, delayed, or out of date. |
| No stop check applies | The sample is a fixed control, has no comparison to monitor, or its acquisition has ended. |

Samples used in several comparisons can show more than one reason.

## Viewing a live report

The report updates when new results are available. When precision monitoring
is disabled, each live update normally requires a new BAM file from every
sample still receiving data.

Keep `report.html`, `report_state.json`, and the
`report_snapshot_revision_*.html` files in the same directory during
analysis. If the browser remains on **Waiting for update** when opening the
report from disk, serve the output directory locally. Replace `output` with
the path to the output directory.

```bash
python3 -m http.server 8000 --bind 127.0.0.1 --directory output
```

Open [the local report](http://127.0.0.1:8000/report.html) in a browser.

The final `report.html` is self-contained and can be shared after the
workflow completes.

## Reusing completed samples

After a run with quantification completes, final counts and quality control
tables are available under `<out_dir>/<group>/<alias>/`. A `FINAL` file
identifies a completed sample.

To add samples to an experiment, update the sample sheet and use the same
`out_dir`. seq_lm reuses completed samples with matching inputs and rebuilds
the combined report. Reuse is independent of Nextflow's `-resume` option.

Use a new output directory to recompute results after changes to the BAM
files, reference files, or strand setting. Run one workflow at a time per
output directory.
