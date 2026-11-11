# Quality control and reports

> **_Author Note:_**  This documentation was generated using OpenAI GPT-6.1 and manually checked for correctness.

The workflow writes `qc_report.html` to the output directory. The **Quality
Control** tab presents read counts, alignment summaries, and distributions
of read length, read quality, and mapping quality for each sample. When
both reference files are available, **Transcript Biotypes** shows the
distribution of assigned reads across transcript biotypes.
[Quantification guide](quantification.md)

Additional tabs appear when results from the enabled analyses are available.

- **Differential Expression** presents an overview of the samples and
  results for each group comparison. [Differential expression guide](differential_expression.md)
- **Gene Set Enrichment** shows GSVA scores and fry enrichment results.
  **Time Course** appears when the sample sheet includes time points.
- **Independent Components** shows iModulon activities and comparisons with
  the control group. [iModulon guide](imodulon_analysis.md)
- **DEA Termination** or **ICA Termination** appears when stopping checks are
  enabled. It shows how much collecting more reads is expected to improve the
  precision of the results, whether the stopping rule was met, and what action
  was taken. The tab shows a waiting message until the first check is available.
  [Automated stopping guide](dea_termination.md)

## Viewing a live report

The report updates when new results are available. When precision monitoring
is disabled, each live update normally requires a new BAM file from every
sample still receiving data.

Keep `qc_report.html`, `qc_report_state.json`, and the
`qc_report_snapshot_revision_*.html` files in the same directory during
analysis. If the browser remains on **Waiting for update** when opening the
report from disk, serve the output directory locally. Replace `output` with
the path to the output directory.

```bash
python3 -m http.server 8000 --bind 127.0.0.1 --directory output
```

Open [the local report](http://127.0.0.1:8000/qc_report.html) in a browser.

The final `qc_report.html` is self-contained and can be shared after the
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
