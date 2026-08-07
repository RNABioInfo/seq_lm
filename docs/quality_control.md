## Analysis Report and Quality Control Workflow

The quality control workflow runs once per sample chunk. A chunk is the set of
new BAM files observed for one sample at a synchronized live-analysis batch
index. Multiple startup BAMs are block-concatenated, while a one-BAM chunk is
reused directly. Neither path sorts or indexes the chunk, so each QC result
describes one sample chunk rather than a final accumulated sample BAM.

The workflow currently runs two tabular QC steps:

* `bam-qc-table` performs a sequential scan using NanoGet's per-alignment metric
  calculations and writes the NanoPlot-compatible `NanoPlot-data.tsv.gz`
  without requiring a coordinate index or rendering unused plots.
* `samtools flagstat -O tsv` reads the same chunk BAM and writes a TSV
  alignment summary.

The main workflow publishes an EPI2ME-displayable HTML analysis report after
both the QC and edgeR outputs for a synchronized batch are complete. Raw
NanoPlot, flagstat, and edgeR tables remain report inputs in Nextflow work
directories. The stable HTML shell is the user-facing entry point; a small
state file points it to the immutable snapshot for the newest complete batch.

Expected published output layout:

```text
qc_report.html
qc_report_state.json
qc_report_snapshot_<batch_index>.html
```

`qc_report.html` retains its historical filename for compatibility, but its
display title is `seq_lm analysis report`. It is published at the workflow
output root with overwrite enabled so EPI2ME can display one stable report path
while live chunks arrive. The shell polls `qc_report_state.json` and refreshes
only its internal report frame when a complete snapshot becomes available.
Active Bootstrap tabs, dropdown selections, and scroll position are restored
after the frame update, so the browser page itself is not reloaded. The report
has exactly three top-level sections:

* **Quality Control** contains the existing read-flow, metrics, read-length,
  read-quality, mapping-quality, and sample views.
* **Differential Analysis** contains a global PCA and one tab per edgeR
  contrast with logFC-versus-logCPM, volcano, and top-gene heatmap plots.
* **Gene Set Enrichment** contains one tab per edgeR contrast with a signed
  directional-fry summary and the same Bootstrap dropdown-tab selector used by
  read-flow samples. Concise labels appear in the selector, while each barcode
  view exposes its complete label and statistics in **Gene-set details**.

The report sample TSV uses this shape:

```text
name	group	chunks_seen	latest_batch_index	qc_dir
rep_1	control	2	2	qc_results
```

Inside the report work directory, `qc_results` contains the raw QC tables for
all chunks currently known to the report:

```text
qc_results/
  nanoplot/
    <group>/
      <name>/
        nanoplot_data_chunk_<batch_index>.tsv.gz
  flagstat/
    <group>/
      <name>/
        flagstat_data_chunk_<batch_index>.tsv
```

The BAM sample sheet uses `alias`, `group`, and `bam_dir` columns, plus optional
`is_live` and `order` columns. `is_live` accepts case-insensitive `true` or
`false`; a missing or blank value defaults to `true`. A sample is watched only
when both `--live_analysis` and its row-level `is_live` value are true.

Final samples are processed once from all BAMs present at startup and must have
at least one BAM. Later QC checkpoints wait for one new BAM from every live
sample, retain the final-sample QC inputs, and refresh the cumulative report.
Only live samples require `STOP` files. The workflow treats `alias` as the
sample name and requires at least two rows where `group` is `control`.

QC-only startup batches are not published as incomplete reports. The first
report appears once edgeR has a current quantification for every experiment
sample; later reports are refreshed once per complete QC/edgeR batch pair.

Required local Docker images:

* `seq_lm/quality_control`: provides NanoGet for the per-read QC table.
* `seq_lm/samtools`: provides samtools for sorting, merging, indexing, and
  flagstat.
* `seq_lm/dea_r`: provides edgeR and writes the differential report inputs.
* `seq_lm/report`: provides ezCharts, Bokeh, pandas, and scikit-learn for the
  combined HTML report.

The QC workflow intentionally does not consume the reference annotation. The
previous `bamstats`-style annotation-aware QC path has been removed from the
active chunk QC workflow.
