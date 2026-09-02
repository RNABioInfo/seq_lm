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
all enabled outputs for a synchronized batch are complete. With paired
reference inputs, this includes the matching cumulative Oarfish biotype summary
and, when enabled, edgeR readiness or results. During
live analysis, raw NanoPlot and flagstat tables remain report inputs in
Nextflow work directories. When the stream completes, they are also persisted
with the sample's `FINAL` checkpoint under `--out_dir` so later invocations can
restore QC without rescanning BAMs. The stable HTML shell is the user-facing
entry point; a small state file points it to the immutable snapshot for the
newest complete batch.

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
after the frame update using hierarchy-based tab keys, so the browser page
itself is not reloaded. Hidden Bokeh and ECharts views are resized when their
tabs open. The report has three primary tabs:

* **Quality Control** contains the existing read-flow, metrics, read-length,
  read-quality, mapping-quality, and sample views. When both reference inputs
  are present, Read flow also contains a cross-sample 100% stacked horizontal
  bar chart of Oarfish abundance by transcript biotype beneath the sample
  Sankey selector.
* **Differential Analysis** separates its overview from contrasts, then uses
  contrast and plot-type subtabs for logFC-versus-logCPM, volcano, and top-gene
  heatmap plots.
* **Gene Set Enrichment** separates raw GSVA scores, GSVA limma testing, and
  directional fry. GSVA includes score heatmaps, raw score distributions,
  coverage, a multi-contrast dot plot, and per-contrast volcano and heatmap
  views. GSVA heatmap rows retain their selection order and sample columns
  retain metadata order; their scales map low values to blue and high values
  to red.
  GSVA views show gene-set identifiers rather than descriptions. fry retains
  the gene-set dropdown used for barcode plots. Concise
  identifiers appear in long dropdowns, while barcode views expose complete
  labels and statistics in **Gene-set details**.

Nested tab content uses reduced padding, while the primary tab bar uses the
workflow brand color to remain visually distinct from analysis, contrast, and
plot-type subtabs.

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
at least one BAM. Later QC checkpoints wait for one new BAM from every currently
active live sample, retain final or stopped samples' QC inputs, and refresh the
cumulative report. A sample leaves later batch barriers after its `STOP` is
drained, so other samples may continue independently. Only live samples require
`STOP` files. The workflow treats `alias` as the
sample name and requires at least two rows where `group` is `control`.

QC-only startup batches are not published as incomplete reports. Without
references, reports follow the QC batch stream. With references, a report waits
for the matching complete cumulative quantification/biotype batch; if
differential expression is enabled it also waits for the matching readiness or
result record. This prevents report figures from mixing live batches.

Required local Docker images:

* `seq_lm/quality_control`: provides NanoGet for the per-read QC table.
* `seq_lm/samtools`: provides samtools for sorting, merging, indexing, and
  flagstat.
* `seq_lm/dea_r`: provides edgeR and writes the differential report inputs.
* `seq_lm/report`: provides ezCharts, Bokeh, pandas, and scikit-learn for the
  combined HTML report.

The chunk-level QC workflow intentionally does not consume the reference
annotation. Transcript-biotype composition is produced by the standalone
quantification workflow and routed into the Quality Control → Read flow report
panel. See [Quantification and transcript-biotype QC](quantification.md).
