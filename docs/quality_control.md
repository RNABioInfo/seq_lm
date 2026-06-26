## Quality Control Workflow

The quality control workflow runs once per sample chunk. A chunk is the set of
new BAM files observed for one sample at a synchronized live-analysis batch
index. The chunk BAMs are sorted, merged, and indexed before QC starts, so each
QC result describes one merged sample chunk rather than a final accumulated
sample BAM.

The workflow currently runs two tabular QC steps:

* `NanoPlot` reads the merged chunk BAM and writes `NanoPlot-data.tsv.gz`.
* `samtools flagstat -O tsv` reads the same merged chunk BAM and writes a TSV
  alignment summary.

The permanent report generator is not implemented yet. TEMPORARY: the main
workflow currently accumulates the NanoPlot and flagstat TSVs per sample over
time, writes placeholder manifests, and also publishes a throwaway
EPI2ME-displayable HTML report after every sample chunk QC update. That
temporary report only lists samples with QC results seen so far and must be
removed when the real live QC report is implemented.

Expected published output layout:

```text
<group>/<alias>/qc/
  chunk_<batch_index>/
    nanoplot/
      NanoPlot-data.tsv.gz
    samtools_flagstat/
      <sample_id>_<batch_index>.flagstat.tsv
  report_inputs/
    <sample_id>_chunk_<batch_index>_qc_report_inputs.tsv
temporary_qc_report.html
```

The `report_inputs` manifest contains one row per QC TSV currently known for
the sample:

```text
sample_id	batch_index	metric	qc_tsv
sample_1	0	nanoplot	...
sample_1	0	samtools_flagstat	...
```

TEMPORARY: `temporary_qc_report.html` is published at the workflow output root
with overwrite enabled so EPI2ME can display one stable report path while live
chunks arrive. The file intentionally uses snake_case naming and contains only a
temporary notice plus the current QC-result sample list.

Required local Docker images:

* `seqlm/quality_control`: provides NanoPlot.
* `seqlm/samtools`: provides samtools for sorting, merging, indexing, and
  flagstat.

The QC workflow intentionally does not consume the reference annotation. The
previous `bamstats`-style annotation-aware QC path has been removed from the
active chunk QC workflow.
