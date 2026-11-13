# Differential expression and gene sets

> **_Author Note:_**  This documentation was generated using OpenAI GPT-6.1 and manually checked for correctness.

seq_lm uses edgeR for differential expression analysis (DEA) of transcript
abundance between sample groups. Supply a reference genome and annotation
and add the following settings to the `params` block of the run configuration.
[Reference file guide](quantification.md)

```groovy
differential_expression = true
gene_set_enrichment = false
```

Both options default to `true`. Set `gene_set_enrichment = false` to run
DEA without a gene-set file.

## Group comparisons

The default analysis compares each target group with the group named `control`.
It requires at least two samples in `control` and at least one target group.

For custom comparisons, set `gene_contrast_file = "contrasts.csv"`. Use a
headerless CSV with one comparison per row. Enter the target group in the
first column and the reference group in the second. Use the group names
from the sample sheet.

```csv
treated,baseline
recovery,treated
```

A positive log₂ fold change indicates higher expression in the target group.
DEA and gene-set analysis use the groups listed in the comparisons.
The quality control report includes all samples.

By default, each group in a comparison must contain at least two samples
with at least 10,000 assigned reads each. All samples in the comparison must
have positive assigned abundance. Use `min_read_count` and
`min_replicate_sample_count` to change these thresholds. During live analysis,
comparisons may become available at different times and estimates may change
as additional reads are included.

## Results and interpretation

Open **Differential Expression** in the HTML report for an overview of the
samples and a volcano plot, MA plot, and heatmap for each comparison.
`de_lfc_cutoff` sets the minimum absolute log₂ fold change (default `1.0`).
`de_padj_cutoff` sets the adjusted p-value cutoff (default `0.05`).

The latest published results are in `<out_dir>/differential_expression/latest/`.
`contrasts.tsv` links comparison directories to group names, and
`contrast_status.tsv` records which comparisons met the analysis criteria.
For each completed comparison, the result table is `edgeR_results.tsv`.
Earlier results are in the adjacent `batch_*` directories.

## Gene-set enrichment

Set `gene_set_enrichment = true` and `gene_sets = "pathways.gmt"` with
differential expression enabled. A GMT file has one gene set per line,
with the set name, description, and member IDs separated by tabs.

```text
set_A	Example gene set	gene_1	gene_2	gene_3
```

Use member IDs that match the quantification output or the reference
annotation, such as gene IDs or locus tags. Each set needs at least two
mapped features that pass expression filtering. For GSVA, at least two of
those features must vary across samples. Consult the coverage tables if an
expected set is missing.
[Example gene-set files](../data/gene_sets/README.md)

Open **Gene Set Enrichment** in the report.

- **fry enrichment** tests for coordinated changes in expression within a
  gene set. Positive bars indicate higher expression in the target group,
  and negative bars indicate higher expression in the reference group.
  Full results are in `fry_results.tsv`.
- **GSVA scores** presents sample-level scores from gene set variation
  analysis (GSVA). Scores are available in `gsva_scores.tsv`.
- **GSVA differential** compares GSVA scores between groups. Positive
  differences indicate higher scores in the target group. Test results
  are in `gsva_limma_results.tsv`.

## Time course

Add an `order` column to the sample sheet to plot GSVA scores and expression
over time. Record elapsed time in whole minutes for every sample in the
selected groups. Replicates in a group share the same time point. Assign
one group per time point and include at least two time points. Open
**Gene Set Enrichment → GSVA scores → Time Course** to see the plots.

DEA can also inform when to end a live sequencing run.
[Automated stopping guide](dea_termination.md)
