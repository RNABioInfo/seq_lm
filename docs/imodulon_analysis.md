# iModulon analysis

> **_Author Note:_**  This documentation was generated using OpenAI GPT-6.1 and manually checked for correctness.

seq_lm estimates iModulon activity for each sample using an existing
independent component analysis (ICA) model. Supply a reference genome and
annotation and add the following settings to the `params` block of the run
configuration.
[Reference file guide](quantification.md)

```groovy
ica_analysis = true
ica_matrix = "/data/weights.tsv"
```

To run iModulon analysis with quality control and quantification, set
`differential_expression` and `gene_set_enrichment` to `false`.

## Sample groups

The default control group is `control`. To use a different group, set
`ica_control_group`. The control group must contain at least two samples.
To select specific groups, set `ica_condition_groups` to a comma-separated
list of group names from the sample sheet, such as `"early,late"`.
The analysis always includes the control group and defaults to all groups
when `ica_condition_groups` is unset.

By default, each selected group must contain at least two samples with at
least 10,000 assigned reads each. All selected samples must have positive
assigned abundance. Use `min_read_count` and `min_replicate_sample_count` to
change these thresholds. The report indicates when selected groups have
insufficient assigned reads.

## Model files

Supply a CSV or TSV of gene weights with a header row. Use the first column
for gene IDs and each remaining column for one component. Use unique gene
and component IDs and retain the model's original weights.

```text
model_gene_id	component_1	component_2
gene_A	0.35	-0.08
gene_B	-0.20	0.42
```

Choose a model compatible with the reference annotation and expression scale.
seq_lm sums transcript counts by gene and scales abundance to counts per
million assigned reads. By default, it adds a pseudocount of 1 before taking
the base-2 logarithm. Set `ica_pseudocount` and `ica_log_base` to change
these values.

seq_lm matches model genes to the reference annotation by gene ID or locus
tag. If the IDs differ, set `ica_gene_map` to a CSV or TSV with the columns
`gene_id` and `model_gene_id`. Use the annotation's gene ID in the first
column and the model's gene ID in the second. Each ID can appear only once.

The default requires a match for every model gene. To allow partial coverage,
set `ica_min_gene_coverage` to a lower fraction, such as `0.9`. Review gene
matches in `gene_mapping.tsv` and coverage by component in
`component_coverage.tsv`.

To display component names and functions, set `ica_imodulon_table` to the
CSV or TSV annotation table from the same model. Provide one row per
component with a `name` and its exact model ID in `component_id`, `k`, or
the first column. You can also include `function`, `regulator`, and `category`.

## Results and interpretation

Open **Independent Components** in the HTML report. **Sample component
activities** shows activities relative to the control group. **Component
details** shows sample distributions and group comparisons. The estimated
effect is the target group mean minus the control group mean, in the units
of the supplied model. Positive differences indicate higher activity in the
target group. `ica_padj_cutoff` sets the adjusted p-value cutoff (default `0.05`).

For time-course plots, add an `order` column with elapsed time in whole
minutes for every selected sample. Replicates in a group share the same time
point. Assign one group per time point and include at least two time points.
The plots appear under **Time course**.

`<out_dir>/ica/latest.json` records the path and status of the latest
analysis (`ready` or `deferred`). Results marked `ready` include
`activities.tsv` for sample activities and `differential_activity.tsv` for
group differences, confidence intervals, and adjusted p-values.

## Sequencing termination

Set `stopping_method = "ica"` to base stopping on the expected reduction in
uncertainty of activity differences. This mode requires at least one
target group and the ICA control group. By default, the stopping rule uses
all components. To use a subset, set `ica_stopping_components` to a text
file with one exact component ID per line.

See the [automated stopping guide](dea_termination.md) for monitoring and MinKNOW setup.
