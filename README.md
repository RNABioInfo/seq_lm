# seq_lm

Performs live or one-time transcriptome analysis of aligned Oxford Nanopore BAM files.
The results are presented as an interactive, self-updating HTML report.
The report includes read quality control, Oarfish quantification, transcript biotype quantification, edgeR differential expression analysis, fry and GSVA gene-set analysis, and iModulon independent component analysis.
Analysis modes may be flexibly combined as needed.
Additionally, live analysis automates the workflow by integrating directly with MinKNOW.

[Quickstart](#quickstart) · [Analysis modes](#analysis-modes) ·
[Parameters](#analysis-and-live-run-parameters) · [Guides](#guides)

## Quickstart

Requires **Nextflow ≥26.04.3**, **Java 17**, and **Docker** or **Singularity/Apptainer**.
Nextflow fetches published images from Docker Hub.
Windows users must run inside **WSL2** with Docker Desktop's WSL integration;
EPI2ME Desktop on Windows is not supported unless the Nextflow version is updated.

### 1. Setup

Install Docker separately.
If Nextflow is not already installed, use the included environment.yml file to create a conda environment.

```bash
git clone https://github.com/RNABioInfo/seq_lm.git
cd seq_lm

# If Nextflow is not already set up
conda env create --file environment.yml
conda activate seq-lm
nextflow -version
```

### 2. Prepare sample sheet

Save a `samples.csv` file containing the paths to the aligned BAM directories:

```csv
alias,group,bam_dir,is_live
control_1,control,/data/bams/control_1,true
control_2,control,/data/bams/control_2,true
treated_1,treated,/data/bams/treated_1,true
treated_2,treated,/data/bams/treated_2,true
```

| Field | Required | Meaning |
| --- | --- | --- |
| `alias` | Yes | Sample name; unique within its group. |
| `group` | Yes | Condition label. Without a gene contrast file, DEA requires at least two `control` samples. ICA uses its independently configured control. |
| `bam_dir` | Yes | Existing directory searched recursively for BAMs; accessible to the container. |
| `is_live` | No | Watch this sample when `live_analysis` is also true. Blank or omitted means true. |
| `order` | No | Signed integer elapsed minutes for temporal analysis. |

BAMs must match the provided reference genome.
One-shot samples require at least one BAM file at startup. Live directories may initially be empty, but the directory must exist.

Including the `order` column automatically enables time-course views for active gene-set enrichment and independent component analysis (ICA).
For each analysis, assign every selected sample an `order` value. Each selected group must have a single timepoint across all its replicates,
Each timepoint must be unique across the selected groups, and at least two timepoints are required.

### 3. Configure and run

Create a `run.config` file. This example runs QC, quantification,
biotypes, and differential expression analysis once on the available BAM files:

```groovy
params {
    sample_sheet = "/path/to/samples.csv"
    reference_genome = "/path/to/reference.fa"
    reference_annotation = "/path/to/annotation.gtf"
    strand_specificity = "forward"
    out_dir = "/path/to/results"
    live_analysis = false
    differential_expression = true
    gene_set_enrichment = false
    ica_analysis = false
}
```

```bash
nextflow run . -profile standard -c run.config -w /path/to/work
```

Use `-profile singularity` for Singularity or `-profile apptainer` for Apptainer.

For NCBI prokaryotic GTFs missing transcript/exon records, convert the annotation before running:

```bash
bin/oarfish-gtf-convert genomic.gtf genomic.oarfish.gtf
```

The converted file should be used as `reference_annotation`. See the
[annotation guide](docs/quantification.md#prokaryotic-annotations).

### 4. View results or run live

`qc_report.html` can be opened in the results directory. For live analysis, set `live_analysis = true` in `run.config` before launching.
The report updates automatically as BAM file segments are received.

In `disabled` or `log` monitoring mode, live sample analysis automatically finish when a validated MinKNOW `final_summary.txt` or `final_summary_*.txt` is present and the summary and BAM directory have remained unchanged for 30 seconds. This also works when analysis starts after sequencing has finished.
Standard dated run directories, `bam_pass`/`bam_fail`, and barcode directories under those folders are automatically recognized.
For other layouts, add the optional `minknow_run_dir` sample-sheet column naming the dedicated directory for that run.
See [completion detection and its limitations](docs/minknow_completion.md).

For a manual completion signal, create a `STOP` file in the BAM directory (automatically in `terminate` monitoring mode):

```bash
touch /data/bams/treated_1/STOP
```

The workflow finalizes any pending data and saves the final report as a self-contained HTML file.

## Analysis modes

All modes include QC. Supply the FASTA and GTF/GFF3 reference files for any of the non-QC-only modes.

| Mode | Additional inputs | `differential_expression` | `gene_set_enrichment` | `ica_analysis` |
| --- | --- | --- | --- | --- |
| QC only | None | false | false | false |
| Quantification and biotypes | Both references | false | false | false |
| Differential expression | Both references | true | false | false |
| Gene-set analysis | Both references, `gene_sets` GMT | true | true | false |
| iModulon analysis | Both references, `ica_matrix` CSV/TSV | false | false | true |
| Combined analysis | Both references, GMT, ICA matrix | true | true | true |

For iModulon analysis, seqLM projects expression data onto a supplied, compatible gene-by-component weight matrix. The workflow does not train a new ICA model.
`ica_imodulon_table = "/path/to/iM_table.csv"` can be specified to provide names and functions from the same model, with one annotation row per component.
The **Independent Components** tab shows **Time course** when time metadata is present, otherwise **Component details**.
Diagnostics contains the component annotations. See the
[iModulon guide](docs/imodulon_analysis.md) for model formats and interpretation.

## Outputs and restarting

| Path under `out_dir` | Contents |
| --- | --- |
| `qc_report.html` | Main report. Live updates during sequencing |
| `qc_report_state.json`, `qc_report_snapshot_revision_*.html` | Live-report state and immutable snapshots. |
| `<group>/<alias>/` | Sample QC, quantification, and `FINAL` completion marker. |
| `differential_expression/` | Differential-expression and gene-set results. |
| `ica/` | iModulon snapshots; `latest.json` identifies the newest snapshot and status. |
| `termination/<method>/` | Precision monitoring state and per-analysis audits for `dea` or `ica`. |
| `execution/` | Nextflow report, timeline, and trace. |

`out_dir` may be reused to restore finalized samples and extend an experiment. Changes to inputs or references prevent checkpoint reuse.
This is separate from Nextflow `-resume` functionality.

## EPI2ME Desktop (Linux/macOS)

EPI2ME should be fully closed, the `seq-lm` environment activated if applicable, and the application pointed to the installed Nextflow:

```bash
# Linux
LABS_NXF_PATH="$(command -v nextflow)" /usr/lib/epi2me/EPI2ME

# macOS
LABS_NXF_PATH="$(command -v nextflow)" /Applications/EPI2ME.app/Contents/MacOS/EPI2ME
```

In **Workflows → Import workflow**, enter `https://github.com/rnabioinfo/seq_lm`, then select **Run this workflow** and provide the required inputs.

## Sequencing termination

Set `monitoring_behavior = "log"` to record stop eligibility but not terminate the acquisition, or `"terminate"` to stop them automatically.
Both methods forecast the expected enhancement in precision over the next 30 minutes by default.
Choose `stopping_method = "dea"` or `"ica"` to base the stopping on either of these analyses.
DEA monitoring requires enabled differential expression. ICA monitoring requires enabled ICA analysis, its model, and at least one condition group besides the ICA control.
Automatic termination also requires live analysis and MinKNOW credentials.
See [precision termination](docs/dea_termination.md) and [ICA monitoring](docs/imodulon_analysis.md#sequencing-termination).

With live monitoring enabled, stable new BAMs trigger quantification for each sample independently.
All live samples sharing a MinKNOW acquisition must qualify before the workflow can stop it.
The report shows the corresponding **DEA Termination** or **ICA Termination** tab.

Generate credentials in the activated `seq-lm` environment.
This step must be performed on the computer that is conducting the sequencing.
Manual certificate setup is also possible.

```bash
bin/seq-run-manager cert --minknow-client-certs-directory /path/to/minknow/conf/rpc-client-certs
```

Add these settings inside the existing `params` block in `run.config`:

```groovy
live_analysis = true
monitoring_behavior = "terminate"
minknow_client_certificate = "/path/to/minknow_cert.pem"
minknow_client_private_key = "/path/to/minknow_key.pem"
minknow_ca_certificate = "/path/to/minknow_cert.crt"
```

Each BAM directory's parent must contain one matching MinKNOW `*sample_sheet*.csv`.
The workflow creates `STOP` only after MinKNOW confirms termination. 
See the [precision termination guide](docs/dea_termination.md)
for uncertainty calculation, interruption guards, outputs, and container requirements.

## Analysis and live-run parameters

Values belong in the config's `params` block. Defaults can be found in `nextflow.config`.

### Inputs and references

| Parameter | Default | Description |
| --- | --- | --- |
| `sample_sheet` | Unset | CSV with sample names, groups, and BAM directories. |
| `reference_genome` | Unset | Reference FASTA; required with the annotation. |
| `reference_annotation` | Unset | Transcript GTF/GFF3; required with the FASTA. |

### Live BAM input and run completion

| Parameter | Default | Description |
| --- | --- | --- |
| `live_analysis` | `true` | Watch samples marked live in the sample sheet. |
| `bam_poll_interval_seconds` | `5` | Seconds between BAM-directory scans. |
| `bam_stability_polls` | `3` | Unchanged observations required before accepting a BAM. |
| `detect_minknow_completion` | `true` | Detect final MinKNOW summaries for live samples in disabled/log monitoring mode. |
| `minknow_completion_settle_seconds` | `30` | Positive integer unchanged file period for the final summary and BAM inventory before closing a sample. |

### Shared DEA and ICA readiness

| Parameter | Default | Description |
| --- | --- | --- |
| `min_read_count` | `10000` | Assigned-read threshold for counting ready replicates in DEA and ICA. |
| `min_replicate_sample_count` | `2` | Minimum libraries per group that must meet `min_read_count` for DEA and ICA. Every included library must have positive depth. Libraries with positive depth below the threshold still enter the analysis. |

### Differential expression (DEA)

| Parameter | Default | Description |
| --- | --- | --- |
| `differential_expression` | `true` | Enable edgeR; requires both references. |
| `gene_contrast_file` | Unset | Headerless CSV with one `target_group,reference_group` pair per row. DEA compares the target group against the reference group. |
| `de_lfc_cutoff` | `1.0` | Absolute log₂ fold-change threshold for edgeR glmTreat. |
| `de_padj_cutoff` | `0.05` | Maximum adjusted p-value for differential-expression calls. |

### Gene-set enrichment

| Parameter | Default | Description |
| --- | --- | --- |
| `gene_set_enrichment` | `true` | Enable fry and GSVA; requires edgeR and GMT gene sets. |
| `gene_sets` | Unset | GMT gene sets for fry and GSVA. |

### iModulon analysis (ICA)

| Parameter | Default | Description |
| --- | --- | --- |
| `ica_analysis` | `false` | Enable iModulon projection; requires a matrix and both references. |
| `ica_matrix` | Unset | CSV/TSV gene-by-component weight matrix. |
| `ica_gene_map` | Unset | Optional CSV/TSV with gene_id and model_gene_id columns. |
| `ica_min_gene_coverage` | `1.0` | Minimum fraction of model genes covered by annotation targets. |
| `ica_imodulon_table` | Unset | Optional component names and functions for the exact model. |
| `ica_control_group` | `"control"` | ICA reference group; requires at least two samples. |
| `ica_condition_groups` | Unset | Comma-separated ICA conditions; control is included automatically. Unset selects all samples. |
| `ica_log_base` | `2.0` | Expression logarithm base. Must exceed one. |
| `ica_pseudocount` | `1.0` | Positive value added before logging gene abundance. |
| `ica_padj_cutoff` | `0.05` | Maximum BH-adjusted p-value for component activity tests. |

### Sequencing termination parameters

#### Shared stopping settings

| Parameter | Default | Description |
| --- | --- | --- |
| `monitoring_behavior` | `"disabled"` | Monitoring action: disabled, log, or terminate. |
| `stopping_method` | `"dea"` | Source of contrast uncertainty intervals for the precision forecast, either dea or ica. |
| `stopping_forecast_minutes` | `30` | Forecast horizon. Evaluation intervals must be no longer than this for automatic stopping. |
| `stopping_gain_cutoff` | `0.01` | The forecast and the measured reduction in aggregate U must both be below this fraction. 0.01 means 1%. |
| `stopping_information_radius_scale` | `1.0` | Declared aggregation scale in log2 units. |
| `stopping_min_points` / `stopping_min_seconds` | `6` / `3600` | Minimum informative evaluations and elapsed time. |
| `stopping_confirmations` | `3` | Consecutive qualifying stopping checkpoints. |
| `stopping_max_processing_lag_seconds` | `600` | Maximum processing lag, independent of file-output interval. |
| `stopping_interruption_fraction` / `stopping_recovery_fraction` | `0.25` / `0.5` | Drop and recovery fractions relative to healthy production. |
| `stopping_recovery_intervals` | `2` | Consecutive productive arrival intervals, also spanning at least one third of the horizon. |

#### Method-specific stopping inputs

| Parameter | Default | Description |
| --- | --- | --- |
| `dea_stopping_genes` | Unset | Optional GMT file that limits the DEA stopping score to the features. Omit to use all features. |
| `ica_stopping_components` | Unset | Optional file with one exact model component ID per line that limits the ICA stopping score. Omit to use all components. |

#### MinKNOW connection

| Parameter | Default | Description |
| --- | --- | --- |
| `minknow_host` | `"host.docker.internal"` | MinKNOW manager host as seen from the selected container runtime. |
| `minknow_port` | `9501` | MinKNOW manager port. |
| `minknow_client_certificate` | Unset | PEM client certificate; required for termination. |
| `minknow_client_private_key` | Unset | PEM client private key; required for termination. |
| `minknow_ca_certificate` | Unset | PEM root CA certificate; required for termination. |

### Outputs and logging

| Parameter | Default | Description |
| --- | --- | --- |
| `out_dir` | `"output"` | Published results directory. |
| `monochrome_logs` | `false` | Disable colored workflow log messages. |
| `disable_ping` | `false` | Disable workflow start/completion telemetry. |

### Help and validation

| Parameter | Default | Description |
| --- | --- | --- |
| `help` | `false` | Display workflow help and exit. |
| `version` | `false` | Display the workflow version and exit. |
| `show_hidden_params` | `false` | Include hidden parameters in help output. |
| `validate_params` | `true` | Validate parameters against the workflow schema. |

## Guides

- [QC and report](docs/quality_control.md)
- [Quantification and annotations](docs/quantification.md)
- [Differential expression, gene sets, and checkpoints](docs/differential_expression.md)
- [iModulon analysis](docs/imodulon_analysis.md)

## License

See [LICENSE](LICENSE).