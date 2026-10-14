# seq_lm

Live or one-shot analysis of aligned Oxford Nanopore transcriptome BAMs, with an
interactive HTML report. Includes read QC, Oarfish quantification, transcript
biotypes, edgeR differential expression, fry/GSVA gene-set analysis, and optional
iModulon projection and MinKNOW stop control.

[Quickstart](#quickstart) · [Analysis modes](#analysis-modes) ·
[Parameters](#analysis-and-live-run-parameters) · [Guides](#guides)

## Quickstart

Requires **Nextflow ≥26.04.3**, **Java 17**, and **Docker** or
**Singularity/Apptainer**. Nextflow fetches published images from Docker Hub.
Windows users must run inside **WSL2** with Docker Desktop's WSL integration;
EPI2ME Desktop on Windows is not supported.

### 1. Install

Install Docker separately, then use the supplied Conda environment for Nextflow,
Java, and host utilities:

```bash
git clone https://github.com/RNABioInfo/seq_lm.git
cd seq_lm
conda env create --file environment.yml
conda activate seq-lm
nextflow -version
```

### 2. Prepare samples

Save `samples.csv` with paths to your aligned BAM directories:

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
| `group` | Yes | Condition label. Without a gene contrast file, DEA requires at least two `control` samples. ICA uses its independent configured control. |
| `bam_dir` | Yes | Existing directory searched recursively for BAMs; accessible to the container. |
| `is_live` | No | Watch this sample when `live_analysis` is also true. Blank or omitted means true. |
| `order` | No | Signed integer elapsed minutes for temporal analysis. |

Headers are case-sensitive. BAMs must match the reference genome, and read IDs
must be unique across chunks of each sample. One-shot samples need at least one
BAM at startup; live directories may start empty.

For temporal analysis, add `order` to every row: each group must have one time,
each time one group, and at least two times are required. Samples represent
independent replicates. This enables ICA time-course views automatically;
gene-set time-course views additionally require `timeline_analysis = true`.

### 3. Configure and run

Save `run.config`, replacing the paths. This example runs QC, quantification,
biotypes, and differential expression once on the available BAMs:

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
    timeline_analysis = false
    ica_analysis = false
}
```

```bash
nextflow run . -profile standard -c run.config -w /path/to/work
```

Choose `strand_specificity = "reverse"` for reads in the opposite transcript
orientation or `"unstranded"` to accept both orientations, such as unstranded
cDNA libraries. The default `"forward"` preserves direct-RNA behavior. This
setting applies to all samples; see [strand specificity](docs/quantification.md#strand-specificity)
for CLI usage and checkpoint compatibility.

Use `-profile singularity` for Singularity/Apptainer. Set parameters in the config:
quote strings, but leave booleans and numbers unquoted. Explicit boolean/numeric
CLI values reach this workflow as strings and fail schema validation.

For NCBI prokaryotic GTFs missing transcript/exon records, convert the annotation
before running:

```bash
bin/oarfish-gtf-convert genomic.gtf genomic.oarfish.gtf
```

Use the converted file as `reference_annotation`; see the
[annotation guide](docs/quantification.md#prokaryotic-annotations).

### 4. View results or run live

Open `qc_report.html` in your results directory. For live analysis, set
`live_analysis = true` in `run.config` before launching. The report updates as
stable BAM chunks arrive.

In `disabled` or `log` monitoring mode, live samples automatically finish when a
validated MinKNOW `final_summary.txt` or `final_summary_*.txt` is present and the
summary and BAM inventory have remained unchanged for 30 seconds, with no pending
unstable BAMs. This also works when analysis starts after sequencing has finished.
Standard dated run directories, `bam_pass`/`bam_fail`, and barcode directories
beneath those folders are recognized. For other layouts, add the optional
`minknow_run_dir` sample-sheet column naming the dedicated directory for that run.
See [completion detection and its limitations](docs/minknow_completion.md).

For a manual completion signal, create a `STOP` file in the BAM directory
(automatically in `terminate` monitoring mode):

```bash
touch /data/bams/treated_1/STOP
```

Repeat for every live sample. The workflow drains pending data and saves the
final report as a self-contained HTML file.

## Analysis modes

Edit the switches below in `run.config`. All modes include QC. Supply the FASTA
and GTF/GFF3 together; remove both reference settings for QC-only analysis.

| Mode | Additional inputs | `differential_expression` | `gene_set_enrichment` | `ica_analysis` |
| --- | --- | --- | --- | --- |
| QC only | None | false | false | false |
| Quantification and biotypes | Both references | false | false | false |
| Differential expression | Both references | true | false | false |
| Gene-set analysis | Both references, `gene_sets` GMT | true | true | false |
| iModulon analysis | Both references, `ica_matrix` CSV/TSV | false | false | true |
| Combined analysis | Both references, GMT, ICA matrix | true | true | true |

Keep `timeline_analysis = false` unless gene-set enrichment is enabled and the
sample sheet contains complete time metadata.

For iModulons, seqLM projects expression onto a supplied, compatible
gene-by-component weight matrix; it does not train a new ICA model. The repository
implementation uses NumPy/SciPy/pandas rather than PyModulon and tests
target-minus-control mean activities with Welch tests and BH correction, not
PyModulon DiMA. Add
`ica_imodulon_table = "/path/to/iM_table.csv"` for names and functions from the
**same model**, with one annotation row per component. The report shows
**Time course details** when time metadata is present, otherwise **Component details**.
Diagnostics contains the component annotations. See the
[iModulon guide](docs/imodulon_analysis.md) for model formats and interpretation.

## Outputs and restarting

| Path under `out_dir` | Contents |
| --- | --- |
| `qc_report.html` | Main report; live updates during sequencing, self-contained on completion. |
| `qc_report_state.json`, `qc_report_snapshot_revision_*.html` | Live-report state and immutable snapshots. |
| `<group>/<alias>/` | Sample QC, quantification, and `FINAL` completion marker. |
| `differential_expression/` | Differential-expression and gene-set results. |
| `ica/` | iModulon snapshots; `latest.json` identifies the newest snapshot and readiness. |
| `stability/` | Stability decisions and stop-control records. |
| `execution/` | Nextflow report, timeline, and trace. |

Reuse `out_dir` to restore finalized samples and extend an experiment. Changed
inputs or references prevent checkpoint reuse. This is separate from Nextflow
`-resume`; use a new output directory for a fresh analysis. Pin a workflow
revision with `-r` when launching `RNABioInfo/seq_lm` directly.

## EPI2ME Desktop (Linux/macOS)

Fully quit EPI2ME, activate `seq-lm`, and point it to the installed Nextflow:

```bash
# Linux
LABS_NXF_PATH="$(command -v nextflow)" /usr/lib/epi2me/EPI2ME

# macOS
LABS_NXF_PATH="$(command -v nextflow)" /Applications/EPI2ME.app/Contents/MacOS/EPI2ME
```

Keep the terminal open. In **Workflows → Import workflow**, enter
`https://github.com/rnabioinfo/seq_lm`, then select **Run this workflow** and fill
in the inputs. Repeat the launch command after restarting EPI2ME.

## Stability monitoring and MinKNOW termination

Set `monitoring_behavior = "log"` to record stop eligibility, or `"terminate"`
to stop eligible acquisitions. Select `stopping_method = "dea"` (the default)
for differential-expression stability, or `"ica"` for iModulon convergence.
ICA monitoring requires enabled ICA analysis and its model, and checks every
component for every required sample. Termination also requires live analysis and
MinKNOW credentials. See [ICA monitoring](docs/imodulon_analysis.md#sequencing-termination)
for the fold-change tolerance and informative checkpoint window.

Generate credentials in the activated `seq-lm` environment:

```bash
bin/seq-run-manager cert \
    --minknow-client-certs-directory /path/to/minknow/conf/rpc-client-certs
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
The workflow creates `STOP` only after MinKNOW confirms termination. See the
[monitoring guide](docs/differential_expression.md#differential-expression-stability)
for sample matching, stability criteria, and connection setup.

## Analysis and live-run parameters

Values belong in the config's `params` block. Defaults are from `nextflow.config`.
`-c`, `-profile`, `-w`, and `-resume` are Nextflow launcher options.

**Inputs and references**

| Parameter | Default | Description |
| --- | --- | --- |
| `sample_sheet` | Unset | CSV with sample names, groups, and BAM directories. |
| `reference_genome` | Unset | Reference FASTA; required with the annotation. |
| `reference_annotation` | Unset | Transcript GTF/GFF3; required with the FASTA. |
| `gene_contrast_file` | Unset | Headerless CSV of target,reference pairs; selects DEA comparisons and groups. |
| `gene_sets` | Unset | GMT gene sets for fry and GSVA. |

**Live acquisition**

| Parameter | Default | Description |
| --- | --- | --- |
| `live_analysis` | `true` | Watch samples marked live in the sample sheet. |
| `bam_poll_interval_seconds` | `5` | Seconds between BAM-directory scans. |
| `bam_stability_polls` | `3` | Unchanged observations required before accepting a BAM. |
| `detect_minknow_completion` | `true` | Detect final MinKNOW summaries for live samples in disabled/log monitoring mode. |
| `minknow_completion_settle_seconds` | `30` | Positive integer quiet period for the final summary and BAM inventory before closing a sample. |

**Differential expression and gene sets**

| Parameter | Default | Description |
| --- | --- | --- |
| `differential_expression` | `true` | Enable edgeR; requires both references. |
| `gene_set_enrichment` | `true` | Enable fry and GSVA; requires edgeR and GMT gene sets. |
| `timeline_analysis` | `false` | Enable temporal gene-set views; requires enrichment and sample order. |
| `min_read_count` | `10000` | Minimum assigned reads per sample for edgeR readiness. |
| `min_replicate_sample_count` | `2` | Minimum ready samples per group. |
| `de_lfc_cutoff` | `1.0` | Absolute log₂ fold-change threshold for edgeR glmTreat. |
| `de_padj_cutoff` | `0.05` | Maximum adjusted p-value for differential-expression calls. |

**iModulon analysis**

| Parameter | Default | Description |
| --- | --- | --- |
| `ica_analysis` | `false` | Enable iModulon projection; requires a matrix and both references. |
| `ica_control_group` | `"control"` | ICA reference group; requires at least two samples. |
| `ica_condition_groups` | Unset | Comma-separated ICA conditions; control is included automatically. Unset selects all samples. |
| `ica_matrix` | Unset | CSV/TSV gene-by-component weight matrix. |
| `ica_imodulon_table` | Unset | Optional component names and functions for the exact model. |
| `ica_gene_map` | Unset | Optional CSV/TSV with gene_id and model_gene_id columns. |
| `ica_log_base` | `2.0` | Expression logarithm base; must exceed one. |
| `ica_pseudocount` | `1.0` | Positive value added before logging gene abundance. |
| `ica_min_gene_coverage` | `1.0` | Minimum fraction of model genes covered by annotation targets. |
| `ica_min_read_count` | `10000` | Minimum assigned abundance in every sample for ICA readiness. |
| `ica_padj_cutoff` | `0.05` | Maximum BH-adjusted p-value for component activity tests. |

**Stability monitoring**

| Parameter | Default | Description |
| --- | --- | --- |
| `monitoring_behavior` | `"disabled"` | Stability action: disabled, log, or terminate. |
| `stopping_method` | `"dea"` | Analysis controlling stopping: dea or ica. |
| `ica_max_fold_change` | `1.05` | Largest allowed component-wise reconstructed expression change across the ICA checkpoint window. |
| `ica_min_read_growth_fraction` | `0.1` | Minimum assigned-abundance increase for a new ICA checkpoint; 0.1 means 10%. |
| `num_stable_batches` | `3` | DEA: consecutive stable comparisons. ICA: informative updates in the window, requiring four checkpoints at the default. |
| `stability_max_feature_diff_fraction` | `0.05` | Maximum fraction of filtered features added or removed. |
| `stability_max_median_abs_lfc_delta` | `0.05` | Maximum median absolute change in log₂ fold change. |
| `stability_min_jaccard_similarity` | `0.95` | Minimum overlap of sufficiently large DE-call sets. |
| `stability_min_de_calls_for_fraction_metrics` | `20` | DE-call union size at which fractional metrics apply. |
| `stability_max_small_set_call_changes` | `2` | Maximum DE-call changes for smaller call sets. |

**MinKNOW connection**

| Parameter | Default | Description |
| --- | --- | --- |
| `minknow_host` | `"host.docker.internal"` | MinKNOW manager host. |
| `minknow_port` | `9501` | MinKNOW manager port. |
| `minknow_client_certificate` | Unset | PEM client certificate; required for termination. |
| `minknow_client_private_key` | Unset | PEM client private key; required for termination. |
| `minknow_ca_certificate` | Unset | PEM root CA certificate; required for termination. |

**Outputs and logging**

| Parameter | Default | Description |
| --- | --- | --- |
| `out_dir` | `"output"` | Published results directory. |
| `monochrome_logs` | `false` | Disable colored workflow log messages. |
| `disable_ping` | `false` | Disable workflow start/completion telemetry. |

**Help and validation**

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
