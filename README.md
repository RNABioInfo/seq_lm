# seq_lm

`seq_lm` is a Nextflow workflow for live or one-shot analysis of aligned Oxford
Nanopore transcriptome sequencing data. It watches per-sample BAM directories,
updates quality-control and expression analyses as stable BAM files arrive, and
publishes a single interactive HTML report throughout the experiment.

The workflow supports:

- chunk-level read and alignment quality control with NanoGet and `samtools`;
- cumulative transcript quantification with Oarfish;
- annotation-based transcript-biotype composition;
- differential-expression analysis with edgeR;
- directional gene-set testing with edgeR `fry` and sample-level scoring with
  GSVA; and
- optional monitoring that records when results stabilize or asks MinKNOW to
  stop eligible acquisitions.

The input BAMs must already be aligned to the genome represented by the supplied
reference FASTA and annotation. Live results are provisional: expression and
enrichment statistics are recalculated as additional reads arrive.

## Analysis modes

Reference inputs and analysis switches determine which parts of the workflow
run. `--differential_expression` and `--gene_set_enrichment` both default to
`true`.

| Mode | Required inputs and options | Main results |
| --- | --- | --- |
| Quality control only | Sample sheet; set `--differential_expression=false --gene_set_enrichment=false` and omit both references | Per-chunk QC and the integrated report |
| QC, quantification, and biotypes | Sample sheet, reference FASTA, and GTF/GFF3; set `--differential_expression=false --gene_set_enrichment=false` | QC, cumulative Oarfish quantification, and transcript-biotype composition |
| Differential expression | Sample sheet, reference FASTA, and GTF/GFF3; set `--gene_set_enrichment=false` | QC, quantification, biotypes, and edgeR contrasts |
| Full analysis | Sample sheet, reference FASTA, GTF/GFF3, and GMT gene sets | QC, quantification, edgeR, fry, and GSVA |

`--reference_genome` and `--reference_annotation` must always be supplied
together. Gene-set enrichment requires differential expression and a
`--gene_sets` GMT file.

## Quickstart

### Requirements and platform support

- [Nextflow](https://www.nextflow.io/) **26.04.3 or newer** and Java 17.
- [Docker](https://www.docker.com/products/docker-desktop/) for the default
  `standard` profile, or
  [Singularity/Apptainer](https://docs.sylabs.io/guides/latest/user-guide/) for
  the `singularity` profile.
- Linux, macOS, or Windows through WSL2.

Nextflow 26 is required because the workflow uses Nextflow's static type system
and strict syntax. On Windows, this workflow must be launched from a Nextflow 26
installation inside **WSL2**; launching it from the EPI2ME Desktop UI on Windows
is not supported. Use Linux paths in WSL, for example `/mnt/c/data/experiment`,
and enable Docker Desktop's WSL integration when using Docker.

On Linux and macOS, the workflow can be launched either from the command line or
through EPI2ME Desktop. The repository's Conda environment installs the required
host-side Nextflow, Java, OpenSSL, MinKNOW API, and utility dependencies. The
analysis programs themselves remain isolated in versioned Docker or Singularity
containers.

### 1. Install the host environment

Clone the canonical repository and create the environment with Conda or Mamba:

```bash
git clone https://github.com/RNABioInfo/seq_lm.git
cd seq_lm
conda env create --file environment.yml
conda activate seq-lm
```

Update an existing environment after `environment.yml` changes with:

```bash
conda env update --file environment.yml --prune
conda activate seq-lm
```

Confirm the runtime and container engine before starting an analysis:

```bash
nextflow -version
docker run --rm hello-world
nextflow run . --help
```

The reported Nextflow version must be at least `26.04.3`. If using
Singularity/Apptainer, verify that runtime instead of Docker and select
`-profile singularity` in the commands below.

### 2. Use Nextflow 26 in EPI2ME Desktop on Linux or macOS

EPI2ME Desktop may bundle an older Nextflow release. Fully quit EPI2ME, activate
the `seq-lm` environment, and launch the application with `LABS_NXF_PATH`
pointing to the Nextflow 26 executable.

Linux:

```bash
conda activate seq-lm
LABS_NXF_PATH="$(command -v nextflow)" /usr/lib/epi2me/EPI2ME
```

macOS:

```bash
conda activate seq-lm
LABS_NXF_PATH="$(command -v nextflow)" /Applications/EPI2ME.app/Contents/MacOS/EPI2ME
```

Keep the terminal open while EPI2ME is running. Repeat this launch method after
restarting the application so that it continues to use Nextflow 26.

In EPI2ME:

1. Open **Workflows**, choose **Import workflow**, and enter
   `https://github.com/RNABioInfo/seq_lm`.
2. Open `seq_lm` from **Installed workflows** and select **Run this workflow**.
3. Complete the required input, reference, analysis, and output fields.
4. Select **Launch workflow**, then monitor the run and open `qc_report.html`
   from the results view.

See the official
[EPI2ME workflow import and launch guide](https://epi2me.nanoporetech.com/epi2me-docs/quickstart/)
for general UI instructions.

### 3. Prepare the sample sheet

Pass a comma-separated sample sheet through `--sample_sheet`. Header names are
case-sensitive.

```csv
alias,group,bam_dir,is_live,order
control_1,control,/data/bams/control_1,false,0
control_2,control,/data/bams/control_2,false,0
treated_1,treated,/data/bams/treated_1,true,1
treated_2,treated,/data/bams/treated_2,true,1
```

| Field | Required | Meaning |
| --- | --- | --- |
| `alias` | Yes | Sample name used in reports and output paths. It must be unique within its group. |
| `group` | Yes | Experimental condition. At least two rows must use `control`, matched case-insensitively. Each other group is contrasted separately with the control group. |
| `bam_dir` | Yes | Existing directory searched recursively for `.bam` files. The path must be visible to the selected container runtime. |
| `is_live` | No | `true`, `false`, or blank, case-insensitively. Blank or missing values default to `true`. A row is watched only when this value and `--live_analysis` are both true. |
| `order` | No | Signed integer elapsed time in minutes. Every row must provide it when `--timeline_analysis` is enabled. With active ICA, including this column automatically enables ICA time-course analysis. In temporal mode, each group is one independent-replicate time point: a group must have one order value and an order value must identify one group. |

Effectively live samples may start with an empty `bam_dir`. Samples that are not
live are processed once and must contain at least one BAM when the workflow
starts. A new live BAM is accepted only after its size and modification time are
unchanged for `--bam_stability_polls` consecutive scans; the default is three
polls at five-second intervals.

Read identifiers must be globally unique across BAM chunks belonging to the
same sample. Standard Nanopore UUID read identifiers normally satisfy this
requirement.

### 4. Run the workflow

The following commands use Docker through the default `standard` profile and
place Nextflow intermediates and published results in explicit directories.

#### Full live analysis

```bash
nextflow run . \
    -profile standard \
    -w /path/to/seq_lm_work \
    --sample_sheet /path/to/samples.csv \
    --reference_genome /path/to/reference.fa \
    --reference_annotation /path/to/annotation.gtf \
    --gene_sets /path/to/pathways.gmt \
    --out_dir /path/to/seq_lm_output
```

#### One-shot full analysis

Process all BAMs present at startup and exit without watching for new files:

```bash
nextflow run . \
    -profile standard \
    -w /path/to/seq_lm_work \
    --sample_sheet /path/to/samples.csv \
    --live_analysis=false \
    --reference_genome /path/to/reference.fa \
    --reference_annotation /path/to/annotation.gtf \
    --gene_sets /path/to/pathways.gmt \
    --out_dir /path/to/seq_lm_output
```

#### Quality control only

```bash
nextflow run . \
    -profile standard \
    -w /path/to/seq_lm_qc_work \
    --sample_sheet /path/to/samples.csv \
    --differential_expression=false \
    --gene_set_enrichment=false \
    --out_dir /path/to/seq_lm_qc_output
```

#### Quantification and transcript-biotype QC

```bash
nextflow run . \
    -profile standard \
    -w /path/to/seq_lm_quant_work \
    --sample_sheet /path/to/samples.csv \
    --reference_genome /path/to/reference.fa \
    --reference_annotation /path/to/annotation.gtf \
    --differential_expression=false \
    --gene_set_enrichment=false \
    --out_dir /path/to/seq_lm_quant_output
```

#### Differential expression without gene-set enrichment

```bash
nextflow run . \
    -profile standard \
    -w /path/to/seq_lm_de_work \
    --sample_sheet /path/to/samples.csv \
    --reference_genome /path/to/reference.fa \
    --reference_annotation /path/to/annotation.gtf \
    --gene_set_enrichment=false \
    --out_dir /path/to/seq_lm_de_output
```

The GitHub-hosted workflow can also be launched without a local clone by
replacing `.` with `RNABioInfo/seq_lm`. Pin a release or revision with
Nextflow's `-r` option when reproducibility across future workflow updates is
required.

### 5. Finish a live analysis

When no more BAMs will arrive for a live sample, create a `STOP` marker inside
that sample's `bam_dir`:

```bash
touch /data/bams/treated_1/STOP
```

The workflow drains any pending stable BAMs, finalizes the sample, and removes
it from later synchronization barriers. Create a marker for every remaining
live sample so the workflow can complete normally and replace the live report
shell with the final self-contained report.

## Stability monitoring and MinKNOW termination

`--monitoring_behavior` controls actions based on differential-expression
stability:

- `disabled` performs no stability action and is the default;
- `log` records when each sample would become eligible to stop; and
- `terminate` sends a stop request to the corresponding MinKNOW acquisition and
  creates the local `STOP` marker only after MinKNOW confirms the request.

Stability is assessed independently for every non-control-versus-control
contrast. With the default `--num_stable_batches 3`, the first successful edgeR
snapshot establishes a baseline, so at least four successful snapshots are
needed before a sample can become eligible.

Termination requires a client certificate, private key, and MinKNOW root CA.
After activating the host environment, generate credentials with:

```bash
seq-run-manager cert \
    --minknow-client-certs-directory /path/to/minknow/conf/rpc-client-certs
```

Then provide the generated files when launching the workflow:

```bash
nextflow run . \
    -profile standard \
    -w /path/to/seq_lm_work \
    --sample_sheet /path/to/samples.csv \
    --reference_genome /path/to/reference.fa \
    --reference_annotation /path/to/annotation.gtf \
    --gene_sets /path/to/pathways.gmt \
    --monitoring_behavior terminate \
    --minknow_host host.docker.internal \
    --minknow_port 9501 \
    --minknow_client_certificate ~/.config/seq-run-manager/minknow/minknow_cert.pem \
    --minknow_client_private_key ~/.config/seq-run-manager/minknow/minknow_key.pem \
    --minknow_ca_certificate ~/.config/seq-run-manager/minknow/minknow_cert.crt \
    --out_dir /path/to/seq_lm_output
```

For termination, the direct parent of each `bam_dir` should contain exactly one
file matching `*sample_sheet*.csv`. The workflow matches its `alias` to the
MinKNOW sheet's `sample_id` to discover the `protocol_run_id`. Missing or
ambiguous metadata disables termination for that sample without stopping the
analysis. Failed stop requests are warned about and retried after the next
stable batch.

See [Differential Expression Workflow](docs/differential_expression.md#differential-expression-stability)
for stability metrics, audit fields, and certificate behavior.

## Outputs and continuing an experiment

The main published results under `--out_dir` include:

```text
qc_report.html
qc_report_state.json
qc_report_snapshot_<batch_index>.html
<group>/<alias>/
  FINAL
  quantification/
  qc/
differential_expression/
  batch_<batch_index>/
  latest/
stability/
  batch_<analysis_index>/
execution/
  report.html
  timeline.html
  trace.txt
```

`qc_report.html` is the stable EPI2ME entry point. During live analysis it loads
the newest complete immutable snapshot without reloading the outer page. After
successful completion it is atomically replaced by the latest self-contained
snapshot. Depending on the selected analysis mode, the report contains Quality
Control, Transcript biotypes, Differential Analysis, Result Stability, Gene Set
Enrichment, and Temporal Analysis views. Temporal Analysis is enabled with
`--timeline_analysis`, requires gene-set enrichment, and summarizes a single
trajectory over the sample-sheet `order` values. For a selected gene set, the
report shows its raw GSVA-score trajectory and a gene-by-time heatmap. Heatmap
cells are gene-wise z-scores of the replicate-mean
`log2(TMM-normalized CPM + 1)` values, and average-linkage clustering groups
genes with similar temporal profiles. Hovering retains each cell's absolute
mean logCPM, SD, group, and replicate count. This view is descriptive: row
standardization removes between-gene abundance differences, and clustering
does not establish co-regulation or a statistically significant time effect.

When a sample finishes, `seq_lm` persists its quantification and raw QC inputs
and writes `FINAL` last. A later invocation using the same `--out_dir` validates
and restores finalized samples, allowing new samples to extend an experiment
without rerunning completed sample-level work. This mechanism is independent of
Nextflow `-resume`. Changed BAMs, references, manifests, or derived artifacts
cause validation to stop rather than silently reuse stale results. To recompute
a finalized sample, remove its complete `<out_dir>/<group>/<alias>/` directory
before launching a new run.

See the detailed documentation for output contracts and interpretation:

- [Analysis report and quality control](docs/quality_control.md)
- [Quantification and transcript-biotype QC](docs/quantification.md)
- [Differential expression, fry, GSVA, stability, and checkpoints](docs/differential_expression.md)
- [Fixed-matrix iModulon projection and differential activity](docs/imodulon_analysis.md)

## Prokaryotic annotations

NCBI prokaryotic GTF files often omit the transcript and exon records Oarfish
needs for protein-coding targets. Convert such an annotation before using it:

```bash
bin/oarfish-gtf-convert genomic.gtf genomic.oarfish.gtf
```

The converter preserves declared RNA models and creates one gene-sized,
single-exon transcript for each eligible protein-coding gene. It does not infer
operons, transcript boundaries, or untranslated regions. Prefer an
experimentally curated transcript annotation when those units matter.

## Interpretation

Differential-expression, fry, and GSVA results describe associations in the
tested organism, conditions, and sampling design. They do not establish causal
drivers, biochemical pathway activity, or biological mechanism. Results from
live batches may change as sequencing depth and statistical power increase;
use the final checkpoint for downstream interpretation.

Temporal figures are descriptive summaries of independent biological
replicates, not time-course significance tests. Their lines connect measured
minutes but do not estimate unobserved intermediate states. Time-associated
bulk-expression patterns may also reflect composition, batch, or other
variables confounded with sampling time.

## License and links

This project is distributed under the terms in [LICENSE](LICENSE).

- [Workflow repository](https://github.com/RNABioInfo/seq_lm)
- [Nextflow documentation](https://www.nextflow.io/docs/latest/)
- [EPI2ME Desktop documentation](https://epi2me.nanoporetech.com/epi2me-docs/)
- [Docker documentation](https://docs.docker.com/)
- [SingularityCE user guide](https://docs.sylabs.io/guides/latest/user-guide/)
