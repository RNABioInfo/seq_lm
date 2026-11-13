## Quickstart

> **_Author Note:_**  This documentation was generated using OpenAI GPT-6.1 and manually checked for correctness.

### 1. Installation

seq_lm requires Nextflow 26.04.3 or later, Java 17, and a container runtime
(Docker, Singularity, or Apptainer). On Windows, use WSL2 with Docker
Desktop's WSL integration enabled.

If Nextflow is unavailable, create and activate the supplied Conda environment
from the repository directory.

```bash
conda env create --file environment.yml
conda activate seq-lm
```

### 2. Sample sheet

Save one row per sample in `samples.csv`. The `alias` must be unique within
each `group`. Set `bam_dir` to the directory containing the sample's aligned
BAM files. seq_lm searches this directory and its subdirectories.

```csv
alias,group,bam_dir
control_1,control,/data/control_1
control_2,control,/data/control_2
treated_1,treated,/data/treated_1
treated_2,treated,/data/treated_2
```

### 3. Quality control for existing data

Save the following settings in `run.config`. This configuration analyzes
existing BAM files once and produces a quality control report. At least one
BAM file per sample must be available at startup.

```groovy
params {
    sample_sheet = "samples.csv"
    out_dir = "output"
    live_analysis = false
    differential_expression = false
    gene_set_enrichment = false
}
```

```bash
nextflow run . -profile standard -c run.config
```

Use the `standard` profile for Docker, `singularity` for Singularity, or
`apptainer` for Apptainer. When the run completes, open
`output/report.html`.
[Guide to the report](quality_control.md)

### Additional analysis options

- Set `live_analysis = true` to analyze BAM files as they become available.
  A live sample directory may be empty at startup. In the sample sheet, set
  `is_live` to `false` for completed samples. Blank or omitted values default
  to `true`. [Live run completion](minknow_completion.md)
- Supply `reference_genome` and `reference_annotation` to estimate transcript
  abundance. [Quantification guide](quantification.md)
- Enable differential expression, gene-set enrichment, or iModulon analysis
  according to the study design. [Differential expression and gene sets](differential_expression.md), [iModulon analysis](imodulon_analysis.md)
- Run `nextflow run . --help` to list the workflow parameters. [README and EPI2ME setup](../README.md)
