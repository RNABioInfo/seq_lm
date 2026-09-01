# Workflow template

This repository contains a [nextflow](https://www.nextflow.io/) workflow
template that can be used as the basis for creating new workflows.

> This workflow is not intended to be used by end users.

## Introduction

This section of documentation typically contains an overview of the workflow in terms of motivation
and bioinformatics methods, listing any key tools or algorithms employed, whilst also describing its
range of use-cases and what a suitable input dataset should look like.

## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[Singularity](https://docs.sylabs.io/guides/latest/user-guide/) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either Docker or Singularity is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit out website](https://labs.epi2me.io/wfindex).

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-template --help
```

to see the options for the workflow.

### Sample sheet

Pass a comma-separated sample sheet with `--sample_sheet`. The header names are
case-sensitive. A sheet containing every supported field looks like this:

```csv
alias,group,bam_dir,is_live,order
control_1,control,/data/bams/control_1,false,0
control_2,control,/data/bams/control_2,false,0
treated_1,treated,/data/bams/treated_1,true,1
treated_2,treated,/data/bams/treated_2,true,1
```

| Field | Required | Accepted values | Meaning |
| --- | --- | --- | --- |
| `alias` | Yes | Non-empty text | Sample name used in reports and output paths. An alias must be unique within its `group`; the same alias may be used in different groups. |
| `group` | Yes | Non-empty text | Experimental condition used for comparisons. At least two rows must belong to the `control` group (matched case-insensitively). Each non-control group is tested separately against the control group during differential-expression analysis. |
| `bam_dir` | Yes | Path to an existing directory | Directory containing the sample's BAM files. BAMs ending in `.bam` are discovered recursively. The directory must exist when the workflow starts. |
| `is_live` | No | `true`, `false`, or blank (case-insensitive) | Controls whether the directory is eligible to be watched for new BAMs. A missing or blank value defaults to `true`. A sample is watched only when both this value and `--live_analysis` are true. |
| `order` | No | Integer | Timeline position for the sample. It may be omitted during normal analysis, but every row must provide it when `--timeline_analysis` is enabled. |

Samples that are watched may start with an empty `bam_dir`. Samples that are
not effectively live, either because `is_live` is `false` or because
`--live_analysis=false`, are processed once and must have at least one BAM at
startup. Live directories are rescanned every `--bam_poll_interval_seconds`
seconds. A new BAM is accepted only after its size and modification time remain
unchanged for `--bam_stability_polls` consecutive scans (three by default).
This polling approach also detects files written by Windows applications below
WSL-mounted paths such as `/mnt/c`. The sample sheet must contain at least one
row.

For example, run the workflow with:

```bash
nextflow run epi2me-labs/wf-template --sample_sheet samples.csv
```

**Workflow outputs**

The primary outputs of the workflow include:

* a simple text file providing a summary of sequencing reads,
* an HTML report document detailing the primary findings of the workflow.

See [Quality Control Workflow](docs/quality_control.md) for the live QC report.
See [Differential Expression Workflow](docs/differential_expression.md) for
cumulative BAM quantification and per-batch edgeR, fry, and GSVA refreshes.
Differential expression and gene-set enrichment can be disabled independently;
gene-set enrichment requires differential expression. QC-only runs do not
require reference files, and differential-only runs do not require a GMT file.
Optional DE-stability monitoring can log sample-level stopping decisions or
create per-sample `STOP` markers after every relevant contrast remains stable.




## Useful links

* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [singularity](https://docs.sylabs.io/guides/latest/user-guide/)
