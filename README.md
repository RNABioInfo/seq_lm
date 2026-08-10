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




## Useful links

* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [singularity](https://docs.sylabs.io/guides/latest/user-guide/)
