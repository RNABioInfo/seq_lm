## Quickstart

Run from the repository root with Nextflow 26.04 or later and Docker or
Singularity installed. The configured runtime fetches versioned container images.

Show the available parameters:

```bash
nextflow run . --help
```

Analyze a sample sheet containing `alias`, `group`, and `bam_dir`:

```bash
nextflow run . -profile standard --sample_sheet samples.csv --out_dir output
```

Use `-profile singularity` for Singularity. Live analysis watches eligible BAM
directories; use `--live_analysis false` to analyze a completed batch.

Supply both `--reference_genome` and `--reference_annotation` to enable transcript
quantification. Enable DEA with `--differential_expression true`, or ICA with
`--ica_analysis true --ica_matrix weights.csv`. See [the README](../README.md)
for analysis modes and sequencing termination settings.

The output directory contains `qc_report.html`, per-sample QC, and the enabled
quantification and analysis results. Finalized sample results include validated
`FINAL` checkpoints for restoration on subsequent runs.
