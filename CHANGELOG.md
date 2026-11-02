# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

- **Breaking:** Use only the uncertainty-versus-depth curve for DEA termination
  forecasts and calibrate its own empirical error margin. Remove the gene-model
  read-depth projection and its fitted-mean/dispersion TSV exports. Publish
  `curve_gain` and curve-only monitoring state version 3; earlier state is rejected.

- Validate stopping settings only for the enabled DEA or ICA method; remove
  obsolete-name migration checks and skip inactive stopping configuration.

- Reuse the existing stable BAM trigger for post-DEA R checks; remove the separate
  timed watcher and global R batch barriers. Forecast from duration-weighted
  arrival rates, with configurable horizon and explicit cadence/lag guards.

- **Breaking:** Replace DEA result-stability stopping with 95% edgeR QL-profile radii and a
  causal configurable-horizon marginal-improvement forecast (30-minute default); default cutoff 1%.
- Use independent two-group fits, timestamped BAM-triggered snapshots, explicit
  pause/recovery protection, and current-format persistent monitoring state.
- Retain ICA stopping as an alternative and preserve GSVA/fry/ICA analysis.
- Use `stopping_method=dea` and the `dea_` parameter, method and package prefix
  for precision monitoring. Remove old DEA stability thresholds, executable and
  audit reader without aliases. Rename readiness checking to `check-de-readiness`.
- Add precision/model exports, monitoring plots and saved-checkpoint evaluation.
- Existing R/report containers are sufficient; no image publication is required.

- Run cumulative Oarfish quantification whenever paired genome and GTF/GFF3
  references are available, independently of differential expression, and add
  annotation-based transcript-biotype composition to Quality Control → Read
  flow for every matching report batch.
- Discover optional MinKNOW protocol run IDs beside sample BAM directories and
  stop stable acquisitions through `seq-run-manager`, with nonfatal retryable
  failures and conditional credential validation.
- Persist final per-sample Oarfish and QC artifacts with a validated `FINAL`
  manifest and restore them on later invocations without relying on `-resume`.
- Continue differential snapshot numbering across invocations and publish both
  immutable `batch_<N>` trees and the newest result under
  `<out_dir>/differential_expression/latest`.

## [v4.1.0]
### Added
- Optional `required_sample_types` field added to fastqingress. The sample sheet must contain at least one of each sample type provided to be deemed valid.
- Configuration for running demo data in AWS

### Changed
- Updated GitHub issue templates to force capture of more information.
- Removed glibc hack from post-test script
- Updated LICENSE to BSD-4-Clause
- Enum choices are enumerated in the `--help` output
- Enum choices are enumerated as part of the error message when a user has selected an invalid choice
- Bumped minimum required Nextflow version to 22.10.8

### Fixed
- Replaced `--threads` option in fastqingress with hardcoded values to remove warning about undefined `param.threads`

## [v4.0.0]
### Added
- GitHub issues template.
- Return of metadata with fastqingress.
- Check of number of samples and barcoded directories.
- Example of how to use the metadata from `fastqingress`.
- Implemented `--version`
- `fastcat_extra_args` option to `fastq_ingress` to pass arbitrary arguments to `fastcat` (defaults to empty string).
- `fastcat_stats` option to `fastq_ingress` to force generation of `fastcat` stats even when the input is only a single file (default is false).

### Changed
- Use `bgzip` for compression instead of `pigz`.
- pre-commit now uses `flake8` v5.0.4.
- Report is now created with [ezCharts](https://github.com/epi2me-labs/ezcharts).
- The workflow now also publishes the metadata emitted by `fastq_ingress` in `$params.out_dir`.
- `fastq_ingress` now returns `[metamap, path-to-fastcat-seqs | null, path-to-fastcat-stats | null]`.
- Bumped base container to v0.2.0.
- Use groovy script to ping after workflow has run.
- Removed sanitize fastq option.
- fastq_ingress now removes unclassified read folders by default.
- Workflow name and version is now more prominently displayed on start

### Fixed
- Output argument in Fastqingress homogenised.
- Sanitize fastq intermittent null object error.
- Add `*.pyc` and `*.pyo` ignores to wf-template .gitignore

### Note
- Bumped version to `v4` to align versioning with Launcher v4

## [v0.2.0]
### Added
- default process label parameter
- Added `params.wf.example_cmd` list to populate `--help`

### Changed
- Update WorkflowMain.groovy to provide better `--help`

## [v0.1.0]
### Changed
- `sample_name` to `sample_id` throughout to mathc MinKNOW samplesheet.

### Added
- Singularity profile include in base config.
- Numerous other changes that have been lost to the mists of time.

## [v0.0.7]
### Added
- Fastqingress module for common handling of (possibly multiplexed) inputs.
- Optimized container size through removal of various conda cruft.

### Changed
- Use mamba by default for building conda environments.
- Cut down README to items specific to workflow.

### Fixed
- Incorrect specification of conda environment file in Nextflow config.

## [v0.0.6]
### Changed
- Explicitely install into base conda env

## [v0.0.5]
### Added
- Software versioning report example.

## [v0.0.4]
### Changed
- Version bump to test CI.

## [v0.0.3]
### Changed
- Moved all CI to templates.
- Use canned aplanat report components.

## [v0.0.2]
### Added
- CI release checks.
- Create pre-releases in CI from dev branch.

## [v0.0.1]
* First release.
