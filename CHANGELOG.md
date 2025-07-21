# nf-core/veba: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of nf-core/veba, created with the [nf-core](https://nf-co.re/) template.

### `Added`
- Added `fastqpreprocessor` subworkflow
  - Includes official nf-core `fastp` module for quality trimming and filtering
  - Integrates `fastqc` and `multiqc` for quality control and reporting
  - Outputs trimmed/merged/failed reads and generates QC reports

### `Fixed`

### `Dependencies`

### `Deprecated`
