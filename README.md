# Integrative Transcriptomics of Grapevine Bud Dormancy Release

Version 1.0.0 — public submission release.

This repository provides analysis code, parameter files, computational provenance documentation, selected frozen processed resources, and figure-source data supporting the manuscript:

“Integrative Transcriptomics Identifies Reproducible Coexpression Programs across Grapevine Genotypes during Bud Dormancy Release.”

## Overview

The repository contains the frozen materials selected for public submission. Scientific results are frozen. Exact replay is not claimed for historical stages whose stochastic or command-level provenance was not fully recovered.

## Data availability

Original public transcriptomic datasets are available from NCBI GEO/SRA under:

- GSE124820
- GSE273240
- GSE184114
- GSE277812
- GSE337039

Raw FASTQ files are not redistributed here. GSE337039 raw reads were reprocessed from the public archive as described in the manuscript.

Accession links, study roles, and metadata are provided in [metadata/dataset_manifest.tsv](metadata/dataset_manifest.tsv). The discovery workflow retains the documented QC subsets and frozen QC table. Distinct external tissue and experimental contexts remain identified in the provenance records.

No Zenodo archive is used for the initial submission. Additional archival deposition may be added later if required by the journal or reviewers.

## Repository structure

- `scripts/` contains base, legacy discovery, extended, final method-revision, figure, and utility implementations.
- `config/` contains the parameter catalog.
- `metadata/` contains design and mapping resources.
- `environment/` separates historical stage records from current observations.
- `figure_source_data/` contains the public figure-source tables and panel maps.
- `docs/` documents source roles, dependencies, workflow boundaries, and accepted reproducibility limitations.

## Authoritative analysis workflow

See [docs/AUTHORITATIVE_ANALYSIS_MAP.md](docs/AUTHORITATIVE_ANALYSIS_MAP.md). Only `scripts/method_revision_v2/15_whole_module_validation_v2.py` supplies the final B = 10,000 random-set/Stouffer inference, with P = (b + 1)/(10,000 + 1) and square-root effective-sample-size weights. The older whole-module script remains `HISTORICAL_REQUIRED` for contextual projection. Main WGCNA powers 12 and 16 are automatic-selection results; the preservation seed remains unknown.

## Reproducibility levels

[Level 1A](docs/REPRODUCIBILITY_LEVELS.md) covers public matrix acquisition and quality control. Level 1B starts from the retained frozen processed inputs and results. These levels are not presented as one fully demonstrated base-to-legacy conversion chain. Level 2 covers the separate GSE337039 raw-read branch.

## Public-data acquisition and frozen inputs

The four processed-matrix datasets feed the [parallel base and legacy branches](docs/MATRIX_PROCESSING_BRANCHES.md). Historical execution records confirm that the legacy parser was run; direct base-to-legacy conversion was not established. Frozen GO mappings, rather than current remote responses, are authoritative for the reported analysis.

Follow [docs/workflow.md](docs/workflow.md) for the documented workflow and retained-resource layout. The `materialize_data` utility only copies and verifies files; it does not start scientific analysis. No implicit one-command scientific workflow is provided.

## GSE337039 raw-read validation

See [docs/GSE337039_RAW_PROCESSING.md](docs/GSE337039_RAW_PROCESSING.md) for the documented processing tools and reference assembly: fastp 0.24.3, STAR 2.7.11b, featureCounts 2.0.8, GCF_000003745.3, and forward-stranded counts. Raw FASTQ, SRA, BAM, and SAM files and genome indexes are not included.

## Historical and final analysis scripts

[docs/HISTORICAL_VS_FINAL_SCRIPTS.md](docs/HISTORICAL_VS_FINAL_SCRIPTS.md) distinguishes retained implementation and context from final inference. Public-copy comments flag historical defaults. No scientific script defaults, functions, or arguments were changed during release metadata cleanup.

## Software environments

See [environment/README.md](environment/README.md). Historical R and Python evidence is recorded by stage, while current inventories are identified separately. No single fully locked historical environment is claimed. The scikit-learn 1.8.0 observation for the saved forest estimator is supported only within the recorded model-metadata scope.

## Known reproducibility limitations

The following limitations remain in force and are detailed in [docs/REPRODUCIBILITY_LIMITATIONS.md](docs/REPRODUCIBILITY_LIMITATIONS.md):

- The dynamic exact command line is not fully recovered.
- Target-per-module selection provenance is incomplete.
- The preservation stochastic seed was not recorded.
- No fully locked single historical environment is available.
- Exact end-to-end replay is not claimed for all historical stages.

The previously recorded Figure 2, Supplementary Figure S1, and Supplementary Figure S2 terminology blocker is **RESOLVED** through verified label-only corrections. Source values and scientific geometry are unchanged. See [docs/FIGURE_TERMINOLOGY_CLEANUP.md](docs/FIGURE_TERMINOLOGY_CLEANUP.md).

## Citation

See [CITATION.cff](CITATION.cff). Please cite the associated manuscript when available. Author metadata and the manuscript DOI will be added only after confirmation by the author team.

## License

Repository code is released under the [MIT License](LICENSE). Original public transcriptomic datasets remain subject to their source-database or provider terms. This repository does not relicense third-party source data.

## Associated manuscript

“Integrative Transcriptomics Identifies Reproducible Coexpression Programs across Grapevine Genotypes during Bud Dormancy Release.”

Manuscript drafts and internal discussions are not distributed in this repository.
