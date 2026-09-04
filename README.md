# Integrative Transcriptomics of Grapevine Bud Dormancy Release

Version 0.1.0-rc2. **UPLOAD-READY BUT NOT PUBLICLY RELEASED. PUBLIC_RELEASE = BLOCKED.** Licenses, creators and final author/asset review remain pending. Suggested repository name: `grapevine-bud-dormancy-transcriptomics` (no repository has been created).

## Overview

The repository provides the analysis code, frozen processed inputs, formal outputs, figure source data, and documented computational provenance supporting the manuscript. Large frozen inputs/results are supplied in the companion Dataset directory; small figure source tables are retained here for inspection. Scientific results are frozen. Exact replay is not claimed for several historical stochastic or command-level steps documented in the reproducibility limitations.

## Public datasets

| Dataset | Role | Original samples | Starting source |
| --- | --- | ---: | --- |
| GSE124820 | Discovery; cross-genotype meta-analysis; reference WGCNA |189|Public processed/source counts|
| GSE273240 | Preservation and contextual projection |90|Public processed/source counts|
| GSE184114 | Acclimation/deacclimation context |74|Public processed/source counts|
| GSE277812 | Developmental context |27|Public processed/source counts|
| GSE337039 | Independent cold validation |60|Public SRA raw reads|

Accession links and metadata are in [dataset_manifest.tsv](metadata/dataset_manifest.tsv). Discovery QC-B retains172 samples and QC-C147; the frozen189-row QC table is preserved. Distinct external tissue/experimental contexts remain documented.

## Repository structure

scripts/ contains base, legacy discovery, extended, final method revision, figure and utility implementations. config/ is a parameter catalog. metadata/ contains small design/mapping resources. environment/ separates historical stage records from CURRENT observations. figure_source_data/Figure1–6 and FigureS1/S2 contain small source tables and panel maps. docs/ explains source roles, dependencies, limitations and the Dataset restoration layout.

## Authoritative analysis workflow

See [AUTHORITATIVE_ANALYSIS_MAP.md](docs/AUTHORITATIVE_ANALYSIS_MAP.md). Only `scripts/method_revision_v2/15_whole_module_validation_v2.py` supplies final B=10,000 random-set/Stouffer inference, with P=(b+1)/(10,000+1) and sqrt(effective n) weights. The older whole-module script remains HISTORICAL_REQUIRED for contextual projection. Main WGCNA powers12/16 are automatic-selection results; preservation seed remains UNKNOWN.

## Reproducibility levels

[Level1A](docs/REPRODUCIBILITY_LEVELS.md) covers public matrix acquisition/QC; Level1B starts from retained frozen processed inputs/results. They are not one fully demonstrated base→legacy conversion chain. Level2 covers the separate GSE337039 raw-read branch.

## Public-data acquisition

The four processed-matrix datasets feed [parallel base/legacy branches](docs/MATRIX_PROCESSING_BRANCHES.md). Historical execution records confirm the legacy parser was run; direct base→legacy conversion was not established. Frozen GO mappings, rather than latest remote responses, are authoritative.

## Frozen processed inputs

Obtain the companion Dataset directory and follow [workflow.md](docs/workflow.md). The explicit materialize_data utility only copies/verifies files; it never begins scientific analysis. No implicit one-command scientific workflow is provided.

## GSE337039 raw-read validation

See [GSE337039_RAW_PROCESSING.md](docs/GSE337039_RAW_PROCESSING.md): fastp0.24.3, STAR2.7.11b, featureCounts2.0.8, GCF_000003745.3, forward-stranded counts. Raw FASTQ/SRA/BAM/SAM and genome indexes are not included.

## Historical vs final analysis scripts

[HISTORICAL_VS_FINAL_SCRIPTS.md](docs/HISTORICAL_VS_FINAL_SCRIPTS.md) distinguishes retained implementation/context from final inference. Public-copy comments flag old defaults. No scientific script defaults, functions or arguments were changed.

## Software environments

See [environment/README.md](environment/README.md). Historical R/Python evidence is recorded by stage; CURRENT inventories are separate. No single fully locked historical environment is claimed. Forest sklearn1.8.0 is supported only within the saved-estimator metadata scope.

## Processed-data archive

The local `zenodo_upload_ready_v0.1.0-rc2` supplies large processed inputs, final outputs, models, full frozen figure sources/exports and supplementary material. It has not been uploaded and has no assigned DOI. DOI and repository links remain TODO. Physical duplicate data are avoided using documented layout aliases.

## Known reproducibility limitations

See [REPRODUCIBILITY_LIMITATIONS.md](docs/REPRODUCIBILITY_LIMITATIONS.md). Complete dynamic CLI/cwd/wrapper were NOT RECOVERED; target-per-module remains UNKNOWN. Preservation external seed remains UNKNOWN. These accepted provenance limitations do not call for invented parameters. The initially recorded Figure2/S1/S2 graphic-terminology blocker is RESOLVED by verified label-only corrections; source values and scientific geometry are unchanged. See [figure correction record](docs/FIGURE_TERMINOLOGY_CLEANUP.md).

## Citation

CITATION.cff is a syntactically valid YAML draft. Confirmed creators, optional ORCID, repository URL and archive/manuscript identifiers must be supplied by authors. No identities or identifiers are invented.

## License

See [LICENSE_DECISION_REQUIRED.md](LICENSE_DECISION_REQUIRED.md). Code and scientific-data licensing must be approved separately. LICENSE is an explicit no-release placeholder, not a grant.

## Associated manuscript

“Integrative Transcriptomics Identifies Reproducible Coexpression Programs across Grapevine Genotypes during Bud Dormancy Release”. Manuscript DOI: TODO if available. Manuscript drafts and internal discussions are not distributed.

SHA256SUMS.tsv covers every other repository file; its own hash is in the external GITHUB_UPLOAD_MANIFEST.tsv. docs/FILE_PROVENANCE.tsv provides source identities and transformations without private paths.
