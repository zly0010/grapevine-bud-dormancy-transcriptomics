# Reproducibility levels

## Level 1A — Public-data acquisition and matrix/QC reconstruction

Uses public processed/source matrices for GSE124820, GSE273240, GSE184114 and GSE277812. The base harmonization branch prepares GSM-named matrices and QC/design tables; the separately recorded legacy parser branch prepares legacy matrix representations from shared public sources. These are not a fully demonstrated unified base-to-legacy conversion chain.

## Level 1B — Frozen processed-input reproduction / inspection

Starts from Dataset-archived legacy matrices, fixed VSTs, design files, module membership and formal result tables. Restore the documented logical layout using the copy-only utility before selecting a downstream stage. Exact replay is not claimed for the historical stochastic/CLI gaps described in REPRODUCIBILITY_LIMITATIONS.md.

## Level 2 — GSE337039 raw-read validation

Starts from public SRA FASTQ, with fastp0.24.3, STAR2.7.11b, featureCounts2.0.8, reference GCF_000003745.3 and forward-stranded counts. The archive retains processed validation inputs/results; raw FASTQ/SRA/BAM/SAM and STAR indexes are not redistributed. This branch is distinct from the old GSE337039 processed matrix route.
