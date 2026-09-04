# Matrix processing branches

```text
Public processed/source data
  +-- base harmonization: scripts/base_pipeline/04_prepare_count_matrices.py
  |     -> GSM-named base matrices, sample maps, design and QC
  +-- legacy discovery: scripts/legacy_discovery_pipeline/07_parse_all_datasets.py
        -> legacy data/processed/GSE*/counts_matrix.txt
```

Historical execution records confirm that the legacy parser was run during the early legacy-processing stage. Its current public copy changes only two private path literals plus a comment header; original-source and public-copy hashes are recorded in FILE_PROVENANCE. The historical execution did not capture a separate source-code hash, so byte identity with the executed script is not asserted.

Five inspected compressed public-matrix pairs shared identical source bytes. GSE124820 uses189 individual count files. Base and legacy column naming differ; GSE184114 legacy processing merges the40/34-sample contexts with Accl_/Deaccl_ prefixes. Static gene/sample-ID relations do not prove a historical base→legacy operation. The direct base→legacy historical command was not found, and no direct conversion is claimed.

For a future Level1A legacy reconstruction, populate `data/legacy/data/raw/GSE124820/counts/` with the original per-sample count files, and `data/legacy/data/raw/GSE273240/`, `GSE184114/`, `GSE277812/` with the corresponding public source matrices in their saved uncompressed names. The archived parser also contains the historic GSE337039 processed route; that output must not enter final validation. No parsing logic or default job list was changed during packaging. Level1B instead begins from the retained archive inputs and does not require this acquisition branch.

Fixed VSTs are samples×genes, including172×16529 discovery values. They include normalization/filtering steps and are not merely transposed raw counts. Their numerical values and design files are retained unchanged.
