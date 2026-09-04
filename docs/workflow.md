# Workflow and data layout

Read AUTHORITATIVE_ANALYSIS_MAP.md before selecting any script. Filename order is not a dependency graph. Public matrices feed parallel base and legacy branches; reference discovery results then support meta-analysis, robustness, contextual projection, GO/evidence integration and dynamic-network prioritization. The GSE337039 raw branch supplies an additional processed validation input to the final external synthesis.

Only method_revision_v2 supplies final random-set/Stouffer inference; historical contextual projection remains separately documented. Source tables and finalized figures can be inspected without running analysis.

Place the Dataset archive next to this repository, then explicitly invoke the copy-only utility:

```bash
python scripts/utilities/materialize_data.py --archive ../zenodo_upload_ready_v0.1.0-rc2 --destination .
```

This verifies hashes and restores logical `data/legacy` and `outputs/advanced` paths from docs/DATA_LAYOUT.tsv. It never starts scientific analysis. Aliases may restore one archived physical file to several documented logical paths. Some historic scripts overwrite their outputs; a future scientific rerun should use a separate working copy after reviewing the stage and limitations. No exact historical dynamic CLI is supplied.

Frozen GO mapping tables are authoritative inputs; do not fetch a latest remote replacement to claim equality. Raw GSE337039 instructions are in GSE337039_RAW_PROCESSING.md. Figure sources and supplementary sources have their own maps.
