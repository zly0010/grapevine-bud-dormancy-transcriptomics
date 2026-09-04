# Supplementary source mapping

**Table S1 — FOUND.** The archived `supplementary/TableS1/qc_sensitivity_comparison.txt` is an existing table with variety, qc_scheme, n_lrt_sig and n_deg_day10. `scripts/legacy_discovery_pipeline/02b_cross_variety_consensus.py` creates it from existing per-variety DESeq2 outputs. The manuscript comparison uses the B_no_fail and C_pass_only rows; A_all rows remain in the existing source for provenance. No new statistics or manuscript table extraction was used as a substitute.

**Table S2 — FOUND.** The final workbook `Supplementary_Table_S2_GSE124820_sample_QC.xlsx` was copied unchanged. Source: `GSE124820_raw_counts_design_qc.tsv`, SHA-256 `a2e9605b55b6413da346c3be5ab289203845ea261e5423c1fe86254a7c1ea8fa`. Keep all 189 rows and their original order. Select these ten source columns: sample_accession, genotype, short_genotype, time, time_unit, replicate, library_size, detected_genes, qc_status, qc_reason. Preserve their values exactly. Store time, replicate, library_size and detected_genes numerically; other fields as text.

Append included_in_QC_B and included_in_QC_C with this deterministic lookup:

| Existing qc_status | included_in_QC_B | included_in_QC_C |
|---|---|---|
| pass | YES | YES |
| watch | YES | NO |
| fail | NO | NO |

Do not recompute QC status or invent ratios. The final workbook contains a Sample_QC sheet and README rule table, with VLOOKUP formulas and cached YES/NO values. The frozen input has 147 pass, 25 watch, 17 fail; QC-B has 172 included, QC-C has 147. The workbook hash is recorded in DATA_MANIFEST.tsv. This specification replaces dependence on the original internal spreadsheet-generation runtime; no new workbook is generated here.

**Figures S1 and S2 — FOUND.** See figure_source_data/FigureS1 and FigureS2. S1 panel c is the fold-enrichment versus nominal-P scatter plot, using the same GO-test table as panel b (including ontology and module labels); no artificial panel-c table was created.
