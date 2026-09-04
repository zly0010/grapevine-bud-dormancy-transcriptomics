# Historical versus final scripts

Only `scripts/method_revision_v2/15_whole_module_validation_v2.py` is authoritative for final external random-set and Stouffer statistics: B=10,000, P=(b+1)/(10,000+1), weights=sqrt(effective n).

`scripts/extended_pipeline/15_whole_module_validation.py` is HISTORICAL_REQUIRED for contextual projection. Its old coherence/null/integration outputs are HISTORICAL_ONLY and do not replace the final CSVs. They live in `final_results/historical_context/15_whole_module_validation/` in the Dataset archive. Figure4c/d retain existing contextual score/effect sources; Figure4a/b/e use method_revision_v2.

The dynamic script remains the archived implementation, while its defaults do not reproduce all known frozen settings. Comment headers flag this distinction; defaults, arguments and functions are unchanged. The complete historical CLI and target selection remain unknown.

Historical remote mapping scripts explain the frozen mapping origin. Their current network responses would not establish identical inputs. The frozen mapping layer supplies authoritative GO inputs.

The file-level manifest and AUTHORITATIVE_ANALYSIS_MAP identify other historic and utility scripts. Original data column/file names containing variety are retained to preserve scientific schemas; publication graphic terminology was corrected separately without changing those schemas or running these scripts.
