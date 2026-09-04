# Stage-specific software evidence

historical_R_versions.tsv and historical_Python_versions.tsv consolidate directly supported historical facts from the completed provenance audit. software_versions.tsv covers GSE337039 raw processing. No package versions were obtained by new historical searches or inferred from current installations.

CURRENT observations from the rc1 preparation are named current_R_environment.* and current_Python_environment.*. They are not historical locks. Preserved sessionInfo files under historical/ are supporting stage records; matching/duplicated text is not an independent execution.

WGCNA1.74 is confirmed only for stability analysis. The formal main WGCNA runtime was R4.4.1; its package version and preservation external seed remain unknown. Early OpenScience numpy2.4.4 must not replace frozen dynamic numpy2.5.1. sklearn1.8.0 is supported by the three retained forest estimators, not every process in the pipeline.

No unified fully locked historical environment, conda solve or renv lock was fabricated. Missing versions remain UNKNOWN. No packages or scientific workflows were installed/executed during rc2 packaging.
