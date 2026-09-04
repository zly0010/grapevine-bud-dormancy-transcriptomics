# Documented reproducibility limitations

## WGCNA stochastic provenance

The formal WGCNA run and selected powers are supported by preserved execution logs. However, no explicit external RNG seed was recovered for the historical modulePreservation run. Exact stochastic bit-for-bit reproduction of the preservation permutations is therefore not claimed.

The recorded runtime was R4.4.1. The retained source expression is normalized publicly as `source("scripts/legacy_discovery_pipeline/wgcna_r_fixed.R")`. Historical execution used a local absolute path; the public path shown here is repository-relative. The exact outer shell invocation/cwd was not independently recovered.

Selected power12 for GSE124820 and16 for GSE273240 are automatic-selection results, not manually fixed historical CLI flags. The saved rule searches1:30 for signed-network negative-slope/signedR2 candidates, using tolerances0.80 and0.699 respectively, with the recorded connectivity preference. Preservation used200 permutations. WGCNA1.74 is confirmed for the later stability-analysis stage; the main-run WGCNA package version was not directly recovered.

## Dynamic command-level provenance

The frozen inputs, major run parameters, model artifacts, and final outputs are preserved. However, the complete historical command-line invocation, including the target-per-module selection argument, was not recovered. The archived script and frozen outputs therefore provide computational provenance and reusable implementation, but an exact command-level replay of the historical perturbation run is not claimed.

Recorded facts:282genes,49TFs,57 genotype×time centroids,53 adjacent transitions,60trees per target,seeds1/2/3,CV30trees and2,000 module-label permutations. Final saved flags identify20 perturbation combinations and11TFs. Complete historical CLI:NOT RECOVERED. target-per-module:UNKNOWN / NOT RECOVERED. cwd/wrapper:NOT RECOVERED. Resume use is also unknown. No partial command is presented as an exact reproduction command.

RUN_INFO records Python3.14.4,numpy2.5.1,pandas3.0.2. scikit-learn1.8.0 is supported by metadata in the three saved forest estimators; this does not fill every dependency of every process or the sparse artifact. All six retained RUN_INFO input hashes were matched in the completed provenance audit.

## Matrices, GO and environments

Public source data feed parallel base-harmonization and legacy-parser branches. A historical legacy-parser execution was recovered; a direct base→legacy command was not found, so no such direct conversion is claimed. Source bytes at the historical execution were not all independently hashed.

GO uses frozen mapping tables associated with Ensembl Plants release62, ASM3070453v1, UniProtKB organism29760 and query/freeze date2026-07-16. Remote raw responses were not fully preserved; the frozen mapping layer is the authoritative reproducibility input. Re-querying the latest UniProt is not an instruction for obtaining identical mappings.

Historical software evidence is stage-specific. No single fully locked historical environment has been generated. readr/tidyverse versions and the formal main WGCNA package version remain unknown. Current environment snapshots are labeled CURRENT, not used to fill historical gaps.

Several plotting scripts recompute descriptive correlations/PCA/summaries when run; they were not executed in release preparation. Existing source values and frozen exports were copied. Figure1c preserves gene order and existing vector/image geometry, not a newly reconstructed dendrogram object.

The parallel STAR wrapper can attempt simultaneous index creation; construct the reference index before parallel workers. Historic/raw/projection routes are distinguished from final inference in the authoritative map.

These are documented provenance and replay limitations, not claims that the frozen scientific findings are unreliable, and not requirements to invent/recover missing parameters before submission. Licenses, creators and unresolved graphic terminology are separate publication gates.
