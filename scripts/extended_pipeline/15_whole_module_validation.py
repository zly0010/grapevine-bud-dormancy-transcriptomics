# WARNING: This script belongs to a historical validation stage.
# Its default random-set size does not represent the final manuscript inference.
# Final random-set and Stouffer statistics were generated with
# scripts/method_revision_v2/15_whole_module_validation_v2.py using B = 10,000.
# Retained contextual projection outputs are documented in the authoritative map.

#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import math
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.stats as st
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.multitest import multipletests


TARGET_MODULES = ("blue", "turquoise", "brown")
SEED = 20260718


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, na_rep="")


def read_vst(path: Path, genes: list[str]) -> pd.DataFrame:
    header = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
    first = header[0]
    selected = [first] + [gene for gene in genes if gene in header]
    frame = pd.read_csv(path, sep="\t", usecols=selected).set_index(first)
    frame.index = frame.index.astype(str)
    return frame


def standardize(frame: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    matrix = frame.to_numpy(dtype=float)
    sd = matrix.std(axis=0, ddof=1)
    keep = np.isfinite(sd) & (sd > 0)
    matrix = matrix[:, keep]
    genes = frame.columns.to_numpy(dtype=object)[keep].astype(str).tolist()
    z = (matrix - matrix.mean(axis=0)) / matrix.std(axis=0, ddof=1)
    return z, genes


def mean_pairwise_correlation(z: np.ndarray) -> float:
    n, p = z.shape
    if n < 3 or p < 2:
        return math.nan
    total_correlation = float(np.square(z.sum(axis=1)).sum() / (n - 1))
    return (total_correlation - p) / (p * (p - 1))


def reference_loadings(
    reference: pd.DataFrame, modules: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    loading_rows: list[pd.DataFrame] = []
    qc_rows: list[dict[str, object]] = []
    for module in TARGET_MODULES:
        genes = modules.loc[modules.module_color == module, "gene_id"].astype(str).tolist()
        z, usable = standardize(reference[genes])
        _, singular, vt = np.linalg.svd(z, full_matrices=False)
        loadings = vt[0, :]
        score = z @ loadings
        mean_score = z.mean(axis=1)
        if np.corrcoef(score, mean_score)[0, 1] < 0:
            loadings = -loadings
            score = -score
        loading_rows.append(
            pd.DataFrame(
                {"module_color": module, "gene_id": usable, "discovery_pc1_loading": loadings}
            )
        )
        qc_rows.append(
            {
                "module_color": module,
                "genes": len(usable),
                "discovery_pc1_variance": float(singular[0] ** 2 / np.square(z).sum()),
                "discovery_mean_pairwise_correlation": mean_pairwise_correlation(z),
            }
        )
    return pd.concat(loading_rows, ignore_index=True), pd.DataFrame(qc_rows)


def project_and_test_coherence(
    dataset: str,
    vst: pd.DataFrame,
    modules: pd.DataFrame,
    loadings: pd.DataFrame,
    permutations: int,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    z_all, genes_all = standardize(vst)
    index = {gene: i for i, gene in enumerate(genes_all)}
    scores: list[pd.DataFrame] = []
    qc: list[dict[str, object]] = []
    null_rows: list[pd.DataFrame] = []

    for module in TARGET_MODULES:
        reference_genes = modules.loc[modules.module_color == module, "gene_id"].astype(str).tolist()
        module_loadings = loadings.loc[loadings.module_color == module].set_index("gene_id")
        usable = [gene for gene in reference_genes if gene in index and gene in module_loadings.index]
        if len(usable) < 20:
            raise RuntimeError(f"{dataset}/{module}: fewer than 20 usable genes")
        indices = np.array([index[gene] for gene in usable], dtype=int)
        z = z_all[:, indices]
        weight = module_loadings.loc[usable, "discovery_pc1_loading"].to_numpy(dtype=float)
        weight = weight / np.sqrt(np.square(weight).sum())
        projected = z @ weight
        projected = (projected - projected.mean()) / projected.std(ddof=1)
        scores.append(
            pd.DataFrame(
                {
                    "dataset": dataset,
                    "sample_id": vst.index,
                    "module_color": module,
                    "projected_module_score": projected,
                }
            )
        )

        observed = mean_pairwise_correlation(z)
        null = np.empty(permutations, dtype=float)
        for b in range(permutations):
            sampled = rng.choice(z_all.shape[1], size=len(usable), replace=False)
            null[b] = mean_pairwise_correlation(z_all[:, sampled])
        null_mean = float(null.mean())
        null_sd = float(null.std(ddof=1))
        coherence_z = (observed - null_mean) / null_sd if null_sd > 0 else math.nan
        empirical_p = (1 + int(np.sum(null >= observed))) / (permutations + 1)
        empirical_z = float(st.norm.isf(empirical_p))
        if observed < null_mean:
            empirical_z = -empirical_z
        coverage = len(usable) / len(reference_genes)
        qc.append(
            {
                "dataset": dataset,
                "module_color": module,
                "samples": len(vst),
                "reference_genes": len(reference_genes),
                "usable_genes": len(usable),
                "gene_coverage": coverage,
                "coverage_grade": "full" if coverage >= 0.60 else "partial" if coverage >= 0.40 else "inadequate",
                "observed_mean_pairwise_correlation": observed,
                "null_mean": null_mean,
                "null_sd": null_sd,
                "coherence_z": coherence_z,
                "empirical_p": empirical_p,
                "calibrated_empirical_z": empirical_z,
                "permutations": permutations,
            }
        )
        null_rows.append(
            pd.DataFrame(
                {
                    "dataset": dataset,
                    "module_color": module,
                    "permutation": np.arange(1, permutations + 1),
                    "null_mean_pairwise_correlation": null,
                }
            )
        )
    return pd.concat(scores, ignore_index=True), pd.DataFrame(qc), pd.concat(null_rows, ignore_index=True)


def type3_effects(data: pd.DataFrame, formula: str, dataset: str, family: str) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for module in TARGET_MODULES:
        subset = data.loc[data.module_color == module].copy()
        model = smf.ols(formula, data=subset).fit()
        design = np.asarray(model.model.exog)
        if np.linalg.matrix_rank(design) != design.shape[1]:
            raise RuntimeError(f"Rank-deficient model: {dataset}/{family}/{module}")
        table = anova_lm(model, typ=3).reset_index(names="effect")
        table = table.loc[~table.effect.isin(["Intercept", "Residual"])].copy()
        table = table.rename(columns={"F": "f_value", "PR(>F)": "p_value"})
        residual_ss = float(np.square(model.resid).sum())
        table["partial_eta_squared"] = table.sum_sq / (table.sum_sq + residual_ss)
        table["dataset"] = dataset
        table["analysis_family"] = family
        table["module_color"] = module
        rows.append(table)
    return pd.concat(rows, ignore_index=True)[
        ["dataset", "analysis_family", "module_color", "effect", "df", "f_value", "p_value", "partial_eta_squared"]
    ]


def leave_one_dataset_out(coherence: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    datasets = sorted(coherence.dataset.unique())
    for module in TARGET_MODULES:
        module_data = coherence.loc[coherence.module_color == module]
        for omitted in ["NONE"] + datasets:
            kept = module_data if omitted == "NONE" else module_data.loc[module_data.dataset != omitted]
            # Use the finite-permutation calibrated Z for evidence synthesis.
            # The raw null-SD Z is retained as an effect-size diagnostic but can
            # become extremely large for modules containing thousands of genes.
            z = kept.calibrated_empirical_z.to_numpy(dtype=float)
            weights = np.sqrt(kept.samples.to_numpy(dtype=float))
            combined = float(np.sum(weights * z) / np.sqrt(np.square(weights).sum()))
            rows.append(
                {
                    "module_color": module,
                    "omitted_dataset": omitted,
                    "datasets_retained": len(kept),
                    "weighted_stouffer_z": combined,
                    "one_sided_p": float(st.norm.sf(combined)),
                    "all_dataset_z_positive": bool(np.all(z > 0)),
                }
            )
    result = pd.DataFrame(rows)
    result["fdr"] = multipletests(result.one_sided_p, method="fdr_bh")[1]
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("source_project", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--permutations", type=int, default=200)
    parser.add_argument("--gse337039-vst", type=Path)
    parser.add_argument("--gse337039-design", type=Path)
    args = parser.parse_args()
    if args.permutations < 100:
        raise ValueError("Use at least 100 permutations")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    root = args.source_project
    module_path = root / "results_corrected/07_wgcna_fixed/GSE124820/module_assignments.txt"
    ref_path = root / "results_corrected/07_wgcna_fixed/GSE124820_vst_fixed.txt"
    modules = pd.read_csv(module_path, sep="\t")
    reference = read_vst(ref_path, modules.gene_id.astype(str).tolist())
    loadings, discovery_qc = reference_loadings(reference, modules)

    datasets = {
        "GSE273240": (
            root / "results_corrected/07_wgcna_fixed/GSE273240_vst_fixed.txt",
            root / "results_corrected/03_deseq2_gse273240/sample_design_gse273240.txt",
        ),
        "GSE184114": (
            root / "results_corrected/07_wgcna_fixed/GSE184114_vst_fixed.txt",
            root / "results_corrected/01_sample_design_and_qc/sample_design_gse184114.txt",
        ),
        "GSE277812": (
            root / "results_corrected/07_wgcna_fixed/GSE277812_vst_fixed.txt",
            root / "results_corrected/05_deseq2_gse277812/sample_design_gse277812.txt",
        ),
    }
    if (args.gse337039_vst is None) != (args.gse337039_design is None):
        raise ValueError("Provide both --gse337039-vst and --gse337039-design")
    if args.gse337039_vst is not None:
        datasets["GSE337039"] = (args.gse337039_vst, args.gse337039_design)
    for paths in datasets.values():
        if not all(path.exists() for path in paths):
            raise FileNotFoundError(paths)

    rng = np.random.default_rng(SEED)
    score_frames: list[pd.DataFrame] = []
    qc_frames: list[pd.DataFrame] = []
    null_frames: list[pd.DataFrame] = []
    design_frames: dict[str, pd.DataFrame] = {}
    for dataset, (vst_path, design_path) in datasets.items():
        print(f"Projecting fixed discovery loadings: {dataset}", flush=True)
        vst = read_vst(vst_path, modules.gene_id.astype(str).tolist())
        design = pd.read_csv(design_path, sep="\t")
        if dataset == "GSE337039" and "run_accession" in design.columns:
            if "sample_id" in design.columns:
                design = design.rename(columns={"sample_id": "source_sample_id"})
            design = design.rename(columns={"run_accession": "sample_id"})
        elif "sample_id" not in design.columns and "run_accession" in design.columns:
            design = design.rename(columns={"run_accession": "sample_id"})
        design.sample_id = design.sample_id.astype(str)
        if set(vst.index) != set(design.sample_id):
            raise RuntimeError(f"{dataset}: VST/design sample mismatch")
        vst = vst.loc[design.sample_id]
        scores, qc, null = project_and_test_coherence(
            dataset, vst, modules, loadings, args.permutations, rng
        )
        score_frames.append(scores)
        qc_frames.append(qc)
        null_frames.append(null)
        design_frames[dataset] = design

    scores = pd.concat(score_frames, ignore_index=True)
    coherence = pd.concat(qc_frames, ignore_index=True)
    nulls = pd.concat(null_frames, ignore_index=True)

    effects: list[pd.DataFrame] = []
    gse337_timecourse: pd.DataFrame | None = None
    s273 = scores.loc[scores.dataset == "GSE273240"].merge(design_frames["GSE273240"], on="sample_id", validate="many_to_one")
    effects.append(type3_effects(
        s273,
        "projected_module_score ~ C(deac_phase, Sum) * C(day, Sum) * C(treatment, Sum)",
        "GSE273240",
        "ABA_deacclimation",
    ))

    s184 = scores.loc[scores.dataset == "GSE184114"].merge(design_frames["GSE184114"], on="sample_id", validate="many_to_one")
    for phase in sorted(s184.phase.unique()):
        # Time zero contains duplicated baseline labels that do not form a full
        # treatment-by-time factorial; the estimable arm-level model starts > 0 h.
        arm = s184.loc[(s184.phase == phase) & (s184.time_h > 0)].copy()
        effects.append(type3_effects(
            arm,
            "projected_module_score ~ C(time_h, Sum) * C(treatment, Sum)",
            "GSE184114",
            f"{phase}_within_tissue_arm",
        ))

    s277 = scores.loc[scores.dataset == "GSE277812"].merge(design_frames["GSE277812"], on="sample_id", validate="many_to_one")
    effects.append(type3_effects(
        s277,
        "projected_module_score ~ C(stage_factor, Sum) * C(node_factor, Sum)",
        "GSE277812",
        "stage_node",
    ))
    if "GSE337039" in design_frames:
        s337 = scores.loc[scores.dataset == "GSE337039"].merge(
            design_frames["GSE337039"], on="sample_id", validate="many_to_one"
        )
        gse337_timecourse = (
            s337.groupby(["cultivar", "time_index", "module_color", "condition"], as_index=False)
            .projected_module_score.agg(["mean", "std", "count"])
            .reset_index()
        )
        for cultivar in sorted(s337.cultivar.unique()):
            arm = s337.loc[s337.cultivar == cultivar].copy()
            effects.append(type3_effects(
                arm,
                "projected_module_score ~ C(condition, Sum) * C(time_index, Sum)",
                "GSE337039",
                f"{cultivar}_chilling_time_course",
            ))
    effect_table = pd.concat(effects, ignore_index=True)
    effect_table["fdr"] = multipletests(effect_table.p_value, method="fdr_bh")[1]

    loo = leave_one_dataset_out(coherence)
    write_tsv(loadings, args.output_dir / "01_discovery_pc1_loadings.tsv")
    write_tsv(discovery_qc, args.output_dir / "02_discovery_module_qc.tsv")
    write_tsv(scores, args.output_dir / "03_independent_projected_scores.tsv")
    write_tsv(coherence, args.output_dir / "04_whole_module_coherence.tsv")
    write_tsv(nulls, args.output_dir / "05_coherence_null_permutations.tsv")
    write_tsv(effect_table, args.output_dir / "06_projected_score_effects.tsv")
    write_tsv(loo, args.output_dir / "07_leave_one_dataset_out_coherence.tsv")
    if gse337_timecourse is not None:
        write_tsv(gse337_timecourse, args.output_dir / "08_gse337039_module_timecourse.tsv")
        means = gse337_timecourse.pivot_table(
            index=["cultivar", "time_index", "module_color"],
            columns="condition",
            values="mean",
        ).reset_index()
        means.columns.name = None
        means["controlled_minus_natural"] = means["Controlled_4C"] - means["Natural"]
        write_tsv(means, args.output_dir / "09_gse337039_controlled_minus_natural.tsv")

    lines = [
        "# Whole-module independent validation",
        "",
        f"Discovery PC1 loadings were fixed in GSE124820 and projected without refitting into {len(datasets)} independent datasets.",
        f"Each module's coexpression was compared with {args.permutations} size-matched random gene sets per dataset.",
        "GSE184114 tissue arms were modeled separately.",
        "",
    ]
    for module in TARGET_MODULES:
        full = loo.loc[(loo.module_color == module) & (loo.omitted_dataset == "NONE")].iloc[0]
        sensitivity = loo.loc[(loo.module_color == module) & (loo.omitted_dataset != "NONE")]
        lines.append(
            f"- {module}: full coherence Z={full.weighted_stouffer_z:.2f}, FDR={full.fdr:.3g}; "
            f"minimum leave-one-dataset-out Z={sensitivity.weighted_stouffer_z.min():.2f}."
        )
    if "GSE337039" in design_frames:
        lines += ["", "GSE337039 condition-by-time effects:"]
        interaction = effect_table.loc[
            (effect_table.dataset == "GSE337039") &
            effect_table.effect.str.contains(":", regex=False)
        ]
        for row in interaction.sort_values(["analysis_family", "module_color"]).itertuples(index=False):
            lines.append(
                f"- {row.analysis_family}, {row.module_color}: interaction FDR={row.fdr:.3g}, "
                f"partial eta-squared={row.partial_eta_squared:.3f}."
            )
    (args.output_dir / "SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (args.output_dir / "INPUT_SHA256.tsv").write_text(
        "file\tsha256\n" + "\n".join(
            f"{path.as_posix()}\t{sha256(path)}"
            for path in [module_path, ref_path] + [p for pair in datasets.values() for p in pair]
        ) + "\n",
        encoding="utf-8",
    )
    print("WHOLE_MODULE_VALIDATION_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
