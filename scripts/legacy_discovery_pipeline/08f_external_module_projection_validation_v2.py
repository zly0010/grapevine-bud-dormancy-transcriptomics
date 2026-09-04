#!/usr/bin/env python
from __future__ import annotations

import hashlib
import math
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import statsmodels
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.multitest import multipletests


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, na_rep="")


def read_vst_selected(path: Path, target_genes: set[str]) -> pd.DataFrame:
    header = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
    if not header:
        raise RuntimeError(f"Empty VST header: {path}")
    first_column = header[0]
    selected = [first_column] + [column for column in header[1:] if column in target_genes]
    frame = pd.read_csv(path, sep="\t", usecols=selected)
    frame = frame.set_index(first_column)
    frame.index = frame.index.astype(str)
    return frame


def project_scores(
    vst: pd.DataFrame,
    modules: pd.DataFrame,
    target_modules: list[str],
    dataset_name: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    score_frames: list[pd.DataFrame] = []
    qc_rows: list[dict] = []
    gene_frames: list[pd.DataFrame] = []
    for module_name in target_modules:
        reference_genes = modules.loc[modules["module_color"] == module_name, "gene_id"].astype(str).tolist()
        common_genes = [gene for gene in reference_genes if gene in vst.columns]
        if len(common_genes) < 20:
            raise RuntimeError(f"{dataset_name} has fewer than 20 genes for {module_name}")
        matrix = vst[common_genes].to_numpy(dtype=float)
        gene_sd = np.nanstd(matrix, axis=0, ddof=1)
        keep = np.isfinite(gene_sd) & (gene_sd > 0)
        usable_genes = np.asarray(common_genes, dtype=object)[keep].tolist()
        matrix = matrix[:, keep]
        if len(usable_genes) < 20:
            raise RuntimeError(f"{dataset_name} has fewer than 20 non-zero-variance genes for {module_name}")
        gene_mean = matrix.mean(axis=0)
        gene_sd = matrix.std(axis=0, ddof=1)
        z = (matrix - gene_mean) / gene_sd
        score = z.mean(axis=1)
        u, singular_values, _ = np.linalg.svd(z, full_matrices=False)
        pc1 = u[:, 0] * singular_values[0]
        score_pc1_correlation = float(np.corrcoef(score, pc1)[0, 1])
        if score_pc1_correlation < 0:
            pc1 = -pc1
            score_pc1_correlation = -score_pc1_correlation
        pc1_variance = float(singular_values[0] ** 2 / np.sum(singular_values**2))
        score_frames.append(
            pd.DataFrame(
                {
                    "sample_id": vst.index,
                    "module_color": module_name,
                    "module_score": score,
                    "oriented_pc1": pc1,
                }
            )
        )
        qc_rows.append(
            {
                "dataset": dataset_name,
                "module_color": module_name,
                "reference_module_genes": len(reference_genes),
                "common_genes": len(common_genes),
                "usable_nonzero_variance_genes": len(usable_genes),
                "reference_gene_coverage_percent": 100 * len(common_genes) / len(reference_genes),
                "pc1_variance_explained_percent": 100 * pc1_variance,
                "score_pc1_correlation": score_pc1_correlation,
            }
        )
        gene_frames.append(
            pd.DataFrame(
                {"dataset": dataset_name, "module_color": module_name, "gene_id": usable_genes}
            )
        )
    return pd.concat(score_frames, ignore_index=True), pd.DataFrame(qc_rows), pd.concat(gene_frames, ignore_index=True)


def normalize_effect_name(effect: str) -> str:
    replacements = {
        "C(time_h, Sum)": "time_factor",
        "C(treatment, Sum)": "treatment",
        "C(phase, Sum)": "phase",
        "C(stage_factor, Sum)": "stage_factor",
        "C(node_factor, Sum)": "node_factor",
    }
    for old, new in replacements.items():
        effect = effect.replace(old, new)
    return effect


def fit_type3(
    data: pd.DataFrame,
    formula: str,
    dataset: str,
    analysis_family: str,
    module_color: str,
) -> pd.DataFrame:
    model = smf.ols(formula, data=data).fit()
    design = np.asarray(model.model.exog)
    if np.linalg.matrix_rank(design) != design.shape[1]:
        raise RuntimeError(f"Rank-deficient model: {dataset}/{analysis_family}/{module_color}")
    table = anova_lm(model, typ=3).reset_index(names="effect")
    table = table.loc[~table["effect"].isin(["Intercept", "Residual"]), :].copy()
    table = table.rename(columns={"sum_sq": "sum_sq", "df": "df", "F": "f_value", "PR(>F)": "p_value"})
    table["effect"] = table["effect"].map(normalize_effect_name)
    residual_ss = float(np.sum(model.resid**2))
    table["partial_eta_squared"] = table["sum_sq"] / (table["sum_sq"] + residual_ss)
    table["dataset"] = dataset
    table["analysis_family"] = analysis_family
    table["module_color"] = module_color
    table["residual_df"] = float(model.df_resid)
    return table[
        [
            "dataset",
            "analysis_family",
            "module_color",
            "effect",
            "sum_sq",
            "df",
            "f_value",
            "p_value",
            "partial_eta_squared",
            "residual_df",
        ]
    ]


def group_summary(data: pd.DataFrame, group_columns: list[str]) -> pd.DataFrame:
    result = (
        data.groupby(group_columns, observed=True)["module_score"]
        .agg(n="size", mean="mean", sd="std")
        .reset_index()
    )
    result["se"] = result["sd"] / np.sqrt(result["n"])
    return result


def main() -> int:
    if len(sys.argv) != 3:
        raise SystemExit("Usage: 08f_external_module_projection_validation_v2.py <project_root> <output_dir>")
    project_root = Path(sys.argv[1]).resolve()
    output_dir = Path(sys.argv[2]).resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise RuntimeError(f"Refusing to overwrite non-empty output directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    files = {
        "modules": project_root / "results_corrected/07_wgcna_fixed/GSE124820/module_assignments.txt",
        "vst184": project_root / "results_corrected/07_wgcna_fixed/GSE184114_vst_fixed.txt",
        "design184": project_root / "results_corrected/01_sample_design_and_qc/sample_design_gse184114.txt",
        "vst277": project_root / "results_corrected/07_wgcna_fixed/GSE277812_vst_fixed.txt",
        "design277": project_root / "results_corrected/05_deseq2_gse277812/sample_design_gse277812.txt",
    }
    missing = [str(path) for path in files.values() if not path.exists()]
    if missing:
        raise RuntimeError("Missing required files: " + "; ".join(missing))

    print("Reading formal module assignment and selecting target VST columns", flush=True)
    modules = pd.read_csv(files["modules"], sep="\t")
    target_modules = ["blue", "turquoise", "brown"]
    expected_sizes = {"blue": 1045, "turquoise": 2333, "brown": 680}
    if len(modules) != 5000 or modules["gene_id"].duplicated().any():
        raise RuntimeError("Formal module assignment fails 5,000-gene uniqueness check")
    observed_sizes = modules["module_color"].value_counts().to_dict()
    if any(observed_sizes.get(module) != size for module, size in expected_sizes.items()):
        raise RuntimeError("Reference module sizes do not match formal values")
    target_genes = set(modules.loc[modules["module_color"].isin(target_modules), "gene_id"].astype(str))
    vst184 = read_vst_selected(files["vst184"], target_genes)
    vst277 = read_vst_selected(files["vst277"], target_genes)
    design184 = pd.read_csv(files["design184"], sep="\t")
    design277 = pd.read_csv(files["design277"], sep="\t")
    design184["sample_id"] = design184["sample_id"].astype(str)
    design277["sample_id"] = design277["sample_id"].astype(str)
    if len(vst184) != 74 or len(design184) != 74 or set(vst184.index) != set(design184["sample_id"]):
        raise RuntimeError("GSE184114 VST/design mismatch")
    if len(vst277) != 27 or len(design277) != 27 or set(vst277.index) != set(design277["sample_id"]):
        raise RuntimeError("GSE277812 VST/design mismatch")
    vst184 = vst184.loc[design184["sample_id"], :]
    vst277 = vst277.loc[design277["sample_id"], :]
    if vst184.isna().any().any() or vst277.isna().any().any():
        raise RuntimeError("External VST matrix contains missing values")

    print("Computing deterministic projected module scores and SVD coherence diagnostics", flush=True)
    score184, qc184, genes184 = project_scores(vst184, modules, target_modules, "GSE184114")
    score277, qc277, genes277 = project_scores(vst277, modules, target_modules, "GSE277812")
    projection_qc = pd.concat([qc184, qc277], ignore_index=True)
    projection_genes = pd.concat([genes184, genes277], ignore_index=True)
    scores184 = score184.merge(design184, on="sample_id", how="left", validate="many_to_one")
    scores277 = score277.merge(design277, on="sample_id", how="left", validate="many_to_one")

    effect_frames: list[pd.DataFrame] = []
    cell_frames: list[pd.DataFrame] = []
    print("Fitting estimable Type III external validation models", flush=True)
    for module_name in target_modules:
        module_data = scores184.loc[scores184["module_color"] == module_name, :].copy()
        for phase_name in ["Acclimation", "Deacclimation"]:
            data = module_data.loc[(module_data["phase"] == phase_name) & (module_data["time_h"] > 0), :].copy()
            family = f"phase_stratified_{phase_name}"
            effect_frames.append(
                fit_type3(
                    data,
                    "module_score ~ C(time_h, Sum) * C(treatment, Sum)",
                    "GSE184114",
                    family,
                    module_name,
                )
            )
            means = group_summary(data, ["time_h", "treatment"])
            means.insert(0, "phase", phase_name)
            means.insert(0, "module_color", module_name)
            means.insert(0, "analysis_family", family)
            means.insert(0, "dataset", "GSE184114")
            cell_frames.append(means)

    for module_name in target_modules:
        data = scores184.loc[
            (scores184["module_color"] == module_name) & (scores184["time_h"].isin([24, 48])), :
        ].copy()
        effect_frames.append(
            fit_type3(
                data,
                "module_score ~ C(phase, Sum) * C(time_h, Sum) * C(treatment, Sum)",
                "GSE184114",
                "common_24_48h_cross_phase",
                module_name,
            )
        )

    for module_name in target_modules:
        data = scores277.loc[scores277["module_color"] == module_name, :].copy()
        effect_frames.append(
            fit_type3(
                data,
                "module_score ~ C(stage_factor, Sum) * C(node_factor, Sum)",
                "GSE277812",
                "stage_by_node",
                module_name,
            )
        )
        means = group_summary(data, ["stage_factor", "node_factor"])
        means.insert(0, "module_color", module_name)
        means.insert(0, "analysis_family", "stage_by_node")
        means.insert(0, "dataset", "GSE277812")
        cell_frames.append(means)

    effects = pd.concat(effect_frames, ignore_index=True)
    if len(effects) != 48:
        raise RuntimeError(f"Expected 48 external effect tests but observed {len(effects)}")
    effects["fdr_bh_global"] = multipletests(effects["p_value"].to_numpy(), method="fdr_bh")[1]
    effects["fdr_bh_within_family"] = np.nan
    for _, indices in effects.groupby(["dataset", "analysis_family"]).groups.items():
        effects.loc[indices, "fdr_bh_within_family"] = multipletests(
            effects.loc[indices, "p_value"].to_numpy(), method="fdr_bh"
        )[1]
    effects["significant_global_fdr_0_05"] = effects["fdr_bh_global"] < 0.05
    effects["significant_family_fdr_0_05"] = effects["fdr_bh_within_family"] < 0.05
    effects = effects.sort_values(
        ["dataset", "analysis_family", "effect", "fdr_bh_global", "module_color"], kind="stable"
    ).reset_index(drop=True)
    cell_means = pd.concat(cell_frames, ignore_index=True, sort=False)
    significant = effects.loc[effects["significant_global_fdr_0_05"], :].copy()

    write_tsv(projection_qc, output_dir / "module_projection_coverage_and_qc.tsv")
    write_tsv(projection_genes, output_dir / "projected_module_gene_lists.tsv")
    write_tsv(scores184, output_dir / "GSE184114_projected_module_scores.tsv")
    write_tsv(scores277, output_dir / "GSE277812_projected_module_scores.tsv")
    write_tsv(effects, output_dir / "external_projection_type3_effect_tests.tsv")
    write_tsv(significant, output_dir / "external_projection_significant_global_fdr.tsv")
    write_tsv(cell_means, output_dir / "external_projection_cell_means.tsv")

    parameters = pd.DataFrame(
        {
            "parameter": [
                "reference_modules",
                "projection_method",
                "gene_standardization",
                "external_networks",
                "GSE184114_phase_models",
                "GSE184114_cross_phase_model",
                "GSE184114_zero_hour",
                "GSE277812_model",
                "anova_type",
                "global_multiple_testing",
                "family_multiple_testing",
                "primary_significance",
            ],
            "value": [
                ",".join(target_modules),
                "Mean of projected reference-module gene z-scores per sample",
                "Within each external dataset, center and scale each usable gene across samples",
                "No external WGCNA network was constructed",
                "Within phase and post-treatment only: score ~ time_factor * treatment",
                "Common 24/48 h subset: score ~ phase * time_factor * treatment",
                "Excluded from factorial tests because ABA zero-hour cells do not exist; retained in score output",
                "score ~ stage_factor * node_factor",
                "Type III F tests with sum-to-zero contrasts via statsmodels",
                "BH across all 48 external projection effect tests",
                "BH within each dataset x analysis-family block",
                "Global FDR < 0.05",
            ],
        }
    )
    write_tsv(parameters, output_dir / "parameters.tsv")

    manifest_rows = []
    for role, path in files.items():
        stat = path.stat()
        manifest_rows.append(
            {
                "role": role,
                "path": str(path),
                "bytes": stat.st_size,
                "modified_time": pd.Timestamp(stat.st_mtime, unit="s", tz="UTC").isoformat(),
                "md5": md5sum(path),
            }
        )
    write_tsv(pd.DataFrame(manifest_rows), output_dir / "source_file_manifest.tsv")
    (output_dir / "sessionInfo.txt").write_text(
        "\n".join(
            [
                f"Python: {sys.version}",
                f"Platform: {platform.platform()}",
                f"pandas: {pd.__version__}",
                f"numpy: {np.__version__}",
                f"scipy: {scipy.__version__}",
                f"statsmodels: {statsmodels.__version__}",
            ]
        ),
        encoding="utf-8",
    )

    report = [
        "# External reference-module projection validation",
        "",
        "- Status: **COMPLETE**",
        "- Reference: formal GSE124820 blue, turquoise and brown module membership.",
        "- Projection: sample-wise mean of within-dataset gene z-scores for common reference-module genes.",
        "- GSE184114 and GSE277812 were not subjected to new WGCNA network construction.",
        "- Primary multiplicity threshold: BH across all 48 external effect tests; global FDR < 0.05.",
        "",
        "## Projection coverage and coherence",
        "",
        "| Dataset | Module | Reference genes | Common genes | Coverage | PC1 variance | Score-PC1 correlation |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in projection_qc.itertuples(index=False):
        report.append(
            f"| {row.dataset} | {row.module_color} | {row.reference_module_genes} | {row.common_genes} | "
            f"{row.reference_gene_coverage_percent:.2f}% | {row.pc1_variance_explained_percent:.2f}% | "
            f"{row.score_pc1_correlation:.3f} |"
        )
    report.extend(
        [
            "",
            "## Model strategy",
            "",
            "- GSE184114 phase-specific post-treatment models: `score ~ time_factor * treatment`.",
            "- GSE184114 common 24/48 h cross-phase model: `score ~ phase * time_factor * treatment`.",
            "- GSE277812 complete model: `score ~ stage_factor * node_factor`.",
            "- The full GSE184114 time grid was not forced into one factorial model because phases use different time points and zero-hour ABA cells are structurally absent.",
            "",
            "## Globally significant effects",
            "",
        ]
    )
    if significant.empty:
        report.extend(["No projection effect passed global FDR < 0.05.", ""])
    else:
        report.extend(
            [
                "| Dataset | Analysis | Module | Effect | F | Partial eta-squared | P | Global FDR |",
                "|---|---|---|---|---:|---:|---:|---:|",
            ]
        )
        for row in significant.itertuples(index=False):
            report.append(
                f"| {row.dataset} | {row.analysis_family} | {row.module_color} | {row.effect} | "
                f"{row.f_value:.4f} | {row.partial_eta_squared:.4f} | {row.p_value:.3g} | "
                f"{row.fdr_bh_global:.3g} |"
            )
        report.append("")
    report.extend(
        [
            "## Guardrails",
            "",
            "- Projection supports reproducibility of coordinated module expression, not preservation of external network topology.",
            "- Module scores have a fixed expression direction: higher score means higher average standardized expression of reference-module genes.",
            "- PC1 metrics are coherence diagnostics only and were not used for hypothesis testing.",
            "- Cross-phase GSE184114 inference is limited to the shared 24/48 h window.",
        ]
    )
    (output_dir / "EXTERNAL_MODULE_PROJECTION_REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    print("External module projection validation completed", flush=True)
    print(projection_qc.to_string(index=False), flush=True)
    print(f"Globally significant effects: {len(significant)} / {len(effects)}", flush=True)
    if not significant.empty:
        print(
            significant[
                ["dataset", "analysis_family", "module_color", "effect", "f_value", "p_value", "fdr_bh_global"]
            ].to_string(index=False),
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
