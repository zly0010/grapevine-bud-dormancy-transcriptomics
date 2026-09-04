# WARNING: The script defaults do not fully represent the frozen manuscript run.
# The frozen run used 60 trees per target, three seeds, and 2,000 module-label
# permutations. The complete historical command and target-per-module invocation
# were not recovered. See docs/REPRODUCIBILITY_LIMITATIONS.md.

#!/usr/bin/env python3
"""Directed core-module network inference and in-silico perturbation.

The analysis is intentionally restricted to robust blue/turquoise module genes.
It uses the dynGENIE3 ODE formulation with deterministic Extra-Trees models and
an independent sparse linear ODE model.  PlantTFDB is used only to identify
Vitis vinifera transcription factors; archived cross-species motif edges are
not used.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import pickle
import shutil
import subprocess
import sys
import time
import warnings
from dataclasses import dataclass
from pathlib import Path

PROJECT = Path(os.environ.get("GRAPE_ADVANCED_PROJECT", Path.cwd())).resolve()
LEGACY = Path(os.environ.get("GRAPE_DISCOVERY_PROJECT", PROJECT)).resolve()
PROJECT_RESULTS = PROJECT / "results" if (PROJECT / "results").is_dir() else PROJECT / "04_分析结果"
LEGACY_RESULTS = LEGACY / "results_corrected"
PACKAGE_METADATA = PROJECT / "05_数据与元数据" / "validation_metadata" / "metadata"
PACKAGE_DYNGENIE3 = PROJECT / "06_源代码" / "third_party" / "dynGENIE3_official"
OFFICIAL_REPO = Path(
    os.environ.get(
        "DYNGENIE3_REPOSITORY",
        PACKAGE_DYNGENIE3 if PACKAGE_DYNGENIE3.is_dir() else PROJECT / "env" / "dynGENIE3_official",
    )
).resolve()
OFFICIAL = OFFICIAL_REPO / "dynGENIE3_python"
if OFFICIAL.exists() and str(OFFICIAL) not in sys.path:
    sys.path.insert(0, str(OFFICIAL))

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import ElasticNetCV


OUT = Path(
    os.environ.get(
        "GRAPE_PERTURB_OUTPUT",
        PROJECT / "07_复现工作区" / "reproduced_results" / "26_virtual_perturbation",
    )
).resolve()
GIT_EXE = os.environ.get("GIT_EXE", shutil.which("git") or "git")
VST_PATH = LEGACY_RESULTS / "07_wgcna_fixed" / "GSE124820_vst_fixed.txt"
DESIGN_PATH = LEGACY_RESULTS / "01_sample_design_and_qc" / "sample_design.tsv"
MODULE_PATH = LEGACY_RESULTS / "07_wgcna_fixed" / "GSE124820" / "module_assignments.txt"
KME_PATH = LEGACY_RESULTS / "07_wgcna_fixed" / "GSE124820" / "kME_table.txt"
META_PATH = PROJECT_RESULTS / "10_cross_variety_meta" / "04_robust_meta_genes.tsv"
TF_PATH = Path(
    os.environ.get(
        "VITIS_TF_BRIDGE",
        PACKAGE_METADATA / "02_tf_id_bridge.tsv"
        if PACKAGE_METADATA.is_dir()
        else PROJECT / "data" / "metadata" / "02_tf_id_bridge.tsv",
    )
)

MODULES = ("blue", "turquoise")
VARIETIES = ("Vamu", "Vvcs", "Vvri", "Vrip")
SCENARIOS = ("KO", "KD50", "OE1SD")


@dataclass
class ForestModel:
    models: list
    alphas: np.ndarray
    vim: np.ndarray
    genes: list[str]
    regulators: list[str]
    regulator_idx: np.ndarray
    seed: int


@dataclass
class SparseModel:
    intercepts: np.ndarray
    coefs: np.ndarray
    alphas: np.ndarray
    x_mean: np.ndarray
    x_sd: np.ndarray
    genes: list[str]
    regulators: list[str]
    regulator_idx: np.ndarray


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_inputs(target_per_module: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str], list[str]]:
    vst = pd.read_csv(VST_PATH, sep="\t", index_col=0)
    design = pd.read_csv(DESIGN_PATH, sep="\t")
    design = design.loc[design["qc_B"].astype(str).str.lower().isin(["true", "1"])].copy()
    design = design.loc[design.sample_id.isin(vst.index)].copy()
    design["time"] = pd.to_numeric(design["time"], errors="raise")

    modules = pd.read_csv(MODULE_PATH, sep="\t")
    kme = pd.read_csv(KME_PATH, sep="\t")
    meta = pd.read_csv(META_PATH, sep="\t")
    tf = pd.read_csv(TF_PATH, sep="\t")
    tf = tf.loc[tf.high_confidence.astype(str).str.lower().isin(["true", "1"])]
    tf = tf.dropna(subset=["vitvi_id"]).sort_values(["vitvi_id", "bitscore"], ascending=[True, False])
    tf = tf.drop_duplicates("vitvi_id")[["vitvi_id", "tf_family", "pident", "query_coverage", "best_second_ratio"]]

    ann = modules.merge(kme, on=["gene_id", "module_color"], validate="one_to_one")
    ann = ann.merge(
        meta[["gene_id", "estimate", "padj", "robust_meta_gene", "all_varieties_same_sign", "all_loo_same_sign"]],
        on="gene_id", how="left", validate="one_to_one",
    )
    ann = ann.merge(tf, left_on="gene_id", right_on="vitvi_id", how="left", validate="one_to_one")
    ann["kME"] = [abs(row[f"ME{row.module_color}"]) for _, row in ann.iterrows()]
    for col in ("robust_meta_gene", "all_varieties_same_sign", "all_loo_same_sign"):
        ann[col] = ann[col].fillna(False).astype(bool)
    ann["is_tf"] = ann.tf_family.notna()
    ann["selection_score"] = ann.kME * np.log1p(ann.estimate.abs().fillna(0))

    eligible = ann.loc[
        ann.module_color.isin(MODULES)
        & ann.robust_meta_gene
        & ann.all_varieties_same_sign
        & ann.all_loo_same_sign
        & (ann.kME >= 0.75)
    ].copy()
    selected_parts = []
    for module in MODULES:
        block = eligible.loc[eligible.module_color.eq(module)].sort_values(
            ["kME", "selection_score"], ascending=False
        )
        selected_parts.append(block.head(target_per_module))
    selected = pd.concat(selected_parts, ignore_index=True)

    regulators_df = ann.loc[
        ann.module_color.isin(MODULES)
        & ann.robust_meta_gene
        & ann.all_varieties_same_sign
        & ann.all_loo_same_sign
        & ann.is_tf
        & (ann.kME >= 0.70)
    ].sort_values(["module_color", "kME"], ascending=[True, False])
    genes = sorted(set(selected.gene_id).union(regulators_df.gene_id))
    regulators = sorted(set(regulators_df.gene_id).intersection(genes))
    ann = ann.loc[ann.gene_id.isin(genes)].copy().sort_values("gene_id")
    if len(regulators) < 20 or len(genes) < 100:
        raise RuntimeError(f"Insufficient directed-network input: {len(regulators)} regulators, {len(genes)} genes")
    if not set(genes).issubset(vst.columns):
        raise RuntimeError("Selected genes are missing from the discovery VST matrix")

    agg = (
        design[["sample_id", "variety", "time"]]
        .merge(vst[genes], left_on="sample_id", right_index=True, validate="one_to_one")
        .groupby(["variety", "time"], as_index=False)[genes]
        .median()
        .sort_values(["variety", "time"])
    )
    return vst, design, agg, genes, regulators, ann


def trajectories(agg: pd.DataFrame, genes: list[str], include: tuple[str, ...]) -> tuple[list[np.ndarray], list[np.ndarray]]:
    ts, tp = [], []
    for variety in include:
        block = agg.loc[agg.variety.eq(variety)].sort_values("time")
        if len(block) < 5:
            continue
        ts.append(block[genes].to_numpy(dtype=float))
        tp.append(block.time.to_numpy(dtype=float))
    return ts, tp


def estimate_degradation(ts: list[np.ndarray], tp: list[np.ndarray]) -> np.ndarray:
    ngenes = ts[0].shape[1]
    values = np.zeros((len(ts), ngenes), dtype=float)
    for i, (x, t) in enumerate(zip(ts, tp)):
        tmin, tmax = float(np.min(t)), float(np.max(t))
        xmin, xmax = np.maximum(np.min(x, axis=0), 1e-6), np.maximum(np.max(x, axis=0), 1e-6)
        values[i] = np.maximum((np.log(xmax) - np.log(xmin)) / max(abs(tmax - tmin), 1e-6), 1e-6)
    return np.maximum(values.max(axis=0), 1e-4)


def dynamic_training_matrix(
    ts: list[np.ndarray], tp: list[np.ndarray], alphas: np.ndarray, regulator_idx: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    xrows, yrows, present, future, dts = [], [], [], [], []
    for x, t in zip(ts, tp):
        dt = np.diff(t)
        ok = dt > 0
        x0, x1, dt = x[:-1][ok], x[1:][ok], dt[ok]
        xrows.append(x0[:, regulator_idx])
        yrows.append((x1 - x0) / dt[:, None] + x0 * alphas[None, :])
        present.append(x0)
        future.append(x1)
        dts.append(dt)
    return (
        np.vstack(xrows), np.vstack(yrows), np.vstack(present), np.vstack(future), np.concatenate(dts)
    )


def fit_forest(
    ts: list[np.ndarray], tp: list[np.ndarray], genes: list[str], regulators: list[str],
    ntrees: int, seed: int,
) -> ForestModel:
    regulator_idx = np.array([genes.index(g) for g in regulators], dtype=int)
    alphas = estimate_degradation(ts, tp)
    X, Y, _, _, _ = dynamic_training_matrix(ts, tp, alphas, regulator_idx)
    models, vim = [], np.zeros((len(regulators), len(genes)), dtype=float)
    for j, gene in enumerate(genes):
        model = ExtraTreesRegressor(
            n_estimators=ntrees, max_features="sqrt", min_samples_leaf=2,
            bootstrap=False, random_state=seed * 10000 + j, n_jobs=-1,
        )
        model.fit(X, Y[:, j])
        # Single-row prediction with joblib parallelism is much slower than
        # serial prediction; training remains parallel, prediction does not.
        model.n_jobs = 1
        models.append(model)
        imp = model.feature_importances_.astype(float)
        if gene in regulators:
            imp[regulators.index(gene)] = 0.0
        if imp.sum() > 0:
            imp /= imp.sum()
        vim[:, j] = imp
    return ForestModel(models, alphas, vim, genes, regulators, regulator_idx, seed)


def fit_sparse(
    ts: list[np.ndarray], tp: list[np.ndarray], genes: list[str], regulators: list[str], seed: int
) -> SparseModel:
    regulator_idx = np.array([genes.index(g) for g in regulators], dtype=int)
    alphas = estimate_degradation(ts, tp)
    X, Y, _, _, _ = dynamic_training_matrix(ts, tp, alphas, regulator_idx)
    x_mean = X.mean(axis=0)
    x_sd = X.std(axis=0, ddof=1)
    x_sd[x_sd < 1e-8] = 1.0
    Z = (X - x_mean) / x_sd
    coefs = np.zeros((len(regulators), len(genes)), dtype=float)
    intercepts = np.zeros(len(genes), dtype=float)
    for j, gene in enumerate(genes):
        model = ElasticNetCV(
            l1_ratio=[0.5, 0.8, 0.95, 1.0],
            alphas=np.logspace(-4, 0.7, 24), cv=3, max_iter=30000,
            random_state=seed * 10000 + j, selection="cyclic", n_jobs=-1,
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            model.fit(Z, Y[:, j])
        beta = model.coef_ / x_sd
        if gene in regulators:
            beta[regulators.index(gene)] = 0.0
        coefs[:, j] = beta
        intercepts[j] = model.intercept_ - np.dot(model.coef_, x_mean / x_sd)
    return SparseModel(intercepts, coefs, alphas, x_mean, x_sd, genes, regulators, regulator_idx)


def one_step_forest(model: ForestModel, x: np.ndarray, dt: float) -> np.ndarray:
    xr = x[model.regulator_idx].reshape(1, -1)
    production = np.array([m.predict(xr)[0] for m in model.models])
    return x + (production - model.alphas * x) * dt


def one_step_sparse(model: SparseModel, x: np.ndarray, dt: float) -> np.ndarray:
    xr = x[model.regulator_idx]
    production = model.intercepts + xr @ model.coefs
    return x + (production - model.alphas * x) * dt


def one_step_batch(model, x: np.ndarray, dt: float) -> np.ndarray:
    """Advance many perturbation states together without changing the ODE."""
    xr = x[:, model.regulator_idx]
    if isinstance(model, ForestModel):
        production = np.column_stack([m.predict(xr) for m in model.models])
    else:
        production = model.intercepts[None, :] + xr @ model.coefs
    return x + (production - model.alphas[None, :] * x) * dt


def predict_fold(model, agg: pd.DataFrame, genes: list[str], held_out: str, method: str) -> dict:
    block = agg.loc[agg.variety.eq(held_out)].sort_values("time")
    actual_all, pred_all, null_all = [], [], []
    for i in range(len(block) - 1):
        x0 = block.iloc[i][genes].to_numpy(dtype=float)
        x1 = block.iloc[i + 1][genes].to_numpy(dtype=float)
        dt = float(block.iloc[i + 1].time - block.iloc[i].time)
        pred = one_step_forest(model, x0, dt) if method == "forest" else one_step_sparse(model, x0, dt)
        actual_all.append(x1 - x0)
        pred_all.append(pred - x0)
        null_all.append(np.zeros_like(x0))
    actual = np.concatenate(actual_all)
    pred = np.concatenate(pred_all)
    null = np.concatenate(null_all)
    pearson = float(np.corrcoef(actual, pred)[0, 1]) if np.std(pred) > 0 else np.nan
    spear = float(spearmanr(actual, pred, nan_policy="omit").statistic)
    rmse = float(np.sqrt(np.mean((actual - pred) ** 2)))
    null_rmse = float(np.sqrt(np.mean((actual - null) ** 2)))
    return {
        "method": method, "held_out_variety": held_out, "n_transitions": len(block) - 1,
        "n_gene_transition_pairs": int(actual.size), "delta_pearson": pearson,
        "delta_spearman": spear, "rmse": rmse, "persistence_rmse": null_rmse,
        "rmse_improvement_fraction": 1.0 - rmse / null_rmse if null_rmse > 0 else np.nan,
    }


def simulate(model, x0: np.ndarray, regulator: str | None, scenario: str, steps: int, dt: float, bounds):
    path = np.zeros((steps + 1, len(x0)), dtype=float)
    path[0] = x0
    reg_idx = model.genes.index(regulator) if regulator is not None else None
    reg_sd = bounds[2][reg_idx] if reg_idx is not None else 0.0
    clamp_value = None
    if reg_idx is not None:
        if scenario == "KO":
            clamp_value = 0.0
        elif scenario == "KD50":
            clamp_value = 0.5 * x0[reg_idx]
        elif scenario == "OE1SD":
            clamp_value = min(x0[reg_idx] + reg_sd, bounds[1][reg_idx])
        path[0, reg_idx] = clamp_value
    for step in range(1, steps + 1):
        prev = path[step - 1]
        nxt = one_step_forest(model, prev, dt) if isinstance(model, ForestModel) else one_step_sparse(model, prev, dt)
        nxt = np.clip(nxt, bounds[0], bounds[1])
        if reg_idx is not None:
            nxt[reg_idx] = clamp_value
        path[step] = nxt
    return path


def module_scores(path: np.ndarray, ann: pd.DataFrame, genes: list[str], center, scale) -> dict[str, np.ndarray]:
    z = (path - center[None, :]) / scale[None, :]
    out = {}
    for module in MODULES:
        block = ann.loc[ann.module_color.eq(module)]
        idx = np.array([genes.index(g) for g in block.gene_id], dtype=int)
        weights = np.array([float(block.loc[block.gene_id.eq(genes[i]), "kME"].iloc[0]) for i in idx])
        out[module] = np.average(z[:, idx], axis=1, weights=np.abs(weights))
    return out


def perturb_model(
    model, agg: pd.DataFrame, ann: pd.DataFrame, genes: list[str], regulators: list[str],
    model_label: str, seed: int, steps: int = 8, dt: float = 0.25,
) -> tuple[pd.DataFrame, dict[tuple[str, str], np.ndarray]]:
    values = agg[genes].to_numpy(dtype=float)
    center, scale = values.mean(axis=0), values.std(axis=0, ddof=1)
    scale[scale < 1e-8] = 1.0
    low = np.zeros(len(genes), dtype=float)
    high = values.max(axis=0) + 2.0 * values.std(axis=0, ddof=1)
    high = np.maximum(high, values.max(axis=0) + 1e-3)
    bounds = (low, high, scale)
    day0 = agg.loc[agg.time.eq(0), genes].median(axis=0).to_numpy(dtype=float)
    conditions = [(None, "baseline")] + [(r, s) for r in regulators for s in SCENARIOS]
    states = np.repeat(day0[None, :], len(conditions), axis=0)
    clamp_idx = np.full(len(conditions), -1, dtype=int)
    clamp_val = np.full(len(conditions), np.nan, dtype=float)
    for row_i, (regulator, scenario) in enumerate(conditions[1:], 1):
        idx = genes.index(regulator)
        clamp_idx[row_i] = idx
        if scenario == "KO":
            clamp_val[row_i] = 0.0
        elif scenario == "KD50":
            clamp_val[row_i] = 0.5 * day0[idx]
        else:
            clamp_val[row_i] = min(day0[idx] + scale[idx], high[idx])
        states[row_i, idx] = clamp_val[row_i]
    for _ in range(steps):
        states = np.clip(one_step_batch(model, states, dt), low[None, :], high[None, :])
        active = np.flatnonzero(clamp_idx >= 0)
        states[active, clamp_idx[active]] = clamp_val[active]

    baseline_final = states[0]
    rows, gene_effects = [], {}
    for row_i, (regulator, scenario) in enumerate(conditions[1:], 1):
        reg_ann = ann.loc[ann.gene_id.eq(regulator)].iloc[0]
        delta_z = (states[row_i] - baseline_final) / scale
        # Exclude the directly clamped regulator from every outcome statistic;
        # otherwise KO/KD can create a tautological module shift without any
        # propagated network response.
        delta_z[genes.index(regulator)] = 0.0
        module_delta = {}
        for module in MODULES:
            block = ann.loc[ann.module_color.eq(module)]
            idx = np.array([genes.index(g) for g in block.gene_id], dtype=int)
            weights = block.set_index("gene_id").loc[
                [genes[i] for i in idx], "kME"
            ].abs().to_numpy(dtype=float)
            module_delta[module] = float(np.average(delta_z[idx], weights=weights))
        delta_blue = module_delta["blue"]
        delta_turq = module_delta["turquoise"]
        alignment = float((-delta_blue + delta_turq) / math.sqrt(2.0))
        rows.append({
            "model": model_label, "seed": seed, "gene_id": regulator,
            "module_color": reg_ann.module_color, "tf_family": reg_ann.tf_family,
            "scenario": scenario, "delta_blue": delta_blue,
            "delta_turquoise": delta_turq, "trajectory_alignment": alignment,
            "global_effect_rms_z": float(np.sqrt(np.mean(delta_z ** 2))),
            "affected_abs_z_ge_0_25": int(np.sum(np.abs(delta_z) >= 0.25)),
        })
        gene_effects[(regulator, scenario)] = delta_z
    return pd.DataFrame(rows), gene_effects


def permutation_fdr(summary: pd.DataFrame, effects: dict, ann: pd.DataFrame, genes: list[str], nperm: int, seed: int):
    rng = np.random.default_rng(seed)
    module_labels = ann.set_index("gene_id").loc[genes, "module_color"].to_numpy()
    weights = ann.set_index("gene_id").loc[genes, "kME"].abs().to_numpy(dtype=float)
    valid = np.isin(module_labels, MODULES)
    module_labels = module_labels[valid]
    weights = weights[valid]
    pvals = []
    for _, row in summary.iterrows():
        dz = effects[(row.gene_id, row.scenario)][valid]
        observed = abs(float(row.forest_alignment_mean))
        null = np.zeros(nperm)
        for b in range(nperm):
            perm = rng.permutation(module_labels)
            db = float(np.average(dz[perm == "blue"], weights=weights[perm == "blue"]))
            dt = float(np.average(dz[perm == "turquoise"], weights=weights[perm == "turquoise"]))
            null[b] = abs((-db + dt) / math.sqrt(2.0))
        pvals.append((1.0 + float(np.sum(null >= observed))) / (nperm + 1.0))
    summary = summary.copy()
    summary["permutation_p"] = pvals
    order = np.argsort(summary.permutation_p.to_numpy())
    ranked = summary.permutation_p.to_numpy()[order] * len(summary) / np.arange(1, len(summary) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    fdr = np.empty_like(ranked)
    fdr[order] = np.minimum(ranked, 1.0)
    summary["permutation_fdr"] = fdr
    return summary


def edge_tables(main_models: list[ForestModel], sparse: SparseModel, ann: pd.DataFrame) -> pd.DataFrame:
    mean_vim = np.mean([m.vim for m in main_models], axis=0)
    top_rank = np.argsort(np.argsort(-mean_vim, axis=0), axis=0) + 1
    rows = []
    for r, regulator in enumerate(main_models[0].regulators):
        for j, target in enumerate(main_models[0].genes):
            if regulator == target:
                continue
            beta = float(sparse.coefs[r, j])
            if top_rank[r, j] <= 5 or abs(beta) > 1e-10:
                ranks = []
                for model in main_models:
                    ranks.append(int(np.argsort(np.argsort(-model.vim[:, j]))[r] + 1))
                rows.append({
                    "regulator": regulator, "target": target,
                    "regulator_module": ann.set_index("gene_id").loc[regulator, "module_color"],
                    "target_module": ann.set_index("gene_id").loc[target, "module_color"],
                    "mean_forest_importance": float(mean_vim[r, j]),
                    "forest_rank": int(top_rank[r, j]),
                    "top5_seed_count": int(sum(x <= 5 for x in ranks)),
                    "sparse_beta": beta, "sparse_nonzero": abs(beta) > 1e-10,
                    "consensus_edge": bool(top_rank[r, j] <= 5 and sum(x <= 5 for x in ranks) >= 2 and abs(beta) > 1e-10),
                })
    return pd.DataFrame(rows).sort_values(["consensus_edge", "mean_forest_importance"], ascending=[False, False])


def benchmark(target_per_module: int, ntrees: int):
    OUT.mkdir(parents=True, exist_ok=True)
    _, _, agg, genes, regulators, ann = read_inputs(target_per_module)
    ts, tp = trajectories(agg, genes, VARIETIES)
    start = time.time()
    model = fit_forest(ts, tp, genes, regulators, ntrees, seed=1)
    elapsed = time.time() - start
    result = {
        "mode": "benchmark", "genes": len(genes), "regulators": len(regulators),
        "transitions": int(sum(len(x) - 1 for x in ts)), "ntrees": ntrees,
        "elapsed_seconds": elapsed, "estimated_full_seconds": elapsed * 7,
        "finite_vim": bool(np.isfinite(model.vim).all()),
    }
    (OUT / "benchmark.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


def fold_perturb(target_per_module: int, cv_ntrees: int):
    """Test whether perturbation directions survive leave-one-variety refitting."""
    OUT.mkdir(parents=True, exist_ok=True)
    src = OUT / "source_data"
    src.mkdir(exist_ok=True)
    _, _, agg, genes, regulators, ann = read_inputs(target_per_module)
    frames = []
    for fold, held in enumerate(VARIETIES, 1):
        print(f"Leave-one perturbation fold {fold}/4: {held}", flush=True)
        include = tuple(x for x in VARIETIES if x != held)
        fts, ftp = trajectories(agg, genes, include)
        model = fit_forest(fts, ftp, genes, regulators, cv_ntrees, seed=200 + fold)
        frame, _ = perturb_model(model, agg, ann, genes, regulators, "leave_one_forest", 200 + fold)
        frame.insert(0, "held_out_variety", held)
        frames.append(frame)
    all_folds = pd.concat(frames, ignore_index=True)
    all_folds.to_csv(src / "08_leave_one_variety_perturbation.tsv", sep="\t", index=False)
    fold_summary = (
        all_folds.groupby(["gene_id", "module_color", "tf_family", "scenario"], as_index=False)
        .agg(
            fold_alignment_mean=("trajectory_alignment", "mean"),
            fold_alignment_min=("trajectory_alignment", "min"),
            fold_alignment_max=("trajectory_alignment", "max"),
        )
    )
    main = pd.read_csv(src / "04_perturbation_consensus_summary.tsv", sep="\t")
    main = main[["gene_id", "scenario", "forest_alignment_mean", "robust_computational_hit"]]
    fold_summary = fold_summary.merge(main, on=["gene_id", "scenario"], validate="one_to_one")
    fold_summary["fold_direction_consistent"] = (
        np.sign(fold_summary.fold_alignment_min) == np.sign(fold_summary.fold_alignment_max)
    )
    fold_summary["fold_main_direction_concordant"] = (
        np.sign(fold_summary.fold_alignment_mean) == np.sign(fold_summary.forest_alignment_mean)
    )
    fold_summary["leave_one_direction_stable"] = (
        fold_summary.fold_direction_consistent & fold_summary.fold_main_direction_concordant
    )
    fold_summary.to_csv(src / "09_leave_one_variety_perturbation_summary.tsv", sep="\t", index=False)
    print(json.dumps({
        "fold_tests": int(len(fold_summary)),
        "direction_stable": int(fold_summary.leave_one_direction_stable.sum()),
        "robust_hits_direction_stable": int(
            (fold_summary.robust_computational_hit & fold_summary.leave_one_direction_stable).sum()
        ),
    }, indent=2), flush=True)


def full(target_per_module: int, ntrees: int, cv_ntrees: int, nperm: int, resume: bool):
    OUT.mkdir(parents=True, exist_ok=True)
    src = OUT / "source_data"
    models_dir = OUT / "models"
    src.mkdir(exist_ok=True)
    models_dir.mkdir(exist_ok=True)
    _, design, agg, genes, regulators, ann = read_inputs(target_per_module)
    ann.to_csv(src / "01_network_gene_inventory.tsv", sep="\t", index=False)
    agg[["variety", "time"] + genes].to_csv(src / "02_aggregated_time_series.tsv", sep="\t", index=False)

    ts, tp = trajectories(agg, genes, VARIETIES)
    forest_models = []
    perturb_frames, effects_by_seed = [], []
    for seed in (1, 2, 3):
        model_path = models_dir / f"dynGENIE3_ET_seed{seed}.pkl"
        if resume and model_path.exists():
            print(f"Loading forest seed {seed}/3", flush=True)
            with model_path.open("rb") as handle:
                model = pickle.load(handle)
        else:
            print(f"Fitting forest seed {seed}/3", flush=True)
            model = fit_forest(ts, tp, genes, regulators, ntrees, seed)
        forest_models.append(model)
        frame, effects = perturb_model(model, agg, ann, genes, regulators, "dynGENIE3_ET", seed)
        perturb_frames.append(frame)
        effects_by_seed.append(effects)
        with model_path.open("wb") as handle:
            pickle.dump(model, handle, protocol=pickle.HIGHEST_PROTOCOL)

    sparse_path = models_dir / "sparse_dynamic.pkl"
    if resume and sparse_path.exists():
        print("Loading main sparse dynamic model", flush=True)
        with sparse_path.open("rb") as handle:
            sparse = pickle.load(handle)
    else:
        print("Fitting main sparse dynamic model", flush=True)
        sparse = fit_sparse(ts, tp, genes, regulators, seed=1)
    sparse_frame, sparse_effects = perturb_model(sparse, agg, ann, genes, regulators, "sparse_dynamic", 1)
    with sparse_path.open("wb") as handle:
        pickle.dump(sparse, handle, protocol=pickle.HIGHEST_PROTOCOL)

    perturb_all = pd.concat(perturb_frames + [sparse_frame], ignore_index=True)
    perturb_all.to_csv(src / "03_all_perturbation_runs.tsv", sep="\t", index=False)

    forest_summary = (
        pd.concat(perturb_frames, ignore_index=True)
        .groupby(["gene_id", "module_color", "tf_family", "scenario"], as_index=False)
        .agg(
            forest_alignment_mean=("trajectory_alignment", "mean"),
            forest_alignment_min=("trajectory_alignment", "min"),
            forest_alignment_max=("trajectory_alignment", "max"),
            delta_blue_mean=("delta_blue", "mean"),
            delta_turquoise_mean=("delta_turquoise", "mean"),
            global_effect_rms_z=("global_effect_rms_z", "mean"),
            affected_genes=("affected_abs_z_ge_0_25", "mean"),
        )
    )
    sparse_small = sparse_frame[["gene_id", "scenario", "trajectory_alignment", "delta_blue", "delta_turquoise"]].rename(
        columns={
            "trajectory_alignment": "sparse_alignment", "delta_blue": "sparse_delta_blue",
            "delta_turquoise": "sparse_delta_turquoise",
        }
    )
    summary = forest_summary.merge(sparse_small, on=["gene_id", "scenario"], validate="one_to_one")
    summary["seed_direction_consistent"] = np.sign(summary.forest_alignment_min) == np.sign(summary.forest_alignment_max)
    summary["model_direction_concordant"] = np.sign(summary.forest_alignment_mean) == np.sign(summary.sparse_alignment)
    mean_effects = {
        key: np.mean([effect[key] for effect in effects_by_seed], axis=0)
        for key in effects_by_seed[0]
    }
    summary = permutation_fdr(summary, mean_effects, ann, genes, nperm, seed=20260806)
    summary["robust_computational_hit"] = (
        summary.seed_direction_consistent & summary.model_direction_concordant
        & (summary.permutation_fdr < 0.05)
        & (summary.global_effect_rms_z >= summary.global_effect_rms_z.quantile(0.75))
    )
    summary = summary.sort_values(
        ["robust_computational_hit", "forest_alignment_mean"], ascending=[False, False]
    )
    summary.to_csv(src / "04_perturbation_consensus_summary.tsv", sep="\t", index=False)

    edges = edge_tables(forest_models, sparse, ann)
    edges.to_csv(src / "05_directed_consensus_edges.tsv", sep="\t", index=False)

    validation_path = src / "06_leave_one_variety_prediction.tsv"
    strength_path = src / "07_leave_one_variety_regulator_strength.tsv"
    if resume and validation_path.exists() and strength_path.exists():
        print("Reusing completed leave-one-variety predictive validation", flush=True)
        validation_df = pd.read_csv(validation_path, sep="\t")
    else:
        validation = []
        fold_edge_rows = []
        for fold, held in enumerate(VARIETIES, 1):
            print(f"Leave-one-variety fold {fold}/4: {held}", flush=True)
            include = tuple(x for x in VARIETIES if x != held)
            fts, ftp = trajectories(agg, genes, include)
            fmodel = fit_forest(fts, ftp, genes, regulators, cv_ntrees, seed=100 + fold)
            smodel = fit_sparse(fts, ftp, genes, regulators, seed=100 + fold)
            validation.append(predict_fold(fmodel, agg, genes, held, "forest"))
            validation.append(predict_fold(smodel, agg, genes, held, "sparse"))
            for r, regulator in enumerate(regulators):
                fold_edge_rows.append({
                    "held_out_variety": held, "regulator": regulator,
                    "outgoing_strength": float(fmodel.vim[r].sum()),
                })
        validation_df = pd.DataFrame(validation)
        validation_df.to_csv(validation_path, sep="\t", index=False)
        pd.DataFrame(fold_edge_rows).to_csv(strength_path, sep="\t", index=False)

    method_corr = float(spearmanr(summary.forest_alignment_mean, summary.sparse_alignment).statistic)
    run_info = {
        "analysis_date": "2026-08-06", "genes": len(genes), "regulators": len(regulators),
        "module_gene_counts": ann.groupby("module_color").size().to_dict(),
        "module_regulator_counts": ann.loc[ann.gene_id.isin(regulators)].groupby("module_color").size().to_dict(),
        "aggregated_time_points": int(len(agg)), "transitions": int(sum(len(x) - 1 for x in ts)),
        "forest_seeds": [1, 2, 3], "ntrees": ntrees, "cv_ntrees": cv_ntrees,
        "permutations_per_intervention": nperm,
        "forest_sparse_alignment_spearman": method_corr,
        "robust_hit_count": int(summary.robust_computational_hit.sum()),
        "consensus_edge_count": int(edges.consensus_edge.sum()),
        "input_sha256": {str(p): sha256(p) for p in [VST_PATH, DESIGN_PATH, MODULE_PATH, KME_PATH, META_PATH, TF_PATH]},
        "software": {"python": sys.version, "numpy": np.__version__, "pandas": pd.__version__},
        "dynGENIE3_official_commit": subprocess.run(
            [
                str(GIT_EXE), "-c", f"safe.directory={OFFICIAL_REPO}",
                "-C", str(OFFICIAL_REPO), "rev-parse", "HEAD",
            ],
            capture_output=True, text=True, check=False,
        ).stdout.strip(),
    }
    (OUT / "RUN_INFO.json").write_text(json.dumps(run_info, ensure_ascii=False, indent=2), encoding="utf-8")
    report = [
        "# Virtual perturbation analysis", "",
        f"- Directed core network: {len(genes)} genes and {len(regulators)} high-confidence Vitis TF regulators.",
        f"- Discovery trajectories: {len(agg)} variety-time centroids and {run_info['transitions']} adjacent transitions.",
        f"- Robust computational hits: {run_info['robust_hit_count']} of {len(summary)} TF-scenario tests.",
        f"- Consensus directed edges: {run_info['consensus_edge_count']}.",
        f"- Forest/sparse perturbation-rank Spearman correlation: {method_corr:.3f}.",
        "", "These are model-based perturbation results, not experimental causal validation.",
    ]
    (OUT / "SUMMARY.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(json.dumps(run_info, ensure_ascii=False, indent=2))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["benchmark", "full", "fold-perturb"], default="benchmark")
    parser.add_argument("--target-per-module", type=int, default=160)
    parser.add_argument("--ntrees", type=int, default=100)
    parser.add_argument("--cv-ntrees", type=int, default=30)
    parser.add_argument("--nperm", type=int, default=1000)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if args.mode == "benchmark":
        benchmark(min(args.target_per_module, 40), min(args.ntrees, 40))
    elif args.mode == "fold-perturb":
        fold_perturb(args.target_per_module, args.cv_ntrees)
    else:
        full(args.target_per_module, args.ntrees, args.cv_ntrees, args.nperm, args.resume)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
