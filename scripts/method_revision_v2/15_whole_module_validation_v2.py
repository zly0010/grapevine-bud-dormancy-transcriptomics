#!/usr/bin/env python3
"""Versioned external whole-module random-set validation.

This script preserves the original fixed reference-module membership and
dataset-specific available-gene background, increases the random-set count,
and uses only pre-defined sqrt(n) weights for signed Stouffer synthesis.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.stats as st
from statsmodels.stats.multitest import multipletests


TARGET_MODULES = ("blue", "turquoise", "brown")
SEED = 20260718


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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


def random_null(
    z_all: np.ndarray,
    set_size: int,
    permutations: int,
    rng: np.random.Generator,
    batch_size: int = 64,
) -> np.ndarray:
    """Generate size-matched null sets without replacement, in memory-safe batches."""
    n_samples, n_genes = z_all.shape
    if set_size > n_genes:
        raise ValueError("random-set size exceeds available-gene background")
    null = np.empty(permutations, dtype=float)
    for start in range(0, permutations, batch_size):
        stop = min(start + batch_size, permutations)
        draws = np.vstack(
            [rng.choice(n_genes, size=set_size, replace=False) for _ in range(stop - start)]
        )
        row_sums = np.empty((n_samples, stop - start), dtype=float)
        for sample_index in range(n_samples):
            row_sums[sample_index, :] = z_all[sample_index, draws].sum(axis=1)
        total_correlation = np.square(row_sums).sum(axis=0) / (n_samples - 1)
        null[start:stop] = (total_correlation - set_size) / (set_size * (set_size - 1))
    return null


def signed_empirical_z(empirical_p: float, direction: int) -> float:
    magnitude = abs(float(st.norm.isf(empirical_p)))
    return float(direction * magnitude)


def build_stouffer(summary: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    full_rows: list[dict[str, object]] = []
    loo_rows: list[dict[str, object]] = []
    datasets = sorted(summary.dataset.unique())
    for module in TARGET_MODULES:
        module_data = summary.loc[summary.module_color == module].copy()
        if sorted(module_data.dataset.tolist()) != datasets:
            raise RuntimeError(f"{module}: incomplete dataset coverage")

        def combine(kept: pd.DataFrame) -> float:
            z = kept.signed_empirical_z.to_numpy(dtype=float)
            w = kept.stouffer_weight.to_numpy(dtype=float)
            return float(np.sum(w * z) / np.sqrt(np.square(w).sum()))

        full_z = combine(module_data)
        full_rows.append(
            {
                "module_color": module,
                "datasets_combined": len(module_data),
                "dataset_list": ";".join(module_data.dataset.astype(str)),
                "weight_formula": "sqrt(effective_sample_size)",
                "direction_in_weight": False,
                "combined_z": full_z,
                "combined_one_sided_p": float(st.norm.sf(full_z)),
                "all_dataset_z_positive": bool(np.all(module_data.signed_empirical_z > 0)),
            }
        )
        for omitted in datasets:
            kept = module_data.loc[module_data.dataset != omitted]
            z = combine(kept)
            loo_rows.append(
                {
                    "module_color": module,
                    "omitted_dataset": omitted,
                    "datasets_retained": len(kept),
                    "retained_dataset_list": ";".join(kept.dataset.astype(str)),
                    "weight_formula": "sqrt(effective_sample_size)",
                    "direction_in_weight": False,
                    "combined_z": z,
                    "combined_one_sided_p": float(st.norm.sf(z)),
                    "all_retained_z_positive": bool(np.all(kept.signed_empirical_z > 0)),
                }
            )
    full = pd.DataFrame(full_rows)
    full["fdr_bh_across_modules"] = multipletests(
        full.combined_one_sided_p, method="fdr_bh"
    )[1]
    loo = pd.DataFrame(loo_rows)
    return full, loo


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("source_project", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--permutations", type=int, default=10_000)
    parser.add_argument("--gse337039-vst", type=Path, required=True)
    parser.add_argument("--gse337039-design", type=Path, required=True)
    args = parser.parse_args()
    if args.permutations < 5_000:
        raise ValueError("Use at least 5,000 random sets")
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise RuntimeError(f"Refusing to overwrite non-empty output directory: {args.output_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    root = args.source_project
    module_path = root / "results_corrected/07_wgcna_fixed/GSE124820/module_assignments.txt"
    reference_path = root / "results_corrected/07_wgcna_fixed/GSE124820_vst_fixed.txt"
    modules = pd.read_csv(module_path, sep="\t", dtype=str)
    if len(modules) != 5_000 or modules.gene_id.duplicated().any():
        raise RuntimeError("reference module assignment is not the expected unique 5,000-gene set")

    datasets: dict[str, tuple[Path, Path]] = {
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
        "GSE337039": (args.gse337039_vst, args.gse337039_design),
    }
    all_input_paths = [module_path, reference_path] + [p for pair in datasets.values() for p in pair]
    missing = [p for p in all_input_paths if not p.exists()]
    if missing:
        raise FileNotFoundError(missing)

    rng = np.random.default_rng(SEED)
    summary_rows: list[dict[str, object]] = []
    null_rows: list[pd.DataFrame] = []
    all_genes = modules.gene_id.astype(str).tolist()
    for dataset, (vst_path, design_path) in datasets.items():
        print(f"Random-set validation: {dataset}", flush=True)
        vst = read_vst(vst_path, all_genes)
        design = pd.read_csv(design_path, sep="\t", dtype=str)
        if dataset == "GSE337039" and "run_accession" in design.columns:
            sample_ids = design.run_accession.astype(str)
        elif "sample_id" in design.columns:
            sample_ids = design.sample_id.astype(str)
        elif "run_accession" in design.columns:
            sample_ids = design.run_accession.astype(str)
        else:
            raise RuntimeError(f"{dataset}: no sample identifier column")
        if set(vst.index) != set(sample_ids):
            raise RuntimeError(f"{dataset}: VST/design sample mismatch")
        vst = vst.loc[sample_ids]
        z_all, genes_all = standardize(vst)
        index = {gene: i for i, gene in enumerate(genes_all)}
        for module in TARGET_MODULES:
            reference_genes = modules.loc[modules.module_color == module, "gene_id"].tolist()
            usable = [gene for gene in reference_genes if gene in index]
            indices = np.array([index[gene] for gene in usable], dtype=int)
            observed = mean_pairwise_correlation(z_all[:, indices])
            null = random_null(z_all, len(usable), args.permutations, rng)
            b = int(np.sum(null >= observed))
            empirical_p = (b + 1) / (args.permutations + 1)
            null_mean = float(null.mean())
            direction = 1 if observed >= null_mean else -1
            z_value = signed_empirical_z(empirical_p, direction)
            summary_rows.append(
                {
                    "dataset": dataset,
                    "module_color": module,
                    "effective_sample_size": len(vst),
                    "reference_module_genes": len(reference_genes),
                    "usable_module_genes": len(usable),
                    "gene_coverage": len(usable) / len(reference_genes),
                    "available_gene_background_size": len(genes_all),
                    "matching_rule": "same_size_without_replacement_from_dataset_available_gene_background",
                    "reference_membership_refit": False,
                    "observed_mean_pairwise_correlation": observed,
                    "null_mean": null_mean,
                    "null_sd": float(null.std(ddof=1)),
                    "null_median": float(np.median(null)),
                    "null_q025": float(np.quantile(null, 0.025)),
                    "null_q975": float(np.quantile(null, 0.975)),
                    "null_min": float(null.min()),
                    "null_max": float(null.max()),
                    "extreme_random_sets_b": b,
                    "B": args.permutations,
                    "empirical_p": empirical_p,
                    "effect_direction": "observed_above_null" if direction > 0 else "observed_below_null",
                    "signed_empirical_z": z_value,
                    "stouffer_weight": math.sqrt(len(vst)),
                    "random_seed": SEED,
                }
            )
            null_rows.append(
                pd.DataFrame(
                    {
                        "dataset": dataset,
                        "module_color": module,
                        "random_set_id": np.arange(1, args.permutations + 1),
                        "null_mean_pairwise_correlation": null,
                    }
                )
            )
            print(
                f"  {module}: observed={observed:.6f}, b={b}, p={empirical_p:.8g}, signed Z={z_value:.4f}",
                flush=True,
            )

    summary = pd.DataFrame(summary_rows)
    null_full = pd.concat(null_rows, ignore_index=True)
    meta, loo = build_stouffer(summary)
    summary.to_csv(args.output_dir / "module_randomset_validation_v2.csv", index=False)
    null_full.to_csv(args.output_dir / "module_randomset_null_full_v2.csv", index=False)
    meta.to_csv(args.output_dir / "stouffer_meta_v2.csv", index=False)
    loo.to_csv(args.output_dir / "stouffer_leave_one_dataset_out_v2.csv", index=False)

    provenance = {
        "script": Path(__file__).name,
        "seed": SEED,
        "B": args.permutations,
        "target_modules": list(TARGET_MODULES),
        "weight_formula": "sqrt(effective_sample_size)",
        "direction_in_weight": False,
        "empirical_p_formula": "(b + 1) / (B + 1)",
        "input_sha256": {str(path): sha256(path) for path in all_input_paths},
    }
    (args.output_dir / "RUN_INFO.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print("WHOLE_MODULE_VALIDATION_V2_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
