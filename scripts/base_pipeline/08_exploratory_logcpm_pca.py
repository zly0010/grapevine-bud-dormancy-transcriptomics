from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def read_report(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def logcpm(counts: pd.DataFrame) -> pd.DataFrame:
    lib_size = counts.sum(axis=0)
    cpm = counts.div(lib_size, axis=1) * 1_000_000
    return np.log2(cpm + 1)


def pca_from_matrix(matrix: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
    values = matrix.T.to_numpy(dtype=float)
    values = values - values.mean(axis=0, keepdims=True)
    values = values / np.where(values.std(axis=0, keepdims=True) == 0, 1, values.std(axis=0, keepdims=True))
    _, s, vt = np.linalg.svd(values, full_matrices=False)
    coords = values @ vt[:2].T
    variance = (s**2) / np.sum(s**2)
    pca = pd.DataFrame(coords, index=matrix.columns, columns=["PC1", "PC2"])
    return pca, variance[:2]


def make_plot(data: pd.DataFrame, accession: str, label: str, variance: np.ndarray, out_path: Path) -> None:
    plt.figure(figsize=(7, 5))
    hue = "genotype" if data["genotype"].nunique() > 1 else "treatment"
    style = "treatment" if data["treatment"].nunique() > 1 else None
    sns.scatterplot(data=data, x="PC1", y="PC2", hue=hue, style=style, s=55)
    plt.xlabel(f"PC1 ({variance[0] * 100:.1f}%)")
    plt.ylabel(f"PC2 ({variance[1] * 100:.1f}%)")
    plt.title(f"{accession} {label} logCPM PCA")
    plt.tight_layout()
    plt.savefig(out_path, dpi=220)
    plt.close()


def main() -> int:
    parser = argparse.ArgumentParser(description="Run exploratory logCPM PCA for each count matrix.")
    parser.add_argument("--report", type=Path, default=Path("data/processed/count_matrix_report.tsv"))
    parser.add_argument("--design-qc-dir", type=Path, default=Path("data/processed/design_qc"))
    parser.add_argument("--out-dir", type=Path, default=Path("results/exploratory_pca"))
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    summaries = []
    for _, record in read_report(args.report).iterrows():
        accession = record["accession"]
        label = record["label"]
        counts = pd.read_csv(record["matrix_path"], sep="\t", index_col=0, compression="gzip")
        design = pd.read_csv(args.design_qc_dir / f"{accession}_{label}_design_qc.tsv", sep="\t", dtype=str).fillna("")
        design = design[design["qc_status"] != "fail"].copy()
        keep_samples = [sample for sample in design["sample_accession"] if sample in counts.columns]
        counts = counts[keep_samples]
        keep_genes = (counts >= 10).sum(axis=1) >= max(3, int(counts.shape[1] * 0.10))
        filtered = counts.loc[keep_genes]
        transformed = logcpm(filtered)
        pca, variance = pca_from_matrix(transformed)
        coords = design.merge(pca, left_on="sample_accession", right_index=True, how="inner")
        coords_path = args.out_dir / f"{accession}_{label}_pca_coordinates.tsv"
        coords.to_csv(coords_path, sep="\t", index=False)
        make_plot(coords, accession, label, variance, args.out_dir / f"{accession}_{label}_pca.png")
        summaries.append(
            {
                "accession": accession,
                "label": label,
                "samples_after_qc": str(len(keep_samples)),
                "genes_after_filter": str(filtered.shape[0]),
                "pc1_variance_pct": f"{variance[0] * 100:.2f}",
                "pc2_variance_pct": f"{variance[1] * 100:.2f}",
            }
        )
        print(f"[pca] {accession} {label}: {filtered.shape[0]} genes x {len(keep_samples)} samples")

    pd.DataFrame(summaries).to_csv(args.out_dir / "pca_summary.tsv", sep="\t", index=False)
    print(f"[ok] wrote {args.out_dir / 'pca_summary.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
