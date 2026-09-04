from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def read_designs(design_dir: Path) -> pd.DataFrame:
    tables = []
    for path in sorted(design_dir.glob("*_design.tsv")):
        if path.name == "all_core_design.tsv":
            continue
        tables.append(pd.read_csv(path, sep="\t", dtype=str).fillna(""))
    return pd.concat(tables, ignore_index=True)


def add_flags(design: pd.DataFrame) -> pd.DataFrame:
    design = design.copy()
    design["library_size"] = pd.to_numeric(design["library_size"], errors="coerce")
    design["detected_genes"] = pd.to_numeric(design["detected_genes"], errors="coerce")
    design["qc_status"] = "pass"
    design["qc_reason"] = ""

    for (accession, label), idx in design.groupby(["accession", "label"]).groups.items():
        group = design.loc[idx]
        lib_median = group["library_size"].median()
        detected_median = group["detected_genes"].median()
        fail_mask = (group["library_size"] < lib_median * 0.25) | (group["detected_genes"] < detected_median * 0.70)
        watch_mask = (group["library_size"] < lib_median * 0.50) | (group["detected_genes"] < detected_median * 0.85)

        for row_idx in group.index[watch_mask]:
            design.loc[row_idx, "qc_status"] = "watch"
        for row_idx in group.index[fail_mask]:
            design.loc[row_idx, "qc_status"] = "fail"

        for row_idx in group.index[watch_mask | fail_mask]:
            reasons = []
            if design.loc[row_idx, "library_size"] < lib_median * 0.25:
                reasons.append("library_size_lt_25pct_median")
            elif design.loc[row_idx, "library_size"] < lib_median * 0.50:
                reasons.append("library_size_lt_50pct_median")
            if design.loc[row_idx, "detected_genes"] < detected_median * 0.70:
                reasons.append("detected_genes_lt_70pct_median")
            elif design.loc[row_idx, "detected_genes"] < detected_median * 0.85:
                reasons.append("detected_genes_lt_85pct_median")
            design.loc[row_idx, "qc_reason"] = ";".join(reasons)

    return design


def plot_metric(data: pd.DataFrame, accession: str, label: str, metric: str, out_dir: Path) -> None:
    group = data[(data["accession"] == accession) & (data["label"] == label)].copy()
    group["sample_order"] = range(1, len(group) + 1)
    plt.figure(figsize=(max(8, len(group) * 0.08), 4))
    sns.scatterplot(
        data=group,
        x="sample_order",
        y=metric,
        hue="qc_status",
        palette={"pass": "#2f6f4e", "watch": "#c58b00", "fail": "#b94040"},
        s=28,
    )
    plt.yscale("log" if metric == "library_size" else "linear")
    plt.xlabel("sample order")
    plt.ylabel(metric.replace("_", " "))
    plt.title(f"{accession} {label} {metric}")
    plt.tight_layout()
    safe_label = label.replace("/", "_")
    plt.savefig(out_dir / f"{accession}_{safe_label}_{metric}.png", dpi=180)
    plt.close()


def main() -> int:
    parser = argparse.ArgumentParser(description="Flag low-quality samples and draw QC plots.")
    parser.add_argument("--design-dir", type=Path, default=Path("data/processed/design"))
    parser.add_argument("--out-dir", type=Path, default=Path("results/qc"))
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    flagged = add_flags(read_designs(args.design_dir))
    flagged_path = args.out_dir / "sample_qc_flags.tsv"
    flagged.to_csv(flagged_path, sep="\t", index=False)
    design_qc_dir = Path("data/processed/design_qc")
    design_qc_dir.mkdir(parents=True, exist_ok=True)
    for (accession, label), group in flagged.groupby(["accession", "label"]):
        group.to_csv(design_qc_dir / f"{accession}_{label}_design_qc.tsv", sep="\t", index=False)
    flagged.to_csv(design_qc_dir / "all_core_design_qc.tsv", sep="\t", index=False)

    summary = (
        flagged.groupby(["accession", "label", "qc_status"], dropna=False)
        .size()
        .reset_index(name="n_samples")
        .sort_values(["accession", "label", "qc_status"])
    )
    summary_path = args.out_dir / "sample_qc_flag_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)

    for accession, label in flagged[["accession", "label"]].drop_duplicates().itertuples(index=False):
        plot_metric(flagged, accession, label, "library_size", args.out_dir)
        plot_metric(flagged, accession, label, "detected_genes", args.out_dir)

    xlsx_path = args.out_dir / "qc_design_summary.xlsx"
    with pd.ExcelWriter(xlsx_path) as writer:
        flagged.to_excel(writer, sheet_name="sample_qc_flags", index=False)
        summary.to_excel(writer, sheet_name="flag_summary", index=False)

    print(f"[ok] wrote {flagged_path}")
    print(f"[ok] wrote {summary_path}")
    print(f"[ok] wrote {xlsx_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
