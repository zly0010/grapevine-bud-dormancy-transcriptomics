from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pandas as pd


def read_report(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def summarize_series(values: pd.Series) -> dict[str, str]:
    return {
        "min": f"{values.min():.0f}",
        "median": f"{values.median():.0f}",
        "max": f"{values.max():.0f}",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Compute basic count-matrix QC metrics.")
    parser.add_argument("--report", type=Path, default=Path("data/processed/count_matrix_report.tsv"))
    parser.add_argument("--out-dir", type=Path, default=Path("results/qc"))
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    summary_records: list[dict[str, str]] = []
    sample_records: list[dict[str, str]] = []

    for record in read_report(args.report):
        accession = record["accession"]
        label = record["label"]
        counts = pd.read_csv(Path(record["matrix_path"]), sep="\t", index_col=0, compression="gzip")
        library_size = counts.sum(axis=0)
        detected_genes = (counts > 0).sum(axis=0)
        lib_summary = summarize_series(library_size)
        detected_summary = summarize_series(detected_genes)

        summary_records.append(
            {
                "accession": accession,
                "label": label,
                "genes": str(counts.shape[0]),
                "samples": str(counts.shape[1]),
                "library_size_min": lib_summary["min"],
                "library_size_median": lib_summary["median"],
                "library_size_max": lib_summary["max"],
                "detected_genes_min": detected_summary["min"],
                "detected_genes_median": detected_summary["median"],
                "detected_genes_max": detected_summary["max"],
            }
        )

        for sample in counts.columns:
            sample_records.append(
                {
                    "accession": accession,
                    "label": label,
                    "sample_accession": sample,
                    "library_size": str(int(library_size[sample])),
                    "detected_genes": str(int(detected_genes[sample])),
                }
            )

    summary_path = args.out_dir / "count_qc_summary.tsv"
    sample_path = args.out_dir / "per_sample_count_qc.tsv"
    pd.DataFrame(summary_records).to_csv(summary_path, sep="\t", index=False)
    pd.DataFrame(sample_records).to_csv(sample_path, sep="\t", index=False)
    print(f"[ok] wrote {summary_path}")
    print(f"[ok] wrote {sample_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
