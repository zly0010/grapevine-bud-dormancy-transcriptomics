# WARNING: This acquisition script also contains a historical GSE337039
# processed-matrix branch. It is not the final manuscript validation input.
# Use the documented GSE337039 raw-read branch for final independent validation.

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import defaultdict, deque
from pathlib import Path

import pandas as pd


RAW_SUPPL = Path("data/raw/geo_supplementary")
METADATA_DIR = Path("data/processed/metadata")
COUNTS_DIR = Path("data/processed/counts")
MAP_DIR = Path("data/processed/sample_maps")


def read_metadata(accession: str) -> list[dict[str, str]]:
    path = METADATA_DIR / f"{accession}_sample_metadata.tsv"
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def normalized_title(title: str) -> str:
    return re.sub(r"^Sample\s+\d+_", "", title or "")


def description_library_name(description: str) -> str:
    return re.sub(r"^Library name:\s*", "", description or "")


def build_alias_queues(rows: list[dict[str, str]], experiment: str | None = None) -> dict[str, deque[str]]:
    queues: dict[str, deque[str]] = defaultdict(deque)
    for row in rows:
        if experiment and row.get("characteristic_experiment") != experiment:
            continue
        sample = row["sample_accession"]
        aliases = {
            sample,
            row.get("sample_title", ""),
            normalized_title(row.get("sample_title", "")),
            row.get("sample_description", ""),
            description_library_name(row.get("sample_description", "")),
        }
        for alias in aliases:
            if alias:
                queues[alias].append(sample)
    return queues


def map_columns(columns: list[str], queues: dict[str, deque[str]]) -> tuple[list[str], list[dict[str, str]]]:
    mapped: list[str] = []
    records: list[dict[str, str]] = []
    used_counts: dict[str, int] = defaultdict(int)
    for column in columns:
        if column in queues and queues[column]:
            sample = queues[column].popleft()
            matched = "true"
        else:
            sample = column
            matched = "false"

        used_counts[sample] += 1
        unique_sample = sample if used_counts[sample] == 1 else f"{sample}__dup{used_counts[sample]}"
        mapped.append(unique_sample)
        records.append({"matrix_column": column, "sample_accession": unique_sample, "matched": matched})
    return mapped, records


def write_sample_map(accession: str, label: str, records: list[dict[str, str]]) -> Path:
    MAP_DIR.mkdir(parents=True, exist_ok=True)
    path = MAP_DIR / f"{accession}_{label}_sample_map.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["matrix_column", "sample_accession", "matched"])
        writer.writeheader()
        writer.writerows(records)
    return path


def combine_gse124820() -> tuple[Path, dict[str, str]]:
    accession = "GSE124820"
    rows = read_metadata(accession)
    expected_samples = {row["sample_accession"] for row in rows}
    raw_dir = RAW_SUPPL / accession / "GSE124820_RAW"
    count_files = sorted(raw_dir.glob("*.counts.txt.gz"))

    series: list[pd.Series] = []
    sample_records: list[dict[str, str]] = []
    for path in count_files:
        sample = path.name.split("_", 1)[0]
        frame = pd.read_csv(path, sep="\t", header=None, names=["gene_id", sample], compression="gzip")
        series.append(frame.set_index("gene_id")[sample])
        sample_records.append({"matrix_column": path.name, "sample_accession": sample, "matched": str(sample in expected_samples).lower()})

    matrix = pd.concat(series, axis=1).fillna(0).astype("int64")
    COUNTS_DIR.mkdir(parents=True, exist_ok=True)
    out_path = COUNTS_DIR / f"{accession}_raw_counts.tsv.gz"
    matrix.to_csv(out_path, sep="\t", compression="gzip")
    map_path = write_sample_map(accession, "raw_counts", sample_records)
    return out_path, {
        "accession": accession,
        "label": "raw_counts",
        "genes": str(matrix.shape[0]),
        "samples": str(matrix.shape[1]),
        "matched_samples": str(sum(record["matched"] == "true" for record in sample_records)),
        "unmatched_samples": str(sum(record["matched"] != "true" for record in sample_records)),
        "sample_map": str(map_path),
        "matrix_path": str(out_path),
    }


def read_matrix(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n")
        second = handle.readline().rstrip("\n")
    delimiter = "," if header.count(",") > header.count("\t") else "\t"
    header_fields = next(csv.reader([header], delimiter=delimiter))
    second_fields = next(csv.reader([second], delimiter=delimiter))

    if len(second_fields) == len(header_fields) + 1:
        frame = pd.read_csv(path, sep=delimiter, compression="gzip", header=None, skiprows=1, dtype={0: str})
        frame.columns = ["gene_id"] + header_fields
        frame = frame.set_index("gene_id")
        frame.index.name = "gene_id"
        return frame

    frame = pd.read_csv(path, sep=delimiter, compression="gzip", dtype={0: str})
    first_col = str(frame.columns[0])
    if first_col.startswith("Unnamed") or first_col in {"", "Data.Matrix"}:
        frame = frame.rename(columns={frame.columns[0]: "gene_id"})
    else:
        first_value = str(frame.iloc[0, 0])
        if first_value.startswith("Vitvi"):
            frame = pd.read_csv(path, sep=delimiter, compression="gzip", header=None)
            frame = frame.rename(columns={0: "gene_id"})
        elif str(frame.index[0]).startswith(("Vitvi", "ENS")):
            frame.index.name = "gene_id"
            return frame
    frame = frame.set_index("gene_id")
    frame.index.name = "gene_id"
    return frame


def prepare_matrix(accession: str, label: str, path: Path, experiment: str | None = None) -> tuple[Path, dict[str, str]]:
    rows = read_metadata(accession)
    matrix = read_matrix(path)
    original_columns = [str(column) for column in matrix.columns]
    mapped_columns, sample_records = map_columns(original_columns, build_alias_queues(rows, experiment))
    matrix.columns = mapped_columns

    COUNTS_DIR.mkdir(parents=True, exist_ok=True)
    out_path = COUNTS_DIR / f"{accession}_{label}_raw_counts.tsv.gz"
    matrix.to_csv(out_path, sep="\t", compression="gzip")
    map_path = write_sample_map(accession, label, sample_records)
    return out_path, {
        "accession": accession,
        "label": label,
        "genes": str(matrix.shape[0]),
        "samples": str(matrix.shape[1]),
        "matched_samples": str(sum(record["matched"] == "true" for record in sample_records)),
        "unmatched_samples": str(sum(record["matched"] != "true" for record in sample_records)),
        "sample_map": str(map_path),
        "matrix_path": str(out_path),
    }


def write_report(records: list[dict[str, str]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["accession", "label", "genes", "samples", "matched_samples", "unmatched_samples", "sample_map", "matrix_path"]
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare raw count matrices with GSM sample IDs.")
    parser.add_argument("--report", type=Path, default=Path("data/processed/count_matrix_report.tsv"))
    args = parser.parse_args()

    records: list[dict[str, str]] = []
    _, record = combine_gse124820()
    records.append(record)

    jobs = [
        (
            "GSE273240",
            "tetralone_ABA",
            RAW_SUPPL / "GSE273240" / "GSE273240_tetralone_ABA_gene_count.csv.gz",
            None,
        ),
        (
            "GSE184114",
            "Acclimation",
            RAW_SUPPL / "GSE184114" / "GSE184114_Acclimation_raw_gene_count_matrix.txt.gz",
            "Acclimation",
        ),
        (
            "GSE184114",
            "Deacclimation",
            RAW_SUPPL / "GSE184114" / "GSE184114_Deacclimation_raw_gene_count_matrix.txt.gz",
            "Deacclimation",
        ),
        (
            "GSE337039",
            "processed_read_counts",
            RAW_SUPPL / "GSE337039" / "GSE337039_geo_processed_read_counts.csv.gz",
            None,
        ),
        (
            "GSE277812",
            "raw_counts",
            RAW_SUPPL / "GSE277812" / "GSE277812_raw_counts.txt.gz",
            None,
        ),
    ]
    for accession, label, path, experiment in jobs:
        _, record = prepare_matrix(accession, label, path, experiment)
        records.append(record)

    write_report(records, args.report)
    for record in records:
        print(
            f"[counts] {record['accession']} {record['label']}: "
            f"{record['genes']} genes x {record['samples']} samples; "
            f"matched {record['matched_samples']}, unmatched {record['unmatched_samples']}"
        )
    print(f"[ok] wrote {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
