from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import pandas as pd


METADATA_DIR = Path("data/processed/metadata")
SAMPLE_MAP_DIR = Path("data/processed/sample_maps")
QC_PATH = Path("results/qc/per_sample_count_qc.tsv")
REPORT_PATH = Path("data/processed/count_matrix_report.tsv")


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def numeric_text(value: str) -> str:
    match = re.search(r"\d+", value or "")
    return match.group(0) if match else ""


def parse_gse124820(row: pd.Series) -> dict[str, str]:
    title = row.get("sample_title", "")
    match = re.match(r"Day(?P<day>\d+)_Rep(?P<rep>\d+)_(?P<short_genotype>.+)", title)
    return {
        "tissue": row.get("characteristic_tissue", ""),
        "organism": row.get("sample_organism_ch1", ""),
        "genotype": row.get("characteristic_genotype", ""),
        "short_genotype": match.group("short_genotype") if match else "",
        "time": match.group("day") if match else "",
        "time_unit": "day",
        "replicate": match.group("rep") if match else numeric_text(row.get("characteristic_replicate", "")),
        "treatment": "room_temperature_deacclimation",
        "experiment": "deacclimation_budbreak",
    }


def parse_gse273240(row: pd.Series) -> dict[str, str]:
    title = row.get("sample_title", "")
    match = re.match(r"(?P<assay>Deac\d+)_(?P<days>\d+)d_(?P<treatment>.+)_rep(?P<rep>\d+)", title)
    return {
        "tissue": row.get("characteristic_tissue", ""),
        "organism": row.get("sample_organism_ch1", ""),
        "genotype": row.get("characteristic_genotype", ""),
        "short_genotype": row.get("characteristic_genotype", ""),
        "time": row.get("characteristic_days", "") or (match.group("days") if match else ""),
        "time_unit": "day",
        "replicate": match.group("rep") if match else "",
        "treatment": row.get("characteristic_treatment", "") or (match.group("treatment") if match else ""),
        "experiment": row.get("characteristic_deacclimation_assay", "") or (match.group("assay") if match else ""),
    }


def parse_gse184114(row: pd.Series) -> dict[str, str]:
    time_text = row.get("characteristic_time", "")
    time = "0" if time_text == "pre-treatment" else numeric_text(time_text)
    return {
        "tissue": row.get("characteristic_tissue", ""),
        "organism": row.get("sample_organism_ch1", ""),
        "genotype": row.get("characteristic_genotype", ""),
        "short_genotype": row.get("characteristic_genotype", ""),
        "time": time,
        "time_unit": "hour",
        "replicate": row.get("characteristic_replicate", ""),
        "treatment": row.get("characteristic_treatment", ""),
        "experiment": row.get("characteristic_experiment", ""),
    }


def parse_gse337039(row: pd.Series) -> dict[str, str]:
    title = row.get("sample_title", "")
    date_match = re.match(r"(?P<month>\d+)_(?P<day>\d+)_(?P<year>\d+)_(?P<rep>TWM\d+)_", title)
    date = ""
    replicate = ""
    if date_match:
        year = "20" + date_match.group("year")
        date = f"{year}-{int(date_match.group('month')):02d}-{int(date_match.group('day')):02d}"
        replicate = date_match.group("rep")
    return {
        "tissue": row.get("characteristic_tissue", ""),
        "organism": row.get("sample_organism_ch1", ""),
        "genotype": "Vitis hybrid cultivar",
        "short_genotype": row.get("sample_organism_ch1", ""),
        "time": date,
        "time_unit": "sampling_date",
        "replicate": replicate,
        "treatment": row.get("characteristic_cell_line", ""),
        "experiment": "chilling_fulfillment",
    }


def parse_gse277812(row: pd.Series) -> dict[str, str]:
    title = row.get("sample_title", "")
    match = re.match(r"T(?P<stage>\d+)_(?P<node>\d+)_(?P<rep>\d+)", title)
    stage_map = {"T1": "BBCH63", "T2": "BBCH77", "T3": "BBCH90"}
    stage_code = f"T{match.group('stage')}" if match else ""
    return {
        "tissue": row.get("characteristic_tissue", ""),
        "organism": row.get("sample_organism_ch1", ""),
        "genotype": row.get("characteristic_cultivar", ""),
        "short_genotype": row.get("characteristic_cultivar", ""),
        "time": row.get("characteristic_phenological_stages", "") or stage_map.get(stage_code, ""),
        "time_unit": "BBCH",
        "replicate": match.group("rep") if match else "",
        "treatment": f"node_{row.get('characteristic_node', '')}",
        "experiment": "bud_fertility",
    }


PARSERS = {
    "GSE124820": parse_gse124820,
    "GSE273240": parse_gse273240,
    "GSE184114": parse_gse184114,
    "GSE337039": parse_gse337039,
    "GSE277812": parse_gse277812,
}


def build_design(accession: str, label: str, sample_map: Path, qc: pd.DataFrame) -> pd.DataFrame:
    metadata = read_tsv(METADATA_DIR / f"{accession}_sample_metadata.tsv")
    mapping = read_tsv(sample_map)
    mapping["sample_accession"] = mapping["sample_accession"].str.replace(r"__dup\d+$", "", regex=True)
    mapping["accession"] = accession
    mapping["label"] = label
    merged = mapping.merge(metadata, on="sample_accession", how="left")
    merged = merged.merge(qc, on=["accession", "label", "sample_accession"], how="left")
    parser = PARSERS[accession]

    parsed_rows: list[dict[str, str]] = []
    for _, row in merged.iterrows():
        parsed = parser(row)
        parsed_rows.append(
            {
                "accession": accession,
                "label": label,
                "sample_accession": row["sample_accession"],
                "matrix_column": row.get("matrix_column", ""),
                "sample_title": row.get("sample_title", ""),
                **parsed,
                "library_size": row.get("library_size", ""),
                "detected_genes": row.get("detected_genes", ""),
            }
        )
    return pd.DataFrame(parsed_rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Build cleaned design tables for count matrices.")
    parser.add_argument("--report", type=Path, default=REPORT_PATH)
    parser.add_argument("--qc", type=Path, default=QC_PATH)
    parser.add_argument("--out-dir", type=Path, default=Path("data/processed/design"))
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    report = read_tsv(args.report)
    qc = read_tsv(args.qc)
    all_tables: list[pd.DataFrame] = []
    for _, record in report.iterrows():
        accession = record["accession"]
        label = record["label"]
        design = build_design(accession, label, Path(record["sample_map"]), qc)
        out_path = args.out_dir / f"{accession}_{label}_design.tsv"
        design.to_csv(out_path, sep="\t", index=False)
        all_tables.append(design)
        print(f"[design] {accession} {label}: {len(design)} samples -> {out_path}")

    combined = pd.concat(all_tables, ignore_index=True)
    combined_path = args.out_dir / "all_core_design.tsv"
    combined.to_csv(combined_path, sep="\t", index=False)
    print(f"[ok] wrote {combined_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
