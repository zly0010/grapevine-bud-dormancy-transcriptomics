#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
from datetime import datetime
from pathlib import Path


def read_rows(path: Path, delimiter: str) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_tsv(rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_date(value: str) -> str:
    return datetime.strptime(value, "%m_%d_%y").date().isoformat()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runinfo", type=Path, required=True)
    parser.add_argument("--design", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    runinfo = read_rows(args.runinfo, ",")
    design = read_rows(args.design, "\t")
    if len(runinfo) != 60 or len(design) != 60:
        raise RuntimeError(f"Expected 60 runinfo and design rows; got {len(runinfo)} and {len(design)}")

    design_by_sample = {row["sample_id"]: row for row in design}
    if len(design_by_sample) != 60:
        raise RuntimeError("Design sample_id values are not unique")

    missing = sorted({row["SampleName"] for row in runinfo} - set(design_by_sample))
    if missing:
        raise RuntimeError("RunInfo samples missing from design: " + ", ".join(missing[:5]))

    cultivar_name = {"M": "Marquette", "B": "Brianna"}
    merged: list[dict[str, object]] = []
    for run in runinfo:
        d = design_by_sample[run["SampleName"]]
        date_iso = parse_date(d["date"])
        merged.append(
            {
                "run_accession": run["Run"],
                "experiment_accession": run["Experiment"],
                "biosample": run["BioSample"],
                "sample_id": d["sample_id"],
                "cultivar": cultivar_name[d["cultivar"]],
                "cultivar_code": d["cultivar"],
                "condition": d["treatment"],
                "condition_code": d["chilling"],
                "collection_date": date_iso,
                "replicate": int(d["replicate"]),
                "library_layout": run["LibraryLayout"],
                "read_length": round(float(run["avgLength"])),
                "spots": int(run["spots"]),
                "bases": int(run["bases"]),
                "normalized_sra_size_mb": float(run["size_MB"]),
                "download_url": run["download_path"],
            }
        )

    cultivar_conditions = sorted({(str(row["cultivar"]), str(row["condition"])) for row in merged})
    for cultivar, condition in cultivar_conditions:
        dates = sorted(
            {
                str(row["collection_date"])
                for row in merged
                if row["cultivar"] == cultivar and row["condition"] == condition
            }
        )
        date_to_index = {date: index for index, date in enumerate(dates)}
        for row in merged:
            if row["cultivar"] == cultivar and row["condition"] == condition:
                row["time_index"] = date_to_index[str(row["collection_date"])]

    merged.sort(key=lambda row: (str(row["cultivar"]), str(row["condition"]), int(row["time_index"]), int(row["replicate"])))

    cells: dict[tuple[object, ...], int] = {}
    for row in merged:
        key = (row["cultivar"], row["condition"], row["time_index"])
        cells[key] = cells.get(key, 0) + 1
    if len(cells) != 20 or set(cells.values()) != {3}:
        raise RuntimeError(f"Expected balanced 2 x 2 x 5 x 3 design; observed cells: {cells}")
    if {row["library_layout"] for row in merged} != {"SINGLE"}:
        raise RuntimeError("Not all runs are single-end")

    output = args.output_dir / "GSE337039_sample_manifest.tsv"
    write_tsv(merged, output)

    total_mb = sum(float(row["normalized_sra_size_mb"]) for row in merged)
    total_spots = sum(int(row["spots"]) for row in merged)
    report = [
        "# GSE337039 run-manifest validation",
        "",
        "- Status: MANIFEST_OK",
        f"- Runs: {len(merged)}",
        "- Design: 2 cultivars x 2 conditions x 5 time points x 3 replicates",
        "- Layout: single-end",
        f"- Read length: {sorted({row['read_length'] for row in merged})}",
        f"- Total spots: {total_spots:,}",
        f"- Normalized SRA download size: {total_mb / 1024:.2f} GiB",
        f"- RunInfo SHA256: `{sha256(args.runinfo)}`",
        f"- Manifest SHA256: `{sha256(output)}`",
        "",
        "Converted FASTQ files will be substantially larger than the normalized SRA objects. Reserve at least 120 GB, preferably 200 GB, for raw rescue and temporary files.",
    ]
    (args.output_dir / "GSE337039_MANIFEST_REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print("MANIFEST_OK")
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
