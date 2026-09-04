from __future__ import annotations

import argparse
import csv
import gzip
from html.parser import HTMLParser
import re
import sys
import urllib.request
from pathlib import Path


GEO_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"


def series_bucket(accession: str) -> str:
    match = re.fullmatch(r"GSE(\d+)", accession)
    if not match:
        raise ValueError(f"Invalid GEO series accession: {accession}")
    digits = match.group(1)
    return f"GSE{digits[:-3]}nnn"


def matrix_url(accession: str) -> str:
    bucket = series_bucket(accession)
    return f"{GEO_FTP_BASE}/{bucket}/{accession}/matrix/{accession}_series_matrix.txt.gz"


def matrix_dir_url(accession: str) -> str:
    bucket = series_bucket(accession)
    return f"{GEO_FTP_BASE}/{bucket}/{accession}/matrix/"


class HrefParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.hrefs: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag != "a":
            return
        for key, value in attrs:
            if key == "href" and value:
                self.hrefs.append(value)


def discover_matrix_urls(accession: str) -> list[tuple[str, str]]:
    url = matrix_dir_url(accession)
    with urllib.request.urlopen(url, timeout=120) as response:
        html = response.read().decode("utf-8", errors="replace")

    parser = HrefParser()
    parser.feed(html)

    files = sorted(href for href in parser.hrefs if href.endswith("_series_matrix.txt.gz"))
    if not files:
        files = [Path(matrix_url(accession)).name]

    discovered: list[tuple[str, str]] = []
    for filename in files:
        platform_match = re.search(r"-(GPL\d+)_series_matrix", filename)
        label = platform_match.group(1) if platform_match else "default"
        discovered.append((label, url + filename))
    return discovered


def read_dataset_table(path: Path, core_only: bool) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if core_only:
        rows = [row for row in rows if row.get("priority") == "core" and row.get("use_in_analysis") == "yes"]
    return rows


def download_file(url: str, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and target.stat().st_size > 0:
        print(f"[skip] {target} already exists")
        return
    print(f"[download] {url}")
    with urllib.request.urlopen(url, timeout=120) as response:
        target.write_bytes(response.read())
    print(f"[ok] {target}")


def parse_series_matrix(matrix_gz: Path, platform: str) -> tuple[list[str], list[dict[str, str]]]:
    sample_ids: list[str] = []
    metadata_by_sample: dict[str, dict[str, str]] = {}
    pending_fields: list[tuple[str, list[str]]] = []

    with gzip.open(matrix_gz, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue

            parts = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
            field = parts[0].lstrip("!")
            values = [part.strip('"') for part in parts[1:]]

            if field == "Sample_geo_accession":
                sample_ids = values
                metadata_by_sample = {sample: {"sample_accession": sample, "platform": platform} for sample in sample_ids}
                for pending_field, pending_values in pending_fields:
                    add_sample_field(metadata_by_sample, sample_ids, pending_field, pending_values)
                continue

            if not sample_ids:
                pending_fields.append((field, values))
                continue

            add_sample_field(metadata_by_sample, sample_ids, field, values)

    rows = [metadata_by_sample[sample] for sample in sample_ids]
    return sample_ids, rows


def add_sample_field(
    metadata_by_sample: dict[str, dict[str, str]],
    sample_ids: list[str],
    field: str,
    values: list[str],
) -> None:
    for sample, value in zip(sample_ids, values):
        metadata_by_sample.setdefault(sample, {"sample_accession": sample})
        if field == "Sample_characteristics_ch1" and ": " in value:
            key, parsed_value = value.split(": ", 1)
            safe_key = "characteristic_" + re.sub(r"[^0-9A-Za-z_]+", "_", key.strip().lower()).strip("_")
            metadata_by_sample[sample][safe_key] = parsed_value
        else:
            key = re.sub(r"[^0-9A-Za-z_]+", "_", field.lower()).strip("_")
            if key in metadata_by_sample[sample] and metadata_by_sample[sample][key] != value:
                key = f"{key}_{len([k for k in metadata_by_sample[sample] if k.startswith(key)]) + 1}"
            metadata_by_sample[sample][key] = value


def write_metadata(accession: str, rows: list[dict[str, str]], output_dir: Path, platform: str | None = None) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    suffix = f"_{platform}" if platform else ""
    out_path = output_dir / f"{accession}{suffix}_sample_metadata.tsv"
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)

    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return out_path


def main() -> int:
    parser = argparse.ArgumentParser(description="Download GEO series matrix files and parse sample metadata.")
    parser.add_argument("--datasets", type=Path, default=Path("config/datasets.tsv"))
    parser.add_argument("--raw-dir", type=Path, default=Path("data/raw/geo"))
    parser.add_argument("--metadata-dir", type=Path, default=Path("data/processed/metadata"))
    parser.add_argument("--core-only", action="store_true", help="Download only core datasets marked for analysis.")
    args = parser.parse_args()

    datasets = read_dataset_table(args.datasets, core_only=args.core_only)
    if not datasets:
        print("No datasets selected.", file=sys.stderr)
        return 1

    for row in datasets:
        accession = row["accession"]
        combined_rows: list[dict[str, str]] = []
        for platform, url in discover_matrix_urls(accession):
            target = args.raw_dir / accession / Path(url).name
            download_file(url, target)
            sample_ids, metadata_rows = parse_series_matrix(target, platform)
            combined_rows.extend(metadata_rows)
            out_path = write_metadata(accession, metadata_rows, args.metadata_dir, platform)
            print(f"[metadata] {accession} {platform}: {len(sample_ids)} samples -> {out_path}")

        combined_path = write_metadata(accession, combined_rows, args.metadata_dir)
        print(f"[metadata] {accession} combined: {len(combined_rows)} samples -> {combined_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
