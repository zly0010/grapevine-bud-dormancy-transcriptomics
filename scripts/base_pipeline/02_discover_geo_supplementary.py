from __future__ import annotations

import argparse
import csv
import re
import urllib.request
from html.parser import HTMLParser
from pathlib import Path


GEO_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"


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


def series_bucket(accession: str) -> str:
    match = re.fullmatch(r"GSE(\d+)", accession)
    if not match:
        raise ValueError(f"Invalid GEO series accession: {accession}")
    digits = match.group(1)
    return f"GSE{digits[:-3]}nnn"


def read_dataset_table(path: Path, core_only: bool) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if core_only:
        rows = [row for row in rows if row.get("priority") == "core" and row.get("use_in_analysis") == "yes"]
    return rows


def supplementary_dir_url(accession: str) -> str:
    return f"{GEO_FTP_BASE}/{series_bucket(accession)}/{accession}/suppl/"


def parse_listing(accession: str) -> list[dict[str, str]]:
    url = supplementary_dir_url(accession)
    with urllib.request.urlopen(url, timeout=120) as response:
        html = response.read().decode("utf-8", errors="replace")

    parser = HrefParser()
    parser.feed(html)
    records: list[dict[str, str]] = []
    for href in parser.hrefs:
        if href.startswith("/") or href.startswith("?") or href.startswith("http") or href == "../":
            continue
        if href == "Parent Directory":
            continue
        file_url = url + href
        records.append(
            {
                "accession": accession,
                "filename": href,
                "url": file_url,
                "likely_count_matrix": str(is_likely_count_matrix(href)).lower(),
            }
        )
    return records


def is_likely_count_matrix(filename: str) -> bool:
    lowered = filename.lower()
    if lowered == "filelist.txt":
        return False
    matrix_terms = ["count", "counts", "matrix", "raw", "expression", "featurecounts"]
    return any(term in lowered for term in matrix_terms) and lowered.endswith((".txt", ".txt.gz", ".tsv", ".tsv.gz", ".csv", ".csv.gz", ".tar"))


def write_records(records: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["accession", "filename", "url", "likely_count_matrix"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)


def main() -> int:
    parser = argparse.ArgumentParser(description="Discover GEO supplementary files for selected datasets.")
    parser.add_argument("--datasets", type=Path, default=Path("config/datasets.tsv"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/geo_supplementary_files.tsv"))
    parser.add_argument("--core-only", action="store_true")
    args = parser.parse_args()

    all_records: list[dict[str, str]] = []
    for row in read_dataset_table(args.datasets, args.core_only):
        accession = row["accession"]
        records = parse_listing(accession)
        all_records.extend(records)
        print(f"[supplementary] {accession}: {len(records)} files")
        for record in records:
            tag = "count?" if record["likely_count_matrix"] == "true" else "other "
            print(f"  [{tag}] {record['filename']}")

    write_records(all_records, args.output)
    print(f"[ok] wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
