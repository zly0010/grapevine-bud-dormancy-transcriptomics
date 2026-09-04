from __future__ import annotations

import argparse
import csv
import tarfile
import urllib.request
from pathlib import Path


def read_records(path: Path, likely_only: bool) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        records = list(csv.DictReader(handle, delimiter="\t"))
    if likely_only:
        records = [record for record in records if record.get("likely_count_matrix") == "true"]
    return records


def download_file(url: str, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and target.stat().st_size > 0:
        print(f"[skip] {target} already exists")
        return
    print(f"[download] {url}")
    with urllib.request.urlopen(url, timeout=300) as response:
        target.write_bytes(response.read())
    print(f"[ok] {target}")


def is_safe_tar_member(member_name: str) -> bool:
    path = Path(member_name)
    return not path.is_absolute() and ".." not in path.parts


def extract_tar(tar_path: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path) as archive:
        safe_members = [member for member in archive.getmembers() if is_safe_tar_member(member.name)]
        archive.extractall(output_dir, members=safe_members)
    print(f"[extract] {tar_path} -> {output_dir}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Download selected GEO supplementary files.")
    parser.add_argument("--manifest", type=Path, default=Path("data/processed/geo_supplementary_files.tsv"))
    parser.add_argument("--output-dir", type=Path, default=Path("data/raw/geo_supplementary"))
    parser.add_argument("--likely-only", action="store_true")
    parser.add_argument("--extract-tar", action="store_true")
    args = parser.parse_args()

    for record in read_records(args.manifest, args.likely_only):
        accession = record["accession"]
        filename = record["filename"]
        target = args.output_dir / accession / filename
        download_file(record["url"], target)
        if args.extract_tar and filename.lower().endswith(".tar"):
            extract_tar(target, args.output_dir / accession / target.stem)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
