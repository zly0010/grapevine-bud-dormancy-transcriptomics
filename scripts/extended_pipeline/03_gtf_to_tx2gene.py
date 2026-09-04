#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import re
from pathlib import Path


ATTR = re.compile(r'(\S+)\s+"([^"]+)"')


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("gtf", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    pairs: dict[str, str] = {}
    with open_text(args.gtf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"transcript", "exon", "CDS"}:
                continue
            attrs = dict(ATTR.findall(fields[8]))
            transcript = attrs.get("transcript_id")
            gene = attrs.get("gene_id")
            if not transcript or not gene:
                continue
            previous = pairs.setdefault(transcript, gene)
            if previous != gene:
                raise RuntimeError(f"Transcript {transcript} maps to multiple genes")

    if not pairs:
        raise RuntimeError("No transcript_id/gene_id pairs parsed from GTF")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["transcript_id", "gene_id"])
        writer.writerows(sorted(pairs.items()))
    print(f"TX2GENE_OK transcripts={len(pairs)} output={args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
