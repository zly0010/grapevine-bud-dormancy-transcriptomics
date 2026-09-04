#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path


REFSEQ_COLUMNS = (15, 35, 39, 43)
ID_PATTERN = re.compile(r"LOC\d+|[A-Za-z][A-Za-z0-9_.-]*")
MISSING = {"", "na", "n/a", "none", "null", "-"}


def refseq_ids(value: str) -> set[str]:
    tokens = set(ID_PATTERN.findall(value.strip()))
    return {token for token in tokens if token.lower() not in MISSING}


def read_gtf_gene_ids(path: Path) -> set[str]:
    genes: set[str] = set()
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            match = re.search(r'gene_id "([^"]+)"', fields[8])
            if match:
                genes.add(match.group(1))
    return genes


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("correspondence_tsv", type=Path)
    parser.add_argument("ncbi_gtf", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()

    ref_to_vitvi: dict[str, set[str]] = defaultdict(set)
    vitvi_to_ref: dict[str, set[str]] = defaultdict(set)
    vitvi_seen: set[str] = set()

    with args.correspondence_tsv.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader)
        for fields in reader:
            if not fields:
                continue
            vitvi = fields[0].strip()
            if not vitvi.startswith("Vitvi"):
                continue
            vitvi_seen.add(vitvi)
            ids: set[str] = set()
            for index in REFSEQ_COLUMNS:
                if index < len(fields):
                    ids.update(refseq_ids(fields[index]))
            for refseq in ids:
                ref_to_vitvi[refseq].add(vitvi)
                vitvi_to_ref[vitvi].add(refseq)

    gtf_genes = read_gtf_gene_ids(args.ncbi_gtf)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    mapping_path = args.output_dir / "refseq_to_vitvi_mapping.tsv"

    rows: list[dict[str, object]] = []
    for refseq in sorted(ref_to_vitvi):
        for vitvi in sorted(ref_to_vitvi[refseq]):
            ref_n = len(ref_to_vitvi[refseq])
            vitvi_n = len(vitvi_to_ref[vitvi])
            rows.append(
                {
                    "refseq_gene_id": refseq,
                    "vitvi_id": vitvi,
                    "refseq_to_vitvi_n": ref_n,
                    "vitvi_to_refseq_n": vitvi_n,
                    "one_to_one": ref_n == 1 and vitvi_n == 1,
                    "present_in_ncbi_gtf": refseq in gtf_genes,
                }
            )

    with mapping_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    mapped_gtf = gtf_genes.intersection(ref_to_vitvi)
    one_to_one_gtf = {
        refseq
        for refseq in mapped_gtf
        if len(ref_to_vitvi[refseq]) == 1
        and len(vitvi_to_ref[next(iter(ref_to_vitvi[refseq]))]) == 1
    }
    report = [
        "# RefSeq to Vitvi mapping audit",
        "",
        "- Status: REFSEQ_VITVI_MAPPING_OK",
        f"- VCost/Vitvi genes in correspondence table: {len(vitvi_seen):,}",
        f"- RefSeq IDs in correspondence table: {len(ref_to_vitvi):,}",
        f"- Genes in NCBI GTF: {len(gtf_genes):,}",
        f"- NCBI GTF genes with any Vitvi mapping: {len(mapped_gtf):,}",
        f"- NCBI GTF genes with strict one-to-one Vitvi mapping: {len(one_to_one_gtf):,}",
        "- Primary cross-dataset analyses must use strict one-to-one mappings.",
    ]
    (args.output_dir / "REFSEQ_VITVI_MAPPING_REPORT.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )
    print("REFSEQ_VITVI_MAPPING_OK")
    print(mapping_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
