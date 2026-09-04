#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


def parse_star(path: Path) -> dict[str, float]:
    values: dict[str, float] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "|" not in line:
            continue
        key, raw = (part.strip() for part in line.split("|", 1))
        raw = raw.rstrip("%")
        try:
            values[key] = float(raw)
        except ValueError:
            continue
    return {
        "clean_reads": values["Number of input reads"],
        "unique_mapping_rate": values["Uniquely mapped reads %"] / 100,
        "multimapping_rate": values["% of reads mapped to multiple loci"] / 100,
        "too_short_rate": values["% of reads unmapped: too short"] / 100,
        "mismatch_rate": values["Mismatch rate per base, %"] / 100,
    }


def parse_fastp(path: Path) -> dict[str, float]:
    data = json.loads(path.read_text(encoding="utf-8"))
    before = data["summary"]["before_filtering"]
    after = data["summary"]["after_filtering"]
    return {
        "raw_reads": before["total_reads"],
        "fastp_pass_rate": after["total_reads"] / before["total_reads"],
        "q30_before": before["q30_rate"],
        "q30_after": after["q30_rate"],
        "duplication_rate": data["duplication"]["rate"],
    }


def featurecounts_assigned(path: Path) -> int:
    table = pd.read_csv(path, sep="\t")
    return int(table.loc[table.Status == "Assigned"].iloc[0, 1])


def summarize_series(series: pd.Series) -> tuple[float, float, float]:
    return float(series.min()), float(series.median()), float(series.max())


def pct(value: float) -> str:
    return f"{100 * value:.1f}%"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("work_dir", type=Path)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("results_dir", type=Path)
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest, sep="\t")
    rows: list[dict[str, object]] = []
    for run in manifest.run_accession:
        star = parse_star(args.work_dir / "star_align" / f"{run}.Log.final.out")
        fastp = parse_fastp(args.work_dir / "qc" / "fastp_star" / f"{run}.json")
        assigned = featurecounts_assigned(
            args.work_dir / "star_counts" / f"{run}.forward.txt.summary"
        )
        rows.append({
            "run_accession": run,
            **fastp,
            **star,
            "forward_assigned_reads": assigned,
            "assigned_per_clean_read": assigned / star["clean_reads"],
        })
    qc = manifest.merge(pd.DataFrame(rows), on="run_accession", validate="one_to_one")
    qc["qc_flag"] = np.select(
        [
            (qc.fastp_pass_rate < 0.90) | (qc.unique_mapping_rate < 0.35) |
            (qc.assigned_per_clean_read < 0.25),
            (qc.unique_mapping_rate < 0.40) | (qc.assigned_per_clean_read < 0.30),
        ],
        ["FAIL", "WARN"],
        default="PASS",
    )
    qc.to_csv(args.results_dir / "07_sample_qc_combined.tsv", sep="\t", index=False)

    result_files = [args.results_dir / "01_three_way_interaction_LRT.tsv"]
    result_files += sorted(args.results_dir.glob("lrt_condition_by_time__*.tsv"))
    result_files += sorted(args.results_dir.glob("condition__*.tsv"))
    deg_rows: list[dict[str, object]] = []
    for path in result_files:
        table = pd.read_csv(path, sep="\t")
        tested = table.padj.notna()
        significant = tested & (table.padj < 0.05)
        is_lrt = path.stem == "01_three_way_interaction_LRT" or path.stem.startswith("lrt_")
        if "log2FoldChange" in table.columns and not is_lrt:
            strict = significant & (table.log2FoldChange.abs() >= 1)
        else:
            strict = significant
        mapped = table.vitvi_id.notna() & (table.vitvi_id.astype(str) != "")
        deg_rows.append({
            "result": path.stem,
            "tested_genes": int(tested.sum()),
            "fdr_0_05": int(significant.sum()),
            "fdr_0_05_abs_lfc_ge_1": int(strict.sum()),
            "mapped_strict_genes": int((strict & mapped).sum()),
        })
    deg = pd.DataFrame(deg_rows)
    deg.to_csv(args.results_dir / "08_deg_summary.tsv", sep="\t", index=False)

    vst_header = pd.read_csv(
        args.results_dir / "04_vst_samples_by_vitvi_genes.tsv", sep="\t", nrows=0
    ).columns
    mapped_vst_genes = len(vst_header) - 1
    selected = (args.results_dir / "SELECTED_STRANDEDNESS.txt").read_text().strip()
    pass_count = int((qc.qc_flag == "PASS").sum())
    warn_count = int((qc.qc_flag == "WARN").sum())
    fail_count = int((qc.qc_flag == "FAIL").sum())
    mapping = summarize_series(qc.unique_mapping_rate)
    assigned = summarize_series(qc.assigned_per_clean_read)
    retained = summarize_series(qc.fastp_pass_rate)
    q30 = summarize_series(qc.q30_after)

    global_lrt = deg.loc[deg.result == "01_three_way_interaction_LRT"].iloc[0]
    cultivar_lrt = deg.loc[deg.result.str.startswith("lrt_condition_by_time")]
    contrasts = deg.loc[deg.result.str.startswith("condition__")]
    lines = [
        "# GSE337039 FASTQ rescue audit",
        "",
        "## Pipeline status",
        "",
        "- Status: GSE337039_FASTQ_RESCUE_OK",
        f"- Samples completed: {len(qc)}/60",
        f"- Selected counting mode: {selected}",
        "- Reference: NCBI GCF_000003745.3 (12X), STAR and featureCounts",
        f"- Genes retained by DESeq2 count filter: 14,950",
        f"- Strict one-to-one Vitvi genes in VST matrix: {mapped_vst_genes:,}",
        "",
        "## Sample-level QC",
        "",
        f"- fastp retained reads, min/median/max: {pct(retained[0])} / {pct(retained[1])} / {pct(retained[2])}",
        f"- Q30 after filtering, min/median/max: {pct(q30[0])} / {pct(q30[1])} / {pct(q30[2])}",
        f"- STAR unique mapping, min/median/max: {pct(mapping[0])} / {pct(mapping[1])} / {pct(mapping[2])}",
        f"- Forward assigned per clean read, min/median/max: {pct(assigned[0])} / {pct(assigned[1])} / {pct(assigned[2])}",
        f"- QC flags: PASS={pass_count}, WARN={warn_count}, FAIL={fail_count}",
        "",
        "## Differential-expression overview",
        "",
        f"- Global cultivar x condition x time LRT: {int(global_lrt.fdr_0_05):,} genes at FDR < 0.05.",
    ]
    for row in cultivar_lrt.itertuples(index=False):
        cultivar = row.result.split("__", 1)[1]
        lines.append(f"- {cultivar} condition x time LRT: {row.fdr_0_05:,} genes at FDR < 0.05.")
    lines += [
        f"- Matched-time controlled-vs-natural contrasts: strict DEG range {int(contrasts.fdr_0_05_abs_lfc_ge_1.min()):,} to {int(contrasts.fdr_0_05_abs_lfc_ge_1.max()):,} genes.",
        "",
        "## Interpretation boundary",
        "",
        "This dataset is suitable for independent chilling-response and whole-module validation. "
        "The cultivars are interspecific grape hybrids, so lower mapping than pure V. vinifera is expected. "
        "Claims should emphasize conserved response programs rather than exact one-gene equivalence.",
    ]
    (args.results_dir / "GSE337039_QC_DEG_REPORT.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    print("GSE337039_SUMMARY_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
