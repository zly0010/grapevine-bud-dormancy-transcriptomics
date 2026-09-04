#!/usr/bin/env python3
"""Build evidence-backed GSE124820 QC and GO mapping audit tables."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def recompute_status(row: pd.Series) -> tuple[str, str]:
    fail = (row.library_size < row.library_size_fail_threshold) or (
        row.detected_genes < row.detected_genes_fail_threshold
    )
    watch = (row.library_size < row.library_size_watch_threshold) or (
        row.detected_genes < row.detected_genes_watch_threshold
    )
    status = "fail" if fail else "watch" if watch else "pass"
    reasons: list[str] = []
    if row.library_size < row.library_size_fail_threshold:
        reasons.append("library_size_lt_25pct_median")
    elif row.library_size < row.library_size_watch_threshold:
        reasons.append("library_size_lt_50pct_median")
    if row.detected_genes < row.detected_genes_fail_threshold:
        reasons.append("detected_genes_lt_70pct_median")
    elif row.detected_genes < row.detected_genes_watch_threshold:
        reasons.append("detected_genes_lt_85pct_median")
    return status, ";".join(reasons)


def build_qc(qc_path: Path, output_dir: Path) -> dict[str, object]:
    qc = pd.read_csv(qc_path, sep="\t", dtype=str).fillna("")
    qc["library_size"] = pd.to_numeric(qc.library_size, errors="raise")
    qc["detected_genes"] = pd.to_numeric(qc.detected_genes, errors="raise")
    medians = (
        qc.groupby(["accession", "label"], as_index=False)
        .agg(library_size_group_median=("library_size", "median"), detected_genes_group_median=("detected_genes", "median"))
    )
    qc = qc.merge(medians, on=["accession", "label"], validate="many_to_one")
    qc["library_size_fail_threshold"] = qc.library_size_group_median * 0.25
    qc["library_size_watch_threshold"] = qc.library_size_group_median * 0.50
    qc["detected_genes_fail_threshold"] = qc.detected_genes_group_median * 0.70
    qc["detected_genes_watch_threshold"] = qc.detected_genes_group_median * 0.85
    calculated = qc.apply(recompute_status, axis=1, result_type="expand")
    qc["status_recomputed"] = calculated[0]
    qc["reason_recomputed"] = calculated[1]
    qc["status_matches_frozen_file"] = qc.status_recomputed.eq(qc.qc_status)
    qc["reason_matches_frozen_file"] = qc.reason_recomputed.eq(qc.qc_reason)
    if not qc.status_matches_frozen_file.all() or not qc.reason_matches_frozen_file.all():
        mismatch = qc.loc[~(qc.status_matches_frozen_file & qc.reason_matches_frozen_file)]
        raise RuntimeError(f"QC rule reconstruction mismatch in {len(mismatch)} samples")
    qc["sample_id"] = qc.matrix_column
    qc["geo_accession"] = qc.sample_accession
    qc["variety"] = qc.short_genotype
    qc["status"] = qc.qc_status
    qc["included_primary_analysis"] = qc.status.isin(["pass", "watch"])
    qc["exclusion_reason"] = qc.qc_reason.where(qc.status.eq("fail"), "")
    qc["retention_reason"] = ""
    qc.loc[qc.status.eq("watch"), "retention_reason"] = (
        "retained_by_predefined_primary_QC_B_rule; sensitivity_QC_C_available"
    )
    qc["library_size_watch_flag"] = qc.library_size < qc.library_size_watch_threshold
    qc["library_size_fail_flag"] = qc.library_size < qc.library_size_fail_threshold
    qc["detected_genes_watch_flag"] = qc.detected_genes < qc.detected_genes_watch_threshold
    qc["detected_genes_fail_flag"] = qc.detected_genes < qc.detected_genes_fail_threshold
    columns = [
        "sample_id", "geo_accession", "sample_title", "variety", "time", "replicate",
        "library_size", "detected_genes", "library_size_group_median", "detected_genes_group_median",
        "library_size_watch_threshold", "library_size_fail_threshold",
        "detected_genes_watch_threshold", "detected_genes_fail_threshold",
        "library_size_watch_flag", "library_size_fail_flag",
        "detected_genes_watch_flag", "detected_genes_fail_flag",
        "status", "included_primary_analysis", "qc_reason", "exclusion_reason", "retention_reason",
        "status_matches_frozen_file", "reason_matches_frozen_file",
    ]
    result = qc[columns].sort_values(["variety", "time", "replicate", "sample_id"])
    result.to_csv(output_dir / "Supplementary_Table_QC_GSE124820.csv", index=False)
    return {
        "total": len(result),
        "status_counts": result.status.value_counts().sort_index().to_dict(),
        "included": int(result.included_primary_analysis.sum()),
        "all_rules_reproduced": bool(result.status_matches_frozen_file.all() and result.reason_matches_frozen_file.all()),
        "library_size_median": float(result.library_size_group_median.iloc[0]),
        "detected_genes_median": float(result.detected_genes_group_median.iloc[0]),
    }


def build_go(
    mapping_path: Path,
    go_path: Path,
    module_path: Path,
    output_dir: Path,
) -> dict[str, object]:
    mapping = pd.read_csv(mapping_path, sep="\t", dtype=str).fillna("")
    go = pd.read_csv(go_path, sep="\t", dtype=str).fillna("")
    modules = pd.read_csv(module_path, sep="\t", dtype=str).fillna("")
    if len(modules) != 5_000 or modules.gene_id.duplicated().any():
        raise RuntimeError("formal module assignment is not the expected unique 5,000-gene set")

    mapping_groups = (
        mapping.groupby("vitvi_id").ensembl_id
        .apply(lambda values: sorted({value.strip() for value in values if value.strip()}))
        .to_dict()
    )
    go_unique = go.loc[(go.ensembl_gene_id != "") & (go.go_id != "")].drop_duplicates()
    go_by_gene = go_unique.groupby("ensembl_gene_id")
    go_ids_by_gene = go_by_gene.go_id.apply(lambda x: sorted(set(x))).to_dict()
    types_by_gene = go_by_gene.go_type.apply(lambda x: set(x)).to_dict()

    rows: list[dict[str, object]] = []
    mapping_id_set = set(mapping.vitvi_id)
    for row in modules.itertuples(index=False):
        ids = mapping_groups.get(row.gene_id, [])
        annotated_ids = [gene_id for gene_id in ids if gene_id in go_ids_by_gene]
        types = set().union(*(types_by_gene.get(gene_id, set()) for gene_id in ids)) if ids else set()
        if row.gene_id not in mapping_id_set:
            reason = "gene_id_not_present_in_mapping_table"
        elif len(ids) == 0:
            reason = "gene_id_unmapped_to_ensembl"
        elif len(ids) > 1:
            reason = "one_to_many_mapping_present_not_excluded"
        elif len(annotated_ids) == 0:
            reason = "mapped_ensembl_gene_without_go"
        else:
            reason = "go_annotated"
        rows.append(
            {
                "gene_id": row.gene_id,
                "module_color": row.module_color,
                "recognized_original_id": row.gene_id in mapping_id_set,
                "mapping_cardinality": len(ids),
                "mapping_class": "unmapped" if len(ids) == 0 else "one_to_one" if len(ids) == 1 else "one_to_many",
                "ensembl_gene_id": ";".join(ids),
                "mapped_current_id": len(ids) > 0,
                "go_annotated": len(annotated_ids) > 0,
                "has_BP": "P" in types,
                "has_MF": "F" in types,
                "has_CC": "C" in types,
                "go_term_count": len(set().union(*(set(go_ids_by_gene.get(gene_id, [])) for gene_id in ids))) if ids else 0,
                "annotation_loss_reason": reason,
                "original_id_source": "VCost.v3 Vitvi gene identifiers",
                "reference_annotation": "Ensembl Plants release 62; Vitis vinifera ASM3070453v1",
                "go_source": "UniProtKB REST organism_id:29760 with GO; parsed 2026-07-16",
                "mapping_handling": "one-to-one retained; one-to-many would be reported and retained for audit; none observed in formal 5,000-gene set",
            }
        )
    audit = pd.DataFrame(rows)
    audit.to_csv(output_dir / "GO_annotation_audit_v2.csv", index=False)

    flow_metrics = [
        ("total_wgcna_genes", len(audit)),
        ("recognized_original_ids", int(audit.recognized_original_id.sum())),
        ("mapped_current_ids", int(audit.mapped_current_id.sum())),
        ("one_to_one_mappings", int(audit.mapping_class.eq("one_to_one").sum())),
        ("one_to_many_mappings", int(audit.mapping_class.eq("one_to_many").sum())),
        ("genes_with_any_GO", int(audit.go_annotated.sum())),
        ("genes_with_BP", int(audit.has_BP.sum())),
        ("genes_with_MF", int(audit.has_MF.sum())),
        ("genes_with_CC", int(audit.has_CC.sum())),
    ]
    summary_rows: list[dict[str, object]] = []
    for metric, count in flow_metrics:
        summary_rows.append(
            {"scope": "all_5000", "module_color": "all", "metric": metric, "count": count, "percent_of_scope": 100 * count / len(audit)}
        )
    for reason, count in audit.annotation_loss_reason.value_counts().sort_index().items():
        summary_rows.append(
            {"scope": "all_5000_reason", "module_color": "all", "metric": reason, "count": int(count), "percent_of_scope": 100 * count / len(audit)}
        )
    for module in ["blue", "turquoise", "brown"]:
        subset = audit.loc[audit.module_color == module]
        for metric, series in [
            ("module_genes", pd.Series([True] * len(subset), index=subset.index)),
            ("mapped_current_ids", subset.mapped_current_id),
            ("genes_with_any_GO", subset.go_annotated),
            ("genes_with_BP", subset.has_BP),
            ("genes_with_MF", subset.has_MF),
            ("genes_with_CC", subset.has_CC),
        ]:
            count = int(series.sum())
            summary_rows.append(
                {"scope": "target_module", "module_color": module, "metric": metric, "count": count, "percent_of_scope": 100 * count / len(subset)}
            )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(output_dir / "GO_annotation_mapping_summary_v2.csv", index=False)
    return {
        "total": len(audit),
        "mapped": int(audit.mapped_current_id.sum()),
        "go_any": int(audit.go_annotated.sum()),
        "bp": int(audit.has_BP.sum()),
        "mf": int(audit.has_MF.sum()),
        "cc": int(audit.has_CC.sum()),
        "reason_counts": audit.annotation_loss_reason.value_counts().sort_index().to_dict(),
        "mapping_cardinality": audit.mapping_class.value_counts().sort_index().to_dict(),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--qc", type=Path, required=True)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--go", type=Path, required=True)
    parser.add_argument("--modules", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    expected_outputs = [
        args.output_dir / "Supplementary_Table_QC_GSE124820.csv",
        args.output_dir / "GO_annotation_audit_v2.csv",
        args.output_dir / "GO_annotation_mapping_summary_v2.csv",
        args.output_dir / "AUDIT_TABLE_RUN_INFO.json",
    ]
    existing = [path for path in expected_outputs if path.exists()]
    if existing:
        raise RuntimeError(f"Refusing to overwrite existing audit outputs: {existing}")
    qc_summary = build_qc(args.qc, args.output_dir)
    go_summary = build_go(args.mapping, args.go, args.modules, args.output_dir)
    provenance = {
        "qc_summary": qc_summary,
        "go_summary": go_summary,
        "inputs": {
            str(path): sha256(path)
            for path in [args.qc, args.mapping, args.go, args.modules]
        },
    }
    (args.output_dir / "AUDIT_TABLE_RUN_INFO.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(provenance, ensure_ascii=False, indent=2))
    print("QC_GO_AUDIT_V2_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
