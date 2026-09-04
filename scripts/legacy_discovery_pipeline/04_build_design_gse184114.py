#!/usr/bin/env python
"""Build sample design for GSE184114 corrected analysis."""
import pandas as pd
import re

df = pd.read_csv("data/legacy/data/processed/GSE184114/counts_matrix.txt",
                 sep="\t", index_col=0, nrows=1)
cols = df.columns.tolist()

rows = []
seen = {}  # base_name -> count for dedup
for c in cols:
    m = re.match(r"(Accl|Deaccl)_(Acclimation|Deacclimation)_(Control|ABA)_(\d+)_(\d+)(?:\.(\d+))?", c)
    if not m:
        print(f"UNMATCHED: {c}")
        continue
    phase_code, phase_name, treatment, time_val, rep, dup = m.groups()
    base_key = f"{phase_code}_{phase_name}_{treatment}_{time_val}_{rep}"
    if dup:
        suffix = f".{dup}"
    else:
        suffix = ""
    rows.append({
        "sample_id": c,
        "base_key": base_key,
        "phase": phase_name,
        "treatment": treatment,
        "time_h": int(time_val),
        "rep": int(rep),
        "dup": int(dup) if dup else 0,
    })

design = pd.DataFrame(rows)
print("Total samples:", len(design))
print("\nBy phase x treatment x time:")
summary = design.groupby(["phase", "treatment", "time_h"]).size().reset_index(name="n")
print(summary.to_string(index=False))

# Mark whether time point has complete design (both Control and ABA)
phase_treat_time = design.groupby(["phase", "time_h", "treatment"]).size().reset_index(name="n")
complete = []
for _, row in phase_treat_time.iterrows():
    complete.append(row)

print("\n\nAll combinations:")
for _, r in phase_treat_time.iterrows():
    print(f"  {r['phase']} | {r['treatment']:8s} | {r['time_h']:3d}h | n={r['n']}")

# Determine which time points have both Control and ABA
phase_times = design.groupby(["phase", "time_h"])["treatment"].apply(set).reset_index()
phase_times["n_treatments"] = phase_times["treatment"].apply(len)
phase_times["has_both"] = phase_times["n_treatments"] == 2

print("\n\nTime points with both treatments:")
for _, r in phase_times[phase_times["has_both"]].iterrows():
    print(f"  {r['phase']} {r['time_h']}h: {r['n_treatments']} treatments")

print("\nTime points with only one treatment (baseline):")
for _, r in phase_times[~phase_times["has_both"]].iterrows():
    print(f"  {r['phase']} {r['time_h']}h: {sorted(r['treatment'])}")

# Build interaction model time points (post-treatment only, complete design)
# Per user instructions:
# Acclimation: 2, 4, 24, 48 h (exclude 0h baseline)
# Deacclimation: 6, 12, 24, 48, 72 h (exclude 0h baseline)
print("\n\n=== Interaction model samples ===")
for phase, times in [("Acclimation", [2, 4, 24, 48]),
                     ("Deacclimation", [6, 12, 24, 48, 72])]:
    subset = design[(design["phase"] == phase) & (design["time_h"].isin(times))]
    print(f"{phase}: {len(subset)} samples across {sorted(subset['time_h'].unique())} hours")
    by_t = subset.groupby("treatment").size()
    for t, n in by_t.items():
        print(f"  {t}: {n}")

# Save design
design.to_csv("data/legacy/results_corrected/01_sample_design_and_qc/sample_design_gse184114.txt",
              sep="\t", index=False)
print("\nDesign saved.")
