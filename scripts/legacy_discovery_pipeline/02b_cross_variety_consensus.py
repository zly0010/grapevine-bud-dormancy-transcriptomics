#!/usr/bin/env python3
"""Cross-variety consensus analysis for GSE124820"""

import pandas as pd
import numpy as np
from itertools import combinations
import os

base = "data/legacy/results_corrected/02_deseq2_gse124820"
consensus_dir = "data/legacy/results_corrected/02_deseq2_gse124820/cross_variety_consensus"
os.makedirs(consensus_dir, exist_ok=True)

varieties = ["Vamu", "Vvcs", "Vvri", "Vrip"]
variety_names = {
    "Vamu": "V. amurensis PI588635",
    "Vvcs": "V. vinifera Cab. Sauvignon",
    "Vvri": "V. vinifera Riesling",
    "Vrip": "V. riparia PI588275"
}

print("=" * 60)
print("GSE124820 Cross-Variety Consensus Analysis (QC scheme B)")
print("=" * 60)

# ============================================
# 1. LRT consensus
# ============================================
print("\n[1] LRT gene overlap (padj < 0.05):")

lrt_genes = {}
for v in varieties:
    # Try both naming conventions
    f1 = f"{base}/{v}_B_no_fail/LRT_{v}_B_no_fail.txt"
    f2 = f"{base}/{v}_B_no_fail/LRT_{v}_B.txt"
    f = f1 if os.path.exists(f1) else f2
    if os.path.exists(f):
        df = pd.read_csv(f, sep='\t')
        sig = set(df[df['padj'] < 0.05]['gene_id'].dropna())
        lrt_genes[v] = sig
        print(f"  {v} ({variety_names[v]}): {len(sig)} genes")

# Pairwise overlaps
print("\n  Pairwise overlaps:")
for v1, v2 in combinations(varieties, 2):
    overlap = lrt_genes[v1] & lrt_genes[v2]
    jaccard = len(overlap) / len(lrt_genes[v1] | lrt_genes[v2]) if (lrt_genes[v1] | lrt_genes[v2]) else 0
    print(f"    {v1} & {v2}: {len(overlap)} ({jaccard:.3f})")

# 3+ and 4/4 overlap
all_genes = set()
for v in varieties:
    all_genes |= lrt_genes[v]

consensus = pd.DataFrame({"gene_id": sorted(all_genes)})
for v in varieties:
    consensus[v] = consensus["gene_id"].isin(lrt_genes[v]).astype(int)
consensus["n_support"] = consensus[varieties].sum(axis=1)

print(f"\n  4/4 varieties: {(consensus['n_support'] == 4).sum()}")
print(f"  >=3/4 varieties: {(consensus['n_support'] >= 3).sum()}")
print(f"  >=2/4 varieties: {(consensus['n_support'] >= 2).sum()}")
print(f"  Any (1+): {(consensus['n_support'] >= 1).sum()}")

consensus.to_csv(f"{consensus_dir}/LRT_consensus_table_B.txt", sep='\t', index=False)

# ============================================
# 2. Day10 DEG consensus (shrunken LFC)
# ============================================
print("\n[2] Day10 vs Day0 DEG overlap (shrunken, |log2FC|>1, padj<0.05):")

day10_genes = {}
day10_lfc = {}
for v in varieties:
    f1 = f"{base}/{v}_B_no_fail/DEG_{v}_Day10_vs_Day0_B_no_fail.txt"
    f2 = f"{base}/{v}_B_no_fail/DEG_{v}_Day10_vs_Day0_B.txt"
    f = f1 if os.path.exists(f1) else f2
    if os.path.exists(f):
        df = pd.read_csv(f, sep='\t')
        sig = df[(df['padj_shrink'] < 0.05) & (df['log2FC_shrink'].abs() > 1)]
        day10_genes[v] = set(sig['gene_id'])
        day10_lfc[v] = dict(zip(sig['gene_id'], sig['log2FC_shrink']))
        print(f"  {v}: {len(day10_genes[v])} DEGs")

# Overlap
if len(day10_genes) == 4:
    common_all = set.intersection(*day10_genes.values())
    print(f"\n  Common to all 4: {len(common_all)}")
    
    # Direction consistency
    if common_all:
        rows = []
        for g in sorted(common_all):
            lfcs = {v: day10_lfc[v].get(g, 0) for v in varieties}
            all_up = all(x > 0 for x in lfcs.values())
            all_down = all(x < 0 for x in lfcs.values())
            direction = "up" if all_up else ("down" if all_down else "mixed")
            mean_abs = np.mean([abs(x) for x in lfcs.values()])
            rows.append({"gene_id": g, **lfcs, "consensus": direction, "mean_abs_lfc": mean_abs})
        
        direction_df = pd.DataFrame(rows)
        direction_df = direction_df.sort_values("mean_abs_lfc", ascending=False)
        
        print(f"  Consistent up: {(direction_df['consensus'] == 'up').sum()}")
        print(f"  Consistent down: {(direction_df['consensus'] == 'down').sum()}")
        print(f"  Mixed: {(direction_df['consensus'] == 'mixed').sum()}")
        
        direction_df.to_csv(f"{consensus_dir}/Day10_consensus_DEGs_B.txt", sep='\t', index=False)
        
        print("\n  Top 20 consensus genes:")
        print(direction_df[['gene_id'] + varieties + ['consensus', 'mean_abs_lfc']].head(20).to_string(index=False))

# ============================================
# 3. QC sensitivity comparison
# ============================================
print("\n[3] QC sensitivity comparison:")

qc_summary = []
for v in varieties:
    for qc in ['A_all', 'B_no_fail', 'C_pass_only']:
        lrt_f1 = f"{base}/{v}_{qc}/LRT_{v}_{qc}.txt"
        lrt_f2 = f"{base}/{v}_{qc}/LRT_{v}_{qc.replace('_no_fail','').replace('_all','').replace('_pass_only','')}.txt"
        lrt_f = lrt_f1 if os.path.exists(lrt_f1) else lrt_f2
        if os.path.exists(lrt_f):
            df = pd.read_csv(lrt_f, sep='\t')
            n_lrt = (df['padj'] < 0.05).sum()
            
            deg_f1 = f"{base}/{v}_{qc}/DEG_{v}_Day10_vs_Day0_{qc}.txt"
            deg_f2 = f"{base}/{v}_{qc}/DEG_{v}_Day10_vs_Day0_{qc.replace('_no_fail','').replace('_all','').replace('_pass_only','')}.txt"
            deg_f = deg_f1 if os.path.exists(deg_f1) else deg_f2
            n_deg10 = 0
            if os.path.exists(deg_f):
                df2 = pd.read_csv(deg_f, sep='\t')
                n_deg10 = ((df2['padj_shrink'] < 0.05) & (df2['log2FC_shrink'].abs() > 1)).sum()
            
            qc_summary.append({
                'variety': v, 'qc_scheme': qc,
                'n_lrt_sig': n_lrt, 'n_deg_day10': n_deg10
            })

qc_df = pd.DataFrame(qc_summary)
print(qc_df.to_string(index=False))
qc_df.to_csv(f"{consensus_dir}/qc_sensitivity_comparison.txt", sep='\t', index=False)

# ============================================
# 4. Summary report
# ============================================
print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"\nCorrected analysis: 4 varieties (Vamu, Vvcs, Vvri, Vrip)")
print(f"QC scheme B (default): exclude fail samples")
print(f"\nLRT consensus (padj<0.05):")
print(f"  4/4: {(consensus['n_support'] == 4).sum()} genes")
print(f"  >=3/4: {(consensus['n_support'] >= 3).sum()} genes")
if len(day10_genes) == 4 and common_all:
    print(f"\nDay10 DEG consensus:")
    print(f"  Common 4/4: {len(common_all)}")
    print(f"  Consistent direction: {((direction_df['consensus'] != 'mixed')).sum()}")

print("\nINVALIDATED旧结果:")
print("  - 3832/5240 共有DEG: INVALIDATED (错误的3品种模型)")
print("  - Top shared genes: INVALIDATED")
print("  - 方向一致率: INVALIDATED")
print("  - 染色体富集: INVALIDATED")
print("\n完成!")
