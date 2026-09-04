#!/usr/bin/env python3
"""Build GSE124820 sample design with correct 4-variety grouping"""

import pandas as pd

df = pd.read_csv('data/processed/design_qc/GSE124820_raw_counts_design_qc.tsv', sep='\t')

design = df[['matrix_column', 'sample_title', 'short_genotype', 'time', 'replicate', 'qc_status', 'qc_reason']].copy()
design.columns = ['sample_id', 'sample_title', 'variety', 'time', 'replicate', 'qc_status', 'qc_reason']

variety_map = {
    'Vamu': 'Vitis amurensis PI588635',
    'Vvcs': 'Vitis vinifera Cabernet Sauvignon',
    'Vvri': 'Vitis vinifera Riesling',
    'Vrip': 'Vitis riparia PI588275'
}
design['variety_full'] = design['variety'].map(variety_map)
design['time_factor'] = 'Day' + design['time'].astype(str)
design['qc_A'] = True
design['qc_B'] = design['qc_status'].isin(['pass', 'watch'])
design['qc_C'] = design['qc_status'] == 'pass'

outfile = 'data/legacy/results_corrected/01_sample_design_and_qc/sample_design.tsv'
design.to_csv(outfile, sep='\t', index=False)

print('Sample design saved')
print()
print('Verification:')
print(f'  Total: {len(design)}')
print(f'  qc_A (all): {design["qc_A"].sum()}')
print(f'  qc_B (no fail): {design["qc_B"].sum()}')
print(f'  qc_C (pass only): {design["qc_C"].sum()}')
print()
print('By variety:')
for v, cnt in design.groupby('variety').size().items():
    print(f'  {v}: {cnt}')
print()
print('QC by variety:')
print(pd.crosstab(design['variety'], design['qc_status']))
