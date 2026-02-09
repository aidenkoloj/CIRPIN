#!/usr/bin/env python3
"""
Script to calculate TM-scores for protein pairs from CPDB dataset using MIN.
"""

import numpy as np
import pandas as pd
import os

def tmscore(q, t, cp=False):
    '''Run TM-align and get back TM-align score'''
    if cp:
        output = os.popen(f'/home/ubuntu/TM_tools/TMalign {q} {t} -cp')
    else:
        output = os.popen(f'/home/ubuntu/TM_tools/TMalign {q} {t}')
    
    tms = {"tms": []}
    parse_float = lambda x: float(x.split("=")[1].split()[0])
    
    for line in output:
        line = line.rstrip()
        if line.startswith("TM-score"): 
            tms["tms"].append(parse_float(line))
    
    if tms['tms']:
        min_tms = min(tms['tms'])
    else:
        print(f"Warning: tms['tms'] is empty, setting min_tms to 0 for {q}, {t}")
        min_tms = 0
    
    return min_tms


def convert_pdb_name(pdb):
    """Convert PDB names from '1fp3A' to '1FP3_A'"""
    pdb_code = pdb[:4].upper()
    chain = pdb[4]
    return f'{pdb_code}_{chain}'


# File paths
input_file = '/home/ubuntu/plmCP/CPDB/tm_scores_cpdb.txt'
pdb_path = '/home/ubuntu/CIRPIN/cpdb_comparison/foldseek_CPDB_comparison/cpdb_pdbs_2/'
output_file = '/home/ubuntu/plmCP/CPDB/tm_scores_cpdb_calculated.txt'

# Read data
print('Reading data...')
df = pd.read_csv(input_file, sep=' ')
print(f'Loaded {len(df)} protein pairs')

# Lists to store results
tm_scores_calculated = []
tm_scores_cp_calculated = []

# Process each row
print('Calculating TM-scores...')
for idx, row in df.iterrows():
    pdb1 = convert_pdb_name(row['Protein1'])
    pdb2 = convert_pdb_name(row['Protein2'])
    
    pdb1_full = f'{pdb_path}{pdb1}.pdb'
    pdb2_full = f'{pdb_path}{pdb2}.pdb'
    
    try:
        tm_score_min = tmscore(pdb1_full, pdb2_full, cp=False)
        tm_score_cp_min = tmscore(pdb1_full, pdb2_full, cp=True)
        
        tm_scores_calculated.append(tm_score_min)
        tm_scores_cp_calculated.append(tm_score_cp_min)
    except:
        print(f'Error on row {idx}: {pdb1} vs {pdb2}')
        tm_scores_calculated.append('Error')
        tm_scores_cp_calculated.append('Error')
    
    if (idx + 1) % 100 == 0:
        print(f'Processed {idx + 1}/{len(df)} rows')

# Save results
df['tm_score_calculated'] = tm_scores_calculated
df['tm_score_cp_calculated'] = tm_scores_cp_calculated

print(f'Saving to {output_file}...')
df.to_csv(output_file, sep='\t', index=False)
print('Done!')