"""
Script to calculate TM-scores for protein pairs from CPDB dataset using MIN.
"""

#!/usr/bin/env python3

import numpy as np
import pandas as pd
import os
from multiprocessing import Process, Queue



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


def calculate_tm_with_timeout(pdb1_full, pdb2_full, cp, result_queue):
    """Helper function to run tmscore in a separate process"""
    try:
        result = tmscore(pdb1_full, pdb2_full, cp=cp)
        result_queue.put(result)
    except Exception as e:
        result_queue.put(None)


def run_with_timeout(pdb1_full, pdb2_full, cp, timeout_seconds=10):
    """Run tmscore with a timeout"""
    result_queue = Queue()
    process = Process(target=calculate_tm_with_timeout, args=(pdb1_full, pdb2_full, cp, result_queue))
    process.start()
    process.join(timeout=timeout_seconds)
    
    if process.is_alive():
        # Timeout occurred
        process.terminate()
        process.join()
        return 'Error'
    
    if not result_queue.empty():
        return result_queue.get()
    return 'Error'


# File paths
input_file = '/home/ubuntu/plmCP/CPDB/tm_scores_cpdb.txt'
pdb_path = '/home/ubuntu/CIRPIN/cpdb_comparison/foldseek_CPDB_comparison/cpdb_pdbs_2/'
output_file = '/home/ubuntu/plmCP/CPDB/tm_scores_cpdb_calculated_FIXED_pdbs.txt'

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
    
    print(f'Row {idx + 1}: {pdb1} vs {pdb2}')
    
    tm_score_min = run_with_timeout(pdb1_full, pdb2_full, cp=False, timeout_seconds=10)
    tm_score_cp_min = run_with_timeout(pdb1_full, pdb2_full, cp=True, timeout_seconds=10)
    
    if tm_score_min == 'Error':
        print(f'  Timeout/Error on regular TM-score')
    if tm_score_cp_min == 'Error':
        print(f'  Timeout/Error on CP TM-score')
    
    tm_scores_calculated.append(tm_score_min)
    tm_scores_cp_calculated.append(tm_score_cp_min)
    
    if (idx + 1) % 100 == 0:
        print(f'=== Processed {idx + 1}/{len(df)} rows ===')

# Save results
df['tm_score_calculated'] = tm_scores_calculated
df['tm_score_cp_calculated'] = tm_scores_cp_calculated

print(f'Saving to {output_file}...')
df.to_csv(output_file, sep='\t', index=False)
print('Done!')