### Script to obtain the CATH labels for all AFDB CPs

import sys
import os 
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Get CATH labels for AFDB CPs')

parser.add_argument('--cluster_rep',
                    default='/home/gridsan/akolodziej/TED/ted_365_chunks/cluster_reps_ted_365m.domain_summary.cath.globularity.taxid.tsv',
                    help='cluster_rep file with TED info')

parser.add_argument('--AFDB_CPS',
                    default='/home/gridsan/akolodziej/TED/ted_365_chunks/unique_CPs_AFDB_cluster_reps_combined.txt',
                    help='unique AFDB CPs')

parser.add_argument('--AFDB_CPS_CATH', default='/home/gridsan/akolodziej/TED/ted_365_chunks/unique_CPs_AFDB_cluster_reps_combined_with_CATH.txt',
                    help='unique AFDB CPs with CATH labels')

args = parser.parse_args()

cluster_rep_fp = args.cluster_rep
afdb_cps = args.AFDB_CPS
output = args.AFDB_CPS_CATH


def create_dict_from_tsv(cluster_rep_file, AFDB_CPs):
    ''' Create a dictionary of the TED domain info for CPs
    inputs: 
    cluster_rep_file: file containing TED domain info
    AFDB_CPs: txt file with each line containing a single AFDB TED structure
    
    Returns:
        dictionary: dictionary containing AFDB TED ID as key, TED info as values
        list: list of all AFDB_CPs
    '''
    
    ### Get a single list from false positives of structure names (afdb_tedids)
    afdb_cps = []
    
    with open(AFDB_CPs, 'r') as f:
        for line in f:
            if not line.startswith('AF'):
                continue
            line = line.strip()  # Remove whitespace/newlines
            afdb_cps.append(line)
        
    cluster_reps = pd.read_csv(cluster_rep_file, sep="\t", header=None)
    # Filter rows where the first column (index 0) is in false_pos_pairs
    cluster_reps_cps = cluster_reps[cluster_reps[0].isin(afdb_cps)]
    
    afdb_cp_dict = {}

    for _, row in cluster_reps_cps.iterrows():
        cols = row.tolist()
        ted_id = cols[0]
        afdb_id = ted_id.split("_TED")[0]

        # Parse domain residues
        domain_res = set()
        chopping = cols[3]


        domain_res = []

        for res_range in chopping.split("_"):
            res_start, res_end = res_range.split("-")
            domain_res.extend(list(range(int(res_start), int(res_end) + 1)))

        domain_res = set(domain_res)


        afdb_cp_dict[ted_id] = {
            'afdb_id': afdb_id,
            'chopping': chopping,
            'nres_dom': cols[4],
            'plddt': cols[6],
            'cath_label': cols[13],
            'tax': cols[20],
            'domain_res': domain_res
        }
    return afdb_cp_dict, afdb_cps


    
def write_output(AFDB_CPS,cluster_rep_fp, output):
    
    TED_dictionary, afdb_cp = create_dict_from_tsv(cluster_rep_fp, AFDB_CPS)
    
    with open(output, "w") as f:
        
        for i, ted_id in enumerate(afdb_cp):
            ted_info = TED_dictionary[ted_id]
            cath = ted_info['cath_label']
            f.write(f"{ted_id}\t{cath}\n")
            if i % 1000 == 0:
                print(f'{i} processed lines!')
    
    f.close()
    print('Done getting CATH labels for all unique AFDB CPs!')
        

print('Getting CATH labels!')
write_output(afdb_cps,cluster_rep_fp,output)
    
    