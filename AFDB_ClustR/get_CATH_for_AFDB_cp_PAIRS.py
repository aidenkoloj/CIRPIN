### Script to obtain the CATH labels for all AFDB CPs PAIRS
### NOT FOR UNIQUE; for PAIRS

import sys
import os 
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Get CATH labels for AFDB CPs')

parser.add_argument('--cluster_rep',
                    default='/home/gridsan/akolodziej/TED/ted_365_chunks/cluster_reps_ted_365m.domain_summary.cath.globularity.taxid.tsv',
                    help='cluster_rep file with TED info')

parser.add_argument('--AFDB_CPS_PAIRS',
                    default='/home/gridsan/akolodziej/TED/ted_365_chunks/verified_pairs_AFDB_cluster_reps_combined.txt',
                    help='AFDB CP PAIRS')

parser.add_argument('--AFDB_CPS_PAIRS_CATH', default='/home/gridsan/akolodziej/TED/ted_365_chunks/verified_pairs_AFDB_cluster_reps_combined_with_CATH.txt',
                    help='unique AFDB CP PAIRS with CATH labels')

args = parser.parse_args()

cluster_rep_fp = args.cluster_rep
afdb_cps = args.AFDB_CPS_PAIRS
output = args.AFDB_CPS_PAIRS_CATH


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
            line = line.strip()
            line_parts = line.split()
            afdb_cps.append(line_parts[0])
            afdb_cps.append(line_parts[1])
    print(f'All afdb_cps: {len(afdb_cps)}')
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
        f.write("query\tcath_q\ttarget\tcath_t\tprog_score\tcirpin_score\ttm_score\ttm_score_cp\tm_diff\n")
        with open(AFDB_CPS, "r") as db:
            next(db) # skip header
            for i, line in enumerate(db):
                if not line.startswith('AF'):
                    continue
                columns = line.strip().split('\t')
            
                q = columns[0]
                t = columns[1]
                
                ted_info_q = TED_dictionary[q]
                cath_q = ted_info_q['cath_label']
                
                ted_info_t = TED_dictionary[t]
                cath_t = ted_info_t['cath_label']
                
                prog = float(columns[2])
                cirpin = float(columns[3])
                tm_score = float(columns[4])
                tm_cp = float(columns[5])
                tm_diff = float(columns[6])
                
                f.write(f"{q}\t{cath_q}\t{t}\t{cath_t}\t{prog:.2f}\t{cirpin:.2f}\t{tm_score:.2f}\t{tm_cp:.4f}\t{tm_diff:.4f}\n")
        
       
                if i % 1000 == 0:
                    print(f'{i} processed lines!')
    
    f.close()
    print('Done getting CATH labels for all unique AFDB CPs!')
        

print('Getting CATH labels!')
write_output(afdb_cps,cluster_rep_fp,output)
    
    