print('Loading libraries...')
import sys
import foldcomp
import os
import time
import glob
import argparse
import pandas as pd
import pickle
print('Loaded all libraries')
#------------------------------------

# Arg Parser------------------------------------
# ------------------------------------ Inputs ------------------------------------

parser = argparse.ArgumentParser(description='Get verified pairs AFDB cluster reps Progres/CIRPIN score differences.')
parser.add_argument('--cluster_rep', default='/home/gridsan/akolodziej/TED/ted_365_chunks/cluster_reps_ted_365m.domain_summary.cath.globularity.taxid.tsv', help='cluster_rep file with TED info')
parser.add_argument('--pairs', default='/home/gridsan/akolodziej/CIRPIN_revisions/AFDB/single_cutoffs/putative_cps_0.70_all_pairs.txt', help='putative cp pair file')


# ------------------------------------ Outputs ------------------------------------

parser.add_argument('--log', default='/home/gridsan/akolodziej/CIRPIN_revisions/AFDB/single_cutoffs/verify_pairs_05/verify_putative_pairs.log', help = 'Log')
parser.add_argument('--output_dict', default='/home/gridsan/akolodziej/CIRPIN_revisions/AFDB/single_cutoffs/verify_pairs_05/putative_pairs_ted_info_dict.pkl', help='Dictionary of TED info for putative pairs')
parser.add_argument('--output_dir', default='/home/gridsan/akolodziej/CIRPIN_revisions/AFDB/single_cutoffs/verify_pairs_05/logs', help='Output directory for logs')

args = parser.parse_args()

cluster_rep_file = args.cluster_rep
putative_pairs_fp = args.pairs
log_file = args.log



# Outputs
out_dir = args.output_dir
output_putative_pair_dict = args.output_dict

os.makedirs(out_dir, exist_ok=True)
os.makedirs(os.path.dirname(log_file), exist_ok=True)
sys.stdout = open(log_file, 'w', buffering=1)
sys.stderr = sys.stdout

# ----------------------------------------

def create_dict_from_tsv(file):
    
    putative_pairs_dict = {}

    for _, row in file.iterrows():
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


        putative_pairs_dict[ted_id] = {
            'afdb_id': afdb_id,
            'chopping': chopping,
            'nres_dom': cols[4],
            'plddt': cols[6],
            'cath_label': cols[13],
            'tax': cols[20],
            'domain_res': domain_res
        }
    return putative_pairs_dict


def main():
    tstart =  time.time()
    putative_pairs = []
    with open(putative_pairs_fp, "r") as f:
        next(f)  # skip header
        for line in f:
            parts = line.split()
            putative_pairs.append([
                parts[0],              # label_0
                parts[1],              # label_1
                float(parts[2]),       # progres_score
                float(parts[3])        # cirpin_score
            ])
    
    #  flattening 
    putative_pairs_single_list = []
    for sublist in putative_pairs:
        for item in sublist:
            putative_pairs_single_list.append(item)

    # Filter the cluster reps file for only the putative pairs
    # Read the TSV file
    ### putative_pair_cutoff_list_basename = os.path.splitext(os.path.basename(putative_pairs_fp))[0]
    ### put_pair_tsv_name = os.path.join('/home/gridsan/akolodziej/CIRPIN_revisions/AFDB/single_cutoffs/verify_pairs_05/', f'filtered_putative_pairs_TED_info_{putative_pair_cutoff_list_basename}.tsv')
    
    cluster_reps = pd.read_csv(cluster_rep_file, sep="\t", header=None)
    # Filter rows where the first column (index 0) is in putative_pairs
    
    cluster_reps_put_pairs = cluster_reps[cluster_reps[0].isin(putative_pairs_single_list)]
    
    # Write filtered results to a new file
    ### cluster_reps_put_pairs.to_csv(put_pair_tsv_name, sep="\t", index=False, header=False)


    # Create dictionary of putative pairs info for easy lookup
    putative_pairs_dict = create_dict_from_tsv(cluster_reps_put_pairs)
    with open(output_putative_pair_dict, 'wb') as f:
        pickle.dump(putative_pairs_dict, f)
    print(f'Saved dictionary to {output_putative_pair_dict}', flush=True)
    tend = time.time()
    print(f'Total time taken {(tend-tstart)/60} minutes', flush=True)
if __name__ == "__main__":
    main()