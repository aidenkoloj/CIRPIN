print('Loading libraries...')

import sys
import foldcomp
#import torch
import os
import time
import glob
import argparse
import pandas as pd
import pickle
print('Loaded all libraries')

##### Arg Parser


parser = argparse.ArgumentParser(description='Get verified pairs AFDB cluster reps Progres/CIRPIN score differences.')
parser.add_argument('--cluster_rep', default='/home/gridsan/akolodziej/TED/ted_365_chunks/cluster_reps_ted_365m.domain_summary.cath.globularity.taxid.tsv', help='cluster_rep file with TED info')
parser.add_argument('--pairs', default='/home/gridsan/akolodziej/TED/ted_365_chunks/AFDB_cluster_rep_putative_pairs.pkl', help='putative cp pair file')
parser.add_argument('--output_dir', default='/home/gridsan/akolodziej/TED/ted_365_chunks/putative_CP_pairs_from_AFDB_cluster_reps/', help='Output directory for logs')

# Outputs
parser.add_argument('--output_cp_pairs', default='/home/gridsan/akolodziej/TED/ted_365_chunks/verified_pairs_AFDB_cluster_reps.txt', help='Verified pairs')
parser.add_argument('--output_other_homologous_pairs', default='/home/gridsan/akolodziej/TED/ted_365_chunks/other_homologous_pairs_AFDB_cluster_reps.txt', help='Homologous pairs')
parser.add_argument('--output_false_pairs', default='/home/gridsan/akolodziej/TED/ted_365_chunks/false_positive_pairs_AFDB_cluster_reps.txt', help='false pos pairs')
parser.add_argument('--output_pairs_unique', default='/home/gridsan/akolodziej/TED/ted_365_chunks/unqiue_CPs_AFDB_cluster_reps.txt', help='unique verified structures')
parser.add_argument('--temp_dir', default='/home/gridsan/akolodziej/TED/temp_pdbs_verified', help='Temporary directory for PDB files')
parser.add_argument('--log', default='/home/gridsan/akolodziej/TED/ted_365_chunks/putative_CP_pairs_from_AFDB_cluster_reps/verify_putative_pairs_log.log', help = 'Log')

args = parser.parse_args()

cluster_rep_file = args.cluster_rep
out_dir = args.output_dir
temp_pdbs = args.temp_dir
putative_pairs_fp = args.pairs

# Outputs
output_cp_pairs = args.output_cp_pairs
output_other_homologous_pairs = args.output_other_homologous_pairs
output_false_pairs = args.output_false_pairs
output_pairs_unique = args.output_pairs_unique
log_file = args.log


os.makedirs(out_dir, exist_ok=True)
os.makedirs(temp_pdbs, exist_ok=True)

sys.stdout = open(log_file, 'w', buffering=1)
sys.stderr = sys.stdout

#####


def get_pdb_from_foldcomp(afdb_id,domain_res, ted_id):
    ''' Get PDB from foldcomp and parse into its domain '''
    temp_afdb_file = os.path.join(temp_pdbs, f"{afdb_id}.pdb")
    
    with foldcomp.open("/home/gridsan/akolodziej/Foldcomp/afdb_cluster_reps/afdb_rep_v4", ids=[afdb_id]) as db:
        for (name, pdb) in db:
            with open(temp_afdb_file, "w") as f:
                f.write(pdb)
    # Extract domain from structure
    print(f'Extracting domain {ted_id}...', flush=True)
    temp_dom_file = os.path.join(temp_pdbs, f"{ted_id}.pdb")
    with open(temp_afdb_file) as af_struc:
        with open(temp_dom_file, "w") as af_dom:
            for line in af_struc:
                if line.startswith("ATOM"):
                    resnum = int(line[22:26])
                    if resnum in domain_res:
                        af_dom.write(line)
    os.remove(temp_afdb_file)
    return temp_dom_file

def tmscore(q,t, cp=False):
    ''' Run TM-align and get back TM-align score '''
    if cp:
        output = os.popen(f'/home/gridsan/akolodziej/TM_tools/TMalign {q} {t} -cp')
    else:
        output = os.popen(f'/home/gridsan/akolodziej/TM_tools/TMalign {q} {t}')
    tms = {"tms":[]}
    parse_float = lambda x: float(x.split("=")[1].split()[0])
    # print(output, flush=True)
    # print(q)
    # print(t)
    # print(f'{os.path.exists(q)}')
    # print(f'{os.path.exists(t)}')
    for line in output:
        line = line.rstrip()
        if line.startswith("TM-score"): 
            tms["tms"].append(parse_float(line))
    #error handling, set tm scor eto 0 for now. 
    if tms['tms']:
        min_tms = min(tms['tms'])
    else:
        print(f"Warning: tms['tms'] is empty, setting min_tms to 0 for {q}, {t}")
        min_tms = 0
    return min_tms

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

def verify_pairs(pairs_list, TED_dictionary):
    ''' pairs_list: list of lists of putative pairs
    TED_dictionary: dictionary of putative pairs with domain info for each pair '''
    
    # Now process each putative pair
    verified_pairs = []
    # list of unique, non duplicated structures
    unique_cp_structures = []
    # other homolog pairs
    other_homolog_pairs = []
    # false pos pairs
    false_pos_pairs = []
    # Get number of putative pairs
    num_pairs = len(pairs_list)
    num_pairs_processed = 0


    

    for pair in pairs_list:
        q = pair[0]  # query TED ID
        t = pair[1]  # target TED ID
        # Progres and cirpin scores saved to the list of putative cps too. 
        prog_score = pair[2]
        cirpin_score = pair[3]


        # Get info for query
        if q in TED_dictionary:
            q_info = TED_dictionary[q]
            q_afdb_id = q_info['afdb_id']
            q_domain_res = q_info['domain_res']
            q_plddt = q_info['plddt']
            q_cath_label = q_info['cath_label']
            q_tax = q_info['tax']
            #print(f"Query {q}: AFDB={q_afdb_id}, pLDDT={q_plddt}, CATH={q_cath_label}, Tax={q_tax}")
            #print(f"  Domain residues: {len(q_domain_res)} residues")

        q_pdb = get_pdb_from_foldcomp(q_afdb_id,q_domain_res, q)

        # Get info for target
        if t in TED_dictionary:
            t_info = TED_dictionary[t]
            t_afdb_id = t_info['afdb_id']
            t_domain_res = t_info['domain_res']
            t_plddt = t_info['plddt']
            t_cath_label = t_info['cath_label']
            t_tax = t_info['tax']
            #print(f"Target {t}: AFDB={t_afdb_id}, pLDDT={t_plddt}, CATH={t_cath_label}, Tax={t_tax}")
            #print(f"  Domain residues: {len(t_domain_res)} residues")

        t_pdb = get_pdb_from_foldcomp(t_afdb_id,t_domain_res, t)

        tm_score_cp = tmscore(q_pdb,t_pdb, cp=True)

        if tm_score_cp >= 0.5:
            tm_score = tmscore(q_pdb,t_pdb, cp=False)
            tm_diff = tm_score_cp - tm_score
            if tm_diff > 0:
                verified_pairs.append([q, t, prog_score, cirpin_score, tm_score, tm_score_cp, tm_diff])
                print(
                    f'Domains {q}, {t} have a tm_score -cp of {tm_score_cp:.2f},' 
                      f'tm score: {tm_score:.2f}, tm_diff: {tm_diff:.2f}, progres: {prog_score:.2f},'
                      f'cirpin: {cirpin_score:.2f}!', 
                      flush=True
                )
                if q not in unique_cp_structures:
                    unique_cp_structures.append(q)
                if t not in unique_cp_structures:
                    unique_cp_structures.append(t)
            else:
                other_homolog_pairs.append([q, t, prog_score, cirpin_score, tm_score, tm_score_cp, tm_diff])
            
        else:
            false_pos_pairs.append([q, t, prog_score, cirpin_score, tm_score_cp])
            os.remove(q_pdb)
            os.remove(t_pdb)
            # RUN MICAN ON FALSE POSITIVES
            #print(f'Low TM-score: {tm_score}!', flush=True)

        num_pairs_processed += 1
        print(f'Processed {num_pairs_processed} out of {num_pairs}', flush = True)

        # Check if there are 1000 files in the directory and if there are, move all these to a new subdir
        file_count = len(glob.glob(os.path.join(temp_pdbs, "*.pdb")))
        if file_count >= 1000:
            # Create a timestamped subdirectory
            timestamp = time.strftime("%Y%m%d_%H%M%S")
            archive_subdir = os.path.join(temp_pdbs, f"verified_pdbs_{timestamp}")
            os.makedirs(archive_subdir, exist_ok=True)
    
            # Move all PDB files to the archive subdirectory
            for file in glob.glob(os.path.join(temp_pdbs, "*.pdb")):
                filename = os.path.basename(file)
                os.rename(file, os.path.join(archive_subdir, filename))
    
            print(f"Moved {file_count} files to {archive_subdir}", flush=True)

    with open(output_cp_pairs, 'w', encoding='utf-8') as f:
        # Write header once
        f.write("query\ttarget\tprog_score\tcirpin_score\ttm_score\ttm_score_cp\ttm_diff\n")
    
        # Write each cp pair
        for pair in verified_pairs:
            query, target, prog, cirpin, tm_score, tm_score_cp, tm_diff = pair
            f.write(f"{query}\t{target}\t{prog:.2f}\t{cirpin:.2f}\t{tm_score:.2f}\t{tm_score_cp:.2f}\t{tm_diff:.2f}\n")
    
        # Add a blank line before summary section
        f.write("\n")
        num_cp_pairs = len(verified_pairs)
        hits_percentage_cps = (num_cp_pairs / num_pairs) * 100
    
        # Write summary stats, each on its own line
        f.write(f"Number of putative pairs: {num_pairs}\n")
        f.write(f"Number of CP pairs (TM-score difference > 0): {num_cp_pairs}\n")
        f.write(f"Percent of putative CPs that are verified CPs: {hits_percentage_cps:.2f}%\n")
        
    with open(output_other_homologous_pairs, 'w', encoding='utf-8') as f:
        # Write header once
        f.write("query\ttarget\tprog_score\tcirpin_score\ttm_score\ttm_score_cp\ttm_diff\n")
    
        # Write each verified pair
        for pair in other_homolog_pairs:
            query, target, prog, cirpin, tm_score, tm_score_cp, tm_diff = pair
            f.write(f"{query}\t{target}\t{prog:.2f}\t{cirpin:.2f}\t{tm_score:.2f}\t{tm_score_cp:.2f}\t{tm_diff:.2f}\n")
    
        # Add a blank line before summary section
        f.write("\n")
        num_homologous_pairs = len(other_homolog_pairs)
        hits_percentage_homologs = (num_homologous_pairs / num_pairs) * 100
    
        # Write summary stats, each on its own line
        f.write(f"Number of putative pairs: {num_pairs}\n")
        f.write(f"Number of non-CP homologous pairs (TM-score difference == 0): {num_homologous_pairs}\n")
        f.write(f"Percent of putative CPs that are non-CP homologous pairs: {hits_percentage_homologs:.2f}%\n")

    with open(output_false_pairs, 'w', encoding='utf-8') as f:
        # Write header once
        f.write("query\ttarget\tprog_score\tcirpin_score\ttm_score_cp\n")
    
        # Write each verified pair
        for pair in false_pos_pairs:
            query, target, prog, cirpin, tm_score_cp = pair
            f.write(f"{query}\t{target}\t{prog:.2f}\t{cirpin:.2f}\t{tm_score_cp:.2f}\n")
    
        # Add a blank line before summary section
        f.write("\n")
        num_false_pairs = len(false_pos_pairs)
        false_pairs_percentage = (num_false_pairs / num_pairs) * 100
    
        # Write summary stats, each on its own line
        f.write(f"Number of putative pairs: {num_pairs}\n")
        f.write(f"Number of pairs below TM-align -cp score of 0.5: {num_false_pairs}\n")
        f.write(f"Percent of putative CPs that are false positives: {false_pairs_percentage:.2f}%\n")

    with open(output_pairs_unique, 'w', encoding='utf-8') as f:
        for i in unique_cp_structures:
            #f.write(pair + '\n')
            f.write(f"{i}\n")

    print(f'Saved verified pairs of CPS: {output_cp_pairs}, {output_other_homologous_pairs} {output_false_pairs}, {output_pairs_unique}!', flush = True)
    
def main():
    # Load putative pairs
    with open(putative_pairs_fp, 'rb') as f:
        putative_pairs = pickle.load(f)
    
    # make putative pairs into a single list
    putative_pairs_single_list =[]
    for sublist in putative_pairs:
        for item in sublist:
            putative_pairs_single_list.append(item)

    # Filter the cluster reps file for only the putative pairs
    # Read the TSV file
    putative_pair_cutoff_list_basename = os.path.splitext(os.path.basename(putative_pairs_fp))[0]
    put_pair_tsv_name = os.path.join('/home/gridsan/akolodziej/TED/ted_365_chunks/', f'filtered_putative_pairs_TED_info_{putative_pair_cutoff_list_basename}.tsv')
    
    cluster_reps = pd.read_csv(cluster_rep_file, sep="\t", header=None)
    # Filter rows where the first column (index 0) is in putative_pairs
    cluster_reps_put_pairs = cluster_reps[cluster_reps[0].isin(putative_pairs_single_list)]
    # Write filtered results to a new file
    cluster_reps_put_pairs.to_csv(put_pair_tsv_name, sep="\t", index=False, header=False)


    # Create dictionary of putative pairs info for easy lookup
    putative_pairs_dict = create_dict_from_tsv(cluster_reps_put_pairs)

    # verify pairs
    verify_pairs(putative_pairs, putative_pairs_dict)

if __name__ == "__main__":
    main()