print('Loading libraries...')
import sys
import foldcomp
import os
import time
import glob
import argparse
import pickle
print('Loaded all libraries')

parser = argparse.ArgumentParser(description='Get verified pairs AFDB cluster reps Progres/CIRPIN score differences.')
parser.add_argument('--pairs', default='/home/gridsan/akolodziej/CIRPIN_revisions/AFDB/single_cutoffs/verify_pairs_05', help='putative cp pair file')
parser.add_argument('--chunk_id', default='0', help='Chunk ID for output naming')
parser.add_argument('--ted_dict', default='/home/gridsan/akolodziej/CIRPIN_revisions/AFDB/single_cutoffs/verify_pairs_05/putative_pairs_ted_info_dict.pkl', help='TED dictionary')
parser.add_argument('--foldcomp_db', default='/home/gridsan/akolodziej/Foldcomp/afdb_cluster_reps/afdb_rep_v4', help='Path to Foldcomp database')
parser.add_argument('--output_dir', default='/home/gridsan/akolodziej/CIRPIN_revisions/AFDB/single_cutoffs/verify_pairs_05/outputs', help='Output directory for logs')
parser.add_argument('--temp_dir', default='/home/gridsan/akolodziej/CIRPIN_revisions/AFDB/single_cutoffs/verify_pairs_05/temp_pdbs', help='Temporary directory for PDB files')
parser.add_argument('--log', default='//home/gridsan/akolodziej/CIRPIN_revisions/AFDB/single_cutoffs/verify_pairs_05/outputs/logs/verify_putative_pairs.log', help='Log')
parser.add_argument('--output_cp_pairs', default='/home/gridsan/akolodziej/CIRPIN_revisions/AFDB/single_cutoffs/verify_pairs_05/outputs/verified_pairs/verified_pairs_AFDB_cluster_reps.txt', help='Verified pairs')
parser.add_argument('--output_other_homologous_pairs', default='/home/gridsan/akolodziej/CIRPIN_revisions/AFDB/single_cutoffs/verify_pairs_05/outputs/other_homo/other_homologous_pairs_AFDB_cluster_reps.txt', help='Homologous pairs')
parser.add_argument('--output_false_pairs', default='/home/gridsan/akolodziej/CIRPIN_revisions/AFDB/single_cutoffs/verify_pairs_05/outputs/false_pos/false_positive_pairs_AFDB_cluster_reps.txt', help='false pos pairs')
parser.add_argument('--output_pairs_unique', default='/home/gridsan/akolodziej/CIRPIN_revisions/AFDB/single_cutoffs/verify_pairs_05/outputs/unique_cps/unique_CPs_AFDB_cluster_reps.txt', help='unique verified structures')
parser.add_argument('--output_pairs_asym', default='/home/gridsan/akolodziej/CIRPIN_revisions/AFDB/single_cutoffs/verify_pairs_05/outputs/asym_pairs/asym_pairs_AFDB_cluster_reps.txt', help='tm-align asym pairs')

args = parser.parse_args()

out_dir = args.output_dir
temp_pdbs = f'{args.temp_dir}_{args.chunk_id}'
putative_pairs_fp = args.pairs
log_file = args.log.replace('.log', f'_{args.chunk_id}.log')
putative_pairs_ted_dict_fp = args.ted_dict
foldcomp_db = args.foldcomp_db

output_cp_pairs = args.output_cp_pairs.replace('.txt', f'_{args.chunk_id}.txt')
output_other_homologous_pairs = args.output_other_homologous_pairs.replace('.txt', f'_{args.chunk_id}.txt')
output_false_pairs = args.output_false_pairs.replace('.txt', f'_{args.chunk_id}.txt')
output_pairs_unique = args.output_pairs_unique.replace('.txt', f'_{args.chunk_id}.txt')
output_pairs_asym = args.output_pairs_asym.replace('.txt', f'_{args.chunk_id}.txt')

os.makedirs(out_dir, exist_ok=True)
os.makedirs(temp_pdbs, exist_ok=True)
# make dirs for all the outputs/*
for path in [output_cp_pairs, output_other_homologous_pairs, output_false_pairs, 
             output_pairs_unique, output_pairs_asym, log_file]:
    os.makedirs(os.path.dirname(path), exist_ok=True)


sys.stdout = open(log_file, 'w', buffering=1)
sys.stderr = sys.stdout

print(f'Using {foldcomp_db} as foldcomp database', flush=True)


def get_already_processed(filepath):
    """Read pairs already written to an output file, return set of (q, t) tuples."""
    processed = set()
    if not os.path.exists(filepath):
        return processed
    with open(filepath, 'r') as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                processed.add((parts[0], parts[1]))
    return processed


def get_pdb_from_foldcomp(afdb_id, domain_res, ted_id):
    temp_afdb_file = os.path.join(temp_pdbs, f"{afdb_id}.pdb")
    with foldcomp.open(foldcomp_db, ids=[afdb_id]) as db:
        for (name, pdb) in db:
            with open(temp_afdb_file, "w") as f:
                f.write(pdb)
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


def tmscore(q, t, cp=False):
    if cp:
        output = os.popen(f'/home/gridsan/akolodziej/TM_tools/TMalign {q} {t} -cp')
    else:
        output = os.popen(f'/home/gridsan/akolodziej/TM_tools/TMalign {q} {t}')
    tms = {"tms": []}
    parse_float = lambda x: float(x.split("=")[1].split()[0])
    for line in output:
        line = line.rstrip()
        if line.startswith("TM-score"):
            tms["tms"].append(parse_float(line))
    if tms['tms']:
        return min(tms['tms'])
    else:
        print(f"Warning: tms['tms'] is empty, setting min_tms to 0 for {q}, {t}")
        return 0


def write_headers_if_new(filepath, header):
    """Write header only if file doesn't exist yet."""
    if not os.path.exists(filepath):
        with open(filepath, 'w') as f:
            f.write(header + '\n')


def verify_pairs(pairs_list, TED_dictionary):

    num_pairs = len(pairs_list)

    # Write headers to new files
    write_headers_if_new(output_cp_pairs,
        "query\ttarget\tprog_score\tcirpin_score\ttm_score\ttm_score_cp\ttm_diff\tis_asym")
    write_headers_if_new(output_other_homologous_pairs,
        "query\ttarget\tprog_score\tcirpin_score\ttm_score\ttm_score_cp\ttm_diff")
    write_headers_if_new(output_false_pairs,
        "query\ttarget\tprog_score\tcirpin_score\ttm_score_cp\ttm_score_cp_r")
    write_headers_if_new(output_pairs_asym,
        "query\ttarget\tprog_score\tcirpin_score\ttm_score_r\ttm_score_cp_r\ttm_diff_r")

    # Load already-processed pairs for crash recovery
    already_processed = get_already_processed(output_cp_pairs)
    already_processed |= get_already_processed(output_other_homologous_pairs)
    already_processed |= get_already_processed(output_false_pairs)
    # unique structures already written
    unique_cp_structures = set()
    if os.path.exists(output_pairs_unique):
        with open(output_pairs_unique, 'r') as f:
            for line in f:
                unique_cp_structures.add(line.strip())

    num_pairs_processed = 0
    num_skipped = 0

    # Open all output files once in append mode
    cp_f = open(output_cp_pairs, 'a', buffering=1)
    homolog_f = open(output_other_homologous_pairs, 'a', buffering=1)
    false_f = open(output_false_pairs, 'a', buffering=1)
    asym_f = open(output_pairs_asym, 'a', buffering=1)
    unique_f = open(output_pairs_unique, 'a', buffering=1)

    try:
        for pair in pairs_list:
            q = pair[0]
            t = pair[1]
            prog_score = pair[2]
            cirpin_score = pair[3]

            # Skip if already processed (crash recovery)
            if (q, t) in already_processed:
                num_skipped += 1
                num_pairs_processed += 1
                continue

            q_path = os.path.join(temp_pdbs, q + '.pdb')
            t_path = os.path.join(temp_pdbs, t + '.pdb')

            if not os.path.exists(q_path):
                if q in TED_dictionary:
                    q_info = TED_dictionary[q]
                    q_afdb_id = q_info['afdb_id']
                    q_domain_res = q_info['domain_res']
                else:
                    print(f'Warning: {q} not in TED dictionary, skipping pair', flush=True)
                    num_pairs_processed += 1
                    continue
                q_pdb = get_pdb_from_foldcomp(q_afdb_id, q_domain_res, q)
            else:
                print(f'Reusing {q}: path is {q_path}', flush=True)
                q_pdb = q_path

            if not os.path.exists(t_path):
                if t in TED_dictionary:
                    t_info = TED_dictionary[t]
                    t_afdb_id = t_info['afdb_id']
                    t_domain_res = t_info['domain_res']
                else:
                    print(f'Warning: {t} not in TED dictionary, skipping pair', flush=True)
                    num_pairs_processed += 1
                    continue
                t_pdb = get_pdb_from_foldcomp(t_afdb_id, t_domain_res, t)
            else:
                print(f'Reusing {t}: path is {t_path}', flush=True)
                t_pdb = t_path

            tm_score_cp = tmscore(q_pdb, t_pdb, cp=True)

            if tm_score_cp >= 0.5:
                tm_score = tmscore(q_pdb, t_pdb, cp=False)
                tm_diff = tm_score_cp - tm_score
                if tm_diff > 0:
                    cp_f.write(f"{q}\t{t}\t{prog_score:.2f}\t{cirpin_score:.2f}\t{tm_score:.2f}\t{tm_score_cp:.2f}\t{tm_diff:.2f}\tFalse\n")
                    print(f'Domains {q}, {t} have a tm_score -cp of {tm_score_cp:.2f}, tm score: {tm_score:.2f}, tm_diff: {tm_diff:.2f}, progres: {prog_score:.2f}, cirpin: {cirpin_score:.2f}!', flush=True)
                    if q not in unique_cp_structures:
                        unique_cp_structures.add(q)
                        unique_f.write(f"{q}\n")
                    if t not in unique_cp_structures:
                        unique_cp_structures.add(t)
                        unique_f.write(f"{t}\n")
                else:
                    tm_score_cp_r = tmscore(t_pdb, q_pdb, cp=True)
                    if tm_score_cp_r > tm_score_cp:
                        tm_score_r = tmscore(t_pdb, q_pdb, cp=False)
                        tm_diff_r = tm_score_cp_r - tm_score_r
                        if tm_diff_r > 0:
                            cp_f.write(f"{q}\t{t}\t{prog_score:.2f}\t{cirpin_score:.2f}\t{tm_score_r:.2f}\t{tm_score_cp_r:.2f}\t{tm_diff_r:.2f}\tTrue\n")
                            asym_f.write(f"{q}\t{t}\t{prog_score:.2f}\t{cirpin_score:.2f}\t{tm_score_r:.2f}\t{tm_score_cp_r:.2f}\t{tm_diff_r:.2f}\n")
                            print(f'Domains {q}, {t} have a tm_score -cp of {tm_score_cp_r:.2f}, tm score: {tm_score_r:.2f}, tm_diff: {tm_diff_r:.2f}, progres: {prog_score:.2f}, cirpin: {cirpin_score:.2f}!', flush=True)
                            if q not in unique_cp_structures:
                                unique_cp_structures.add(q)
                                unique_f.write(f"{q}\n")
                            if t not in unique_cp_structures:
                                unique_cp_structures.add(t)
                                unique_f.write(f"{t}\n")
                        else:
                            homolog_f.write(f"{q}\t{t}\t{prog_score:.2f}\t{cirpin_score:.2f}\t{tm_score:.2f}\t{tm_score_cp:.2f}\t{tm_diff:.2f}\n")
                    else:
                        homolog_f.write(f"{q}\t{t}\t{prog_score:.2f}\t{cirpin_score:.2f}\t{tm_score:.2f}\t{tm_score_cp:.2f}\t{tm_diff:.2f}\n")
            else:
                tm_score_cp_r = tmscore(t_pdb, q_pdb, cp=True)
                if tm_score_cp_r >= 0.5:
                    tm_score_r = tmscore(t_pdb, q_pdb, cp=False)
                    tm_diff_r = tm_score_cp_r - tm_score_r
                    if tm_diff_r > 0:
                        cp_f.write(f"{q}\t{t}\t{prog_score:.2f}\t{cirpin_score:.2f}\t{tm_score_r:.2f}\t{tm_score_cp_r:.2f}\t{tm_diff_r:.2f}\tTrue\n")
                        asym_f.write(f"{q}\t{t}\t{prog_score:.2f}\t{cirpin_score:.2f}\t{tm_score_r:.2f}\t{tm_score_cp_r:.2f}\t{tm_diff_r:.2f}\n")
                        print(f'Domains {q}, {t} have a tm_score -cp of {tm_score_cp_r:.2f}, tm score: {tm_score_r:.2f}, tm_diff: {tm_diff_r:.2f}, progres: {prog_score:.2f}, cirpin: {cirpin_score:.2f}!', flush=True)
                        if q not in unique_cp_structures:
                            unique_cp_structures.add(q)
                            unique_f.write(f"{q}\n")
                        if t not in unique_cp_structures:
                            unique_cp_structures.add(t)
                            unique_f.write(f"{t}\n")
                    else:
                        homolog_f.write(f"{q}\t{t}\t{prog_score:.2f}\t{cirpin_score:.2f}\t{tm_score_r:.2f}\t{tm_score_cp_r:.2f}\t{tm_diff_r:.2f}\n")
                else:
                    false_f.write(f"{q}\t{t}\t{prog_score:.2f}\t{cirpin_score:.2f}\t{tm_score_cp:.2f}\t{tm_score_cp_r:.2f}\n")

            num_pairs_processed += 1
            print(f'Processed {num_pairs_processed} out of {num_pairs} (skipped {num_skipped} already done)', flush=True)

            file_count = len(glob.glob(os.path.join(temp_pdbs, "*.pdb")))
            if file_count >= 1000:
                for f in glob.glob(os.path.join(temp_pdbs, "*.pdb")):
                    os.remove(f)

    finally:
        # Always close files cleanly even on crash
        cp_f.close()
        homolog_f.close()
        false_f.close()
        asym_f.close()
        unique_f.close()

    print(f'Done. Processed {num_pairs_processed} pairs, skipped {num_skipped} already completed.', flush=True)


def main():
    tstart = time.time()
    putative_pairs = []
    with open(putative_pairs_fp, "r") as f:
        next(f)  # skip header
        for line in f:
            parts = line.split()
            putative_pairs.append([
                parts[0],
                parts[1],
                float(parts[2]),
                float(parts[3])])

    with open(putative_pairs_ted_dict_fp, 'rb') as f:
        putative_pairs_ted_dict = pickle.load(f)

    verify_pairs(putative_pairs, putative_pairs_ted_dict)
    tend = time.time()
    print(f'Total time taken {(tend - tstart) / 60} minutes', flush=True)


if __name__ == "__main__":
    main()