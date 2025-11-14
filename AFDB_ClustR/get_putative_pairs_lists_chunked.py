
import torch
import pickle 
import argparse
import os
import sys

''' Script to calculate all v all of AFDB cluster reps using Progres/CIRPIN and return pairs that have
a score difference above a cutoff value 
The inputs are the Progres/CIRPIN embeddings of AFDB cluster reps
The script calculates all v all similarity in chunks and then saves the ones above the cutoff to a list'''



parser = argparse.ArgumentParser(description='Get putative pairs from AFDB cluster reps using Progres/CIRPIN score differences.')
parser.add_argument('--progres_fp', default='/home/gridsan/akolodziej/TED/ted_365_chunks/combined_embs_3M_progres.pt', help='Embeddings of AFDB')
parser.add_argument('--CIRPIN_fp', default='/home/gridsan/akolodziej/TED/ted_365_chunks/combined_embs_3M_CIRPIN.pt', help='Embeddings of AFDB')
parser.add_argument('--score_cutoff', default=0.7, type=float, help='score difference cutoff between Progres/CIRPIN scores')
parser.add_argument('--chunk_size', default=500, type=int, help='chunks to calculate all v all')
parser.add_argument('--out_list', default = '/home/gridsan/akolodziej/TED/ted_365_chunks/AFDB_cluster_rep_putative_pairs', help = 'Output file')
parser.add_argument('--n_files', default =80, type = int, help = 'Number of files to split the output into')


def load_embs(fp_prog, fp_cirpin):
    ''' Load progres, cirpin embs and labels '''
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    # Load progres
    prog_data = torch.load(fp_prog, map_location=device)
    progres_emb = prog_data['embeddings']
    progres_labels = prog_data['ids']
    print('Loaded progres embs and labels!', flush=True)
    
    # Load cirpin
    cirpin_data = torch.load(fp_cirpin, map_location=device)
    cirpin_emb = cirpin_data['embeddings']
    cirpin_labels = cirpin_data['ids']
    print('Loaded cirpin embs and labels!', flush=True)
    
    return progres_emb, progres_labels, cirpin_emb, cirpin_labels

def get_putative_pairs_list(progres_pt_fp, cirpin_pt_fp, chunk_size, score_diff_cutoff):
    
    progres_emb, progres_labels, cirpin_emb, cirpin_labels = load_embs(progres_pt_fp, cirpin_pt_fp)
    
    putative_cps = []
    n = len(progres_emb)
    num_putative_pairs = 0
    chunk_num = 0
    
    for i_start in range(0, n, chunk_size):

        i_end = min(i_start + chunk_size, n)
        chunk_i_p = progres_emb[i_start:i_end]
        dot_chunk_p = chunk_i_p @ progres_emb.T
        scaled_chunk_p = (dot_chunk_p + 1) / 2

        chunk_i_c = cirpin_emb[i_start:i_end]
        dot_chunk_c = chunk_i_c @ cirpin_emb.T
        scaled_chunk_c = (dot_chunk_c + 1) / 2

        diff_chunk_scores = scaled_chunk_c - scaled_chunk_p
        mask = diff_chunk_scores > score_diff_cutoff


        indices = torch.nonzero(mask, as_tuple=False)
        # Fix indexing so that it matches labels
        indices[:, 0] += i_start
        

        indices_list = indices.tolist()

        # Save pair and save the Progres and CIRPIN scores 
        for pair in indices_list:
            p0 = pair[0]
            p1 = pair[1]
            p0_reindex = p0 - (i_start)
            putative_cps.append([
                progres_labels[p0],
                progres_labels[p1],
                float(scaled_chunk_p[p0_reindex, p1].item()),
                float(scaled_chunk_c[p0_reindex, p1].item())
            ])

        # for pair in indices_list:
        #     putative_cps.append([progres_labels[pair[0]],progres_labels[pair[1]]])

        num_putative_pairs += len(indices_list)


        #print(scaled_chunk_p[0][1])
        #print(scaled_chunk_c[0][1])
        #print(dot_chunk.shape)
        #print(diff_chunk_scores[0][1])

        chunk_num +=1
        if chunk_num % 1000 == 0:
            print(f'Number of putative pairs: {num_putative_pairs}', flush=True)
            print(f'Processed {chunk_num} chunks!', flush=True)

    return putative_cps



def save_list(cp_list, fp):

    with open(fp, 'wb') as f:
        pickle.dump(cp_list, f)

def divide_putative_cps(putative_cps, output_list, n_files):
    ''' Divide the putative cp file into smaller files '''
    num_putative_cps =len(putative_cps)
    num_files = n_files

    len_noremainder = num_putative_cps - (num_putative_cps % num_files)

    block = len_noremainder/num_files

    istart = 0
    iend = int(block)
    for i in range(0, num_files):
        if i == num_files-1:
            iend = num_putative_cps  # Last sublist takes the rest

        sublist = putative_cps[istart:iend]

        fp = f'{output_list}_{i}_.pkl'
        save_list(sublist,fp)

        istart += int(block)
        iend  += int(block)
        
def main():
    
    print('Generating putative CPs from AFDB cluster reps!')
    log_file = os.path.join('/home/gridsan/akolodziej/TED/ted_365_chunks/', 'generate_putative_AFDB_cluster_rep_CP_pairs.log')

    sys.stdout = open(log_file, 'w', buffering=1)
    sys.stderr = sys.stdout
    
    args = parser.parse_args()
    progres_fp = args.progres_fp
    cirpin_fp = args.CIRPIN_fp
    score_cutoff = args.score_cutoff
    chunk_size = args.chunk_size
    output_list = args.out_list
    n_files = args.n_files
    
    putative_cps = get_putative_pairs_list(progres_fp, cirpin_fp, chunk_size=chunk_size, score_diff_cutoff=score_cutoff)
    # save putative cps to (n_files)
    divide_putative_cps(putative_cps,output_list,n_files)
    #save_list(putative_cps, output_list)
    

    
if __name__ == "__main__":
    main()