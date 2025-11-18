''' Download individual PDZ domains and save their Progres/CIRPIN individual embeddings '''

import sys
sys.path.append('/home/gridsan/akolodziej/progres_link')

# Import CIRPIN/progres with foldcomp functionality
import progres as pg
import foldcomp
import torch
import os
import time
import glob
import argparse

print('Loading libraries...', flush=True)

parser = argparse.ArgumentParser(description='Process domain embeddings for a chunk.')
parser.add_argument('--cluster_rep', default='/home/gridsan/akolodziej/TED/ted_365_chunks/AFDB_PDZ/cluster_reps_ted_365m_PDZ_domains.tsv', help='cluster_rep file with TED info')
parser.add_argument('--model', default='Progres', help='Embed domains using either Progres or CIRPIN')
parser.add_argument('--output_dir', default='/home/gridsan/akolodziej/TED/ted_365_chunks/AFDB_PDZ/', help='Output dir')
parser.add_argument('--pdz_pdbs', default='/home/gridsan/akolodziej/TED/ted_365_chunks/AFDB_PDZ/PDZ_pdbs', help='Output dir for PDZ PDB files')

args = parser.parse_args()

cluster_rep_file = args.cluster_rep
pdz_pdbs = args.pdz_pdbs
model_setting = args.model
out_dir = args.output_dir + f'embs_{model_setting}_single'


os.makedirs(out_dir, exist_ok=True)
os.makedirs(pdz_pdbs, exist_ok=True)

log_file = os.path.join(out_dir, f'download_PDZ_{model_setting}_embed_log.log')
sys.stdout = open(log_file, 'w', buffering=1)
sys.stderr = sys.stdout

print('Loaded all libraries', flush=True)

domids, nres, notes = [], [], []

afdb_id_last = ""
last_afdb_struc = None
structure_list = []

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# Load model once

if model_setting == 'CIRPIN':
    cirpin = '/home/gridsan/akolodziej/progres_link/trained_models/CIRPIN/CIRPIN_model/model_5k_cp_epoch301.pt'
    loaded_trained_model = pg.load_trained_model(device=device, trained_model=cirpin)
else:
    loaded_trained_model =pg.load_trained_model(device=device)


#print(f'Loaded {pg.trained_model_fp}!', flush=True)
start_time = time.time()
structure_count = 0

with open(cluster_rep_file) as f:
    for line in f:
        cols = line.rstrip().split()
        
        ted_id = cols[0]
        chopping = cols[3]
        nres_dom = cols[4]
        plddt = cols[6]
        cath_label = cols[13]
        tax = cols[20]
        afdb_id = ted_id.split("_TED")[0]
        
        # Parse domain residues
        domain_res = []
        for res_range in chopping.split("_"):
            res_start, res_end = res_range.split("-")
            domain_res.extend(list(range(int(res_start), int(res_end) + 1)))
        domain_res = set(domain_res)
        
        temp_afdb_file = os.path.join(pdz_pdbs, f"{afdb_id}.pdb")
        
        # Only fetch structure if it's a new AFDB ID or not cached
        if afdb_id != afdb_id_last or not os.path.exists(temp_afdb_file):
            print(f'Fetching structure for {afdb_id}...', flush=True)
            with foldcomp.open("/home/gridsan/akolodziej/Foldcomp/afdb_cluster_reps/afdb_rep_v4", ids=[afdb_id]) as db:
                for (name, pdb) in db:
                    with open(temp_afdb_file, "w") as f:
                        f.write(pdb)
            afdb_id_last = afdb_id
            last_afdb_struc = temp_afdb_file
        else:
            print(f'Reusing cached structure for {afdb_id}', flush=True)
        
        # Extract domain from structure
        print(f'Extracting domain {ted_id}...', flush=True)
        temp_dom_file = os.path.join(pdz_pdbs, f"{ted_id}.pdb")
        
        if not os.path.exists(temp_dom_file):
            with open(temp_afdb_file) as af_struc:
                with open(temp_dom_file, "w") as af_dom:
                    for line2 in af_struc:
                        if line2.startswith("ATOM"):
                            resnum = int(line2[22:26])
                            if resnum in domain_res:
                                af_dom.write(line2)
        
        # Process metadata
        cath_label = "N/A" if cath_label == "-" else cath_label
        tax = "N/A" if tax == "-" else tax
        assert "-" not in [afdb_id, chopping, plddt, cath_label, tax]
        
        # Append data
        domids.append(ted_id)
        nres.append(int(nres_dom))
        notes.append(f"{afdb_id} {chopping} - pLDDT {plddt} - {cath_label} - {tax}")
        
        structure_list.append(f"{temp_dom_file} {ted_id} {notes[-1]}\n")
        
        structure_count += 1
        print(f"Processed {structure_count}", flush=True)

        domain_emb = pg.embed_structure(querystructure=temp_dom_file, device=device,model=loaded_trained_model)
        torch.save(domain_emb, os.path.join(out_dir, f"{ted_id}.pt"))
                        
            
    
    end = time.time()
    time_took = end - start_time
    print(f"Time taken: {time_took} seconds to process {structure_count} domains", flush=True)
    

print(f"Successfully processed {structure_count} domains!", flush=True)
