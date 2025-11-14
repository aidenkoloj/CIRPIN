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
parser.add_argument('--input', required=True, help='Input chunk TSV file')
parser.add_argument('--output', required=True, help='output dir')
parser.add_argument('--temp_dir', required=True, help='Temporary directory for PDB files')
args = parser.parse_args()

cluster_rep_file = args.input
out_dir = args.output
temp_pdbs = args.temp_dir

os.makedirs(out_dir, exist_ok=True)
os.makedirs(temp_pdbs, exist_ok=True)

log_file = os.path.join(out_dir, 'cluster_rep_progres_embed_log.log')
sys.stdout = open(log_file, 'w', buffering=1)
sys.stderr = sys.stdout

print('Loaded all libraries')

domids, nres, notes = [], [], []

afdb_id_last = ""
last_afdb_struc = None
structure_list = []

# Load CIRPIN model once
trained_model_path = 'CIRPIN/trained_models/CIRPIN_model/CIRPIN_model_5k_cp_epoch301.pt'
model_pg = pg.load_trained_model(trained_model=trained_model_path)
#print(f'Loaded {pg.trained_model_fp}!', flush=True)
start_time = time.time()
structure_count = 0
batch_num = 1

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
        
        temp_afdb_file = os.path.join(temp_pdbs, f"{afdb_id}.pdb")
        
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
        temp_dom_file = os.path.join(temp_pdbs, f"{ted_id}.pdb")
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
        
        if len(structure_list) == 128:
            # Process batch
            temp_list_file = os.path.join(temp_pdbs, 'temp_structure_list.txt')
            out_file = os.path.join(out_dir, f'embs_{batch_num}.pt')
            with open(temp_list_file, "w") as f:
                f.writelines(structure_list)
            
            device = "cuda" if torch.cuda.is_available() else "cpu"
            batch_size = 128
            print(f"Running progres_embed on {device} with batch_size={batch_size}", flush=True)
            pg.progres_embed(
                structurelist=temp_list_file,
                outputfile=out_file,
                fileformat="pdb",
                device=device,
                batch_size=batch_size,
                float_type=torch.float16,
                trained_model='',
                loaded_model=model_pg
            )
            
            end = time.time()
            time_took = end - start_time
            print(f"Time taken: {time_took} seconds to process {structure_count} domains", flush=True)
            
            # Cleanup old PDBs except the last afdb structure
            for pdb_file in glob.glob(os.path.join(temp_pdbs, "*.pdb")):
                if os.path.abspath(pdb_file) != os.path.abspath(last_afdb_struc):
                    os.remove(pdb_file)
            
            structure_list = []
            batch_num += 1

# Process any remaining structures
if structure_list:
    temp_list_file = os.path.join(temp_pdbs, 'temp_structure_list.txt')
    out_file = os.path.join(out_dir, f'embs_{batch_num}.pt')
    with open(temp_list_file, "w") as f:
        f.writelines(structure_list)
    
    device = "cuda" if torch.cuda.is_available() else "cpu"
    batch_size = len(structure_list)
    print(f"Running progres_embed on {device} with batch_size={batch_size}", flush=True)
    pg.progres_embed(
        structurelist=temp_list_file,
        outputfile=out_file,
        fileformat="pdb",
        device=device,
        batch_size=batch_size,
        float_type=torch.float16,
        trained_model='',
        loaded_model=model_pg
    )
    
    end = time.time()
    time_took = end - start_time
    print(f"Time taken: {time_took} seconds to process {structure_count} domains", flush=True)
    
    batch_num += 1

# Final cleanup: remove all temp files
for file in glob.glob(os.path.join(temp_pdbs, "*")):
    os.remove(file)
os.rmdir(temp_pdbs)

print(f"Successfully processed {structure_count} domains!", flush=True)