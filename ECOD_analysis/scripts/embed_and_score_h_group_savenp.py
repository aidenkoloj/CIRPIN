import torch
import sys
import time
from pathlib import Path
import progres as pg
import numpy as np

def run_embeddings(structure_list_file, output_dir, device='cuda'):
    """
    Run Progres and CIRPIN embeddings.
    Returns paths to embedding files.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    batch_size = 64
    embeddings = {}
    
    for model_name, model_path in [('Progres', None), ('CIRPIN', '/home/gridsan/akolodziej/progres_link/trained_models/CIRPIN/CIRPIN_model/model_5k_cp_epoch301.pt')]:
        try:
            print(f'\nEmbedding with {model_name}...')
            start = time.time()
            
            loaded_model = pg.load_trained_model(
                device=device,
                **({"trained_model": model_path} if model_path else {})
            )
            
            out_file = str(output_dir / f'embs_{model_name}.pt')
            
            pg.progres_embed(
                structurelist=structure_list_file,
                outputfile=out_file,
                fileformat='pdb',
                device=str(device),
                batch_size=batch_size,
                float_type=torch.float16,
                trained_model='',
                loaded_model=loaded_model,
            )
            
            embeddings[model_name] = out_file
            print(f'{model_name} done in {time.time() - start:.1f}s')
        
        except Exception as e:
            print(f"Error with {model_name}: {e}", file=sys.stderr)
    
    return embeddings

def load_embeddings(emb_file, device='cpu'):
    """
    Load embeddings from .pt file.
    Expected format: {'embeddings': tensor, 'ids': list}
    Returns tuple of (embeddings_tensor, ids_list)
    """
    data = torch.load(emb_file, map_location=device)
    
    if isinstance(data, dict) and 'embeddings' in data and 'ids' in data:
        embeddings = data['embeddings'].float()
        ids = data['ids']
        return embeddings, ids
    else:
        raise ValueError(f"Unexpected embedding file format in {emb_file}")

def all_vs_all_comparison(progres_embs, progres_ids, cirpin_embs, cirpin_ids, device='cpu'):
    """
    Compare all structures pairwise using matrix multiplication.
    Similarity = ((emb @ emb.T) + 1) / 2

    Returns a dict of numpy arrays:
      {
        'struct1':  shape (N_pairs,)  — string IDs
        'struct2':  shape (N_pairs,)  — string IDs
        'progres':  shape (N_pairs,)  — float32 scores
        'cirpin':   shape (N_pairs,)  — float32 scores
      }
    """
    progres_embs = progres_embs.to(device)
    cirpin_embs  = cirpin_embs.to(device)
    
    progres_sim = ((progres_embs @ progres_embs.T) + 1) / 2
    cirpin_sim  = ((cirpin_embs  @ cirpin_embs.T)  + 1) / 2
    
    n = len(progres_embs)
    i_idx, j_idx = torch.triu_indices(n, n, offset=1, device=device)
    
    # Move scores to CPU as float32 numpy arrays — no Python list ever created
    progres_scores = progres_sim[i_idx, j_idx].cpu().numpy().astype(np.float32)
    cirpin_scores  = cirpin_sim[i_idx,  j_idx].cpu().numpy().astype(np.float32)

    ids_np   = np.array(progres_ids)
    i_np     = i_idx.cpu().numpy()
    j_np     = j_idx.cpu().numpy()

    return {
        'struct1': ids_np[i_np],
        'struct2': ids_np[j_np],
        'progres': progres_scores,
        'cirpin':  cirpin_scores,
    }

def save_all_scores(scores, out_file):
    """
    Save all unfiltered scores dict to a .npz file (fast, compact).
    """
    np.savez(out_file, **scores)
    print(f"All scores saved to: {out_file}.npz")

def filter_and_save(scores, results_file, progres_threshold=0.6, cirpin_threshold=0.9):
    """
    Apply (progres < threshold) & (cirpin > threshold) via numpy boolean mask,
    then write the passing rows to a tab-separated .txt.
    No Python list is created at any point.
    """
    print('Filtering results...')
    mask = (scores['progres'] < progres_threshold) & (scores['cirpin'] > cirpin_threshold)
    n_filtered = int(mask.sum())

    with open(results_file, 'w') as f:
        f.write("struct1\tstruct2\tProgres_score\tCIRPIN_score\n")
        # Write directly from masked numpy arrays — vectorised, no intermediate list
        for s1, s2, p, c in zip(
            scores['struct1'][mask],
            scores['struct2'][mask],
            scores['progres'][mask],
            scores['cirpin'][mask],
        ):
            f.write(f"{s1}\t{s2}\t{p:.4f}\t{c:.4f}\n")

    print(f"Done filtering. {n_filtered} pairs passed. Saved to: {results_file}")
    return n_filtered

def process_group(structure_list_file, group_name, output_dir, device='cuda'):
    """
    Process a single group: embed, compare, filter.
    """
    group_dir = Path(output_dir) / group_name
    group_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\n{'='*60}")
    print(f"Processing {group_name}")
    print(f"{'='*60}")
    
    # Step 1: Run embeddings
    print(f"\nRunning embeddings...")
    embeddings = run_embeddings(structure_list_file, group_dir, device=device)
    
    if 'Progres' not in embeddings or 'CIRPIN' not in embeddings:
        print("Skipping: Could not generate all embeddings", file=sys.stderr)
        return None
    
    # Step 2: Load embeddings
    print(f"\nLoading embeddings...")
    progres_embs, progres_ids = load_embeddings(embeddings['Progres'], device=device)
    cirpin_embs, cirpin_ids   = load_embeddings(embeddings['CIRPIN'],  device=device)
    
    # Step 3: All-vs-all comparison — returns dict of numpy arrays, no Python list
    print(f"\nRunning all-vs-all comparison...")
    scores = all_vs_all_comparison(progres_embs, progres_ids, cirpin_embs, cirpin_ids, device=device)
    print(f"Total comparisons: {len(scores['progres'])}")
    
    # Step 4: Save all unfiltered scores as .npz
    print(f"\nSaving all unfiltered scores...")
    all_scores_file = group_dir / f"{group_name}_all_scores"
    save_all_scores(scores, str(all_scores_file))
    
    # Step 5: Filter and save results
    print(f"\nFiltering results (Progres < 0.6 AND CIRPIN > 0.9)...")
    results_file = group_dir / f"{group_name}_filtered_results.txt"
    n_filtered = filter_and_save(scores, results_file)
    
    return n_filtered

def main():
    """
    Main loop: process all structure_list files in structure_lists_greaterthan10/
    """
    
    if len(sys.argv) < 2:
        print("Usage: python embed_and_score_h_group_savenp.py <batch_dir>", file=sys.stderr)
        sys.exit(1)

    structure_lists_dir = Path(sys.argv[1])
    output_dir = Path('/home/gridsan/akolodziej/ECOD_2/ECOD_results')
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    
    print(f"Using device: {device}")
    print(f"Input directory: {structure_lists_dir}")
    print(f"Output directory: {output_dir}")
    
    if not structure_lists_dir.exists():
        print(f"Error: {structure_lists_dir} not found", file=sys.stderr)
        return
    
    structure_list_files = sorted(structure_lists_dir.glob('structure_list_*.txt'))
    print(f"\nFound {len(structure_list_files)} structure_list files\n")
    
    if not structure_list_files:
        print("No structure_list files found!", file=sys.stderr)
        return
    
    results_summary = {}
    
    for i, structure_list_file in enumerate(structure_list_files, 1):
        group_name = structure_list_file.stem.replace('structure_list_', 'group_')
        
        print(f"\n[{i}/{len(structure_list_files)}] Processing {group_name}")
        
        try:
            n_filtered = process_group(
                str(structure_list_file),
                group_name,
                output_dir,
                device=device
            )
            results_summary[group_name] = n_filtered if n_filtered is not None else 'FAILED'
        
        except Exception as e:
            print(f"Error processing {group_name}: {e}", file=sys.stderr)
            results_summary[group_name] = 'ERROR'
        
        status = results_summary[group_name]
        if isinstance(status, int):
            print(f"✓ {group_name} finished: {status} filtered pairs", flush=True)
        else:
            print(f"✗ {group_name} finished: {status}", flush=True)
    
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    for group_name, count in results_summary.items():
        if isinstance(count, int):
            print(f"{group_name}: {count} filtered pairs")
        else:
            print(f"{group_name}: {count}")

if __name__ == '__main__':
    main()