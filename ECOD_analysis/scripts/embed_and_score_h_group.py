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
    
    batch_size = 32
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
    Returns results as list of (struct1, struct2, progres_score, cirpin_score)
    """
    progres_embs = progres_embs.to(device)
    cirpin_embs = cirpin_embs.to(device)
    
    progres_sim = ((progres_embs @ progres_embs.T) + 1) / 2
    cirpin_sim = ((cirpin_embs @ cirpin_embs.T) + 1) / 2
    
    i_idx, j_idx = torch.triu_indices(len(progres_embs), len(progres_embs), offset=1, device=device)
    
    progres_scores = progres_sim[i_idx, j_idx].cpu().numpy()
    cirpin_scores = cirpin_sim[i_idx, j_idx].cpu().numpy()

    print(f"\nScores calculated, zipping up into results.")
    # results = []
    # for idx, (i, j) in enumerate(zip(i_idx, j_idx)):
    #     struct1 = progres_ids[i.item()]
    #     struct2 = progres_ids[j.item()]
    #     results.append((struct1, struct2, progres_scores[idx], cirpin_scores[idx]))
    
    # return results

    progres_ids = np.array(progres_ids)  # Convert to numpy array
    struct1_ids = progres_ids[i_idx.cpu().numpy()]
    struct2_ids = progres_ids[j_idx.cpu().numpy()]
    
    results = list(zip(struct1_ids, struct2_ids, progres_scores, cirpin_scores))
        
    return results
    
    

def filter_results(results, progres_threshold=0.6, cirpin_threshold=0.9):
    """
    Filter results: Progres < threshold AND CIRPIN > threshold
    """
    # filtered = [
    #     r for r in results
    #     if r[2] < progres_threshold and r[3] > cirpin_threshold
    # ]
    print('Filtering results...')
    struct1, struct2, progres, cirpin = zip(*results)  # Unpack once
    progres = np.array(progres)  # Convert to numpy arrays (once)
    cirpin = np.array(cirpin)
    mask = (progres < progres_threshold) & (cirpin > cirpin_threshold)  # ✅ NumPy vectorized comparison
    filtered = [(struct1[i], struct2[i], progres[i], cirpin[i]) for i in np.where(mask)[0]]
    print('Done filtering results.')

    return filtered

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
    cirpin_embs, cirpin_ids = load_embeddings(embeddings['CIRPIN'], device=device)
    
    # Step 3: All-vs-all comparison
    print(f"\nRunning all-vs-all comparison...")
    all_results = all_vs_all_comparison(progres_embs, progres_ids, cirpin_embs, cirpin_ids, device=device)
    print(f"Total comparisons: {len(all_results)}")
    
    # Step 4: Save all unfiltered scores
    print(f"\nSaving all unfiltered scores...")
    all_scores_file = group_dir / f"{group_name}_all_scores.txt"
    with open(all_scores_file, 'w') as f:
        f.write("struct1\tstruct2\tProgres_score\tCIRPIN_score\n")
        for struct1, struct2, progres, cirpin in all_results:
            f.write(f"{struct1}\t{struct2}\t{progres:.4f}\t{cirpin:.4f}\n")
    print(f"All scores saved to: {all_scores_file}")
    
    # Step 5: Filter results
    print(f"\nFiltering results (Progres < 0.6 AND CIRPIN > 0.9)...")
    filtered = filter_results(all_results, progres_threshold=0.6, cirpin_threshold=0.9)
    print(f"Filtered results: {len(filtered)}")
    
    # Step 6: Save filtered results
    results_file = group_dir / f"{group_name}_filtered_results.txt"
    with open(results_file, 'w') as f:
        f.write("struct1\tstruct2\tProgres_score\tCIRPIN_score\n")
        for struct1, struct2, progres, cirpin in filtered:
            f.write(f"{struct1}\t{struct2}\t{progres:.4f}\t{cirpin:.4f}\n")
    
    print(f"Filtered results saved to: {results_file}")
    
    return filtered

def main():
    """
    Main loop: process all structure_list files in structure_lists_greaterthan10/
    """
    structure_lists_dir = Path('/home/gridsan/akolodziej/ECOD_2/structure_lists/structure_lists_greater_than_10')
    output_dir = Path('/home/gridsan/akolodziej/ECOD_2/ECOD_results')
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    
    print(f"Using device: {device}")
    print(f"Input directory: {structure_lists_dir}")
    print(f"Output directory: {output_dir}")
    
    if not structure_lists_dir.exists():
        print(f"Error: {structure_lists_dir} not found", file=sys.stderr)
        return
    
    # Get all structure_list files
    structure_list_files = sorted(structure_lists_dir.glob('structure_list_*.txt'))
    print(f"\nFound {len(structure_list_files)} structure_list files\n")
    
    if not structure_list_files:
        print("No structure_list files found!", file=sys.stderr)
        return
    
    results_summary = {}
    
    # Process each group
    for i, structure_list_file in enumerate(structure_list_files, 1):
        group_name = structure_list_file.stem.replace('structure_list_', 'group_')
        
        print(f"\n[{i}/{len(structure_list_files)}] Processing {group_name}")
        
        try:
            filtered = process_group(
                str(structure_list_file),
                group_name,
                output_dir,
                device=device
            )
            
            if filtered is not None:
                results_summary[group_name] = len(filtered)
            else:
                results_summary[group_name] = 'FAILED'
        
        except Exception as e:
            print(f"Error processing {group_name}: {e}", file=sys.stderr)
            results_summary[group_name] = 'ERROR'
        
        # Print completion status
        status = results_summary[group_name]
        if isinstance(status, int):
            print(f"✓ {group_name} finished: {status} filtered pairs", flush=True)
        else:
            print(f"✗ {group_name} finished: {status}", flush=True)
    
    # Print summary
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