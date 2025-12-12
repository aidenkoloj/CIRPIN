from pathlib import Path
import os
import pickle
from multiprocessing import Pool, cpu_count
from functools import partial
import subprocess

def tmscore(q, t, cp=False):
    """Run TM-align and get back TM-align score"""
    cmd = ['/home/ubuntu/TM_tools/TMalign', str(q), str(t)]
    if cp:
        cmd.append('-cp')
    
    try:
        # Use subprocess instead of os.popen for better performance
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        output = result.stdout
        
        tms_scores = []
        for line in output.split('\n'):
            if line.startswith("TM-score"):
                score = float(line.split("=")[1].split()[0])
                tms_scores.append(score)
        
        return min(tms_scores) if tms_scores else 0.0
    except (subprocess.TimeoutExpired, Exception) as e:
        print(f"Error processing {q} vs {t}: {e}")
        return 0.0

def process_pair(args):
    """Process a single query-target pair"""
    q, t, threshold = args
    tm = tmscore(q, t, cp=True)
    if tm > threshold:
        return (str(q), str(t), tm)
    return None

def main():
    # Load PDB lists
    CIRPIN_pdbs = [
        p for p in Path("/home/ubuntu/CIRPIN_db/ark_CPDB_analysis/CIRPIN_all_structures").iterdir()
        if p.is_file() and p.suffix == ".pdb"
    ]
    CPDB_pdbs = [
        p for p in Path('/home/ubuntu/CIRPIN_db/ark_CPDB_analysis/CPDB_pdb_structures').iterdir()
        if p.is_file() and p.suffix == ".pdb"
    ]
    
    print(f"Found {len(CPDB_pdbs)} CPDB structures and {len(CIRPIN_pdbs)} CIRPIN structures")
    print(f"Total comparisons: {len(CPDB_pdbs) * len(CIRPIN_pdbs):,}")
    
    # Create all pairs
    threshold = 0.5
    pairs = [(q, t, threshold) for q in CPDB_pdbs for t in CIRPIN_pdbs]
    
    # Use multiprocessing to parallelize
    num_processes = max(1, cpu_count() - 2)  # Leave 2 cores free for system
    print(f"Using {num_processes} processes")
    
    hits = []
    with Pool(processes=num_processes) as pool:
        # Process in chunks to show progress
        # Larger chunks for better efficiency with many cores
        chunk_size = 5000
        for i in range(0, len(pairs), chunk_size):
            chunk = pairs[i:i+chunk_size]
            # chunksize parameter helps with load balancing
            results = pool.map(process_pair, chunk, chunksize=10)
            
            # Filter out None results
            chunk_hits = [r for r in results if r is not None]
            hits.extend(chunk_hits)
            
            percent = (min(i+chunk_size, len(pairs)) / len(pairs)) * 100
            print(f"Progress: {percent:.1f}% ({min(i+chunk_size, len(pairs)):,}/{len(pairs):,}) - Found {len(chunk_hits)} hits in this chunk (Total: {len(hits)})")
    
    print(f"\nTotal hits found: {len(hits)}")
    
    # Save results
    with open('overlap_CIRPIN_scope_CPDB.pkl', "wb") as f:
        pickle.dump(hits, f)
    
    # Also save as human-readable text
    with open('overlap_CIRPIN_scope_CPDB.txt', 'w') as f:
        f.write(f"Query\tTarget\tTM-score\n")
        for q, t, tm in sorted(hits, key=lambda x: x[2], reverse=True):
            f.write(f"{q}\t{t}\t{tm:.4f}\n")

if __name__ == "__main__":
    main()
