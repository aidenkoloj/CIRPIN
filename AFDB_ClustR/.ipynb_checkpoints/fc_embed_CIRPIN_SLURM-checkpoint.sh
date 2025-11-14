#!/bin/bash
#SBATCH --array=0-138 # 25K per chunk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=slurm_logs/job_%A_%a.out
#SBATCH --error=slurm_logs/job_%A_%a.err

# Load environment (adjust as needed)
module load anaconda/2023a-pytorch
source activate prog_mod
# Define chunk file
CHUNK_FILE="CIRPIN/AFDB_ClustR/cluster_reps_chunks/chunk_${SLURM_ARRAY_TASK_ID}.tsv"

# Define output directory
OUT_DIR="CIRPIN/AFDB_ClustR/embs_CIRPIN_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$OUT_DIR"

TEMP_PDB_DIR="CIRPIN/AFDB_ClustR/cluster_reps_chunks/temp_pdb_CIRPIN_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TEMP_PDB_DIR"


# Run script
python CIRPIN/AFDB_ClustR/fc_embed_afdb_clustr_reps_CIRPIN.py --input "$CHUNK_FILE" --output "$OUT_DIR" --temp_dir "$TEMP_PDB_DIR"