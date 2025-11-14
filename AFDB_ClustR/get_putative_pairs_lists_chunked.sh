#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=slurm_logs/job_%A_%a.out
#SBATCH --error=slurm_logs/job_%A_%a.err
#SBATCH --gres=gpu:volta:1

# Load environment (adjust as needed)
module load anaconda/2023a-pytorch
source activate prog_mod

python get_putative_pairs_lists_chunked.py --chunk_size=500 --n_files=80