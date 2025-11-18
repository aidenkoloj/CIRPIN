#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=slurm_logs/job_%A_%a.out
#SBATCH --error=slurm_logs/job_%A_%a.err

# Load environment (adjust as needed)
module load anaconda/2023a-pytorch
source activate prog_mod

python /home/gridsan/akolodziej/TED/ted_365_chunks/AFDB_PDZ/download__embed_individual_afdb_PDZ_domains.py --model "Progres"