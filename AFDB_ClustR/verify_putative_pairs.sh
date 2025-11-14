#!/bin/bash
#SBATCH --array=0-79
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=slurm_logs/job_%A_%a.out
#SBATCH --error=slurm_logs/job_%A_%a.err

# Load environment (adjust as needed)
module load anaconda/2023a-pytorch
source activate prog_mod

PAIR_FILE="/home/gridsan/akolodziej/TED/ted_365_chunks/AFDB_cluster_rep_putative_pairs_${SLURM_ARRAY_TASK_ID}_.pkl"

## Output Files
CP_PAIRS="/home/gridsan/akolodziej/TED/ted_365_chunks/verified_pairs_AFDB_cluster_reps_${SLURM_ARRAY_TASK_ID}_.txt"
HOMO_PAIRS="/home/gridsan/akolodziej/TED/ted_365_chunks/other_homologous_pairs_AFDB_cluster_reps_${SLURM_ARRAY_TASK_ID}_.txt"
FP_PAIRS="/home/gridsan/akolodziej/TED/ted_365_chunks/false_positive_pairs_AFDB_cluster_reps_${SLURM_ARRAY_TASK_ID}_.txt"
UNI_PAIRS="/home/gridsan/akolodziej/TED/ted_365_chunks/unqiue_CPs_AFDB_cluster_reps_${SLURM_ARRAY_TASK_ID}_.txt"
TEMP_DIR="/home/gridsan/akolodziej/TED/temp_pdbs_verified_${SLURM_ARRAY_TASK_ID}"
LOG="/home/gridsan/akolodziej/TED/ted_365_chunks/putative_CP_pairs_from_AFDB_cluster_reps/verify_putative_pairs_log_${SLURM_ARRAY_TASK_ID}_.log"

python /home/gridsan/akolodziej/TED/ted_365_chunks/verify_putative_pairs.py \
--pairs "$PAIR_FILE" \
--output_cp_pairs "$CP_PAIRS" \
--output_other_homologous_pairs "$HOMO_PAIRS" \
--output_false_pairs "$FP_PAIRS" \
--output_pairs_unique "$UNI_PAIRS" \
--log "$LOG" \
--temp_dir "$TEMP_DIR"