#!/bin/bash
#SBATCH -a 1-3
#SBATCH --job-name=ECOD_h_group_scoring_savenp
#SBATCH --output=logs/ECOD_h_group_scoring_%a.log
#SBATCH --error=logs/ECOD_h_group_scoring_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:volta:1
#SBATCH --time=24:00:00
#SBATCH --mem=64G

source /etc/profile.d/modules.sh
module load anaconda/2023a
source activate prog_mod

mkdir -p logs

STRUC_LIST_DIR="/home/gridsan/akolodziej/ECOD_2/structure_lists/structure_lists_greater_than_10/batch${SLURM_ARRAY_TASK_ID}"

echo "Job ${SLURM_ARRAY_TASK_ID} started at $(date)"
echo "Processing batch dir: ${STRUC_LIST_DIR}"

python -u embed_and_score_h_group_savenp.py "${STRUC_LIST_DIR}"

echo "Job ${SLURM_ARRAY_TASK_ID} completed at $(date)"