#!/bin/bash
source /etc/profile.d/modules.sh
module load anaconda/2023a
source activate prog_mod

INPUT_PAIRS=putative_cps_0.50_subset_3_of_3.txt

echo "My LLSUB RANK: $LLSUB_RANK"
echo "My LLSUB ID: $LLSUB_ID"
echo "TMPDIR PATH: $TMPDIR"
echo "Verifying $INPUT_PAIRS"

# Need to chunk and build dicts first
if [[ ! -d "./dictionaries" ]]; then
    echo "ERROR: Required directories ./dictionaries"
    echo "Run build_dict.sh first."
    exit 1
fi
start=$(date +%s)
copy_start=$(date +%s)

WORKDIR=/home/gridsan/akolodziej/CIRPIN_revisions/AFDB/single_cutoffs/verify_pairs_05
FOLDCOMP_DB_OG=/home/gridsan/akolodziej/Foldcomp/afdb_cluster_reps/afdb_rep_v4

mkdir -p $TMPDIR/$LLSUB_ID
cp ${FOLDCOMP_DB_OG}* $TMPDIR/$LLSUB_ID

copy_end=$(date +%s)
echo "Copy completed in $(( (copy_end - copy_start) / 60 )) minutes"

python tm_align_putative_cps.py \
    --foldcomp_db $TMPDIR/$LLSUB_ID/afdb_rep_v4 \
    --pairs ${WORKDIR}/pair_chunks/pairs_chunk_${LLSUB_RANK}.txt \
    --ted_dict ${WORKDIR}/dictionaries/putative_pairs_ted_info_dict_${LLSUB_RANK}.pkl \
    --output_dir  ${WORKDIR}/outputs \
    --temp_dir    $TMPDIR/$LLSUB_ID/temp_pdbs \
    --log         ${WORKDIR}/outputs/logs/verify_putative_pairs_0.50_s3.log \
    --output_cp_pairs               ${WORKDIR}/outputs/verified_pairs/verified_pairs.txt \
    --output_other_homologous_pairs ${WORKDIR}/outputs/other_homo/other_homologous_pairs.txt \
    --output_false_pairs            ${WORKDIR}/outputs/false_pos/false_positive_pairs.txt \
    --output_pairs_unique           ${WORKDIR}/outputs/unique_cps/unique_CPs.txt \
    --output_pairs_asym             ${WORKDIR}/outputs/asym_pairs/asym_pairs.txt \
    --chunk_id $LLSUB_RANK

end=$(date +%s)
echo "Task $LLSUB_RANK completed in $(( (end - start) / 60 )) minutes"