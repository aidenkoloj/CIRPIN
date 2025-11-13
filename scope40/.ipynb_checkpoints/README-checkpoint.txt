Protein Pair Data

This directory contains deduplicated datasets of circular permutation (CP) pairs, 
extension/insertion/rewiring pairs, and similar secondary content pairs, all found in SCOPe 40%. 

================================================================================

FILES

putative_pairs_list_score_diff_0.3_progres_cutoff_0.6_cirpin_cutoff_0.9_.pkl
	Pairs of SCOPe proteins with Progres score < 0.6, CIRPIN score > 0.9

cp_pairs_scope40_dedup_asym.tsv
    Contains all 1,968 verified circular permutation pairs.

other_homologous_pairs_scope40_dedup_asym.tsv
    Contains 109 pairs related by insertions, extensions, or rewirings.

false_pos_pairs_scope40_dedup_asym.tsv
    Contains pairs with similar secondary content that do not represent true 
    circular permutations.

================================================================================

SCRIPTS

deduplicate.ipynb
    Notebook to deduplicate pairs of CPs

verify_putative_cps_scope40.py
    Verify pairs which pairs of proteins with HIGH CIRPIN, LOW Progres scores are true CPs


================================================================================

SUBDIRS:

unique_CPs:
    Structures identified to have at least 1 verified CP in SCOPe40. 


================================================================================

DEDUPLICATION NOTES

All datasets have been deduplicated, taking into account the rare cases where TM-align -cp gives asymmetric scores. The 
"asymmetric" designation refers to the property of TM-align's -cp option, which 
can produce different scores depending on the order of query and target 
structures. This asymmetry required special handling during deduplication. For 
implementation details, please refer to the deduplication notebook.


