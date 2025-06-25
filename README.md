# CIRPIN: Learning Circular Permutation-Invariant Representations to Uncover Putative Protein Homologs

This repository contains the code for CIRPIN, submitted to MLCB 2025 workshop. 

## Installation

1. Installation requires Progres ( [paper link](https://academic.oup.com/bioinformaticsadvances/article/5/1/vbaf042/8107707)). The installation for Progres can be found ( [install link](https://academic.oup.com/bioinformaticsadvances/article/5/1/vbaf042/8107707](https://github.com/greener-group/progres?tab=readme-ov-file)))

2. Within `progres.py` update the `trained_model_fp` variable to `os.path.join(/home/ubuntu/CIRPIN_GNN/trained_models/CIRPIN_model, "CIRPIN_model_5k_cp_epoch301.pt")`.
3. Add databases `ark_scope40_CIRPIN_embed_1_7_25` and `ark_scope40_Progres_embed_1_7_25` to `pre_embedded_dbs`
4. Once the model has been updated, searches using the CIRPIN model can be executed following the progres command: `progres search -q query.pdb -t ark_scope40_CIRPIN_embed_1_7_25`


