# TED-Annotated PDZ Domains in AFDB-ClustR

This directory contains all TED-annotated PDZ domains identified in the AlphaFold Clustered Representatives database (AFDB-ClustR), along with tools for their analysis and embedding generation.

## Directory Contents

### Data Files

#### `AFDB_TED_PDZ_topo_dictionary_h_v_z_s_no_match_diff_025.pkl`
Contains TM-align scores for each TED domain against the four PDZ topologies id entified by CIRPIN. Domains with highly similar scores across topologies (difference < 0.025) are assigned as "no match" to avoid ambiguous classification.

#### `TED_PDZ_domains.txt`
Complete list of all TED-annotated PDZ domains in the dataset.

### Scripts

#### `make_pdz_list_from_AFDB.py`
Extracts all TED-annotated PDZ domains from the AFDB-ClustR database.

#### `download__embed_individual_afdb_PDZ_domains.py`
Downloads individual PDZ domain structures and generates their embeddings using both Progres and CIRPIN models.

#### `download_pdz.sh`
Shell script for batch downloading PDZ domain structures.

#### `PDZ_PCA_analysis.ipynb`
Notebook for PCA analysis of PDZ domains


## PDZ Topology Classification

CIRPIN identified **four distinct PDZ topologies** in this dataset:
- **h-topology** (2hga)
- **v-topology** (2vsv)
- **z-topology** (2z9i)
- **s-topology** (AFDB structure, in PDZ_topos)

Domains with TM-score differences < 0.025 between topologies are classified as "no match" due to ambiguity.

