# CPDB and CIRPIN Data Documentation

## Overview

This repository contains representative structures from network clustering of circular permutation (CP) pairs, comparing data from CIRPIN and CPDB databases.

> **Note:** Access to CPDB is no longer maintained. Data from CPDB was obtained from [https://github.com/YueHuLab/plmCP/](https://github.com/YueHuLab/plmCP/)

---

## Representative Structure Folders

These folders contain representative structures from network clustering of CP pairs:

| Folder | Count | Description |
|--------|-------|-------------|
| `CIRPIN_afdb_rep_network_structures` | 455 | CIRPIN structures from AlphaFold DB |
| `CIRPIN_scope_rep_network_structures` | 72 | CIRPIN structures from SCOPe |
| `CIRPIN_all_network_structures` | 527 | Combined AFDB and SCOPe structures |
| `CPDB_pdb_rep_network_structures` | 129 | CPDB structures from PDB |

---

## Network Cluster Lists

Location: `network_cluster_lists/`

### CIRPIN Data Files

- **`CIRPIN_afdb_represntative_structures_single_cluster.pkl`**  
  Single representative structure from each network cluster (455 total)

- **`CIRPIN_afdb_representative_structures.pkl`**  
  Complete list containing all structures across all network clusters (~20,000 structures)

- **`CIRPIN_SCOPe_representative_structures_single_cluster.pkl`**  
  Single representative structure from each network cluster (72 total)

### CPDB Data Files

- **`CPDB_clusters.pkl`**  
  Complete list containing all structures across all network clusters (1,667 structures)

- **`CPDB_representative_structures_single_cluster.pkl`**  
  Single representative structure from each network cluster (129 total)

---

## CPDB/CIRPIN Comparison

- **`overlap_CIRPIN_scope_CPDB.*`**  
  Files documenting similar structures between CIRPIN CP hits and CPDB hits

---

## CPDB Original Pairs

- **`tm_scores_cpdb.txt`**  
  TM-scores between pairs of circular permutations from the CPDB  
  Original data source: [https://github.com/YueHuLab/plmCP/](https://github.com/YueHuLab/plmCP/)

---

## Analysis

- **`CPDB_analyis.ipynb`**  
  Jupyter notebook containing analysis of CPDB data

- **`compare_dbs_optimize.py`**  
  Calculate TM scores between representative structures from CIRPIN clusters/CPDB clusters to identify novel CPs 

