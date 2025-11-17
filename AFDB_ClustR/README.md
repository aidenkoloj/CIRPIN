# CIRPIN AFDB Cluster Representatives Pipeline

This directory contains all the scripts used to run the schematic in Figure 2 of the CIRPIN paper using the AlphaFold Database (AFDB) clustered representatives dataset.

## Overview

The pipeline identifies putative circular permutation pairs by comparing Progres and CIRPIN embedding scores across ~3 million AFDB cluster representatives. Many operations are performed in chunks to manage memory constraints when working with large-scale protein datasets.

## Key Scripts

### (1) `generate_afdb_clustr_reps.py`
### (2) `fc_embed*`
### (3) `test_all_v_all_afdb.ipynb`
### (4) `get_putative_pairs_lists*`
### (5) `verify_putative_pairs*`
### (6) `combine_verify_pairs_output.py`
### (7) `get_CATH_for_AFDB_*`

**Purpose**

## (1)
### - Used to generate a .tsv of all the TED domains contained within AFDB-ClustR structures

## (2)
### - Used to embed all the TED domains using Foldcomp (fc)

## (3)
### - Small scale tests of all v all
### - Combine all chunked embeddings into single .pt --> combined_embs_3M*.pt

## (4)
### - Run all vs all comparisons on AFDB-ClustR

## (5)
### - Verify the putative pairs 

## (6)
### - Combine the outputs from (5) 

## (7)
### - Scripts to obtain CATH labels for outputs