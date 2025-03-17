#!/bin/bash

## Downloading the necessary data files:
## ========================================================================================================
# cd ../../1_inputs/raw_inputs

## Downloading the AI-predicted scores for nsSNPs from the https://www.dbnsfp.org/download database
## It is located on G-drive in the below link: (gdown could be installed using conda/pip)
# gdown --fuzzy https://drive.google.com/file/d/1QoVE-o8YBiXyp2wQjrmIvplBEdC5QYUE/view?usp=sharing

# mkdir -p alphafold
# cd alphafold
## Downloading the Alphafold2 predictions for UP000005640_9606_HUMAN_v4 from EBI ftp
# wget https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v4.tar
# mkdir -p UP000005640_9606_HUMAN_v4/ && tar -xvf UP000005640_9606_HUMAN_v4.tar -C UP000005640_9606_HUMAN_v4/
# wget https://ftp.ebi.ac.uk/pub/databases/alphafold/accession_ids.csv

## Downloading the PROVEAN scores: link doesn't seem to work anymore, downloaded by Georges
# ftp://ftp.jcvi.org/pub/data/provean/precomputed_scores/PROVEAN_scores_ensembl66_human.tsv ##
## ========================================================================================================

## Processing the data files to get the per-residue scores
## ========================================================================================================

## When in scripts/, run:
# python convert_af_accession2parquet.py > convert_af_acc2parq.log

## Extract pLDDT scores from AF2-predicted .pdb structures
# python extract_plddt_scores.py

## Calculate the RSA scores per-residue based on AF2-predicted .pdb structures
# python calculate_rsa_scores.py
## ========================================================================================================

## Adding scores to them
python merge_scores2clinvar.py > add_plddt_rsa_scores.log