#!/bin/bash

## Downloading the necessary data files:
## ========================================================================================================
# mkdir -p raw_inputs

## Downloading the ClinVar data from 2025-01
# mkdir -p raw_inputs/clinvar
# cd raw_inputs/clinvar
# wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/variant_summary_2025-01.txt.gz


## Downloading the AI-predicted scores for nsSNPs from the https://www.dbnsfp.org/download database
## It is located on G-drive in the below link: (gdown could be installed using conda/pip)
## The readme file is located @1_allele_collection/1_inputs/dbNSFP5.0a.readme.txt
# gdown --fuzzy https://drive.google.com/file/d/1QoVE-o8YBiXyp2wQjrmIvplBEdC5QYUE/view?usp=sharing


# mkdir -p alphafold
# cd alphafold
## Downloading the Alphafold2 predictions for UP000005640_9606_HUMAN_v4 from EBI ftp
# wget https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v4.tar
# mkdir -p UP000005640_9606_HUMAN_v4/ && tar -xvf UP000005640_9606_HUMAN_v4.tar -C UP000005640_9606_HUMAN_v4/
# wget https://ftp.ebi.ac.uk/pub/databases/alphafold/accession_ids.csv


## (OPTIONAL) downloading the PROVEAN scores: link doesn't seem to work anymore
## Information already included by dbnsfp, so no need
# wget ftp://ftp.jcvi.org/pub/data/provean/precomputed_scores/PROVEAN_scores_ensembl66_human.tsv ##
## ========================================================================================================