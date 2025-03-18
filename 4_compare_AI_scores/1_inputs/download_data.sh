#!/bin/bash

## Downloading the necessary data files:
## ========================================================================================================
# cd raw_inputs


## Downloading the AI-predicted scores for nsSNPs from the https://www.dbnsfp.org/download database
## It is located on G-drive in the below link: (gdown could be installed using conda/pip)
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
