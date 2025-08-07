#!/bin/bash

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
# python merge_scores2clinvar.py \
#     --input_file "../../3_outputs/processed_data/dbnsfp/dbNSFP5.0a_variant.clin_var_re-annot_dedup_variants_10k.parquet" \
#     --output_file "../../3_outputs/processed_data/dbnsfp/dbNSFP5.0a_variant.clin_var_re-annot_dedup_variants_10k_scores.parquet" \
#     > add_plddt_rsa_scores_100k.log

python merge_scores2clinvar.py \
    --input_file "../../3_outputs/processed_data/dbnsfp/dbNSFP5.0a_variant.clin_var_re-annot_pdb_variants.parquet" \
    --output_file "../../3_outputs/processed_data/dbnsfp/dbNSFP5.0a_variant.clin_var_re-annot_pdb_variants_plddt_rsa.parquet" \
    > pdb_prioritized_variants_add_plddt_rsa_scores.log

## Append all protein features from G2P portal to allele_collection variants
source /home/shenrunx/software/anaconda3/etc/profile.d/conda.sh && conda activate ai4protein && nohup python process_large_variant_collection.py