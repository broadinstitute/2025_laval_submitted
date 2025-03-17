#!/bin/bash

## Download the necessary data and files for imaging snakemake pipeline 

BASEPATH="s3://cellpainting-gallery/cpg0020-varchamp/broad/workspace"
DATAPATH="../2_analysis/1_snakemake_pipeline/2025_varchamp_snakemake/snakemake_pipeline/inputs"
## To download all of the sources in the JUMP dataset
# batches=`aws s3 ls --no-sign-request "$BASEPATH/backend/" | awk '{print substr($2, 1, length($2)-1)}'`

batches="2024_01_23_Batch_7 2024_02_06_Batch_8"

# mkdir -p inputs/well_profiles
mkdir -p $DATAPATH/single_cell_profiles
mkdir -p $DATAPATH/metadata/platemaps

for batch_id in $batches;
do
    aws s3 sync --no-sign-request --exclude "*.csv" "$BASEPATH/backend/$batch_id" $DATAPATH/single_cell_profiles/$batch_id
    aws s3 sync --no-sign-request "$BASEPATH/metadata/platemaps/$batch_id" $DATAPATH/metadata/platemaps/$batch_id

    ## (Optional) If you want to download the data channel by channel, you can also try: 
    # aws s3 cp --no-sign-request \
    #     $BASEPATH/analysis/$batch_id/2025_01_27_B13A7A8P1_T1/analysis/ \
    #     inputs/analysis/2025_01_27_B13A7A8P1_T1/ \
    #     --exclude "*" \
    #     --include "**/Cells.csv" \
    #     --include "**/Cytoplasm.csv" \
    #     --include "**/Nuclei.csv" \
    #     --include "**/Image.csv" \
    #     --recursive
done