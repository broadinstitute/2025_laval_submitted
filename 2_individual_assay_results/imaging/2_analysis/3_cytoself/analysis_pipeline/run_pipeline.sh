#!/bin/bash

## Modify this path to your actual Conda installation directory:
CONDA_PATH="$HOME/software/anaconda3"

## Source the conda.sh script to initialize Conda in the current shell
## This tells bash how to find and use 'conda activate'
source "$CONDA_PATH/etc/profile.d/conda.sh"

conda activate cytoself

## Download the images files and their CellProfiler features from the Cell Painting Gallery server
python 1_download_image.py --batch_id "2024_02_06_Batch_8"

## Rescale the image intensity and zip them to zarr files for easier processing
python 2_convert_tiff_zarr.py --batch_id "2024_02_06_Batch_8"

## Rescale the image intensity and zip them to zarr files for easier processing
python 3_filter_crop_images.py --batch_id "2024_02_06_Batch_8" --exclude_plate "B8A2R2_P1"

# ulimit -n 100000  # Set to a higher number
## Re-train the cytoself model on the variant painting images
# python 4_train_model.py &> 4_train_model.log
