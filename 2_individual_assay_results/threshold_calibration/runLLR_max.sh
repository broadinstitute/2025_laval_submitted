#!/bin/bash

## run this script by: bash runLLR.sh

source /opt/homebrew/Caskroom/miniconda/base/etc/profile.d/conda.sh
conda activate py310

# Reference: 
# van Loggerenberg et al., 2023, The American Journal of Human Genetics 110, 1769–1786
# https://doi.org/10.1016/j.ajhg.2023.08.012

# Example usage: calculate the LLR of each variant's score per each assay for its likelihood of being pathogenic or benign
# Rscript calcLLR.R LLR_AllScores.csv LLR_RefSet.csv

# install necessary R packages in the conda env
conda install -c conda-forge r-irkernel r-devtools r-tidyverse r-ggplot2 r-dplyr r-rcpp
conda install -c conda-forge r-pbmcapply

R -e "IRkernel::installspec(user = TRUE)"
# Install remaining packages from CRAN
R -e "options(repos = c(CRAN = 'https://cran.rstudio.com/')); install.packages('argparser')"
R -e "options(repos = c(CRAN = 'https://cran.rstudio.com/')); remotes::install_github('jweile/yogitools')"

# Install packages from GitHub
R -e "remotes::install_github('jweile/yogiroc')"
R -e "remotes::install_github('jweile/yogilog')"
R -e "remotes::install_github('jweile/maveLLR')"

## acmg scaler method from https://github.com/badonyi/acmgscaler
## https://www.biorxiv.org/content/10.1101/2025.05.16.654507v1
# R -e "remotes::install_github('badonyi/acmgscaler')"

## Run for DUAL-IPA example
# Rscript calcLLR.R ../dual_ipa/3_outputs/DUALIPA_allScores.csv ../dual_ipa/3_outputs/DUALIPA_referenceSet.csv --printTransitions