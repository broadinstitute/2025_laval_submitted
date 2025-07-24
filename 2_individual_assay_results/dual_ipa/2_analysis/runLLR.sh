#!/bin/bash

## run this script by: bash runLLR.sh

## install necessary R packages in the conda env
# conda install -c conda-forge r-irkernel r-devtools r-tidyverse r-ggplot2 r-dplyr r-rcpp
# conda install -c conda-forge r-pbmcapply

# R -e "IRkernel::installspec(user = TRUE)"
# # Install remaining packages from CRAN
# R -e "options(repos = c(CRAN = 'https://cran.rstudio.com/')); install.packages('argparser')"
# R -e "options(repos = c(CRAN = 'https://cran.rstudio.com/')); remotes::install_github('jweile/yogitools')"

# # Install packages from GitHub
# R -e "remotes::install_github('jweile/yogiroc')"
# R -e "remotes::install_github('jweile/yogilog')"
# R -e "remotes::install_github('jweile/maveLLR')"

source ~/software/anaconda3/etc/profile.d/conda.sh

conda activate py310

Rscript ../../threshold_calibration/calcLLR.R ../3_outputs/DUALIPA_allScores.csv ../3_outputs/DUALIPA_referenceSet.csv --outfile "../3_outputs/DUALIPA_LLR_score" --logfile "../3_outputs/DUALIPA_LLR.log" --printTransitions