# Notes for the PPI assay

The PPI assay detects changes in protein-protein interactions between variant and reference alleles from the same gene. There are three key experimental components: the cell type (yeast or HEK293T), the protein that is split (GAL4 - a transcription factor - Y2H method, or nanoluciferase - a luminescent enzyme - N2H method), and the assay readout (growth - Y2H only, or bioluminescence - either Y2H or N2H). 


## Files in 1_inputs

There are three types of data: sqY2H (semi-quantitative yeast 2-hybrid; this is the main PPI profiling assay conducted in yeast), N2H_validation (nanoluc 2-Hybrid; this is a validation step conducted in orthogonal assays in both yeast and mammalian cells), and nanoluc_Y2H (nanoluc yeast 2-hybrid; another validation that switches from yeast growth to a bioluminescence readout but still uses the Y2H transcription factor).

### sqY2H input data
Within this folder:
- **cp_files**: cherry picking file relating plates and wells to allele metadata
- **mapping_384_to_96.rtf**: a text file mapping source (n = 96) well to desitination (4 reps * 96 source wells = 384 desitination) spots.
- **VarChAMP_1percent_Y2H_scores.tsv**: these are the manual scores attributed to the Y2H images using iScore software. Raw images are available on (Zenodo)[https://zenodo.org/records/15020115]. There are 6 conditions (media; 1-6) that yeast grow on.
- **media_control_score_thresholds.csv**: there are six different controls on each mega-plate (6*384) and each has media-specific thresholds. This file outlines the control-media-specific control thresholds and are used for QA/QC.
- **swimseq_confirmation**: both the activation-domain and dna-binding-domain proteins were sequenced. This folder contains binary labels for whether each protein was successfully sequence confirmed.


### N2H_validation input data
Within this folder:
- **mammalian_N2H**: raw protein-protein interactions for the reference alleles (mN2H_WTPWT) and variant alleles (mN2H_APWT), along with a file that maps the plate and wells from the raw data to allele metadata.
- **yeast_N2H**: raw protein-protein interactions for the reference alleles (WT_PWT) and variant alleles (APWT), along with a file that maps the plate and wells from the raw data to allele metadata and a file denoting the empty wells.
- **N2H sequencing confirmation**: these files contain binary scores denoting whether the expected allele was confirmed by the sequencing results.


## Using 2_analysis
Navigate into the `2_analysis` directory and create the ppi Python virtual environment with all dependencies:
```
mamba create --name ppi --file env-reqs.txt
mamba activate ppi
```
Python notebooks are numbered in the order that they should be run to reproduce all of the output files. Each notebook starts with a description of the analysis tasks that it performs, and contains extensive comments throughoupt. 

## Files in 3_outputs

### sqY2H outputs
The sqY2H analysis produces:
- **raw_individual_scores.csv**: manual scores from 0-4 for each spot. Each biological sample (ie. allele) is represented as 4 technical replicate spots.
- **raw_consensus_scores.csv**: aggregated manual scores from 0-4 for each biological sample representing the consensus of the 4 technical replicates. 
- **normalized_consensus_scores.csv**: aggregated scores that have been QA/QCed, adjusted for autoactivation, and normalized between 0 and 1.
- **edgotype_profiles**: matrix of normalized scores (interactors X alleles) for each profiled gene. 

### N2H_validation outputs
The N2H validation produces a final titration plot that shows that the overall sqY2H results are reproduced in an orthogonal assay. The same outputs are produced for both yeast N2H and mammalian N2H:
- **wtpwt_raw.csv**: normalized luminenscence ratios from the wild type pairwise test
- **wtpwt_counts.csv**: number of sqY2H wild types that can be detected in the N2H assay at different thresholds
- **wtpwt_positive_PPIs.csv**: normalized luminenscence ratios from the wild type pairwise test for all wild types that passed the threshold selected based on the _counts results
- **apwt_allele_df.csv**: normalized luminenscence ratios for both wild types and alleles for the positive_PPI wild type genes, along with perturbation z-scores and labels. 