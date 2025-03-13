# Notes for the PPI assay

The PPI assay detects changes in protein-protein interactions between variant and reference alleles from the same gene. There are three key experimental components: the cell type (yeast or HEK293T), the protein that is split (GAL4 - a transcription factor - Y2H method, or nanoluciferase - a luminescent enzyme - N2H method), and the assay readout (growth - Y2H only, or bioluminescence - either Y2H or N2H). 


## Files in 1_inputs

There are three types of data: sqY2H (semi-quantitative yeast 2-hybrid; this is the main PPI profiling assay conducted in yeast), N2H_validation (nanoluc 2-Hybrid; this is a validation step conducted in orthogonal assays in both yeast and mammalian cells), and nanoluc_Y2H (nanoluc yeast 2-hybrid; another validation that switches from yeast growth to a bioluminescence readout but still uses the Y2H transcription factor).

### sqY2H data
Within this folder:
- **cp_files**: cherry picking file relating plates and wells to allele metadata
- **mapping_384_to_96.rtf**: a text file mapping source (n = 96) well to desitination (4 reps * 96 source wells = 384 desitination) spots.
- **VarChAMP_1percent_Y2H_scores.tsv**: these are the manual scores attributed to the Y2H images using iScore software. Raw images are available on (Zenodo)[https://zenodo.org/records/15020115]. There are 6 conditions (media; 1-6) that yeast grow on. Three different people score, so there are 18 different numbers (r01-r18). There are also 4 batches of experimental data VUSAPWT1B1, VUSAPWT1B2, VUSAPWT2B1, and VUSAPWT6B1. Flo scored the first two (r07-r12), Georges scored the last two (r13-r18). Should add another column (Scorer) and convert all numerical scores to 1-6 to reflect comparable media conditions. 


### N2H_validation data
Within this folder:
- **mammalian_N2H**: raw protein-protein interactions for both the reference alleles (mN2H_WTPWT) and the variant alleles (mN2H_APWT). CHECK TO SEE IF INDIVIDUAL FILES CAN BE REMOVED IF ALSO COMPILED.
- **N2H sequencing confirmation**: 
- **yeast_N2H**: 


### nanoluc_Y2H data