# Allele collection notes
Here we collect and process all of the allele metadata that will be used throughout the rest of the repository. Downstream analyses should read metadata from this location instead of re-downloading to ensure that the same filtering & processing logic is used throughout the project.


## Files in 1_inputs
- **VarChamp_Consolidated_allele_collection.tsv**: allele metadata from the Vidal Lab
- **download_clinvar_gnomad.ipynb**: this notebook downloads allele frequency data from gnomAD and pathogenic vs. benign labels from ClinVar
- **raw_inputs**: this is a symbolic link to store the outputs of the **download_clinvar_gnomad.ipynb** notebook. These files are too large to upload to track through Github and so should be regenerated locally by running the notebook if needed. 


## 2_analysis
The following files are used to process the raw data:
- **1_process_clinvar_gnomad.ipynb**: 
- **2_filter_allele_collection.ipynb**: this script 


## Files in 3_outputs
- processed data is a symbolic link (intermediate files)
- clinvar_gnomad is all columns merged to Tong's data