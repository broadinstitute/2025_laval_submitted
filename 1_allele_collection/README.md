# Allele collection notes
Here we collect and process all of the allele metadata that will be used throughout the rest of the repository. Downstream analyses should read metadata from this location instead of re-downloading to ensure that the same filtering & processing logic is used throughout the project.


## Files in 1_inputs
- **VarChamp_Consolidated_allele_collection.tsv**: allele metadata from the Vidal Lab
- **download_clinvar_gnomad.ipynb**: this notebook downloads allele frequency data from gnomAD and pathogenic vs. benign labels from ClinVar
- **raw_inputs**: this is a symbolic link to store the outputs of the **download_clinvar_gnomad.ipynb** notebook. These files are too large to upload to track through Github and so should be regenerated locally by running the notebook if needed. 


## 2_analysis
The following files are used to process the raw data:
- **1_process_clinvar_gnomad.ipynb**: clean up ClinVar annotations and extract metadata for the VarChamp alleles
- **2_filter_allele_collection.ipynb**: create slim versions of the metadata for all alleles and for reference wild types only


## Files in 3_outputs
- **processed data** is a symbolic link where large, intermediate outputs from the ClinVar and gnomAD annotation processing are stored
- **varchamp_consolidated_alleles_clinvar_gnomad.tsv**: all of the ClinVar and gnomad metadata merged with the allele metadata
- **orf_info.csv**: gene-level metadata for all alleles in the VarChamp data