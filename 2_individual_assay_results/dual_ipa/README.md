# Notes for the DUAL-IPA assay

The analysis pipeline downloads raw FACS data from Zenodo, compiles everything into a single data frame, and merges it with the allele metadata. The output file (facs_single_cell.parquet) is >200MB and so it is not tracked in the repository. To generate it for downstream analysis, run the **2_compile_facs.ipynb** notebook in the **2_analysis** directory. 