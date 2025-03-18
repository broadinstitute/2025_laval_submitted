# Notes for the Imaging assay

The imaging assay detects differences in protein localization and abundance using Variant Painting, an adaptation of the Cell Painting assay where two of the usual channels are dropped (here, ER and RNA) and we add a "GFP" channel that allows us to visualize the subcellular localization of the protein for which we are profiling variant impacts. 


## Files in 1_inputs
All raw imaging data and metadata is publicly available on the Cell Painting Gallery, and so we include a shell script (**download_aws_cp_data.sh**) to download all input data needed for the analysis pipeline. 


## Using 2_analysis
The analysis is split into three parts:
-**1_snakemake_pipeline**: pipeline for detecting protein mislocalization between reference and variant alleles based on single-cell CellProfiler features
-**2_abundance**: flagging changes in abundance between reference and variant alleles based on single-cell CellProfiler features
-**3_cytoself**: using a VQ-VAE to learn a cell representation that has a higher subcellular localization signal


## Files in 3_outputs
Each of the main imaging analysis modules produces results:
-**sm_xgb_classifier_metrics**: variant-reference AUROC scores from binary XGBoost classifiers for each batch of data
-**abundance**: t-test statistics and p-values for changes in GFP intensity for each variant-reference pair