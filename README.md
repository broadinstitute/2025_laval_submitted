# VarChAMP: functional profiling of human genetic variants

Some of the analysis results in this repository have been published in [this bioRxiv preprint](https://www.biorxiv.org/content/10.64898/2026.02.14.705848v2). Please cite this article when using the associated analyses or results.

This repository contains the data analysis code, processed results, and manuscript figures and tables for the first VarChAMP publication. The analyses compare human protein variants with their reference alleles across complementary assays of **protein abundance, subcellular localization, and protein–protein interactions**. They integrate these experimental measurements, compare them with computational variant-effect predictions, and explore selected biological mechanisms.

## Analyses

| Directory | Contents |
| --- | --- |
| [1_allele_collection](1_allele_collection/) | Allele curation and shared metadata, including ClinVar annotations and gnomAD allele frequencies. |
| [2_individual_assay_results](2_individual_assay_results/) | Assay-specific processing, quality control, and statistical threshold calibration for DUAL-IPA abundance measurements, Variant Painting imaging, and protein–protein interaction assays. |
| [3_integrated_assay_results](3_integrated_assay_results/) | Cross-assay comparisons, integrated variant profiles, and analyses of clinical annotations and inheritance. |
| [4_compare_ai_scores](4_compare_ai_scores/) | Comparison of computational variant-effect predictions with experimental measurements and ClinVar annotations. |
| [5_mechanistic_case_studies](5_mechanistic_case_studies/) | Focused analyses of selected genes, including MVK and MLH1, and follow-up functional assays. |
| [6_ppi_analysis](6_ppi_analysis/) | Contextualization and prioritization of protein–protein interactions using published single-cell and single-nucleus gene-expression datasets from CELLxGENE. |

The individual assay workflows are described in the [DUAL-IPA](2_individual_assay_results/dual_ipa/README.md), [imaging](2_individual_assay_results/imaging/README.md), and [PPI](2_individual_assay_results/ppi/README.md) READMEs.

## Repository structure

Analysis sections generally follow the same input-to-output organization:

- **`1_inputs/`**: raw assay data, metadata, and external annotations.
- **`2_analysis/`**: notebooks and scripts for processing, analysis, and visualization.
- **`3_outputs/`**: processed data and intermediate results used by downstream analyses.
- **`4_manuscript_figures/`**: high-resolution main-text and supplementary figures.
- **`5_manuscript_tables/`**: formatted main-text and supplementary tables.

## Exploring and reproducing the analyses

Start with the analysis section of interest and its README or notebook descriptions. Shared allele metadata comes from `1_allele_collection/`; individual assay outputs feed into the integrated analyses and downstream comparisons.

Workflows use Python, R, and, for the imaging pipeline, Snakemake. Follow the environment definitions and execution instructions provided within each analysis directory; there is no single environment or command that runs every analysis. Where notebooks and scripts are numbered, follow their documented execution order and check their input dependencies.

## License

See [LICENSE](LICENSE) for the repository's license terms.
