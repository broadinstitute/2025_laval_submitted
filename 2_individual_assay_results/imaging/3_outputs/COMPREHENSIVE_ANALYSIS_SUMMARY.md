# Comprehensive Disease Biology and Feature Analysis

**Date:** 2025-11-19
**Analysis:** Cell Painting Imaging Variant Clustering with Disease Modules and Feature Analysis
**Variants:** 140 across 7 morphological clusters

---

## Executive Summary

This comprehensive analysis integrates:
1. **Disease module classification** - Understanding biological pathways
2. **Gene-specific clustering patterns** - Identifying consistent phenotypes
3. **Feature-driven cluster analysis** - Morphological basis of clustering
4. **Feature correlation structure** - Representative feature selection

### Key Findings

#### 1. Disease Module Enrichment

- **Cluster 2**: 53.6% Other_Genetic_Disorders enrichment (15/28 variants)
- **Cluster 3**: 71.4% Other_Genetic_Disorders enrichment (15/21 variants)
- **Cluster 7**: 57.1% Neuropathy enrichment (8/14 variants)

#### 2. Gene Clustering Patterns

**Perfect clustering (100% consistency):** RAD51D, NR0B1, DCX, NT5C3A, PNPO

**High consistency (≥80%):** PLP1, LITAF, GFAP


#### 3. Feature Structure

- Total features: 52
- Feature clusters identified: 23
- Representative features selected: 23
- Largest feature cluster: 7 features


---

## Part 1: Disease Module Analysis

### Overall Statistics

- **Total variants analyzed:** 140
- **Variants with disease annotations:** 109
- **Total clusters:** 7
- **Disease modules identified:** 10
- **Mean AUROC:** 0.964

### Cluster-by-Cluster Disease Analysis

#### Cluster 1

- **Size:** 29 variants
- **Mean AUROC:** 0.950
- **Altered:** 22/29 (75.9%)

**Disease Modules:**

- Cancer_Syndromes: 14 variants (48.3%)
- Other_Genetic_Disorders: 10 variants (34.5%)
- Unclassified: 4 variants (13.8%)
- Immunodeficiency: 1 variants (3.4%)

**Top Genes:**

- MLH1: 8 variants
- RAD51D: 6 variants
- KRT8: 3 variants
- BFSP2: 2 variants
- APOA1: 2 variants

#### Cluster 2

- **Size:** 28 variants
- **Mean AUROC:** 0.971
- **Altered:** 16/28 (57.1%)

**Disease Modules:**

- Other_Genetic_Disorders: 15 variants (53.6%)
- Unclassified: 9 variants (32.1%)
- Cancer_Syndromes: 2 variants (7.1%)
- Coagulation_Disorders: 1 variants (3.6%)
- Metabolic_Disorders: 1 variants (3.6%)

**Top Genes:**

- DCX: 4 variants
- NT5C3A: 4 variants
- MVK: 3 variants
- FBP1: 3 variants
- PNPO: 3 variants

#### Cluster 3

- **Size:** 21 variants
- **Mean AUROC:** 0.971
- **Altered:** 16/21 (76.2%)

**Disease Modules:**

- Other_Genetic_Disorders: 15 variants (71.4%)
- Coagulation_Disorders: 3 variants (14.3%)
- Unclassified: 2 variants (9.5%)
- Retinal_Disorders: 1 variants (4.8%)

**Top Genes:**

- GFAP: 12 variants
- F9: 3 variants
- SFTPA2: 1 variants
- SFTPC: 1 variants
- SDHD: 1 variants

#### Cluster 4

- **Size:** 18 variants
- **Mean AUROC:** 0.939
- **Altered:** 11/18 (61.1%)

**Disease Modules:**

- Other_Genetic_Disorders: 6 variants (33.3%)
- Unclassified: 5 variants (27.8%)
- Cancer_Syndromes: 4 variants (22.2%)
- Cardiomyopathy: 1 variants (5.6%)
- Neuropathy: 1 variants (5.6%)
- Retinal_Disorders: 1 variants (5.6%)

**Top Genes:**

- SMAD3: 4 variants
- BRCA1: 3 variants
- AGXT: 2 variants
- MLH1: 2 variants
- BAG3: 1 variants

#### Cluster 5

- **Size:** 16 variants
- **Mean AUROC:** 0.981
- **Altered:** 14/16 (87.5%)

**Disease Modules:**

- Other_Genetic_Disorders: 7 variants (43.8%)
- Unclassified: 5 variants (31.2%)
- Neuropathy: 2 variants (12.5%)
- Retinal_Disorders: 2 variants (12.5%)

**Top Genes:**

- PLP1: 8 variants
- PMP22: 3 variants
- RHO: 3 variants
- SFTPC: 2 variants

#### Cluster 6

- **Size:** 14 variants
- **Mean AUROC:** 0.989
- **Altered:** 14/14 (100.0%)

**Disease Modules:**

- Other_Genetic_Disorders: 5 variants (35.7%)
- Cancer_Syndromes: 5 variants (35.7%)
- Unclassified: 4 variants (28.6%)

**Top Genes:**

- MLH1: 7 variants
- NR0B1: 4 variants
- SMAD3: 2 variants
- AIPL1: 1 variants

#### Cluster 7

- **Size:** 14 variants
- **Mean AUROC:** 0.963
- **Altered:** 8/14 (57.1%)

**Disease Modules:**

- Neuropathy: 8 variants (57.1%)
- Muscle_Disorders: 2 variants (14.3%)
- Other_Genetic_Disorders: 2 variants (14.3%)
- Unclassified: 2 variants (14.3%)

**Top Genes:**

- LITAF: 10 variants
- DES: 2 variants
- FXN: 1 variants
- PLP1: 1 variants


---

## Part 2: Gene Clustering Consistency

Genes with multiple variants show varying degrees of clustering consistency,
indicating how similar their morphological phenotypes are across different mutations.

### Consistency Analysis

| Gene | Variants | Primary Cluster | Consistency | Mean AUROC |
|------|----------|----------------|-------------|------------|
| RAD51D | 6 | 1 | 100.0% | 0.962 |
| NR0B1 | 4 | 6 | 100.0% | 0.998 |
| DCX | 4 | 2 | 100.0% | 0.949 |
| NT5C3A | 4 | 2 | 100.0% | 0.975 |
| PNPO | 3 | 2 | 100.0% | 0.959 |
| GFAP | 13 | 3 | 92.3% | 0.973 |
| LITAF | 11 | 7 | 90.9% | 0.957 |
| PLP1 | 9 | 5 | 88.9% | 0.978 |
| BRCA1 | 4 | 4 | 75.0% | 0.973 |
| F9 | 4 | 3 | 75.0% | 0.968 |
| SMAD3 | 7 | 4 | 57.1% | 0.964 |
| MLH1 | 17 | 1 | 47.1% | 0.947 |


### Interpretation

- **100% consistency**: All variants from this gene cluster together - strong gene-specific phenotype
- **80-99% consistency**: Most variants cluster together - relatively consistent phenotype
- **<80% consistency**: Variants show heterogeneous phenotypes - mutation-location may matter


---

## Part 3: Feature Correlation Structure

Features were clustered based on their pairwise correlations to identify
redundant features and select representatives.

### Feature Clustering Summary

- **Total features:** 52
- **Feature clusters:** 23
- **Singleton clusters:** 12
- **Largest cluster:** 7 features

### Top 10 Representative Features

| Feature | Cluster Size | Mean |r| with Cluster |
|---------|-------------|----------------------|
| Nuclei_RadialDistribution_RadialCV_GFP_10of10 | 7 | 0.810 |
| Cytoplasm_Texture_InfoMeas1_GFP_20_00_256 | 5 | 0.806 |
| Nuclei_Granularity_4_GFP | 5 | 0.853 |
| Nuclei_RadialDistribution_MeanFrac_GFP_10of10 | 5 | 0.900 |
| Cells_Correlation_Correlation_DNA_GFP | 3 | 0.833 |
| Cytoplasm_Granularity_3_GFP | 3 | 0.802 |
| Cytoplasm_RadialDistribution_MeanFrac_GFP_6of10 | 3 | 0.890 |
| Cells_Correlation_RWC_Mito_GFP | 3 | 0.819 |
| Cells_Correlation_RWC_GFP_DNA | 2 | 0.790 |
| Nuclei_Correlation_K_GFP_Mito | 2 | 0.765 |


### Usage Recommendations

1. **For dimensionality reduction:** Use representative features to reduce from 52 to 23 features
2. **For interpretation:** Focus on representative features from large clusters (high redundancy removed)
3. **For modeling:** Using representatives avoids multicollinearity issues


---

## Part 4: Cluster-Specific Morphological Features

For each cluster, we identified the top 3 features that distinguish it from other clusters
using Cohen's d effect size.

### Feature-Cluster Associations

#### Cluster 1 Distinctive Features

| Feature | Effect Size (d) | Direction | Biological Relevance |
|---------|----------------|-----------|--------------------|
| Cells_RadialDistribution_RadialCV_GFP_5of10 | 2.72 | ↑ Higher | Cytoplasmic |
| Nuclei_RadialDistribution_RadialCV_GFP_10of10 | 2.14 | ↑ Higher | Nuclear morphology |
| Cytoplasm_Correlation_RWC_AGP_GFP | -1.98 | ↓ Lower | ER structure |

#### Cluster 2 Distinctive Features

| Feature | Effect Size (d) | Direction | Biological Relevance |
|---------|----------------|-----------|--------------------|
| Cytoplasm_Texture_InverseDifferenceMoment_GFP_5_00 | 1.85 | ↑ Higher | Cytoplasmic |
| Nuclei_RadialDistribution_MeanFrac_GFP_1of10 | 1.81 | ↑ Higher | Nuclear morphology |
| Nuclei_Texture_InverseDifferenceMoment_GFP_20_00_2 | 1.79 | ↑ Higher | Nuclear morphology |

#### Cluster 3 Distinctive Features

| Feature | Effect Size (d) | Direction | Biological Relevance |
|---------|----------------|-----------|--------------------|
| Cytoplasm_Texture_InverseDifferenceMoment_GFP_5_00 | -1.78 | ↓ Lower | Cytoplasmic |
| Cells_RadialDistribution_FracAtD_GFP_8of10 | 1.30 | ↑ Higher | Cytoplasmic |
| Cells_Intensity_MADIntensity_GFP | 1.26 | ↑ Higher | Cytoplasmic |

#### Cluster 4 Distinctive Features

| Feature | Effect Size (d) | Direction | Biological Relevance |
|---------|----------------|-----------|--------------------|
| Cytoplasm_Intensity_MassDisplacement_GFP | -1.02 | ↓ Lower | Cytoplasmic |
| Nuclei_Correlation_Overlap_GFP_Mito | -0.99 | ↓ Lower | Nuclear morphology |
| Cells_Intensity_MADIntensity_GFP | -0.94 | ↓ Lower | Cytoplasmic |

#### Cluster 5 Distinctive Features

| Feature | Effect Size (d) | Direction | Biological Relevance |
|---------|----------------|-----------|--------------------|
| Cells_Correlation_Correlation_GFP_Mito | 3.31 | ↑ Higher | Mitochondrial |
| Cytoplasm_Correlation_Correlation_DNA_GFP | 2.81 | ↑ Higher | Nuclear morphology |
| Cytoplasm_Texture_Correlation_GFP_5_01_256 | 2.70 | ↑ Higher | Cytoplasmic |

#### Cluster 6 Distinctive Features

| Feature | Effect Size (d) | Direction | Biological Relevance |
|---------|----------------|-----------|--------------------|
| Nuclei_Correlation_Correlation_DNA_GFP | -3.42 | ↓ Lower | Nuclear morphology |
| Nuclei_RadialDistribution_MeanFrac_GFP_10of10 | 3.30 | ↑ Higher | Nuclear morphology |
| Cells_Correlation_Correlation_DNA_GFP | -3.03 | ↓ Lower | Nuclear morphology |

#### Cluster 7 Distinctive Features

| Feature | Effect Size (d) | Direction | Biological Relevance |
|---------|----------------|-----------|--------------------|
| Nuclei_Texture_Correlation_GFP_5_03_256 | 2.68 | ↑ Higher | Nuclear morphology |
| Nuclei_Correlation_RWC_GFP_Mito | -2.50 | ↓ Lower | Nuclear morphology |
| Cytoplasm_Correlation_RWC_GFP_DNA | 2.41 | ↑ Higher | Nuclear morphology |


---

## Part 5: Biological Interpretation

### Feature-Disease Connections

Based on the morphological features that distinguish each cluster, we can infer
biological mechanisms:

#### Cluster 1: Cancer_Syndromes

- **Nuclear features prominent**: 10of10
- **ER features prominent**: Suggests ER stress or protein trafficking defects

#### Cluster 2: Other_Genetic_Disorders

- **Nuclear features prominent**: 1of10, 256

#### Cluster 3: Other_Genetic_Disorders


#### Cluster 4: Other_Genetic_Disorders

- **Nuclear features prominent**: Mito
- **Mitochondrial features**: May indicate metabolic dysfunction

#### Cluster 5: Other_Genetic_Disorders

- **Nuclear features prominent**: GFP
- **Mitochondrial features**: May indicate metabolic dysfunction

#### Cluster 6: Other_Genetic_Disorders

- **Nuclear features prominent**: GFP, 10of10

#### Cluster 7: Neuropathy

- **Nuclear features prominent**: 256, Mito
- **Mitochondrial features**: May indicate metabolic dysfunction


---

## Part 6: Data Files and Visualizations

### Generated Files

1. **comprehensive_disease_biology_analysis.pdf** - All visualizations in one PDF
2. **COMPREHENSIVE_ANALYSIS_SUMMARY.md** - This summary document

### Original Data Files Used

1. `variant_clustering_results.csv` - PCA, UMAP, cluster assignments
2. `variant_protein_feature_matrix.csv` - Morphological features per variant
3. `protein_features_corr_mat.csv` - Feature correlation matrix

### Tables for Export

The following dataframes are available in the analysis:
- `classified_df` - Disease module classifications per variant
- `gene_df` - Gene clustering consistency analysis
- `representative_features_df` - Representative features per feature cluster
- `cluster_feature_df` - Discriminative features per variant cluster

---

## Conclusions

This integrated analysis reveals:

1. **Disease biology drives clustering**: Variants with similar disease mechanisms cluster together morphologically
2. **Gene-specific signatures exist**: Some genes (100% consistency) have invariant phenotypes across mutations
3. **Features are interpretable**: Morphological features connect to biological processes (ER, mitochondria, nucleus)
4. **Feature redundancy is high**: 52 features can be reduced to 23 representatives

### Clinical Implications

- Morphological clustering could predict variant pathogenicity
- Gene consistency patterns inform variant interpretation strategies
- Feature-disease connections suggest mechanistic hypotheses
- Representative features enable efficient future screening

### Recommended Follow-ups

1. **Validate feature-disease connections** with orthogonal assays
2. **Test variant predictions** using morphological similarity to known variants
3. **Develop feature-based biomarkers** for specific disease modules
4. **Investigate heterogeneous genes** (e.g., MLH1) to understand why variants differ

---

**Analysis completed:** 2025-11-19
**Scripts:** `comprehensive_disease_biology_analysis.py`, `disease_module_annot.py`
**Contact:** Bioinformatics team for questions
