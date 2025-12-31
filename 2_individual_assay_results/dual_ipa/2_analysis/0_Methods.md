# DUAL-IPA Data Analysis

## Data preprocessing and quality control

Raw single-cell FACS data were compiled from all experimental plates. Single-cell GFP values of zero were replaced with the minimum non-zero GFP value to enable log transformation. For each well, we calculated the GFP:mCherry ratio for each cell and computed both mean and median ratios across all cells in the well.

Wells were subjected to the following quality control filters:
- Variant alleles: Required ≥800 cells and successful sequence confirmation (dualip_seq_final = "Y")
- Reference (wild-type) alleles: Required mean GFP signal ≥100 and GFP:mCherry ratio ≥ (mean_empty + 3×SD_empty), where empty wells (A07, A08) contained no ORF constructs
- Reference alleles also required ≥2 replicates per gene per plate for variance estimation

## Z-score calculation

For each gene, we calculated a reference abundance value as the mean GFP:mCherry ratio of all wild-type replicates on the same plate. We estimated protein abundance change for each missense variant as the log2-transformed fold-change relative to its matched reference:

$$
\text{log2FC} = \log_2\left(\frac{\text{variant GFP:mCherry ratio}}{\text{reference GFP:mCherry ratio}}\right)
$$

To determine whether this change in abundance was beyond baseline assay variability, we computed the distribution of log2 fold-changes for all wild-type-to-wild-type comparisons within each gene. Specifically, for each wild-type replicate, we calculated its log2 fold-change relative to the mean of all wild-type replicates for that gene on that plate. The mean and standard deviation of this wild-type log2FC distribution were used to normalize variant scores.

The final abundance z-score for each variant was calculated as:

$$
\text{Z-score} = \frac{\text{log2FC}_{\text{variant}} - \mu_{\text{WT:WT}}}{\sigma_{\text{WT:WT}}}
$$

where $\mu_{\text{WT:WT}}$ is the mean and $\sigma_{\text{WT:WT}}$ is the standard deviation of log2 fold-changes computed from wild-type replicate comparisons. This approach, which normalizes each variant by gene-specific wild-type variability, was found to provide superior discrimination between pathogenic and benign ClinVar variants compared to alternative normalization strategies (balanced accuracy = 0.68, AUROC = 0.70, balanced AUPRC = 0.73). 

## LLR calculation and calibration

### Evaluating strength of evidence provided by functional impact scores for variant classification

To aid in clinical variant interpretation, we determined a quantitative Bayesian evidence weight for each variant within our maps by translating the abundance scores to log likelihood ratios of pathogenicity (LLRps). To this end, we estimated the probability densities underlying the distributions of the scores of damaging and tolerated variants by kernel density estimation to obtain estimated probability density functions for pathogenic and benign variants: $\pi_+(s)$ and $\pi_-(s)$, respectively. The pathogenic:benign log likelihood ratio for a variant with a given functional impact score, $s$, was calculated as the ratio of the estimated probability density functions evaluated at $s$:

$$
\log(\Lambda) = \log\left(\frac{P(s|\text{Pathogenic})}{P(s|\text{Benign})}\right) \approx \log\left(\frac{\pi_+(s)}{\pi_-(s)}\right)
$$

The LLR calculation was performed using the maveLLR R package (https://github.com/jweile/maveLLR) with the following parameters: kernel density estimation with Epanechnikov kernel, bandwidth selection using Scott's rule of thumb ("nrd"), and outlier suppression strength of 0.1. Positive and negative reference sets were defined as ClinVar variants annotated as Pathogenic/Likely Pathogenic and Benign/Likely Benign, respectively. Confidence intervals for LLR values were estimated using 10,000 bootstrap iterations, sampling from normal distributions centered at each variant's score with standard error derived from biological replicates.

We next calibrated the relationship between LLR values and evidence strength categories within the ACMG/AMP framework, using an approach adapted from Tavtigian et al. Briefly, LLRs for descending evidence levels are modeled to decrease by a factor of 2 at each level starting from a fixed LLR, the "pathogenic very strong" (PVSt) level, which is chosen such that all ACMG/AMP evidence combination rulesets result in posterior probability of pathogenicity that satisfy the following constraints: >99% for all pathogenic (P) rules, >90% for all likely pathogenic (LP) rules, <10% for all likely benign (LB) rules, and <1% for all benign (B) rules. We followed Tavtigian et al. in modeling a global prior probability of 0.1, as this enabled a set of LLR values for evidence at different evidence strength levels, which was internally consistent across the set of "path rules" established within the ACMG guidelines.

### Determination of pathogenicity thresholds

Based on the LLR calibration and examination of the score distributions for ClinVar pathogenic and benign variants, we established the following thresholds for protein abundance change classification:

- **Pathogenic-supporting abundance change**: Abundance z-score ≤ -2.522 or ≥ 1.994
- **Benign-supporting abundance**: -0.566 ≤ abundance z-score ≤ 0.932
- **Uncertain**: All other scores

These thresholds correspond to transitions in ACMG/AMP evidence categories and were determined by evaluating the LLR function across the score range. Variants with z-scores ≤ -2.522 or ≥ 1.994 demonstrate protein abundance changes consistent with loss-of-function or gain-of-function mechanisms that align with pathogenic variants in ClinVar. Conversely, variants with z-scores within the benign range (-0.566 to 0.932) show abundance levels similar to wild-type, consistent with benign variants.

The direction of abundance change (increase vs. decrease) was also annotated, with scores < 0 indicating decreased abundance and scores > 0 indicating increased abundance relative to wild-type.

---

## Manuscript Method

### DUAL-IPA protein abundance measurement and analysis

Raw single-cell FACS data measuring GFP and mCherry fluorescence were compiled from all experimental plates. Single-cell GFP values of zero were replaced with the minimum non-zero GFP value to enable subsequent log transformation. For each well, the GFP:mCherry ratio was calculated for individual cells, and both mean and median ratios were computed across all cells in the well. Wells were retained for analysis only if they contained at least 800 cells. Variant alleles additionally required successful sequence confirmation, while reference wild-type alleles were filtered to retain only those with mean GFP signal ≥100 and GFP:mCherry ratio exceeding the mean of empty wells (containing no ORF constructs) by at least three standard deviations. Reference alleles also required at least two replicates per gene per plate to enable variance estimation.

For each gene, we calculated a reference abundance value as the mean GFP:mCherry ratio of all wild-type replicates on the same plate. Protein abundance change for each missense variant was estimated as the log2-transformed fold-change of the variant's mean GFP:mCherry ratio relative to its matched reference. To assess whether observed changes exceeded baseline assay variability, we computed the distribution of log2 fold-changes for wild-type-to-wild-type comparisons. Specifically, for each wild-type replicate, we calculated its log2 fold-change relative to the mean of all wild-type replicates for that gene on that plate. The final abundance z-score for each variant was calculated by subtracting the mean and dividing by the standard deviation of this wild-type log2 fold-change distribution. This gene-specific normalization approach provided superior discrimination between pathogenic and benign ClinVar variants (balanced accuracy = 0.68, AUROC = 0.70, balanced AUPRC = 0.73) compared to alternative normalization strategies.

To translate abundance z-scores into quantitative Bayesian evidence weights for clinical variant interpretation, we calculated log likelihood ratios of pathogenicity (LLRp) using kernel density estimation. Probability density functions for pathogenic and benign variants, π₊(s) and π₋(s), were estimated from ClinVar variants annotated as Pathogenic/Likely Pathogenic and Benign/Likely Benign, respectively. The LLRp for a variant with score s was calculated as log[π₊(s)/π₋(s)] using the maveLLR R package with an Epanechnikov kernel, Scott's rule of thumb bandwidth selection, and outlier suppression strength of 0.1. Confidence intervals were estimated using 10,000 bootstrap iterations. LLR values were calibrated to ACMG/AMP evidence strength categories following Tavtigian et al., modeling evidence levels to decrease by a factor of 2 at each descending level from the pathogenic very strong threshold, with a global prior probability of 0.1.

Based on LLR calibration and examination of ClinVar variant score distributions, we established pathogenicity thresholds corresponding to transitions in ACMG/AMP evidence categories. Variants with abundance z-scores ≤-2.522 or ≥1.994 were classified as having pathogenic-supporting abundance changes, indicating protein levels consistent with loss-of-function or gain-of-function mechanisms. Variants with z-scores between -0.566 and 0.932 were classified as having benign-supporting abundance, indicating protein levels similar to wild-type. All other variants were classified as uncertain. The direction of abundance change was annotated as decreased (z-score < 0) or increased (z-score > 0) relative to wild-type.