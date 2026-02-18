# Figure 2 Panel Generation Plan for VarChAMP

## Overview

Generate Figure 2 panels for the VarChAMP manuscript. Read the master dataset and create panels based on **available features only**.

---

## Data Source

```python
import pandas as pd
import sys
sys.path.insert(0, '../../1_allele_collection/2_analysis')
from plot_utils import palette_dict, clinvar_order_dict

# Load data
DATA_PATH = '../3_outputs/0_all_gene_variants_assayed_summary.tsv'
df = pd.read_csv(DATA_PATH, sep='\t', low_memory=False)
```

---

## Color Palettes

### ClinVar Palette (IMPORT FROM plot_utils.py)

Use the existing palette from `1_allele_collection/2_analysis/plot_utils.py`:

```python
# Import existing palettes - DO NOT REDEFINE
from plot_utils import palette_dict, clinvar_order_dict

clinvar_palette = palette_dict['clinvar_clnsig_clean']
# Contains: Pathogenic, PLP, Benign, BLB, Conflicting, VUS, Others
```

### Assay-Specific Colors (Nature-Friendly)

```python
assay_palette = {
    'Dual-IPA': '#E69F00',       # orange
    'Imaging': '#0072B2',        # blue
    'PPI': '#009E73',            # bluish green
}
```

### Structural/Categorical Palettes

```python
# RSA (Relative Solvent Accessibility) - Binary at threshold 0.25
rsa_palette = {
    'Buried': '#D55E00',         # vermillion (RSA < 0.25)
    'Exposed': '#56B4E9',        # sky blue (RSA >= 0.25)
}

# pLDDT (AlphaFold confidence) - Binary at threshold 70
plddt_palette = {
    'Disordered': '#D55E00',     # vermillion (pLDDT < 70)
    'Structured': '#56B4E9',     # sky blue (pLDDT >= 70)
}

# Secondary structure (3-state DSSP)
secondary_structure_palette = {
    'H (helix)': '#CC79A7',      # reddish purple
    'B (undefined)': '#E69F00',  # orange (strand/undefined)
    'C (loop/coil)': '#56B4E9',  # sky blue
}

# Inheritance mode
inheritance_palette = {
    'autosomal dominant': '#CC79A7',   # reddish purple
    'autosomal recessive': '#0072B2',  # blue
    'x-linked': '#D55E00',             # vermillion
    'multiple': '#708090',             # slate
}
```

---

## Available Columns Summary

Based on `0_all_gene_variants_assayed_summary.tsv`:

| Category | Available Columns |
|----------|-------------------|
| **Core scores** | `dualipa_score`, `imaging_score`, `ppi_score` |
| **Hit flags** | `dualipa_hit`, `imaging_hit`, `ppi_hit` |
| **Clinical** | `ClinicalSignificance`, `ClinSigSimple`, `clinvar_clnsig_clean` |
| **Gene/variant** | `symbol`, `aa_change`, `gene_variant` |
| **RSA** | `rsa` (0-1 continuous) - binarize at 0.25 |
| **pLDDT** | `plddt` (19-99 continuous) - binarize at 70 |
| **Secondary structure** | `g2p_Secondary structure (DSSP 3-state)*` |
| **Disease** | `disease_modules`, `disease_module_count`, `OMIM_disease_names` |
| **Inheritance** | `inheritance_pattern`, `inheritance_source` |
| **PPI details** | `ppi_num_edges` |

**NOT AVAILABLE** (panels requiring these should be skipped):
- Regional constraint (gnomAD_RMC)
- Grantham score (can be computed from `aa_change` if needed)
- Edgotype classification

---

## Panel Specifications

### BLOCK A: Pathogenic vs Benign Discrimination

#### Panel A: Per-Assay Score Distributions (P vs B)

**Required columns**: `ClinicalSignificance`, `dualipa_score`, `imaging_score`, `ppi_score`

**Logic**:
```python
# Create binary classification
def classify_clinvar(cs):
    if pd.isna(cs):
        return 'Other'
    cs = str(cs)
    if 'Benign' in cs and 'Pathogenic' not in cs and 'Conflicting' not in cs:
        return 'B/LB'
    elif 'Pathogenic' in cs and 'Benign' not in cs and 'Conflicting' not in cs:
        return 'P/LP'
    else:
        return 'Other'

df['clin_binary'] = df['ClinicalSignificance'].apply(classify_clinvar)
plot_df = df[df['clin_binary'].isin(['P/LP', 'B/LB'])]
```

**Plot type**: 3 violin plots (one per assay), split by P/LP vs B/LB
**Colors**: Use `clinvar_palette` - P/LP (#CA7682), B/LB (#1D7AAB)
**Add**: Mann-Whitney U test p-values, effect size annotations

---

#### Panel B: Assay Orthogonality (UpSet Plot)

**Required columns**: `dualipa_hit`, `imaging_hit`, `ppi_hit`, `ClinicalSignificance`

**Logic**:
```python
# For pathogenic variants only
patho_df = df[df['clin_binary'] == 'P/LP'].copy()

# Use existing hit columns directly
# dualipa_hit, imaging_hit, ppi_hit are boolean/nan
```

**Plot type**: UpSet plot showing intersection sizes for pathogenic variants
**Colors**: Use `assay_palette`

---

#### Panel C: Gene-Level Heatmap (Assay Responsiveness)

**Required columns**: `symbol`, `ClinicalSignificance`, `dualipa_hit`, `imaging_hit`, `ppi_hit`

**Logic**:
```python
# For each gene with pathogenic variants in >=2 assays:
# Calculate fraction of pathogenic variants that are hits per assay
# Filter to genes with >=3 pathogenic variants tested

gene_stats = df[df['clin_binary'] == 'P/LP'].groupby('symbol').agg({
    'dualipa_hit': ['sum', 'count'],
    'imaging_hit': ['sum', 'count'],
    'ppi_hit': ['sum', 'count']
})
# Calculate fractions: hits / tested
```

**Plot type**: Clustered heatmap
- Rows = genes (filter to >=3 pathogenic variants)
- Columns = 3 assays
- Color = fraction of pathogenic variants that are hits (0 to 1)
**Color scale**: Sequential, white to orange (`#E69F00`)
**Clustering**: Hierarchical on rows

---

### BLOCK B: Structural Context

---

#### Threshold References

**RSA threshold (0.25)**:
- Miller S, Janin J, Lesk AM, Chothia C. (1987) "Interior and surface of monomeric proteins." *J Mol Biol* 196:641-656. [PMID: 3681970]
- Rost B, Sander C. (1994) "Conservation and prediction of solvent accessibility in protein families." *Proteins* 20:216-226. [PMID: 7892171]
- The 25% RSA threshold is the most commonly used cutoff for buried/exposed classification, roughly dividing residues into two equally-sized subsets.

**pLDDT threshold (70)**:
- Jumper J, et al. (2021) "Highly accurate protein structure prediction with AlphaFold." *Nature* 596:583-589. [PMID: 34265844]
- pLDDT > 70 corresponds to "confident" backbone predictions; pLDDT > 90 indicates high accuracy for both backbone and side chains.
- See also: [EBI AlphaFold pLDDT documentation](https://www.ebi.ac.uk/training/online/courses/alphafold/inputs-and-outputs/evaluating-alphafolds-predicted-structures-using-confidence-scores/plddt-understanding-local-confidence/)

---

#### Panel D: RSA (Relative Solvent Accessibility) Stratification

**Required columns**: `rsa`, `dualipa_hit`, `imaging_hit`, `ppi_hit`, `dualipa_score`, `imaging_score`, `ppi_score`, `ClinicalSignificance`

**Binary threshold**: RSA < 0.25 = Buried, RSA >= 0.25 = Exposed

Dataset characteristics:
- ~45%/55% split (612 buried, 755 exposed)
- Pathogenic variants enriched in buried residues (56% vs 44%)
- Benign variants enriched in exposed residues (17% vs 83%)

**Logic**:
```python
from scipy.stats import fisher_exact, mannwhitneyu

RSA_THRESHOLD = 0.25
df['rsa_binary'] = df['rsa'].apply(
    lambda x: 'Buried' if x < RSA_THRESHOLD else 'Exposed'
)

patho_df = df[df['clin_binary'] == 'P/LP']
```

##### Panel D1: Hit Rate Comparison (Fisher's Exact Test)

**Plot type**: Grouped bar chart with hit rates
- X-axis: RSA bins (Buried, Exposed)
- Y-axis: Hit rate (fraction of hits / total tested)
- Grouped by assay

**Statistics** (per assay):
```python
def compute_fisher_stats(df, assay_hit_col, group_col):
    """Compute Fisher's exact test for hit enrichment."""
    results = {}
    for assay in ['dualipa_hit', 'imaging_hit', 'ppi_hit']:
        # Build contingency table: [[buried_hit, buried_nohit], [exposed_hit, exposed_nohit]]
        buried = df[df[group_col] == 'Buried'][assay].dropna()
        exposed = df[df[group_col] == 'Exposed'][assay].dropna()

        table = [
            [buried.sum(), len(buried) - buried.sum()],
            [exposed.sum(), len(exposed) - exposed.sum()]
        ]
        odds_ratio, p_value = fisher_exact(table, alternative='two-sided')

        results[assay] = {
            'odds_ratio': odds_ratio,
            'p_value': p_value,
            'buried_rate': buried.mean(),
            'exposed_rate': exposed.mean()
        }
    return results
```

**Annotation format**:
- Above each pair of bars: `OR = X.XX, p = X.XXe-X`
- Use `*` for p < 0.05, `**` for p < 0.01, `***` for p < 0.001

##### Panel D2: Score Distribution Comparison (Mann-Whitney U Test)

**Plot type**: Violin/box plots with continuous scores
- X-axis: RSA bins (Buried, Exposed)
- Y-axis: Assay score (continuous)
- Separate subplot per assay

**Statistics** (per assay):
```python
def compute_mannwhitney_stats(df, assay_score_col, group_col):
    """Compute Mann-Whitney U test for score differences."""
    results = {}
    for assay in ['dualipa_score', 'imaging_score', 'ppi_score']:
        buried = df[df[group_col] == 'Buried'][assay].dropna()
        exposed = df[df[group_col] == 'Exposed'][assay].dropna()

        stat, p_value = mannwhitneyu(buried, exposed, alternative='two-sided')

        # Effect size: rank-biserial correlation
        n1, n2 = len(buried), len(exposed)
        r = 1 - (2 * stat) / (n1 * n2)  # rank-biserial correlation

        results[assay] = {
            'U_statistic': stat,
            'p_value': p_value,
            'effect_size_r': r,
            'buried_median': buried.median(),
            'exposed_median': exposed.median()
        }
    return results
```

**Annotation format**:
- Above each violin: `U = XXXX, p = X.XXe-X, r = X.XX`
- r = rank-biserial correlation (effect size)

**Colors**: Use `assay_palette` for bars/violins

---

#### Panel E: pLDDT Stratification

**Required columns**: `plddt`, `dualipa_hit`, `imaging_hit`, `ppi_hit`, `dualipa_score`, `imaging_score`, `ppi_score`, `ClinicalSignificance`

**Binary threshold**: pLDDT >= 70 = Structured, pLDDT < 70 = Disordered

Dataset characteristics:
- ~85%/15% split (1156 structured, 211 disordered)
- Pathogenic variants enriched in structured regions (92%)
- Benign variants more represented in disordered regions (42%)

**Logic**:
```python
PLDDT_THRESHOLD = 70
df['plddt_binary'] = df['plddt'].apply(
    lambda x: 'Structured' if x >= PLDDT_THRESHOLD else 'Disordered'
)
```

##### Panel E1: Hit Rate Comparison (Fisher's Exact Test)

**Plot type**: Grouped bar chart with hit rates (same as D1)
- X-axis: pLDDT bins (Disordered, Structured)
- Y-axis: Hit rate (fraction)
- Grouped by assay

**Statistics**: Same Fisher's exact test approach as Panel D1
**Annotation format**: `OR = X.XX, p = X.XXe-X`

##### Panel E2: Score Distribution Comparison (Mann-Whitney U Test)

**Plot type**: Violin/box plots with continuous scores (same as D2)
- X-axis: pLDDT bins (Disordered, Structured)
- Y-axis: Assay score (continuous)
- Separate subplot per assay

**Statistics**: Same Mann-Whitney U approach as Panel D2
**Annotation format**: `U = XXXX, p = X.XXe-X, r = X.XX`

**Colors**: Use `assay_palette`

---

#### Panel F: Secondary Structure Context

**Required columns**: `g2p_Secondary structure (DSSP 3-state)*`, `dualipa_hit`, `imaging_hit`, `ppi_hit`, `ClinicalSignificance`

**Logic**:
```python
# Column: 'g2p_Secondary structure (DSSP 3-state)*'
# Values: 'H (helix)', 'C (loop/coil)', 'B (undefined)'

df['ss_3state'] = df['g2p_Secondary structure (DSSP 3-state)*']
# Filter out NaN values for this analysis
```

**Plot type**: Grouped bar chart
- X-axis: Secondary structure (helix, strand/undefined, coil)
- Y-axis: Hit rate
- Grouped by assay
**Colors**: Use `assay_palette`

---

### BLOCK C: Disease Biology Patterns

#### Panel G: Disease Module Stratification

**Required columns**: `disease_modules`, `disease_module_count`, `symbol`, `dualipa_hit`, `imaging_hit`, `ppi_hit`

**Logic**:
```python
# Parse disease_modules (may be semicolon-separated)
# For each module with >=5 variants, calculate hit rates per assay

df['modules_list'] = df['disease_modules'].str.split(';')
# Explode and aggregate
```

**Plot type**: Heatmap
- Rows = disease modules (filter to modules with >=5 variants)
- Columns = 3 assays
- Color = fraction of pathogenic variants that are hits
**Color scale**: Sequential, white to blue (`#0072B2`)

---

#### Panel H: Inheritance Mode Patterns

**Required columns**: `inheritance_pattern`, `dualipa_hit`, `imaging_hit`, `ppi_hit`, `ClinicalSignificance`

**Logic**:
```python
# Use inheritance_pattern column directly
# Values: 'autosomal dominant', 'autosomal recessive', 'x-linked recessive',
#         'x-linked dominant', 'x-linked', 'multiple', 'unknown', 'somatic'

# Simplify to major categories
inheritance_map = {
    'autosomal dominant': 'Dominant',
    'autosomal recessive': 'Recessive',
    'x-linked recessive': 'X-linked',
    'x-linked dominant': 'X-linked',
    'x-linked': 'X-linked',
    'multiple': 'Multiple',
}
df['inheritance_simple'] = df['inheritance_pattern'].map(inheritance_map)
```

**Plot type**: Grouped bar chart
- X-axis: Inheritance mode (Dominant, Recessive, X-linked)
- Y-axis: Hit rate
- Grouped by assay
**Colors**: Use `assay_palette` for bars

---

### BLOCK D: Mechanistic Diversity

#### Panel I: Within-Gene Mechanistic Heterogeneity

**Required columns**: `symbol`, `aa_change`, `dualipa_hit`, `imaging_hit`, `ppi_hit`

**Logic**:
```python
# Identify genes with:
# - >=5 pathogenic variants tested in multiple assays
# - Hits in >=2 different assays (demonstrating heterogeneity)

# Extract protein position from aa_change
import re
def extract_position(aa):
    match = re.search(r'(\d+)', str(aa))
    return int(match.group(1)) if match else None

df['aa_pos'] = df['aa_change'].apply(extract_position)

# Select exemplar genes with good coverage
candidates = df.groupby('symbol').filter(
    lambda x: (x['clin_binary'] == 'P/LP').sum() >= 5
)['symbol'].unique()
```

**Plot type**: Lollipop plot or scatter for 2-3 selected genes
- X-axis: Protein position
- Y-axis: Categorical (which assay detected hit)
- Color by assay
**Colors**: Use `assay_palette`

---

## Panels NOT Available (Skip)

The following panels from the original plan cannot be generated due to missing columns:

| Panel | Reason | Missing Column |
|-------|--------|----------------|
| ~~Panel G: gnomAD Constraint~~ | No regional constraint data | `gnomAD_RMC` |
| ~~Panel H: Grantham Distance~~ | No pre-computed Grantham score | `Grantham_score` |
| ~~Panel I: Edgetic Analysis~~ | No edgotype classification | `edgotype` |

---

## Execution Instructions

### Step 1: Setup and Load Data

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import sys

# Import existing palettes
sys.path.insert(0, '../../1_allele_collection/2_analysis')
from plot_utils import palette_dict, clinvar_order_dict

# Load data
DATA_PATH = '../3_outputs/0_all_gene_variants_assayed_summary.tsv'
df = pd.read_csv(DATA_PATH, sep='\t', low_memory=False)

# Setup ClinVar binary classification
def classify_clinvar(cs):
    if pd.isna(cs):
        return 'Other'
    cs = str(cs)
    if 'Benign' in cs and 'Pathogenic' not in cs and 'Conflicting' not in cs:
        return 'B/LB'
    elif 'Pathogenic' in cs and 'Benign' not in cs and 'Conflicting' not in cs:
        return 'P/LP'
    else:
        return 'Other'

df['clin_binary'] = df['ClinicalSignificance'].apply(classify_clinvar)

print(f"Total variants: {len(df)}")
print(f"P/LP: {(df['clin_binary'] == 'P/LP').sum()}")
print(f"B/LB: {(df['clin_binary'] == 'B/LB').sum()}")
```

### Step 2: Define Palettes

```python
# ClinVar palette from plot_utils
clinvar_palette = palette_dict['clinvar_clnsig_clean']
clinvar_binary_palette = {'P/LP': '#CA7682', 'B/LB': '#1D7AAB'}

# Assay palette
assay_palette = {
    'Dual-IPA': '#E69F00',
    'Imaging': '#0072B2',
    'PPI': '#009E73',
}
assay_order = ['Dual-IPA', 'Imaging', 'PPI']
```

### Step 3: Generate Each Panel

See individual panel specifications above.

### Step 4: Output

- Save each panel as separate PDF: `figure2_panel_[a,b,c,d1,d2,e1,e2,f,g,h,i].pdf`
- Also save combined multi-panel figure: `figure2_combined.pdf`
- Export source data for each panel as CSV: `figure2_panel_[*]_data.csv`
- Export statistical test results as TSV: `figure2_stats_summary.tsv`

---

## Figure Specifications (Nature Format)

- **Resolution**: 300 dpi minimum
- **Font**: Arial or Helvetica, 6-8pt for labels
- **Panel labels**: 8pt bold lowercase (a, b, c...)
- **Figure width**:
  - Single column: 89mm
  - 1.5 column: 120-136mm
  - Double column: 183mm
- **Axes**: Include units in parentheses
- **Statistical annotations**: Show p-values or significance stars

---

## Summary Checklist

| Panel | Description | Required Columns | Status | Statistical Test |
|-------|-------------|------------------|--------|------------------|
| A | Score distributions (P vs B) | `*_score`, `ClinicalSignificance` | **Available** | Mann-Whitney U |
| B | Assay orthogonality (UpSet) | `*_hit`, `ClinicalSignificance` | **Available** | - |
| C | Gene-level heatmap | `symbol`, `*_hit`, `ClinicalSignificance` | **Available** | - |
| D1 | RSA stratification (hit rates) | `rsa`, `*_hit` | **Available** | Fisher's exact |
| D2 | RSA stratification (scores) | `rsa`, `*_score` | **Available** | Mann-Whitney U |
| E1 | pLDDT stratification (hit rates) | `plddt`, `*_hit` | **Available** | Fisher's exact |
| E2 | pLDDT stratification (scores) | `plddt`, `*_score` | **Available** | Mann-Whitney U |
| F | Secondary structure | `g2p_Secondary structure (DSSP 3-state)*`, `*_hit` | **Available** | Chi-square |
| G | Disease modules | `disease_modules`, `*_hit` | **Available** | - |
| H | Inheritance patterns | `inheritance_pattern`, `*_hit` | **Available** | Chi-square |
| I | Within-gene heterogeneity | `symbol`, `aa_change`, `*_hit` | **Available** | - |

**Panels removed** (missing data):
- ~~gnomAD constraint~~ - no `gnomAD_RMC`
- ~~Grantham distance~~ - no `Grantham_score`
- ~~Edgetic analysis~~ - no `edgotype`
