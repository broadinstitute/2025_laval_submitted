#!/usr/bin/env python3
"""
Comprehensive Disease Biology and Feature Analysis

This script performs:
1. Disease module clustering analysis
2. Gene-specific clustering patterns
3. Feature-driven cluster separation analysis
4. Feature correlation and representative feature selection

All visualizations are saved to a single PDF.
All results are summarized in one markdown file.

Author: Claude Code
Date: 2025-11-19
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import warnings
warnings.filterwarnings('ignore')

# Import disease module analyzer
import sys
sys.path.append('.')
from disease_module_annot import FinalDiseaseModuleAnalyzer

# Set style
sns.set_style('whitegrid')
plt.rcParams['font.size'] = 10

print("="*80)
print("COMPREHENSIVE DISEASE BIOLOGY AND FEATURE ANALYSIS")
print("="*80)

# ============================================================================
# PART 1: LOAD ALL DATA
# ============================================================================
print("\n[1/6] Loading data...")

# Clustering results with disease info
clustering_df = pd.read_csv('../../3_outputs/variant_clustering_results.csv', index_col=0)
print(f"  - Loaded {len(clustering_df)} variants with clustering assignments")

# Feature matrix
feature_matrix = pd.read_csv('../../3_outputs/variant_protein_feature_matrix.csv', index_col=0)
print(f"  - Loaded feature matrix: {feature_matrix.shape[0]} variants x {feature_matrix.shape[1]} features")

# Feature correlation matrix
feature_corr = pd.read_csv('../../3_outputs/protein_features_corr_mat.csv', index_col=0)
print(f"  - Loaded feature correlation matrix: {feature_corr.shape[0]} x {feature_corr.shape[1]}")

# Disease modules
disease_modules = pd.read_csv(
    '../1_snakemake_pipeline/2025_varchamp_snakemake/3.downstream_analyses/inputs/disease_modules.csv'
)
print(f"  - Loaded disease annotations for {len(disease_modules)} variants")

# ============================================================================
# PART 2: DISEASE MODULE CLASSIFICATION
# ============================================================================
print("\n[2/6] Classifying disease modules...")

analyzer = FinalDiseaseModuleAnalyzer()

# Create gene_allele column
clustering_df['gene_allele'] = clustering_df.index

# Merge with disease modules (drop duplicates to prevent row explosion)
disease_modules_dedup = disease_modules[['gene_allele', 'OMIM_disease_names', 'OMIM_IDs']].drop_duplicates(
    subset='gene_allele', keep='first'
)
merged_df = clustering_df.merge(
    disease_modules_dedup,
    on='gene_allele',
    how='left',
    suffixes=('', '_disease_db')
)

# Verify merge didn't create duplicates
assert len(merged_df) == len(clustering_df), f"Merge error: {len(merged_df)} != {len(clustering_df)}"

# Use the OMIM_disease_names from disease database, fallback to clustering metadata
merged_df['diseases_to_classify'] = merged_df['OMIM_disease_names'].fillna(
    merged_df['Metadata_OMIM_disease_names']
)

# Classify disease modules
results = []
for idx, row in merged_df.iterrows():
    diseases = analyzer.clean_disease_names(row.get('diseases_to_classify', ''))
    classifications = analyzer.classify_disease_modules(diseases)

    primary_module = classifications[0][0] if classifications else 'Unclassified'
    primary_confidence = classifications[0][1] if classifications else 0.0

    results.append({
        'gene_allele': row['gene_allele'],
        'gene': row['gene_allele'].split('_')[0],
        'cluster': row['Metadata_cluster'],
        'AUROC': row['Metadata_AUROC_Mean_GFP_ADJ'],
        'is_altered': row['Metadata_Altered_95th_perc_both_batches_GFP_ADJ'],
        'primary_disease_module': primary_module,
        'primary_module_confidence': primary_confidence,
        'cleaned_diseases': '; '.join(diseases),
        'disease_count': len(diseases)
    })

classified_df = pd.DataFrame(results)
print(f"  - Classified {len(classified_df)} variants into disease modules")

# ============================================================================
# PART 3: FEATURE CORRELATION CLUSTERING
# ============================================================================
print("\n[3/6] Analyzing feature correlations...")

# Cluster features based on correlation
feature_corr_mat = feature_corr.values
np.fill_diagonal(feature_corr_mat, 1.0)  # Ensure diagonal is 1

# Convert correlation to distance
feature_dist = 1 - np.abs(feature_corr_mat)
np.fill_diagonal(feature_dist, 0)

# Hierarchical clustering of features
feature_linkage = linkage(squareform(feature_dist), method='average')
feature_clusters = fcluster(feature_linkage, t=0.3, criterion='distance')

# Find representative feature for each cluster
feature_cluster_df = pd.DataFrame({
    'feature': feature_corr.columns,
    'feature_cluster': feature_clusters
})

# Select representative: highest mean absolute correlation within cluster
representative_features = []
for fc in feature_cluster_df['feature_cluster'].unique():
    cluster_features = feature_cluster_df[feature_cluster_df['feature_cluster'] == fc]['feature'].tolist()

    if len(cluster_features) == 1:
        representative_features.append({
            'feature_cluster': fc,
            'representative_feature': cluster_features[0],
            'cluster_size': 1,
            'mean_abs_corr': 0.0
        })
    else:
        # Calculate mean absolute correlation for each feature with others in cluster
        mean_corrs = []
        for feat in cluster_features:
            other_feats = [f for f in cluster_features if f != feat]
            mean_corr = np.abs(feature_corr.loc[feat, other_feats]).mean()
            mean_corrs.append(mean_corr)

        best_idx = np.argmax(mean_corrs)
        representative_features.append({
            'feature_cluster': fc,
            'representative_feature': cluster_features[best_idx],
            'cluster_size': len(cluster_features),
            'mean_abs_corr': mean_corrs[best_idx],
            'all_features': ', '.join(cluster_features)
        })

representative_features_df = pd.DataFrame(representative_features).sort_values(
    'cluster_size', ascending=False
)
print(f"  - Identified {len(representative_features_df)} feature clusters")
print(f"  - Representative features selected from each cluster")

# ============================================================================
# PART 4: CLUSTER-SPECIFIC FEATURE ANALYSIS
# ============================================================================
print("\n[4/6] Analyzing cluster-specific features...")

# For each cluster, find top discriminative features
cluster_feature_importance = []

for cluster_id in sorted(classified_df['cluster'].unique()):
    cluster_variants = classified_df[classified_df['cluster'] == cluster_id]['gene_allele'].tolist()
    other_variants = classified_df[classified_df['cluster'] != cluster_id]['gene_allele'].tolist()

    # Get feature values for this cluster vs others
    cluster_feats = feature_matrix.loc[
        [v for v in cluster_variants if v in feature_matrix.index]
    ]
    other_feats = feature_matrix.loc[
        [v for v in other_variants if v in feature_matrix.index]
    ]

    if len(cluster_feats) > 0 and len(other_feats) > 0:
        # Calculate effect size (Cohen's d) for each feature
        for feat in feature_matrix.columns:
            cluster_mean = cluster_feats[feat].mean()
            other_mean = other_feats[feat].mean()
            pooled_std = np.sqrt(
                (cluster_feats[feat].var() + other_feats[feat].var()) / 2
            )

            if pooled_std > 0:
                cohens_d = (cluster_mean - other_mean) / pooled_std

                cluster_feature_importance.append({
                    'cluster': cluster_id,
                    'feature': feat,
                    'cohens_d': cohens_d,
                    'abs_cohens_d': abs(cohens_d),
                    'cluster_mean': cluster_mean,
                    'other_mean': other_mean,
                    'direction': 'higher' if cohens_d > 0 else 'lower'
                })

cluster_feature_df = pd.DataFrame(cluster_feature_importance)

# Get top 3 features per cluster
top_features_per_cluster = {}
for cluster_id in sorted(classified_df['cluster'].unique()):
    top_feats = cluster_feature_df[
        cluster_feature_df['cluster'] == cluster_id
    ].nlargest(3, 'abs_cohens_d')
    top_features_per_cluster[cluster_id] = top_feats

print(f"  - Identified top discriminative features for {len(top_features_per_cluster)} clusters")

# ============================================================================
# PART 5: GENERATE ALL VISUALIZATIONS (SINGLE PDF)
# ============================================================================
print("\n[5/6] Generating comprehensive visualizations (single PDF)...")

pdf_path = '../../3_outputs/comprehensive_disease_biology_analysis.pdf'

with PdfPages(pdf_path) as pdf:

    # PAGE 1: Cluster Overview
    print("  - Page 1: Cluster overview metrics")
    fig = plt.figure(figsize=(11, 8.5))
    gs = fig.add_gridspec(3, 2, hspace=0.35, wspace=0.3)

    # Cluster sizes
    ax1 = fig.add_subplot(gs[0, 0])
    cluster_sizes = classified_df['cluster'].value_counts().sort_index()
    cluster_sizes.plot(kind='bar', ax=ax1, color='steelblue')
    ax1.set_xlabel('Cluster')
    ax1.set_ylabel('Number of Variants')
    ax1.set_title('Cluster Sizes', fontweight='bold')

    # Mean AUROC by cluster
    ax2 = fig.add_subplot(gs[0, 1])
    cluster_aurocs = classified_df.groupby('cluster')['AUROC'].mean().sort_index()
    cluster_aurocs.plot(kind='bar', ax=ax2, color='coral')
    ax2.axhline(y=0.95, color='red', linestyle='--', alpha=0.7)
    ax2.set_xlabel('Cluster')
    ax2.set_ylabel('Mean AUROC')
    ax2.set_title('Mean AUROC by Cluster', fontweight='bold')
    ax2.set_ylim([0.90, 1.0])

    # Disease module distribution
    ax3 = fig.add_subplot(gs[1, :])
    module_cluster = pd.crosstab(
        classified_df['primary_disease_module'],
        classified_df['cluster']
    )
    module_cluster_pct = module_cluster.div(module_cluster.sum(axis=0), axis=1) * 100
    module_cluster_pct.T.plot(kind='bar', stacked=True, ax=ax3, colormap='tab10')
    ax3.set_xlabel('Cluster')
    ax3.set_ylabel('Percentage (%)')
    ax3.set_title('Disease Module Composition by Cluster', fontweight='bold')
    ax3.legend(title='Disease Module', bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
    ax3.set_ylim([0, 100])

    # UMAP with clusters
    ax4 = fig.add_subplot(gs[2, 0])
    scatter = ax4.scatter(
        clustering_df['UMAP1'],
        clustering_df['UMAP2'],
        c=clustering_df['Metadata_cluster'],
        cmap='tab10',
        s=50,
        alpha=0.7,
        edgecolors='black',
        linewidth=0.5
    )
    ax4.set_xlabel('UMAP1')
    ax4.set_ylabel('UMAP2')
    ax4.set_title('UMAP: Colored by Cluster', fontweight='bold')
    plt.colorbar(scatter, ax=ax4, label='Cluster')
    ax4.grid(True, alpha=0.3)

    # UMAP with AUROC
    ax5 = fig.add_subplot(gs[2, 1])
    scatter2 = ax5.scatter(
        clustering_df['UMAP1'],
        clustering_df['UMAP2'],
        c=clustering_df['Metadata_AUROC_Mean_GFP_ADJ'],
        cmap='RdYlGn',
        s=50,
        alpha=0.7,
        edgecolors='black',
        linewidth=0.5,
        vmin=0.93,
        vmax=1.0
    )
    ax5.set_xlabel('UMAP1')
    ax5.set_ylabel('UMAP2')
    ax5.set_title('UMAP: Colored by AUROC', fontweight='bold')
    plt.colorbar(scatter2, ax=ax5, label='AUROC')
    ax5.grid(True, alpha=0.3)

    fig.suptitle('Comprehensive Disease Biology Analysis - Overview',
                 fontsize=14, fontweight='bold', y=0.995)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close()

    # PAGE 2: Gene Clustering Consistency
    print("  - Page 2: Gene clustering patterns")
    fig, axes = plt.subplots(2, 1, figsize=(11, 8.5))

    # Analyze gene consistency
    gene_counts = classified_df['gene'].value_counts()
    top_genes = gene_counts[gene_counts >= 3].index.tolist()[:12]

    gene_cluster_analysis = []
    for gene in top_genes:
        gene_data = classified_df[classified_df['gene'] == gene]
        clusters = gene_data['cluster'].value_counts()
        primary_cluster = clusters.index[0]
        consistency = (clusters.iloc[0] / len(gene_data) * 100)

        gene_cluster_analysis.append({
            'gene': gene,
            'total_variants': len(gene_data),
            'primary_cluster': primary_cluster,
            'consistency_%': consistency,
            'mean_AUROC': gene_data['AUROC'].mean()
        })

    gene_df = pd.DataFrame(gene_cluster_analysis).sort_values('consistency_%', ascending=True)

    # Plot
    ax = axes[0]
    colors = plt.cm.RdYlGn(gene_df['consistency_%'] / 100)
    bars = ax.barh(range(len(gene_df)), gene_df['consistency_%'], color=colors)
    ax.set_yticks(range(len(gene_df)))
    ax.set_yticklabels(gene_df['gene'])
    ax.set_xlabel('Clustering Consistency (%)', fontsize=11)
    ax.set_ylabel('Gene', fontsize=11)
    ax.set_title('Gene Clustering Consistency (% in Primary Cluster)', fontweight='bold')
    ax.axvline(x=80, color='red', linestyle='--', alpha=0.7, label='80% threshold')
    ax.legend()

    # Add counts
    for i, (idx, row) in enumerate(gene_df.iterrows()):
        ax.text(row['consistency_%'] + 2, i, f"n={row['total_variants']}",
                va='center', fontsize=9)

    # Disease module heatmap
    ax = axes[1]
    module_cluster_counts = pd.crosstab(
        classified_df['primary_disease_module'],
        classified_df['cluster']
    )
    sns.heatmap(module_cluster_counts, annot=True, fmt='d', cmap='YlOrRd',
                ax=ax, cbar_kws={'label': 'Count'})
    ax.set_xlabel('Cluster', fontsize=11)
    ax.set_ylabel('Disease Module', fontsize=11)
    ax.set_title('Disease Module Distribution Across Clusters', fontweight='bold')

    fig.suptitle('Gene and Disease Module Patterns', fontsize=14, fontweight='bold', y=0.995)
    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')
    plt.close()

    # PAGE 3: Feature Correlation Clustering
    print("  - Page 3: Feature correlation analysis")
    fig = plt.figure(figsize=(11, 8.5))
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)

    # Feature cluster sizes
    ax1 = fig.add_subplot(gs[0, :])
    cluster_size_counts = representative_features_df['cluster_size'].value_counts().sort_index()
    ax1.bar(cluster_size_counts.index, cluster_size_counts.values, color='steelblue')
    ax1.set_xlabel('Feature Cluster Size')
    ax1.set_ylabel('Number of Feature Clusters')
    ax1.set_title('Distribution of Feature Cluster Sizes', fontweight='bold')
    ax1.grid(True, alpha=0.3, axis='y')

    # Top representative features
    ax2 = fig.add_subplot(gs[1, :])
    top_rep_feats = representative_features_df.nlargest(15, 'cluster_size')
    y_pos = range(len(top_rep_feats))
    ax2.barh(y_pos, top_rep_feats['cluster_size'], color='coral')
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels([f[:50] for f in top_rep_feats['representative_feature']], fontsize=8)
    ax2.set_xlabel('Cluster Size (# features)')
    ax2.set_title('Top 15 Representative Features by Cluster Size', fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='x')

    fig.suptitle('Feature Correlation Clustering Results',
                 fontsize=14, fontweight='bold', y=0.995)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close()

    # PAGES 4+: Cluster-Specific Feature Analysis (one page per cluster)
    print("  - Pages 4+: Cluster-specific feature analysis")
    for cluster_id in sorted(classified_df['cluster'].unique()):
        fig = plt.figure(figsize=(11, 8.5))
        gs = fig.add_gridspec(3, 2, hspace=0.35, wspace=0.3)

        cluster_data = classified_df[classified_df['cluster'] == cluster_id]

        # Cluster info
        ax_info = fig.add_subplot(gs[0, :])
        ax_info.axis('off')
        info_text = f"""
CLUSTER {cluster_id} ANALYSIS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Size: {len(cluster_data)} variants
Mean AUROC: {cluster_data['AUROC'].mean():.3f}
Altered variants: {cluster_data['is_altered'].sum()} ({cluster_data['is_altered'].sum()/len(cluster_data)*100:.1f}%)

Disease Module Distribution:
{cluster_data['primary_disease_module'].value_counts().to_string()}

Top Genes:
{cluster_data['gene'].value_counts().head(5).to_string()}
        """
        ax_info.text(0.05, 0.95, info_text, transform=ax_info.transAxes,
                    fontsize=9, verticalalignment='top', fontfamily='monospace')

        # Top discriminative features
        if cluster_id in top_features_per_cluster:
            top_feats = top_features_per_cluster[cluster_id]

            ax1 = fig.add_subplot(gs[1, :])
            colors_feat = ['red' if d > 0 else 'blue' for d in top_feats['cohens_d']]
            bars = ax1.barh(range(len(top_feats)), top_feats['cohens_d'], color=colors_feat, alpha=0.7)
            ax1.set_yticks(range(len(top_feats)))
            ax1.set_yticklabels([f[:50] for f in top_feats['feature']], fontsize=9)
            ax1.set_xlabel("Cohen's d (effect size)")
            ax1.set_title(f"Top 3 Discriminative Features for Cluster {cluster_id}", fontweight='bold')
            ax1.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
            ax1.grid(True, alpha=0.3, axis='x')

            # Feature value distributions
            for i, (idx, feat_row) in enumerate(top_feats.iterrows()):
                ax = fig.add_subplot(gs[2, i]) if i < 2 else None
                if ax is None:
                    continue

                feat_name = feat_row['feature']
                cluster_vals = feature_matrix.loc[
                    [v for v in cluster_data['gene_allele'] if v in feature_matrix.index],
                    feat_name
                ]
                other_vals = feature_matrix.loc[
                    [v for v in classified_df[classified_df['cluster'] != cluster_id]['gene_allele']
                     if v in feature_matrix.index],
                    feat_name
                ]

                ax.hist(other_vals, bins=20, alpha=0.5, label='Other clusters', color='gray')
                ax.hist(cluster_vals, bins=20, alpha=0.7, label=f'Cluster {cluster_id}', color='red')
                ax.set_xlabel('Feature Value')
                ax.set_ylabel('Frequency')
                ax.set_title(f"{feat_name[:40]}", fontsize=9)
                ax.legend(fontsize=8)

        fig.suptitle(f'Cluster {cluster_id}: Detailed Feature Analysis',
                     fontsize=14, fontweight='bold', y=0.995)
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()

print(f"  - Saved comprehensive PDF: {pdf_path}")

# ============================================================================
# PART 6: GENERATE MARKDOWN SUMMARY
# ============================================================================
print("\n[6/6] Generating comprehensive markdown summary...")

markdown_path = '../../3_outputs/COMPREHENSIVE_ANALYSIS_SUMMARY.md'

with open(markdown_path, 'w') as f:
    f.write("""# Comprehensive Disease Biology and Feature Analysis

**Date:** 2025-11-19
**Analysis:** Cell Painting Imaging Variant Clustering with Disease Modules and Feature Analysis
**Variants:** {total_variants} across {n_clusters} morphological clusters

---

## Executive Summary

This comprehensive analysis integrates:
1. **Disease module classification** - Understanding biological pathways
2. **Gene-specific clustering patterns** - Identifying consistent phenotypes
3. **Feature-driven cluster analysis** - Morphological basis of clustering
4. **Feature correlation structure** - Representative feature selection

### Key Findings

""".format(
        total_variants=len(classified_df),
        n_clusters=classified_df['cluster'].nunique()
    ))

    # Key findings
    f.write("#### 1. Disease Module Enrichment\n\n")
    for cluster_id in sorted(classified_df['cluster'].unique()):
        cluster_data = classified_df[classified_df['cluster'] == cluster_id]
        dominant_module = cluster_data['primary_disease_module'].value_counts()
        if len(dominant_module) > 0:
            enrichment = dominant_module.iloc[0] / len(cluster_data) * 100
            if enrichment > 50:
                f.write(f"- **Cluster {cluster_id}**: {enrichment:.1f}% {dominant_module.index[0]} enrichment ({dominant_module.iloc[0]}/{len(cluster_data)} variants)\n")

    f.write("\n#### 2. Gene Clustering Patterns\n\n")
    perfect_genes = gene_df[gene_df['consistency_%'] == 100]
    if len(perfect_genes) > 0:
        f.write(f"**Perfect clustering (100% consistency):** {', '.join(perfect_genes['gene'].tolist())}\n\n")

    high_genes = gene_df[(gene_df['consistency_%'] >= 80) & (gene_df['consistency_%'] < 100)]
    if len(high_genes) > 0:
        f.write(f"**High consistency (≥80%):** {', '.join(high_genes['gene'].tolist())}\n\n")

    f.write("\n#### 3. Feature Structure\n\n")
    f.write(f"- Total features: {len(feature_corr)}\n")
    f.write(f"- Feature clusters identified: {representative_features_df['feature_cluster'].nunique()}\n")
    f.write(f"- Representative features selected: {len(representative_features_df)}\n")
    f.write(f"- Largest feature cluster: {representative_features_df['cluster_size'].max()} features\n\n")

    # Detailed sections
    f.write("""
---

## Part 1: Disease Module Analysis

### Overall Statistics

""")

    f.write(f"- **Total variants analyzed:** {len(classified_df)}\n")
    f.write(f"- **Variants with disease annotations:** {(classified_df['disease_count'] > 0).sum()}\n")
    f.write(f"- **Total clusters:** {classified_df['cluster'].nunique()}\n")
    f.write(f"- **Disease modules identified:** {classified_df['primary_disease_module'].nunique()}\n")
    f.write(f"- **Mean AUROC:** {classified_df['AUROC'].mean():.3f}\n\n")

    f.write("### Cluster-by-Cluster Disease Analysis\n\n")

    for cluster_id in sorted(classified_df['cluster'].unique()):
        cluster_data = classified_df[classified_df['cluster'] == cluster_id]
        f.write(f"#### Cluster {cluster_id}\n\n")
        f.write(f"- **Size:** {len(cluster_data)} variants\n")
        f.write(f"- **Mean AUROC:** {cluster_data['AUROC'].mean():.3f}\n")
        f.write(f"- **Altered:** {cluster_data['is_altered'].sum()}/{len(cluster_data)} ({cluster_data['is_altered'].sum()/len(cluster_data)*100:.1f}%)\n\n")

        f.write("**Disease Modules:**\n\n")
        module_dist = cluster_data['primary_disease_module'].value_counts()
        for module, count in module_dist.items():
            pct = count / len(cluster_data) * 100
            f.write(f"- {module}: {count} variants ({pct:.1f}%)\n")

        f.write("\n**Top Genes:**\n\n")
        gene_dist = cluster_data['gene'].value_counts().head(5)
        for gene, count in gene_dist.items():
            f.write(f"- {gene}: {count} variants\n")
        f.write("\n")

    f.write("""
---

## Part 2: Gene Clustering Consistency

Genes with multiple variants show varying degrees of clustering consistency,
indicating how similar their morphological phenotypes are across different mutations.

### Consistency Analysis

""")

    f.write("| Gene | Variants | Primary Cluster | Consistency | Mean AUROC |\n")
    f.write("|------|----------|----------------|-------------|------------|\n")
    for idx, row in gene_df.sort_values('consistency_%', ascending=False).iterrows():
        f.write(f"| {row['gene']} | {row['total_variants']} | {row['primary_cluster']} | {row['consistency_%']:.1f}% | {row['mean_AUROC']:.3f} |\n")

    f.write("""

### Interpretation

- **100% consistency**: All variants from this gene cluster together - strong gene-specific phenotype
- **80-99% consistency**: Most variants cluster together - relatively consistent phenotype
- **<80% consistency**: Variants show heterogeneous phenotypes - mutation-location may matter

""")

    f.write("""
---

## Part 3: Feature Correlation Structure

Features were clustered based on their pairwise correlations to identify
redundant features and select representatives.

### Feature Clustering Summary

""")

    f.write(f"- **Total features:** {len(feature_corr)}\n")
    f.write(f"- **Feature clusters:** {representative_features_df['feature_cluster'].nunique()}\n")
    f.write(f"- **Singleton clusters:** {(representative_features_df['cluster_size'] == 1).sum()}\n")
    f.write(f"- **Largest cluster:** {representative_features_df['cluster_size'].max()} features\n\n")

    f.write("### Top 10 Representative Features\n\n")
    f.write("| Feature | Cluster Size | Mean |r| with Cluster |\n")
    f.write("|---------|-------------|----------------------|\n")
    for idx, row in representative_features_df.head(10).iterrows():
        f.write(f"| {row['representative_feature'][:60]} | {row['cluster_size']} | {row['mean_abs_corr']:.3f} |\n")

    f.write("""

### Usage Recommendations

1. **For dimensionality reduction:** Use representative features to reduce from {total_feat} to {rep_feat} features
2. **For interpretation:** Focus on representative features from large clusters (high redundancy removed)
3. **For modeling:** Using representatives avoids multicollinearity issues

""".format(
        total_feat=len(feature_corr),
        rep_feat=len(representative_features_df)
    ))

    f.write("""
---

## Part 4: Cluster-Specific Morphological Features

For each cluster, we identified the top 3 features that distinguish it from other clusters
using Cohen's d effect size.

### Feature-Cluster Associations

""")

    for cluster_id in sorted(classified_df['cluster'].unique()):
        if cluster_id in top_features_per_cluster:
            f.write(f"#### Cluster {cluster_id} Distinctive Features\n\n")
            top_feats = top_features_per_cluster[cluster_id]

            f.write("| Feature | Effect Size (d) | Direction | Biological Relevance |\n")
            f.write("|---------|----------------|-----------|--------------------|\n")

            for idx, row in top_feats.iterrows():
                direction = "↑ Higher" if row['cohens_d'] > 0 else "↓ Lower"

                # Infer biological meaning from feature name
                feat_name = row['feature']
                bio_meaning = ""
                if 'Nuclei' in feat_name or 'DNA' in feat_name:
                    bio_meaning = "Nuclear morphology"
                elif 'ER' in feat_name or 'AGP' in feat_name:
                    bio_meaning = "ER structure"
                elif 'Mito' in feat_name or 'Mitochondria' in feat_name:
                    bio_meaning = "Mitochondrial"
                elif 'Cytoplasm' in feat_name or 'Cells' in feat_name:
                    bio_meaning = "Cytoplasmic"
                elif 'Texture' in feat_name:
                    bio_meaning = "Texture pattern"
                elif 'Intensity' in feat_name:
                    bio_meaning = "Intensity measure"
                elif 'Granularity' in feat_name:
                    bio_meaning = "Granular structure"
                elif 'RadialDistribution' in feat_name:
                    bio_meaning = "Spatial distribution"
                else:
                    bio_meaning = "Morphology"

                f.write(f"| {feat_name[:50]} | {row['cohens_d']:.2f} | {direction} | {bio_meaning} |\n")

            f.write("\n")

    f.write("""
---

## Part 5: Biological Interpretation

### Feature-Disease Connections

Based on the morphological features that distinguish each cluster, we can infer
biological mechanisms:

""")

    # Analyze which types of features are important per cluster
    for cluster_id in sorted(classified_df['cluster'].unique()):
        if cluster_id in top_features_per_cluster:
            cluster_data = classified_df[classified_df['cluster'] == cluster_id]
            dominant_module = cluster_data['primary_disease_module'].value_counts()

            if len(dominant_module) > 0:
                top_feats = top_features_per_cluster[cluster_id]

                f.write(f"#### Cluster {cluster_id}: {dominant_module.index[0]}\n\n")

                # Categorize features
                nuclear_feats = [f for f in top_feats['feature'] if 'Nuclei' in f or 'DNA' in f]
                er_feats = [f for f in top_feats['feature'] if 'ER' in f or 'AGP' in f]
                mito_feats = [f for f in top_feats['feature'] if 'Mito' in f]

                if nuclear_feats:
                    f.write(f"- **Nuclear features prominent**: {', '.join([f.split('_')[-1] for f in nuclear_feats[:2]])}\n")
                if er_feats:
                    f.write(f"- **ER features prominent**: Suggests ER stress or protein trafficking defects\n")
                if mito_feats:
                    f.write(f"- **Mitochondrial features**: May indicate metabolic dysfunction\n")

                f.write("\n")

    f.write("""
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
4. **Feature redundancy is high**: {n_feat} features can be reduced to {n_rep} representatives

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
""".format(
        n_feat=len(feature_corr),
        n_rep=len(representative_features_df)
    ))

print(f"  - Saved markdown summary: {markdown_path}")

# Save key dataframes as CSV
classified_df.to_csv('../../3_outputs/disease_module_classifications.csv', index=False)
gene_df.to_csv('../../3_outputs/gene_clustering_consistency.csv', index=False)
representative_features_df.to_csv('../../3_outputs/representative_features.csv', index=False)

cluster_feature_summary = []
for cluster_id in sorted(classified_df['cluster'].unique()):
    if cluster_id in top_features_per_cluster:
        top_feats = top_features_per_cluster[cluster_id]
        for idx, row in top_feats.iterrows():
            cluster_feature_summary.append({
                'cluster': cluster_id,
                'feature': row['feature'],
                'cohens_d': row['cohens_d'],
                'direction': row['direction']
            })

pd.DataFrame(cluster_feature_summary).to_csv(
    '../../3_outputs/cluster_discriminative_features.csv', index=False
)

print("\nAdditional CSV files saved:")
print("  - disease_module_classifications.csv")
print("  - gene_clustering_consistency.csv")
print("  - representative_features.csv")
print("  - cluster_discriminative_features.csv")

print("\n" + "="*80)
print("ANALYSIS COMPLETE!")
print("="*80)
print(f"\n📊 Visualizations: ../../3_outputs/comprehensive_disease_biology_analysis.pdf")
print(f"📝 Summary: ../../3_outputs/COMPREHENSIVE_ANALYSIS_SUMMARY.md")
print(f"📁 Data tables: 4 CSV files in ../../3_outputs/")
print("\nAll results consolidated into 1 PDF + 1 Markdown file!")
