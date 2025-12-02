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

Usage:
  python comprehensive_disease_biology_analysis.py --clustering misloc_var_clustering_results.csv --prefix Misloc_Var_

Arguments:
  --clustering: CSV file containing clustering results with Metadata_cluster column and feature columns
                (non-Metadata_ prefixed columns are treated as features)
  --prefix: Output prefix for all generated files (e.g., "Misloc_Var_")
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
import argparse

# Import disease module analyzer
import sys
sys.path.append('.')
from disease_module_annot import FinalDiseaseModuleAnalyzer

# ============================================================================
# ARGUMENT PARSING
# ============================================================================
parser = argparse.ArgumentParser(description='Comprehensive Disease Biology and Feature Analysis')
parser.add_argument('--clustering', type=str, required=True,
                    help='Clustering results CSV filename (in ../../3_outputs/). '
                         'Should contain Metadata_cluster and feature columns (non-Metadata_ prefixed).')
parser.add_argument('--prefix', type=str, required=True,
                    help='Output prefix for all generated files (e.g., "Misloc_Var_")')
args = parser.parse_args()

# Set style
sns.set_style('whitegrid')
plt.rcParams['font.size'] = 10

print("="*80)
print("COMPREHENSIVE DISEASE BIOLOGY AND FEATURE ANALYSIS")
print("="*80)
if args.prefix:
    print(f"Output prefix: {args.prefix}")

# ============================================================================
# PART 1: LOAD ALL DATA
# ============================================================================
print("\n[1/6] Loading data...")

# Clustering results with disease info and features
clustering_df = pd.read_csv(f'../../3_outputs/{args.clustering}', index_col=0)
print(f"  - Loaded {len(clustering_df)} variants with clustering assignments")

# Extract feature matrix (columns without 'Metadata_' prefix)
feature_cols = [c for c in clustering_df.columns if not c.startswith('Metadata_')]
feature_matrix = clustering_df[feature_cols].copy()
print(f"  - Extracted feature matrix: {feature_matrix.shape[0]} variants x {feature_matrix.shape[1]} features")

# Compute feature correlation matrix
feature_corr = feature_matrix.corr()
print(f"  - Computed feature correlation matrix: {feature_corr.shape[0]} x {feature_corr.shape[1]}")

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
# Detect column naming convention (AUROC vs misloc_score)
if 'Metadata_AUROC_Mean_GFP_ADJ' in merged_df.columns:
    score_col = 'Metadata_AUROC_Mean_GFP_ADJ'
    hit_col = 'Metadata_Altered_95th_perc_both_batches_GFP_ADJ'
    score_name = 'AUROC'
elif 'Metadata_img_misloc_score' in merged_df.columns:
    score_col = 'Metadata_img_misloc_score'
    hit_col = 'Metadata_img_misloc_hit'
    score_name = 'score'
else:
    score_col = None
    hit_col = None
    score_name = 'score'

results = []
for idx, row in merged_df.iterrows():
    diseases = analyzer.clean_disease_names(row.get('diseases_to_classify', ''))
    classifications = analyzer.classify_disease_modules(diseases)

    primary_module = classifications[0][0] if classifications else 'Unclassified'
    primary_confidence = classifications[0][1] if classifications else 0.0

    result = {
        'gene_allele': row['gene_allele'],
        'gene': row['gene_allele'].split('_')[0],
        'cluster': row['Metadata_cluster'],
        'primary_disease_module': primary_module,
        'primary_module_confidence': primary_confidence,
        'cleaned_diseases': '; '.join(diseases),
        'disease_count': len(diseases)
    }
    if score_col:
        result['score'] = row[score_col]
    if hit_col:
        result['is_hit'] = row[hit_col]

    results.append(result)

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

# Get top 5 features per cluster
top_features_per_cluster = {}
for cluster_id in sorted(classified_df['cluster'].unique()):
    top_feats = cluster_feature_df[
        cluster_feature_df['cluster'] == cluster_id
    ].nlargest(5, 'abs_cohens_d')
    top_features_per_cluster[cluster_id] = top_feats

print(f"  - Identified top discriminative features for {len(top_features_per_cluster)} clusters")

# ============================================================================
# PART 5: GENERATE ALL VISUALIZATIONS (SINGLE PDF)
# ============================================================================
print("\n[5/6] Generating comprehensive visualizations (single PDF)...")

pdf_path = f'../../3_outputs/{args.prefix}comprehensive_disease_biology_analysis.pdf'

with PdfPages(pdf_path) as pdf:

    # PAGE 1: Cluster Overview
    print("  - Page 1: Cluster overview metrics")
    fig = plt.figure(figsize=(11, 8.5))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3)

    # Cluster sizes
    ax1 = fig.add_subplot(gs[0, 0])
    cluster_sizes = classified_df['cluster'].value_counts().sort_index()
    cluster_sizes.plot(kind='bar', ax=ax1, color='steelblue')
    ax1.set_xlabel('Cluster')
    ax1.set_ylabel('Number of Variants')
    ax1.set_title('Cluster Sizes', fontweight='bold')

    # Mean score by cluster
    ax2 = fig.add_subplot(gs[0, 1])
    if 'score' in classified_df.columns:
        cluster_scores = classified_df.groupby('cluster')['score'].mean().sort_index()
        cluster_scores.plot(kind='bar', ax=ax2, color='coral')
        ax2.axhline(y=0.95, color='red', linestyle='--', alpha=0.7)
        ax2.set_xlabel('Cluster')
        ax2.set_ylabel(f'Mean {score_name}')
        ax2.set_title(f'Mean {score_name} by Cluster', fontweight='bold')
        ax2.set_ylim([0.90, 1.0])
    else:
        ax2.text(0.5, 0.5, 'Score not available', ha='center', va='center', transform=ax2.transAxes)
        ax2.set_title('Mean Score by Cluster', fontweight='bold')

    # Disease module distribution with custom colors (grey for Unclassified)
    ax3 = fig.add_subplot(gs[1, :])
    module_cluster = pd.crosstab(
        classified_df['primary_disease_module'],
        classified_df['cluster']
    )
    module_cluster_pct = module_cluster.div(module_cluster.sum(axis=0), axis=1) * 100

    # Create custom color palette with grey for Unclassified
    modules = module_cluster_pct.index.tolist()
    base_colors = plt.cm.tab10.colors
    custom_colors = []
    color_idx = 0
    for module in modules:
        if module == 'Unclassified':
            custom_colors.append('#808080')  # Grey for Unclassified
        else:
            custom_colors.append(base_colors[color_idx % len(base_colors)])
            color_idx += 1

    module_cluster_pct.T.plot(kind='bar', stacked=True, ax=ax3, color=custom_colors)
    ax3.set_xlabel('Cluster')
    ax3.set_ylabel('Percentage (%)')
    ax3.set_title('Disease Module Composition by Cluster', fontweight='bold')
    ax3.legend(title='Disease Module', bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
    ax3.set_ylim([0, 100])

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

        result = {
            'gene': gene,
            'total_variants': len(gene_data),
            'primary_cluster': primary_cluster,
            'consistency_%': consistency,
        }
        if 'score' in gene_data.columns:
            result['mean_score'] = gene_data['score'].mean()
        gene_cluster_analysis.append(result)

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

    # PAGES 3+: Cluster-Specific Feature Analysis (one page per cluster)
    print("  - Pages 3+: Cluster-specific feature analysis")

    # Define metadata columns to check for enrichment (excluding disease_modules, using primary_disease_module instead)
    metadata_cols_to_check = [
        'Metadata_inheritance_pattern',
        'Metadata_pLDDT_Category', 'Metadata_RSA_Category',
        'Metadata_g2p_Signal (UniProt)', 'Metadata_g2p_Transit peptide (UniProt)',
        'Metadata_g2p_Transmembrane (UniProt)', 'Metadata_g2p_Topological domain (UniProt)',
        'Metadata_g2p_Intramembrane (UniProt)', 'Metadata_g2p_Lipidation (UniProt)'
    ]

    for cluster_id in sorted(classified_df['cluster'].unique()):
        cluster_data = classified_df[classified_df['cluster'] == cluster_id]

        # Use larger figure for more content
        fig = plt.figure(figsize=(16, 20))
        gs = fig.add_gridspec(4, 2, hspace=0.25, wspace=0.2,
                              height_ratios=[1.4, 0.7, 0.5, 0.5])

        # =====================================================================
        # TOP SECTION: Cluster Info + Variants List + Enrichments
        # =====================================================================
        ax_info = fig.add_subplot(gs[0, :])
        ax_info.axis('off')

        # Sort variants by gene name, then by position (extract number from variant name)
        def extract_position(variant):
            import re
            match = re.search(r'(\d+)', variant.split('_')[-1] if '_' in variant else variant)
            return int(match.group(1)) if match else 0

        cluster_data_sorted = cluster_data.copy()
        cluster_data_sorted['_gene'] = cluster_data_sorted['gene_allele'].apply(lambda x: x.split('_')[0])
        cluster_data_sorted['_pos'] = cluster_data_sorted['gene_allele'].apply(extract_position)
        cluster_data_sorted = cluster_data_sorted.sort_values(['_gene', '_pos'])

        # Build comprehensive info text
        info_lines = [
            f"CLUSTER {cluster_id} ANALYSIS",
            "=" * 100,
            "",
            f"Size: {len(cluster_data)} variants",
        ]
        if 'score' in cluster_data.columns:
            info_lines.append(f"Mean {score_name}: {cluster_data['score'].mean():.3f}")
        if 'is_hit' in cluster_data.columns:
            hit_count = cluster_data['is_hit'].sum()
            info_lines.append(f"Hit variants: {hit_count} ({hit_count/len(cluster_data)*100:.1f}%)")

        # All variants in cluster (sorted)
        info_lines.extend([
            "",
            "ALL VARIANTS (sorted by gene, position):",
            "-" * 50,
        ])
        variants_list = cluster_data_sorted['gene_allele'].tolist()
        # Format as multiple columns for readability
        n_cols = 4
        for i in range(0, len(variants_list), n_cols):
            row_variants = variants_list[i:i+n_cols]
            info_lines.append("  " + "  |  ".join(f"{v:<22}" for v in row_variants))

        # Disease modules (detailed from Metadata_disease_modules)
        info_lines.extend([
            "",
            "DISEASE MODULES (detailed):",
            "-" * 50,
        ])
        if 'Metadata_disease_modules' in clustering_df.columns:
            # Get disease modules for this cluster's variants
            cluster_variants = cluster_data['gene_allele'].tolist()
            disease_modules_col = clustering_df.loc[
                clustering_df.index.isin(cluster_variants), 'Metadata_disease_modules'
            ]
            # Build mapping: module -> list of variants
            module_to_variants = {}
            for variant in cluster_variants:
                if variant in clustering_df.index:
                    dm = clustering_df.loc[variant, 'Metadata_disease_modules']
                    if pd.notna(dm) and dm != '':
                        for m in str(dm).split(';'):
                            m = m.strip()
                            if m:
                                if m not in module_to_variants:
                                    module_to_variants[m] = []
                                module_to_variants[m].append(variant)
            if module_to_variants:
                # Sort by count descending
                sorted_modules = sorted(module_to_variants.items(), key=lambda x: -len(x[1]))
                for module, vars_list in sorted_modules[:10]:
                    vars_str = ', '.join(vars_list)
                    info_lines.append(f"  {module}: {len(vars_list)} [{vars_str}]")
            else:
                info_lines.append("  No detailed disease modules available")

        # Annotation enrichments (including cleaned_diseases from classified_df)
        info_lines.extend([
            "",
            "ANNOTATION ENRICHMENTS (values appearing ≥2 times):",
            "-" * 50,
        ])
        enrichments_found = False

        # First check cleaned_diseases from classified_df (split by "; " and count)
        if 'cleaned_diseases' in cluster_data.columns:
            # Build mapping: disease -> list of variants
            disease_to_variants = {}
            for _, row in cluster_data.iterrows():
                diseases_str = row.get('cleaned_diseases', '')
                if diseases_str and str(diseases_str).strip():
                    variant = row['gene_allele']
                    for d in str(diseases_str).split(';'):
                        d = d.strip()
                        if d:
                            if d not in disease_to_variants:
                                disease_to_variants[d] = []
                            disease_to_variants[d].append(variant)
            # Only show diseases with >=2 variants
            enriched_diseases = {k: v for k, v in disease_to_variants.items() if len(v) >= 2}
            if enriched_diseases:
                # Sort by count descending
                sorted_diseases = sorted(enriched_diseases.items(), key=lambda x: -len(x[1]))
                for disease, vars_list in sorted_diseases:
                    vars_str = ', '.join(vars_list)
                    info_lines.append(f"  Disease: {disease} ({len(vars_list)}x) [{vars_str}]")
                    enrichments_found = True

        # Then check other metadata columns
        cluster_variants = cluster_data['gene_allele'].tolist()
        for col in metadata_cols_to_check:
            if col in clustering_df.columns:
                # Build mapping: value -> list of variants
                value_to_variants = {}
                for variant in cluster_variants:
                    if variant in clustering_df.index:
                        val = clustering_df.loc[variant, col]
                        if pd.notna(val) and str(val).strip() != '':
                            val_str = str(val).strip()
                            if val_str not in value_to_variants:
                                value_to_variants[val_str] = []
                            value_to_variants[val_str].append(variant)
                # Only show values with >=2 variants
                enriched = {k: v for k, v in value_to_variants.items() if len(v) >= 2}
                if enriched:
                    col_short = col.replace('Metadata_', '').replace('g2p_', '')
                    # Sort by count descending
                    sorted_vals = sorted(enriched.items(), key=lambda x: -len(x[1]))
                    for val, vars_list in sorted_vals:
                        vars_str = ', '.join(vars_list)
                        info_lines.append(f"  {col_short}: {val} ({len(vars_list)}x) [{vars_str}]")
                        enrichments_found = True
        if not enrichments_found:
            info_lines.append("  No enriched annotations found")

        info_text = "\n".join(info_lines)
        # Position text at absolute top-left with no margins, larger font
        ax_info.text(-0.02, 1.02, info_text, transform=ax_info.transAxes,
                    fontsize=10, verticalalignment='top', horizontalalignment='left',
                    fontfamily='monospace', linespacing=1.1)

        # =====================================================================
        # MIDDLE SECTION: Top 5 Discriminative Features
        # =====================================================================
        if cluster_id in top_features_per_cluster:
            top_feats = top_features_per_cluster[cluster_id]

            ax1 = fig.add_subplot(gs[1, :])
            colors_feat = ['#d62728' if d > 0 else '#1f77b4' for d in top_feats['cohens_d']]
            y_positions = range(len(top_feats))
            bars = ax1.barh(y_positions, top_feats['cohens_d'], color=colors_feat, alpha=0.8)
            ax1.set_yticks(y_positions)
            ax1.set_yticklabels([f[:55] for f in top_feats['feature']], fontsize=9)
            ax1.set_xlabel("Cohen's d (effect size)", fontsize=10)
            ax1.set_title(f"Top 5 Discriminative Features for Cluster {cluster_id}",
                         fontweight='bold', fontsize=11)
            ax1.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
            ax1.grid(True, alpha=0.3, axis='x')
            # Add value labels
            for i, (idx, row) in enumerate(top_feats.iterrows()):
                ax1.text(row['cohens_d'], i, f" {row['cohens_d']:.2f}",
                        va='center', fontsize=8,
                        ha='left' if row['cohens_d'] >= 0 else 'right')

            # =====================================================================
            # BOTTOM SECTION: Feature value distributions (top 2 features)
            # =====================================================================
            for i, (idx, feat_row) in enumerate(top_feats.head(2).iterrows()):
                ax = fig.add_subplot(gs[2 + i//2, i%2])

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
                ax.hist(cluster_vals, bins=20, alpha=0.7, label=f'Cluster {cluster_id}', color='#d62728')
                ax.set_xlabel('Feature Value', fontsize=9)
                ax.set_ylabel('Frequency', fontsize=9)
                ax.set_title(f"{feat_name[:50]}", fontsize=9)
                ax.legend(fontsize=8, loc='upper right')

        fig.suptitle(f'Cluster {cluster_id}: Comprehensive Analysis',
                     fontsize=14, fontweight='bold', y=0.995)
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()

print(f"  - Saved comprehensive PDF: {pdf_path}")

# ============================================================================
# PART 6: SAVE RESULTS CSV
# ============================================================================
print("\n[6/6] Saving analysis results...")

# Save main results CSV (renamed from disease_module_classifications)
results_csv_path = f'../../3_outputs/{args.prefix}analysis_results.csv'
classified_df.to_csv(results_csv_path, index=False)
print(f"  - Saved results CSV: {results_csv_path}")

print("\n" + "="*80)
print("ANALYSIS COMPLETE!")
print("="*80)
print(f"\nOutputs:")
print(f"  - PDF: ../../3_outputs/{args.prefix}comprehensive_disease_biology_analysis.pdf")
print(f"  - CSV: {results_csv_path}")
