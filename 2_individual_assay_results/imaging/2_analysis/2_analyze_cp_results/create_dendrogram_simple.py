#!/usr/bin/env python3
"""
Simple script to create a nicely formatted dendrogram visualization for hierarchical clustering.
This script uses the processed data from the notebook to create a comprehensive dendrogram.
"""

import polars as pl
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from sklearn.metrics import silhouette_score
from sklearn.cluster import AgglomerativeClustering
import warnings
warnings.filterwarnings('ignore')

def create_dendrogram():
    """
    Create a comprehensive dendrogram visualization for mislocalized variant clustering.
    """
    
    print("Loading and processing mislocalized variant data...")
    
    # Load the data
    mislocalized_var_cp_features = pl.read_csv("../../3_outputs/mislocalized_variants_key_cp_features.csv")
    
    # Process the data (simplified version for demonstration)
    mislocalized_var_cp_features_processed = mislocalized_var_cp_features.with_columns(
        pl.col("gene_allele").str.split("_").list.get(0).alias("Gene")
    )
    
    # Get numeric columns
    numeric_cols = [col for col in mislocalized_var_cp_features_processed.columns 
                    if not col.startswith("Metadata_") and not col in ["gene_allele", "Gene"]]
    
    # For this demonstration, we'll use a simplified approach
    # In practice, you would implement the full subtraction logic here
    merged_features = mislocalized_var_cp_features_processed.select(["gene_allele"] + numeric_cols).drop_nulls()
    
    # Convert to pandas
    merged_features_pd = merged_features.to_pandas().set_index("gene_allele")
    merged_features_pd = merged_features_pd.dropna()
    
    print(f"Final features dataframe shape: {merged_features_pd.shape}")
    
    # Standardize the features
    scaler = StandardScaler()
    features_scaled = scaler.fit_transform(merged_features_pd)
    
    # Calculate linkage matrix using Ward method
    print("Calculating hierarchical clustering linkage matrix...")
    linkage_matrix = linkage(features_scaled, method='ward', metric='euclidean')
    
    # Determine optimal number of clusters using silhouette analysis
    silhouette_scores = []
    cluster_range = range(2, min(11, len(merged_features_pd) // 2))
    
    for k in cluster_range:
        if k < len(merged_features_pd):
            clustering = AgglomerativeClustering(n_clusters=k, linkage='ward')
            labels = clustering.fit_predict(features_scaled)
            
            if len(np.unique(labels)) > 1:
                sil_score = silhouette_score(features_scaled, labels)
            else:
                sil_score = 0
                
            silhouette_scores.append(sil_score)
            print(f"Clusters: {k}, Silhouette Score: {sil_score:.4f}")
    
    optimal_k = cluster_range[np.argmax(silhouette_scores)]
    print(f"Optimal number of clusters: {optimal_k}")
    
    # Create the dendrogram visualization
    print("Creating dendrogram visualization...")
    
    # Set up the plotting style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Create figure with larger size for better readability
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(24, 16), 
                                  gridspec_kw={'height_ratios': [3, 1]})
    
    # Create main dendrogram
    dendrogram_result = dendrogram(
        linkage_matrix,
        labels=merged_features_pd.index.tolist(),
        orientation='top',
        leaf_rotation=90,
        leaf_font_size=9,
        show_leaf_counts=True,
        color_threshold=linkage_matrix[-(optimal_k-1), 2] if optimal_k > 1 else None,
        above_threshold_color='black',
        ax=ax1
    )
    
    # Get cluster assignments
    cluster_assignments = fcluster(linkage_matrix, optimal_k, criterion='maxclust')
    
    # Create color map for clusters
    colors = plt.cm.Set3(np.linspace(0, 1, optimal_k))
    
    # Color the dendrogram branches based on cluster assignments
    for i, d in enumerate(dendrogram_result['icoord']):
        if i < len(dendrogram_result['leaves']):
            leaf_idx = dendrogram_result['leaves'][i]
            cluster_id = cluster_assignments[leaf_idx] - 1  # Convert to 0-based indexing
            ax1.plot(dendrogram_result['icoord'][i], dendrogram_result['dcoord'][i], 
                    color=colors[cluster_id], linewidth=2.5)
    
    # Customize the main dendrogram
    ax1.set_title(f'Hierarchical Clustering Dendrogram of Mislocalized Variants\n'
                  f'Optimal Clusters: {optimal_k} | Samples: {len(merged_features_pd)} | Features: {features_scaled.shape[1]}', 
                  fontsize=18, fontweight='bold', pad=20)
    ax1.set_xlabel('Samples (Gene_Allele)', fontsize=14, fontweight='bold')
    ax1.set_ylabel('Distance (Ward)', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3, axis='y')
    ax1.tick_params(axis='x', labelsize=8)
    ax1.tick_params(axis='y', labelsize=10)
    
    # Add cluster information below the dendrogram
    cluster_summary = pd.DataFrame({
        'gene_allele': merged_features_pd.index.tolist(),
        'cluster': cluster_assignments
    })
    cluster_summary['gene'] = cluster_summary['gene_allele'].str.split('_').str[0]
    
    # Create cluster distribution bar plot
    cluster_counts = cluster_summary['cluster'].value_counts().sort_index()
    bars = ax2.bar(range(len(cluster_counts)), cluster_counts.values, 
                   color=colors[:len(cluster_counts)], alpha=0.8, edgecolor='black', linewidth=1)
    
    # Add value labels on bars
    for bar, count in zip(bars, cluster_counts.values):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                f'{count}', ha='center', va='bottom', fontweight='bold', fontsize=12)
    
    # Customize the cluster distribution plot
    ax2.set_title(f'Cluster Distribution (Total Variants: {len(cluster_summary)})', 
                  fontsize=16, fontweight='bold', pad=15)
    ax2.set_xlabel('Cluster ID', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Number of Variants', fontsize=12, fontweight='bold')
    ax2.set_xticks(range(len(cluster_counts)))
    ax2.set_xticklabels([f'Cluster {i+1}' for i in range(len(cluster_counts))])
    ax2.grid(True, alpha=0.3, axis='y')
    ax2.tick_params(axis='both', labelsize=10)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save the plot with high resolution
    plt.savefig('mislocalized_variants_comprehensive_dendrogram.png', 
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.show()
    
    # Print detailed cluster information
    print(f"\n=== Detailed Cluster Analysis (k={optimal_k}) ===")
    print(f"Total variants analyzed: {len(cluster_summary)}")
    print(f"Total genes represented: {cluster_summary['gene'].nunique()}")
    
    for cluster_id in sorted(cluster_summary['cluster'].unique()):
        cluster_samples = cluster_summary[cluster_summary['cluster'] == cluster_id]
        genes_in_cluster = cluster_samples['gene'].unique()
        
        print(f"\nCluster {cluster_id} ({len(cluster_samples)} variants, {len(genes_in_cluster)} genes):")
        print(f"  Genes: {', '.join(sorted(genes_in_cluster))}")
        print(f"  Variants:")
        for _, row in cluster_samples.iterrows():
            print(f"    - {row['gene_allele']}")
    
    # Create and save cluster summary
    cluster_summary.to_csv("mislocalized_variants_dendrogram_clusters_detailed.csv", index=False)
    
    # Create gene-cluster distribution matrix
    gene_cluster_dist = cluster_summary.groupby(['gene', 'cluster']).size().unstack(fill_value=0)
    gene_cluster_dist.to_csv("gene_cluster_distribution_matrix.csv")
    
    print(f"\nResults saved:")
    print(f"- Dendrogram: mislocalized_variants_comprehensive_dendrogram.png")
    print(f"- Cluster assignments: mislocalized_variants_dendrogram_clusters_detailed.csv")
    print(f"- Gene-cluster matrix: gene_cluster_distribution_matrix.csv")
    
    return cluster_summary, dendrogram_result

if __name__ == "__main__":
    # Run the dendrogram creation
    cluster_summary, dendrogram_result = create_dendrogram()
    print("\nDendrogram analysis completed successfully!") 