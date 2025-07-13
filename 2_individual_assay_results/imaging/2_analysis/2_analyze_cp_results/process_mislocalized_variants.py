#!/usr/bin/env python3
"""
Process mislocalized variant CP features for robust clustering analysis.

This script processes the mislocalized_var_cp_features polars dataframe to:
1. Extract gene names from gene_allele column
2. Group by Gene and Metadata_Bio_Batch and subtract reference gene features
3. Merge gene_allele from different batches by taking median
4. Perform robust clustering analysis using multiple methods
"""

import polars as pl
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import AgglomerativeClustering, KMeans, DBSCAN
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score, calinski_harabasz_score
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import pdist
import warnings
warnings.filterwarnings('ignore')

def create_dendrogram_visualization(features_scaled, sample_names, optimal_k, save_path="dendrogram_visualization.png"):
    """
    Create a nicely formatted dendrogram visualization for hierarchical clustering.
    
    Parameters:
    -----------
    features_scaled : array-like
        Standardized feature matrix
    sample_names : list
        List of sample names (gene_allele)
    optimal_k : int
        Optimal number of clusters
    save_path : str
        Path to save the dendrogram plot
    """
    
    # Calculate linkage matrix
    linkage_matrix = linkage(features_scaled, method='ward', metric='euclidean')
    
    # Create figure with larger size for better readability
    plt.figure(figsize=(20, 12))
    
    # Create dendrogram
    dendrogram_result = dendrogram(
        linkage_matrix,
        labels=sample_names,
        orientation='top',
        leaf_rotation=90,
        leaf_font_size=8,
        show_leaf_counts=True,
        color_threshold=linkage_matrix[-(optimal_k-1), 2] if optimal_k > 1 else None,
        above_threshold_color='black'
    )
    
    # Get cluster assignments from dendrogram
    cluster_assignments = fcluster(linkage_matrix, optimal_k, criterion='maxclust')
    
    # Create color map for clusters
    colors = plt.cm.Set3(np.linspace(0, 1, optimal_k))
    
    # Color the leaves based on cluster assignments
    for i, d in enumerate(dendrogram_result['icoord']):
        if i < len(dendrogram_result['leaves']):
            leaf_idx = dendrogram_result['leaves'][i]
            cluster_id = cluster_assignments[leaf_idx] - 1  # Convert to 0-based indexing
            plt.plot(dendrogram_result['icoord'][i], dendrogram_result['dcoord'][i], 
                    color=colors[cluster_id], linewidth=2)
    
    # Customize the plot
    plt.title(f'Hierarchical Clustering Dendrogram of Mislocalized Variants\n(Optimal Clusters: {optimal_k})', 
              fontsize=16, fontweight='bold', pad=20)
    plt.xlabel('Samples (Gene_Allele)', fontsize=12, fontweight='bold')
    plt.ylabel('Distance (Ward)', fontsize=12, fontweight='bold')
    
    # Add grid for better readability
    plt.grid(True, alpha=0.3, axis='y')
    
    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()
    
    # Create cluster summary
    cluster_summary = pd.DataFrame({
        'gene_allele': sample_names,
        'cluster': cluster_assignments
    })
    
    # Add gene information
    cluster_summary['gene'] = cluster_summary['gene_allele'].str.split('_').str[0]
    
    print(f"\n=== Dendrogram Cluster Summary (k={optimal_k}) ===")
    for cluster_id in sorted(cluster_summary['cluster'].unique()):
        cluster_samples = cluster_summary[cluster_summary['cluster'] == cluster_id]
        print(f"\nCluster {cluster_id} ({len(cluster_samples)} variants):")
        for _, row in cluster_samples.iterrows():
            print(f"  - {row['gene_allele']} ({row['gene']})")
    
    return cluster_summary, dendrogram_result

def process_mislocalized_variants():
    """
    Main function to process mislocalized variant CP features.
    """
    
    # Load the data
    print("Loading mislocalized variant CP features...")
    mislocalized_var_cp_features = pl.read_csv("../../3_outputs/mislocalized_variants_key_cp_features.csv")
    print(f"Loaded dataframe with shape: {mislocalized_var_cp_features.shape}")
    
    # Step 1: Get the gene from "gene_allele" column by split by "_" and get the first element
    print("\nStep 1: Extracting gene names from gene_allele column...")
    mislocalized_var_cp_features_processed = mislocalized_var_cp_features.with_columns(
        pl.col("gene_allele").str.split("_").list.get(0).alias("Gene")
    )
    
    # Step 2: Group the rows by Gene and Metadata_Bio_Batch and subtract all the numeric features 
    # from the row with Gene == reference gene (with no "_aaXXXaa" variant) in all rows of variants
    
    # First, let's identify the reference genes (those without "_" in gene_allele)
    reference_genes = mislocalized_var_cp_features_processed.filter(
        ~pl.col("gene_allele").str.contains("_")
    ).select("Gene").unique().to_series().to_list()
    
    print(f"Reference genes found: {reference_genes}")
    
    # Get numeric columns (excluding metadata columns)
    numeric_cols = [col for col in mislocalized_var_cp_features_processed.columns 
                    if not col.startswith("Metadata_") and not col in ["gene_allele", "Gene"]]
    
    print(f"Number of numeric features: {len(numeric_cols)}")
    
    # Create a dictionary to store reference gene profiles for each batch
    reference_profiles = {}
    
    # For each gene and batch combination, get the reference profile
    for gene in reference_genes:
        for batch in mislocalized_var_cp_features_processed["Metadata_Bio_Batch"].unique():
            ref_profile = mislocalized_var_cp_features_processed.filter(
                (pl.col("Gene") == gene) & 
                (pl.col("Metadata_Bio_Batch") == batch) &
                (~pl.col("gene_allele").str.contains("_"))
            )
            
            if ref_profile.height > 0:
                ref_values = ref_profile.select(numeric_cols).row(0)
                reference_profiles[(gene, batch)] = ref_values
    
    print(f"Reference profiles created for {len(reference_profiles)} gene-batch combinations")
    
    # Step 3: Subtract reference gene features from variant features
    print("\nStep 2: Subtracting reference gene features from variant features...")
    variant_profiles = mislocalized_var_cp_features_processed.filter(
        pl.col("gene_allele").str.contains("_")
    )
    
    # Create a new dataframe with subtracted features
    subtracted_features = []
    
    for row in variant_profiles.iter_rows(named=True):
        gene = row["Gene"]
        batch = row["Metadata_Bio_Batch"]
        
        if (gene, batch) in reference_profiles:
            ref_values = reference_profiles[(gene, batch)]
            
            # Subtract reference values from variant values
            subtracted_row = {"gene_allele": row["gene_allele"], "Gene": gene, "Metadata_Bio_Batch": batch}
            
            for i, col in enumerate(numeric_cols):
                try:
                    variant_val = row[col]
                    ref_val = ref_values[i]
                    if variant_val is not None and ref_val is not None:
                        subtracted_row[col] = variant_val - ref_val
                    else:
                        subtracted_row[col] = None
                except (IndexError, TypeError):
                    subtracted_row[col] = None
            
            subtracted_features.append(subtracted_row)
    
    subtracted_df = pl.DataFrame(subtracted_features)
    print(f"Created subtracted features dataframe with shape: {subtracted_df.shape}")
    
    # Step 4: Merge the gene_allele from different Metadata_Bio_Batch by taking the median of the features
    print("\nStep 3: Merging gene_allele from different batches by taking median...")
    merged_features = subtracted_df.group_by("gene_allele").agg(
        pl.col(col).median().alias(col) for col in numeric_cols
    ).drop_nulls()
    
    print(f"Merged features dataframe shape: {merged_features.shape}")
    
    # Convert to pandas for clustering
    merged_features_pd = merged_features.to_pandas().set_index("gene_allele")
    
    # Remove any remaining NaN values
    merged_features_pd = merged_features_pd.dropna()
    
    print(f"Final features dataframe shape after removing NaN: {merged_features_pd.shape}")
    
    # Step 5: Perform robust clustering methods using scipy and scikit-learn
    print("\nStep 4: Performing robust clustering analysis...")
    
    # Standardize the features
    scaler = StandardScaler()
    features_scaled = scaler.fit_transform(merged_features_pd)
    
    # Method 1: Hierarchical Clustering (Agglomerative)
    print("\n=== Hierarchical Clustering Analysis ===")
    
    # Try different numbers of clusters
    silhouette_scores_hier = []
    calinski_scores_hier = []
    cluster_range = range(2, min(11, len(merged_features_pd) // 2))
    
    for k in cluster_range:
        if k < len(merged_features_pd):
            clustering = AgglomerativeClustering(n_clusters=k, linkage='ward')
            labels = clustering.fit_predict(features_scaled)
            
            if len(np.unique(labels)) > 1:
                sil_score = silhouette_score(features_scaled, labels)
                cal_score = calinski_harabasz_score(features_scaled, labels)
            else:
                sil_score = 0
                cal_score = 0
                
            silhouette_scores_hier.append(sil_score)
            calinski_scores_hier.append(cal_score)
            print(f"Clusters: {k}, Silhouette Score: {sil_score:.4f}, Calinski-Harabasz Score: {cal_score:.2f}")
    
    # Find optimal number of clusters
    optimal_k_hier = cluster_range[np.argmax(silhouette_scores_hier)]
    print(f"Optimal number of clusters (Hierarchical): {optimal_k_hier}")
    
    # Method 2: K-Means Clustering
    print("\n=== K-Means Clustering Analysis ===")
    
    silhouette_scores_kmeans = []
    calinski_scores_kmeans = []
    
    for k in cluster_range:
        if k < len(merged_features_pd):
            kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
            labels = kmeans.fit_predict(features_scaled)
            
            if len(np.unique(labels)) > 1:
                sil_score = silhouette_score(features_scaled, labels)
                cal_score = calinski_harabasz_score(features_scaled, labels)
            else:
                sil_score = 0
                cal_score = 0
                
            silhouette_scores_kmeans.append(sil_score)
            calinski_scores_kmeans.append(cal_score)
            print(f"Clusters: {k}, Silhouette Score: {sil_score:.4f}, Calinski-Harabasz Score: {cal_score:.2f}")
    
    optimal_k_kmeans = cluster_range[np.argmax(silhouette_scores_kmeans)]
    print(f"Optimal number of clusters (K-Means): {optimal_k_kmeans}")
    
    # Method 3: Gaussian Mixture Model
    print("\n=== Gaussian Mixture Model Analysis ===")
    
    silhouette_scores_gmm = []
    calinski_scores_gmm = []
    
    for k in cluster_range:
        if k < len(merged_features_pd):
            gmm = GaussianMixture(n_components=k, random_state=42)
            labels = gmm.fit_predict(features_scaled)
            
            if len(np.unique(labels)) > 1:
                sil_score = silhouette_score(features_scaled, labels)
                cal_score = calinski_harabasz_score(features_scaled, labels)
            else:
                sil_score = 0
                cal_score = 0
                
            silhouette_scores_gmm.append(sil_score)
            calinski_scores_gmm.append(cal_score)
            print(f"Clusters: {k}, Silhouette Score: {sil_score:.4f}, Calinski-Harabasz Score: {cal_score:.2f}")
    
    optimal_k_gmm = cluster_range[np.argmax(silhouette_scores_gmm)]
    print(f"Optimal number of clusters (GMM): {optimal_k_gmm}")
    
    # Method 4: DBSCAN (density-based)
    print("\n=== DBSCAN Clustering Analysis ===")
    
    # Try different eps values
    eps_values = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    best_dbscan_score = -1
    best_dbscan_params = None
    
    for eps in eps_values:
        dbscan = DBSCAN(eps=eps, min_samples=3)
        labels = dbscan.fit_predict(features_scaled)
        
        if len(np.unique(labels)) > 1 and -1 not in labels:  # No noise points
            sil_score = silhouette_score(features_scaled, labels)
            cal_score = calinski_harabasz_score(features_scaled, labels)
            
            if sil_score > best_dbscan_score:
                best_dbscan_score = sil_score
                best_dbscan_params = (eps, len(np.unique(labels)))
            
            print(f"DBSCAN eps={eps}, clusters={len(np.unique(labels))}, Silhouette Score: {sil_score:.4f}")
    
    if best_dbscan_params:
        print(f"Best DBSCAN: eps={best_dbscan_params[0]}, clusters={best_dbscan_params[1]}")
    
    # Visualize clustering results
    print("\nGenerating visualization plots...")
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Plot silhouette scores comparison
    axes[0, 0].plot(list(cluster_range), silhouette_scores_hier, 'o-', label='Hierarchical', linewidth=2)
    axes[0, 0].plot(list(cluster_range), silhouette_scores_kmeans, 's-', label='K-Means', linewidth=2)
    axes[0, 0].plot(list(cluster_range), silhouette_scores_gmm, '^-', label='GMM', linewidth=2)
    axes[0, 0].set_xlabel('Number of Clusters')
    axes[0, 0].set_ylabel('Silhouette Score')
    axes[0, 0].set_title('Silhouette Score Comparison')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # Plot Calinski-Harabasz scores comparison
    axes[0, 1].plot(list(cluster_range), calinski_scores_hier, 'o-', label='Hierarchical', linewidth=2)
    axes[0, 1].plot(list(cluster_range), calinski_scores_kmeans, 's-', label='K-Means', linewidth=2)
    axes[0, 1].plot(list(cluster_range), calinski_scores_gmm, '^-', label='GMM', linewidth=2)
    axes[0, 1].set_xlabel('Number of Clusters')
    axes[0, 1].set_ylabel('Calinski-Harabasz Score')
    axes[0, 1].set_title('Calinski-Harabasz Score Comparison')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    # Perform PCA for visualization
    pca = PCA(n_components=2)
    features_pca = pca.fit_transform(features_scaled)
    
    # Apply best clustering method (using hierarchical as default)
    best_clustering = AgglomerativeClustering(n_clusters=optimal_k_hier, linkage='ward')
    cluster_labels = best_clustering.fit_predict(features_scaled)
    
    # Plot PCA with clusters
    scatter = axes[1, 0].scatter(features_pca[:, 0], features_pca[:, 1], c=cluster_labels, 
                               cmap='viridis', alpha=0.7, s=60)
    axes[1, 0].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)')
    axes[1, 0].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)')
    axes[1, 0].set_title(f'Hierarchical Clustering (k={optimal_k_hier}) - PCA Projection')
    plt.colorbar(scatter, ax=axes[1, 0])
    
    # Plot cluster distribution
    cluster_counts = np.bincount(cluster_labels)
    axes[1, 1].bar(range(len(cluster_counts)), cluster_counts, alpha=0.7, color='skyblue')
    axes[1, 1].set_xlabel('Cluster')
    axes[1, 1].set_ylabel('Number of Variants')
    axes[1, 1].set_title(f'Cluster Distribution (k={optimal_k_hier})')
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('mislocalized_variants_clustering_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Create final clustering results dataframe
    clustering_results = merged_features_pd.copy()
    clustering_results['Cluster'] = cluster_labels
    clustering_results['Gene'] = clustering_results.index.str.split('_').str[0]
    
    # Show cluster assignments
    print(f"\n=== Final Clustering Results (k={optimal_k_hier}) ===")
    cluster_summary = clustering_results.groupby('Cluster').agg({
        'Gene': ['count', 'unique']
    }).round(2)
    print(cluster_summary)
    
    # Show variants in each cluster
    for cluster_id in sorted(clustering_results['Cluster'].unique()):
        cluster_variants = clustering_results[clustering_results['Cluster'] == cluster_id]
        print(f"\nCluster {cluster_id} ({len(cluster_variants)} variants):")
        for variant in cluster_variants.index:
            print(f"  - {variant}")
    
    # Create dendrogram visualization
    print("\nGenerating dendrogram visualization...")
    sample_names = merged_features_pd.index.tolist()
    cluster_summary_df, dendrogram_result = create_dendrogram_visualization(
        features_scaled, sample_names, optimal_k_hier, 
        save_path="mislocalized_variants_dendrogram.png"
    )
    
    # Save results
    clustering_results.to_csv("mislocalized_variants_clustering_results.csv")
    cluster_summary_df.to_csv("mislocalized_variants_dendrogram_clusters.csv", index=False)
    print(f"\nClustering results saved to 'mislocalized_variants_clustering_results.csv'")
    print(f"Dendrogram cluster assignments saved to 'mislocalized_variants_dendrogram_clusters.csv'")
    
    # Additional analysis: Show gene distribution across clusters
    print(f"\n=== Gene Distribution Across Clusters ===")
    gene_cluster_dist = clustering_results.groupby(['Gene', 'Cluster']).size().unstack(fill_value=0)
    print(gene_cluster_dist)
    
    # Save gene distribution
    gene_cluster_dist.to_csv("gene_cluster_distribution.csv")
    print("Gene cluster distribution saved to 'gene_cluster_distribution.csv'")
    
    return clustering_results, merged_features_pd, cluster_summary_df

if __name__ == "__main__":
    # Run the analysis
    clustering_results, features_df, dendrogram_clusters = process_mislocalized_variants()
    print("\nAnalysis completed successfully!") 