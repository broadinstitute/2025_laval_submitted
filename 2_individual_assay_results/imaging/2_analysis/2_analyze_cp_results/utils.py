import numpy as np
import polars as pl
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster, optimal_leaf_ordering
from scipy.spatial.distance import squareform, pdist
from sklearn.preprocessing import StandardScaler
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
import warnings
from tqdm.contrib.concurrent import thread_map
import scipy.stats as ss
warnings.filterwarnings('ignore')


def plot_top_n_important_feat(df, feat_cols, n=10, ax=None, title=""):
    # Calculate mean and standard deviation for each column
    column_stats = df[feat_cols].agg(['mean', 'std']).T
    # Sort columns by the sum of mean and std in descending order
    sorted_columns = column_stats.sort_values(by=["mean"], ascending=False)
    sorted_columns = sorted_columns.head(n)
    # display(sorted_columns)

    # Plot mean and std for each column
    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 4))
    ax.barh(y=sorted_columns.index, width=sorted_columns["mean"], color='skyblue') ## yerr=sorted_columns["std"], capsize=5, edgecolor='black'
    ax.set_yticks(range(len(sorted_columns.index)))
    ax.set_yticklabels(sorted_columns.index, rotation=0, fontsize=12)
    # ax.set_xlabel("Mean")
    ax.set_title(f"Top {n} Features Importance Mean {title}")
    ax.axvline(x=0.01, color='r', linestyle='--')
    return sorted_columns


def rank_int_array(array: np.ndarray,
                   c: float = 3.0 / 8,
                   stochastic: bool = True,
                   seed: int = 0):
    '''
    Perform rank-based inverse normal transformation in a 1d numpy array. If
    stochastic is True ties are given rank randomly, otherwise ties will share
    the same value.

    Copied directly from: https://github.com/carpenter-singh-lab/2023_Arevalo_NatComm_BatchCorrection/blob/cd9bcf99240880a5c9f9858debf70e94f5b4c0f7/preprocessing/transform.py#L8
    Adapted from: https://github.com/edm1/rank-based-INT/blob/85cb37bb8e0d9e71bb9e8f801fd7369995b8aee7/rank_based_inverse_normal_transformation.py
    '''
    rng = np.random.default_rng(seed=seed)

    if stochastic:
        # Shuffle
        ix = rng.permutation(len(array))
        rev_ix = np.argsort(ix)
        array = array[ix]
        # Get rank, ties are determined by their position(hence why we shuffle)
        rank = ss.rankdata(array, method="ordinal")
        rank = rank[rev_ix]
    else:
        # Get rank, ties are averaged
        rank = ss.rankdata(array, method="average")

    x = (rank - c) / (len(rank) - 2 * c + 1)
    return ss.norm.ppf(x)
    

def inverse_normal_transform(data, feat_cols=None, c: float = 3.0 / 8, stochastic: bool = True, seed: int = 0):
    """
    Apply inverse normal transformation (rank-based normalization) to data using rank_int_array.
    
    Parameters:
    -----------
    data : polars.DataFrame
       Input dataframe
    feat_cols : list, optional
       Specific columns to transform. If None, auto-detect numeric columns
    c : float
       Blom's constant for rank transformation
    stochastic : bool
       Whether to use stochastic tie-breaking
    seed : int
       Random seed for stochastic tie-breaking
    """
    from tqdm.contrib.concurrent import thread_map
    import scipy.stats as ss
    
    # Determine columns to transform
    if feat_cols is not None:
        numeric_cols = feat_cols
    else:
        numeric_cols = [col for col in data.columns if "Meta" not in col]
    
    if not numeric_cols:
        return data  # No columns to transform
    
    # Convert to numpy for processing
    data_np = data.select(numeric_cols).to_numpy()
    
    def to_normal(i):
        """Transform column i using rank_int_array"""
        col_data = data_np[:, i].copy()  # Work on copy to avoid race conditions
        mask = ~np.isnan(col_data)
        
        if np.sum(mask) > 0:  # Only transform if there are non-null values
            # Transform non-null values
            transformed_values = rank_int_array(
               col_data[mask], 
               c=c, 
               stochastic=stochastic, 
               seed=seed + i  # Different seed per column for stochastic mode
            ).astype(np.float32)
            
            # Put transformed values back
            col_data[mask] = transformed_values
            data_np[:, i] = col_data
    
    # Apply transformation to all columns in parallel
    thread_map(to_normal, range(len(numeric_cols)), leave=False)
    
    # Convert back to polars DataFrame
    transformed_df = pl.DataFrame(
       data_np, 
       schema={col: pl.Float32 for col in numeric_cols}
    )
    
    # Combine with non-numeric columns if they exist
    non_numeric_cols = [col for col in data.columns if col not in numeric_cols]
    
    if non_numeric_cols:
        # Preserve original column order
        result_df = data.select(non_numeric_cols).hstack(transformed_df)
        return result_df.select(data.columns)  # Maintain original column order
    else:
        return transformed_df


def robust_cosine_clustering_with_annotations(feat_df_clean_well_pd, 
                                              feature_rows=None,
                                              annotation_rows=None,
                                              continuous_annotations=None,
                                              n_clusters=None, 
                                              method='average', 
                                              standardize=True, 
                                              min_clusters=2, 
                                              max_clusters=20,
                                              optimal_ordering=True):
    """
    Perform robust clustering using cosine similarities with categorical and continuous annotations.
    
    Parameters:
    -----------
    feat_df_clean_well_pd : pandas.DataFrame
        Feature matrix with annotations (features/annotations x samples)
    feature_rows : list, optional
        List of row names to use for clustering. If None, will exclude annotation rows
    annotation_rows : list, optional
        List of categorical annotation row names (e.g., ["clinvar_clnsig_clean", "Altered_95th_perc_GFP"])
    continuous_annotations : list, optional
        List of continuous annotation row names (e.g., ["AUROC"])
    n_clusters : int, optional
        Number of clusters. If None, will be determined automatically
    method : str
        Linkage method for hierarchical clustering
    standardize : bool
        Whether to standardize features before clustering
    optimal_ordering : bool
        Whether to use optimal leaf ordering
    
    Returns:
    --------
    dict with clustering results and validation metrics
    """
    # Default annotation rows if not specified
    if annotation_rows is None:
        annotation_rows = ["clinvar_clnsig_clean", "Altered_95th_perc_GFP"]
    
    if continuous_annotations is None:
        continuous_annotations = ["AUROC_Mean_GFP"]
    
    # Identify feature rows (non-annotation rows)
    if feature_rows is None:
        all_annotation_rows = annotation_rows + continuous_annotations
        feature_rows = [row for row in feat_df_clean_well_pd.index 
                       if row not in all_annotation_rows]
    
    # Extract feature data for clustering
    data = feat_df_clean_well_pd.loc[feature_rows].astype(float)
    
    # Step 1: Data preprocessing
    if standardize:
        scaler = StandardScaler()
        data_scaled = pd.DataFrame(
            scaler.fit_transform(data.T).T, 
            index=data.index, 
            columns=data.columns
        )
    else:
        data_scaled = data.copy()
    
    # Step 2: Compute cosine distance between columns (samples)
    # Transpose to cluster columns (samples) -> shape (n_samples, n_features)
    data_T = data_scaled.T.values
    
    # Compute cosine distance between samples
    cosine_dist = pdist(data_T, metric='cosine')
    
    # Step 3: Perform hierarchical clustering
    linkage_matrix = linkage(cosine_dist, method=method)
    
    # Step 4: Optimal leaf ordering to minimize distance jumps
    if optimal_ordering:
        opt_linkage = optimal_leaf_ordering(linkage_matrix, cosine_dist)
    else:
        opt_linkage = linkage_matrix
    
    # Step 5: Determine optimal number of clusters if not specified
    if n_clusters is None:
        n_clusters = find_optimal_clusters_pdist(
            data_scaled.T, cosine_dist, linkage_matrix, 
            min_clusters, max_clusters
        )
    
    # Step 6: Get cluster labels
    cluster_labels = fcluster(opt_linkage, n_clusters, criterion='maxclust')
    
    # Create cluster labels for samples
    sample_cluster_df = pd.DataFrame({
        'Sample': data.columns,
        'Cluster': cluster_labels
    }).set_index('Sample')
    
    # Step 7: Create results dictionary
    results = {
        'cluster_labels': cluster_labels,
        'n_clusters': n_clusters,
        'linkage_matrix': opt_linkage,
        'linkage_method': method,
        'data_scaled': data_scaled,
        'feature_rows': feature_rows,
        'annotation_rows': annotation_rows,
        'continuous_annotations': continuous_annotations,
        'sample_clusters': sample_cluster_df,
        'original_data': feat_df_clean_well_pd
    }
    return results
    

def find_optimal_clusters_pdist(data, distance_vector, linkage_matrix, min_k, max_k):
    """Find optimal number of clusters using multiple criteria with precomputed distances."""
    from sklearn.metrics import silhouette_score, calinski_harabasz_score
    scores = {
        'silhouette': [],
        'calinski_harabasz': []
    }
    k_range = range(min_k, min(max_k + 1, len(data)))
    for k in k_range:
        labels = fcluster(linkage_matrix, k, criterion='maxclust')
        
        try:
            # Silhouette score using precomputed distances
            distance_matrix = squareform(distance_vector)
            sil_score = silhouette_score(distance_matrix, labels, metric='precomputed')
            scores['silhouette'].append(sil_score)
            
            # Calinski-Harabasz score (higher is better)
            ch_score = calinski_harabasz_score(data, labels)
            scores['calinski_harabasz'].append(ch_score)
            
        except:
            scores['silhouette'].append(-1)
            scores['calinski_harabasz'].append(0)
    
    # Find optimal k based on silhouette score
    if scores['silhouette'] and max(scores['silhouette']) > -1:
        optimal_k = k_range[np.argmax(scores['silhouette'])]
    else:
        optimal_k = min_k
    return optimal_k


def create_color_palettes(feat_df_clean_well_pd, annotation_rows, continuous_annotations, n_clusters):
    """Create color-blind friendly palettes for all annotations and clusters."""
    palettes = {}
    # Color-blind friendly palettes
    categorical_colors = sns.color_palette("Set2", 10)
    cluster_colors = sns.color_palette("tab10", n_clusters)
    
    # Create palettes for categorical annotations
    color_idx = 0
    for ann_row in annotation_rows:
        if ann_row in feat_df_clean_well_pd.index:
            unique_vals = sorted(list(feat_df_clean_well_pd.loc[ann_row].unique()))
            unique_vals = [v for v in unique_vals if pd.notna(v)]  # Remove NaN
            
            n_colors_needed = len(unique_vals)
            if n_colors_needed <= len(categorical_colors):
                row_colors = categorical_colors[color_idx:color_idx + n_colors_needed]
            else:
                # If more colors needed, use a larger palette
                row_colors = sns.color_palette("husl", n_colors_needed)
            
            palettes[ann_row] = dict(zip(unique_vals, row_colors))
            color_idx += n_colors_needed
    
    # Create palette for clusters
    cluster_vals = list(range(1, n_clusters + 1))
    palettes['Cluster'] = dict(zip(cluster_vals, cluster_colors))
    
    # For continuous annotations, we'll use a colormap
    continuous_palettes = {}
    for cont_ann in continuous_annotations:
        if cont_ann in feat_df_clean_well_pd.index:
            continuous_palettes[cont_ann] = 'viridis_r'  # Color-blind friendly continuous colormap
    return palettes, continuous_palettes


def create_column_colors_dataframe(feat_df_clean_well_pd, clustering_results, palettes, continuous_palettes):
    """Create a DataFrame with column colors for all annotations and clusters."""
    
    annotation_rows = clustering_results['annotation_rows']
    continuous_annotations = clustering_results['continuous_annotations']
    sample_clusters = clustering_results['sample_clusters']
    
    col_colors_dict = {}
    
    # Add categorical annotations
    for ann_row in annotation_rows:
        if ann_row in feat_df_clean_well_pd.index:
            col_colors_dict[ann_row] = feat_df_clean_well_pd.loc[ann_row].map(palettes[ann_row])
    
    # Add cluster colors
    cluster_series = pd.Series(clustering_results['cluster_labels'], 
                              index=feat_df_clean_well_pd.columns)
    col_colors_dict['Cluster'] = cluster_series.map(palettes['Cluster'])
    
    # Add continuous annotations (normalized to [0,1] for colormap)
    for cont_ann in continuous_annotations:
        if cont_ann in feat_df_clean_well_pd.index:
            values = feat_df_clean_well_pd.loc[cont_ann].astype(float)
            # Normalize values to [0,1]
            normalized_values = (values - values.min()) / (values.max() - values.min())
            
            # Apply colormap
            cmap = plt.cm.get_cmap(continuous_palettes[cont_ann])
            col_colors_dict[cont_ann] = [cmap(val) for val in normalized_values]
    
    # Create DataFrame
    col_colors_df = pd.DataFrame(col_colors_dict, index=feat_df_clean_well_pd.columns)
    return col_colors_df


def create_comprehensive_clustermap(feat_df_clean_well_pd, clustering_results, 
                                    figsize=(15, 18), cmap='vlag', 
                                    save_path=None, show_legends=True):
    """
    Create a comprehensive clustermap with all annotations and cluster information.
    """
    
    # Get the clustering data
    data = clustering_results['data_scaled']
    linkage_matrix = clustering_results['linkage_matrix']
    
    # Create color palettes
    palettes, continuous_palettes = create_color_palettes(
        feat_df_clean_well_pd, 
        clustering_results['annotation_rows'],
        clustering_results['continuous_annotations'],
        clustering_results['n_clusters']
    )
    
    # Create column colors DataFrame
    col_colors_df = create_column_colors_dataframe(
        feat_df_clean_well_pd, clustering_results, palettes, continuous_palettes
    )
    
    # Create the clustermap
    g = sns.clustermap(
        data,
        col_linkage=linkage_matrix,
        col_colors=col_colors_df,
        cmap=cmap,
        figsize=figsize,
        cbar_kws={'label': 'Scaled Feature Value'},
        linewidths=0.01,
        row_cluster=True,  # Also cluster features
        xticklabels=True,  # Hide sample labels for cleaner look
        yticklabels=True
    )
    
    # Customize the plot
    g.ax_heatmap.set_xlabel('Samples (Clustered)', fontsize=12)
    g.ax_heatmap.set_ylabel('Features', fontsize=12)
    
    # Add legends if requested
    if show_legends:
        legend_y_start = 0.95
        legend_spacing = 0.15
        # Legend for categorical annotations
        for i, (ann_name, palette) in enumerate(palettes.items()):
            if ann_name in clustering_results['annotation_rows'] or ann_name == 'Cluster':
                handles = [mpatches.Patch(color=color, label=label) 
                          for label, color in palette.items()]
                
                legend = plt.legend(handles=handles, title=ann_name,
                                  bbox_to_anchor=(0, legend_y_start - i * legend_spacing), 
                                  bbox_transform=g.ax_heatmap.transAxes,
                                  loc='upper left')
                legend.set_title(ann_name, prop={'weight': 'bold'})
                plt.gca().add_artist(legend)
        
        # Add continuous annotation info in title or as text
        cont_info = []
        for cont_ann in clustering_results['continuous_annotations']:
            if cont_ann in feat_df_clean_well_pd.index:
                values = feat_df_clean_well_pd.loc[cont_ann].astype(float)
                cont_info.append(f"{cont_ann}: {values.min():.2f}-{values.max():.2f}")
        
        if cont_info:
            plt.figtext(0.02, 0.02, "Continuous annotations: " + ", ".join(cont_info), 
                       fontsize=10, style='italic')
    
    plt.suptitle(f'Comprehensive Variant Clustering Analysis\n'
                f'Method: {clustering_results["linkage_method"]}, '
                f'Clusters: {clustering_results["n_clusters"]}, '
                f'Features: {len(clustering_results["feature_rows"])}', 
                y=0.98, fontsize=14, weight='bold')
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return g, palettes, continuous_palettes

    
def analyze_cluster_characteristics(feat_df_clean_well_pd, clustering_results):
    """Analyze and summarize characteristics of each cluster."""
    sample_clusters = clustering_results['sample_clusters']
    annotation_rows = clustering_results['annotation_rows']
    continuous_annotations = clustering_results['continuous_annotations']
    
    cluster_summary = {}
    for cluster_id in range(1, clustering_results['n_clusters'] + 1):
        cluster_samples = sample_clusters[sample_clusters['Cluster'] == cluster_id].index
        cluster_summary[f'Cluster_{cluster_id}'] = {
            'n_samples': len(cluster_samples),
            'samples': list(cluster_samples)
        }
        
        # Categorical annotation summaries
        for ann_row in annotation_rows:
            if ann_row in feat_df_clean_well_pd.index:
                cluster_values = feat_df_clean_well_pd.loc[ann_row, cluster_samples]
                value_counts = cluster_values.value_counts()
                cluster_summary[f'Cluster_{cluster_id}'][ann_row] = value_counts.to_dict()
        
        # Continuous annotation summaries
        for cont_ann in continuous_annotations:
            if cont_ann in feat_df_clean_well_pd.index:
                cluster_values = feat_df_clean_well_pd.loc[cont_ann, cluster_samples].astype(float)
                cluster_summary[f'Cluster_{cluster_id}'][cont_ann] = {
                    'mean': cluster_values.mean(),
                    'std': cluster_values.std(),
                    'min': cluster_values.min(),
                    'max': cluster_values.max()
                }
    return cluster_summary

    
# Example usage function
def run_comprehensive_analysis(feat_df_clean_well_pd, 
                               feature_rows=None,
                               annotation_rows=None,
                               continuous_annotations=None,
                               n_clusters=None,
                               method='average',
                               save_path=None):
    """
    Run the complete analysis pipeline.
    """
    
    print("Starting comprehensive clustering analysis...")
    
    # Set defaults if not provided
    if annotation_rows is None:
        annotation_rows = ["clinvar_clnsig_clean", "Altered_95th_perc_GFP"]
    
    if continuous_annotations is None:
        continuous_annotations = ["AUROC"]
    
    # Perform clustering
    print("Performing robust cosine similarity clustering...")
    clustering_results = robust_cosine_clustering_with_annotations(
        feat_df_clean_well_pd,
        feature_rows=feature_rows,
        annotation_rows=annotation_rows,
        continuous_annotations=continuous_annotations,
        n_clusters=n_clusters,
        method=method
    )

    print(f"Optimal number of clusters: {clustering_results['n_clusters']}")
    print(f"Using {len(clustering_results['feature_rows'])} features for clustering")
    
    # Create comprehensive clustermap
    print("Creating comprehensive clustermap...")
    g, palettes, continuous_palettes = create_comprehensive_clustermap(
        feat_df_clean_well_pd, 
        clustering_results,
        save_path=save_path
    )

    # Analyze cluster characteristics
    print("Analyzing cluster characteristics...")
    cluster_summary = analyze_cluster_characteristics(feat_df_clean_well_pd, clustering_results)
    
    # Print summary
    for cluster_name, summary in cluster_summary.items():
        print(f"\n{cluster_name}: {summary['n_samples']} samples")
        for ann_row in annotation_rows:
            if ann_row in summary:
                print(f"  {ann_row}: {summary[ann_row]}")
        for cont_ann in continuous_annotations:
            if cont_ann in summary:
                print(f"  {cont_ann}: mean={summary[cont_ann]['mean']:.3f} ± {summary[cont_ann]['std']:.3f}")
    
    return {
        'clustering_results': clustering_results,
        'clustermap': g,
        'palettes': palettes,
        'continuous_palettes': continuous_palettes,
        'cluster_summary': cluster_summary
    }




import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler
from scipy import stats
from scipy.optimize import minimize_scalar
import warnings
warnings.filterwarnings('ignore')

# Set style
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

class DistributionAnalyzer:
    """
    Class to analyze and compare distributions using Gaussian Mixture Models
    """
    
    def __init__(self, data_dict):
        """
        Initialize with data dictionary where keys are category names and values are arrays
        """
        self.data_dict = data_dict
        self.fitted_models = {}
        self.model_info = {}
        
    def fit_gmm(self, category, n_components_range=(1, 5), criterion='bic'):
        """
        Fit Gaussian Mixture Model with optimal number of components
        """
        data = self.data_dict[category].reshape(-1, 1)
        
        # Test different numbers of components
        scores = []
        models = []
        
        for n_comp in range(n_components_range[0], n_components_range[1] + 1):
            gmm = GaussianMixture(n_components=n_comp, random_state=42)
            gmm.fit(data)
            
            if criterion == 'bic':
                score = gmm.bic(data)
            else:  # aic
                score = gmm.aic(data)
                
            scores.append(score)
            models.append(gmm)
        
        # Select best model (lowest BIC/AIC)
        best_idx = np.argmin(scores)
        best_model = models[best_idx]
        best_n_comp = best_idx + n_components_range[0]
        
        self.fitted_models[category] = best_model
        self.model_info[category] = {
            'n_components': best_n_comp,
            'score': scores[best_idx],
            'means': best_model.means_.flatten(),
            'covariances': best_model.covariances_.flatten(),
            'weights': best_model.weights_
        }
        
        return best_model
    
    def fit_all_categories(self):
        """
        Fit GMM for all categories
        """
        for category in self.data_dict.keys():
            self.fit_gmm(category)
    
    def plot_distributions_with_gmm(self, figsize=(15, 10)):
        """
        Plot original distributions with fitted GMM overlays
        """
        n_categories = len(self.data_dict)
        fig, axes = plt.subplots(2, 2, figsize=figsize)
        axes = axes.flatten()
        
        colors = ['red', 'blue', 'green', 'purple']
        
        for i, (category, data) in enumerate(self.data_dict.items()):
            ax = axes[i]
            
            # Plot histogram
            ax.hist(data, bins=50, density=True, alpha=0.6, 
                   color=colors[i], label=f'{category} (Data)')
            
            # Plot fitted GMM
            if category in self.fitted_models:
                x_range = np.linspace(data.min(), data.max(), 1000)
                gmm = self.fitted_models[category]
                
                # Overall GMM
                log_prob = gmm.score_samples(x_range.reshape(-1, 1))
                prob = np.exp(log_prob)
                ax.plot(x_range, prob, color='black', linewidth=2, 
                       label=f'GMM ({self.model_info[category]["n_components"]} components)')
                
                # Individual components
                for j in range(gmm.n_components):
                    component_prob = gmm.weights_[j] * stats.norm.pdf(
                        x_range, 
                        gmm.means_[j, 0], 
                        np.sqrt(gmm.covariances_[j, 0, 0])
                    )
                    ax.plot(x_range, component_prob, '--', alpha=0.7,
                           label=f'Component {j+1}')
            
            ax.set_title(f'{category}')
            ax.set_xlabel('zscore_log2')
            ax.set_ylabel('Density')
            ax.legend()
            ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()
    
    def compare_with_null(self, null_category='2_Benign'):
        """
        Compare all distributions with the null distribution (2_Benign)
        """
        if null_category not in self.fitted_models:
            raise ValueError(f"Null category {null_category} not found in fitted models")
        
        null_model = self.fitted_models[null_category]
        results = {}
        
        for category, model in self.fitted_models.items():
            if category == null_category:
                continue
                
            # Calculate KL divergence
            kl_div = self._calculate_kl_divergence(null_model, model)
            
            # Identify extra peaks
            extra_peaks = self._identify_extra_peaks(null_model, model, category)
            
            results[category] = {
                'kl_divergence': kl_div,
                'extra_peaks': extra_peaks,
                'n_components': self.model_info[category]['n_components'],
                'null_n_components': self.model_info[null_category]['n_components']
            }
        
        return results
    
    def _calculate_kl_divergence(self, model1, model2, n_samples=10000):
        """
        Calculate KL divergence between two GMM models
        """
        # Sample points for comparison
        data_range = np.linspace(-25, 5, n_samples).reshape(-1, 1)
        
        # Get log probabilities
        log_prob1 = model1.score_samples(data_range)
        log_prob2 = model2.score_samples(data_range)
        
        # Convert to probabilities and normalize
        prob1 = np.exp(log_prob1)
        prob2 = np.exp(log_prob2)
        
        prob1 = prob1 / np.sum(prob1)
        prob2 = prob2 / np.sum(prob2)
        
        # Calculate KL divergence (avoid log(0))
        epsilon = 1e-10
        kl_div = np.sum(prob2 * np.log((prob2 + epsilon) / (prob1 + epsilon)))
        
        return kl_div
    
    def _identify_extra_peaks(self, null_model, target_model, category):
        """
        Identify extra peaks in target distribution compared to null
        """
        null_means = null_model.means_.flatten()
        target_means = target_model.means_.flatten()
        target_weights = target_model.weights_
        
        # Find target peaks that are far from null peaks
        extra_peaks = []
        threshold_distance = 2.0  # Adjust based on your data scale
        
        for i, (mean, weight) in enumerate(zip(target_means, target_weights)):
            # Check if this peak is far from all null peaks
            min_distance = np.min(np.abs(mean - null_means))
            
            if min_distance > threshold_distance and weight > 0.1:  # Significant weight
                extra_peaks.append({
                    'mean': mean,
                    'weight': weight,
                    'component': i,
                    'distance_from_null': min_distance
                })
        
        return extra_peaks
    
    def find_optimal_threshold(self, null_category='2_Benign', target_category='1_Pathogenic'):
        """
        Find optimal threshold to separate null from target distribution
        """
        null_data = self.data_dict[null_category]
        target_data = self.data_dict[target_category]
        
        # Define range for threshold search
        min_val = min(null_data.min(), target_data.min())
        max_val = max(null_data.max(), target_data.max())
        
        def objective(threshold):
            # Calculate misclassification rate
            null_misclass = np.sum(null_data < threshold) / len(null_data)
            target_misclass = np.sum(target_data >= threshold) / len(target_data)
            return null_misclass + target_misclass
        
        # Find optimal threshold
        result = minimize_scalar(objective, bounds=(min_val, max_val), method='bounded')
        optimal_threshold = result.x
        
        # Calculate performance metrics
        null_correct = np.sum(null_data >= optimal_threshold) / len(null_data)
        target_correct = np.sum(target_data < optimal_threshold) / len(target_data)
        
        return {
            'threshold': optimal_threshold,
            'null_accuracy': null_correct,
            'target_accuracy': target_correct,
            'overall_accuracy': (null_correct + target_correct) / 2
        }
    
    def plot_comparison_with_threshold(self, null_category='2_Benign', 
                                     target_category='1_Pathogenic', figsize=(12, 8)):
        """
        Plot comparison between null and target with optimal threshold
        """
        threshold_info = self.find_optimal_threshold(null_category, target_category)
        threshold = threshold_info['threshold']
        
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize)
        
        # Plot 1: Overlapping distributions with threshold
        null_data = self.data_dict[null_category]
        target_data = self.data_dict[target_category]
        
        ax1.hist(null_data, bins=50, density=True, alpha=0.6, 
                color='blue', label=f'{null_category} (NULL)')
        ax1.hist(target_data, bins=50, density=True, alpha=0.6, 
                color='red', label=f'{target_category}')
        
        # Add fitted GMM curves
        x_range = np.linspace(min(null_data.min(), target_data.min()), 
                             max(null_data.max(), target_data.max()), 1000)
        
        null_prob = np.exp(self.fitted_models[null_category].score_samples(x_range.reshape(-1, 1)))
        target_prob = np.exp(self.fitted_models[target_category].score_samples(x_range.reshape(-1, 1)))
        
        ax1.plot(x_range, null_prob, color='blue', linewidth=2, linestyle='--')
        ax1.plot(x_range, target_prob, color='red', linewidth=2, linestyle='--')
        
        # Add threshold line
        ax1.axvline(threshold, color='black', linestyle='-', linewidth=2, 
                   label=f'Optimal Threshold: {threshold:.2f}')
        
        ax1.set_xlabel('zscore_log2')
        ax1.set_ylabel('Density')
        ax1.set_title(f'Distribution Comparison: {null_category} vs {target_category}')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Plot 2: Difference in distributions
        diff = target_prob - null_prob
        ax2.plot(x_range, diff, color='purple', linewidth=2)
        ax2.axhline(0, color='black', linestyle='--', alpha=0.5)
        ax2.axvline(threshold, color='black', linestyle='-', linewidth=2, alpha=0.7)
        ax2.fill_between(x_range, diff, 0, where=(diff > 0), alpha=0.3, color='red')
        ax2.fill_between(x_range, diff, 0, where=(diff < 0), alpha=0.3, color='blue')
        
        ax2.set_xlabel('zscore_log2')
        ax2.set_ylabel('Density Difference')
        ax2.set_title(f'Difference: {target_category} - {null_category}')
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()
        
        # # Print results
        # print(f\"\\nThreshold Analysis Results:\")
        # print(f\"Optimal Threshold: {threshold:.3f}\")
        # print(f\"{null_category} Accuracy: {threshold_info['null_accuracy']:.3f}\")
        # print(f\"{target_category} Accuracy: {threshold_info['target_accuracy']:.3f}\")
        # print(f\"Overall Accuracy: {threshold_info['overall_accuracy']:.3f}\"")
    
#     def generate_summary_report(self):
#         \"\"\"
#         Generate a comprehensive summary report
#         \"\"\"
#         print(\"=== GAUSSIAN MIXTURE MODEL ANALYSIS REPORT ===\")
#         print(\"\\n1. MODEL FITTING RESULTS:\")
        
#         for category, info in self.model_info.items():
#             print(f\"\\n{category}:\"")
#             print(f\"  - Optimal components: {info['n_components']}\")
#             print(f\"  - BIC score: {info['score']:.2f}\")
#             print(f\"  - Component means: {info['means']}\")
#             print(f\"  - Component weights: {info['weights']}\")
        
#         print(\"\\n2. COMPARISON WITH NULL (2_Benign):\")
#         comparison_results = self.compare_with_null()
        
#         for category, results in comparison_results.items():
#             print(f\"\\n{category}:\")
#             print(f\"  - KL Divergence from null: {results['kl_divergence']:.4f}\")
#             print(f\"  - Number of extra peaks: {len(results['extra_peaks'])}\")
            
#             if results['extra_peaks']:
#                 for i, peak in enumerate(results['extra_peaks']):
#                     print(f\"    Peak {i+1}: mean={peak['mean']:.2f}, weight={peak['weight']:.3f}\")

# # Example usage (replace with your actual data):
# # Assuming you have your data loaded in a dictionary format
# # data_dict = {
# #     '1_Pathogenic': pathogenic_zscore_data,
# #     '2_Benign': benign_zscore_data, 
# #     '3_Conflicting': conflicting_zscore_data,
# #     '4_VUS': vus_zscore_data
# # }

# # If you need to load data from your existing variables, do something like:
# # data_dict = {
# #     '1_Pathogenic': df[df['category'] == '1_Pathogenic']['zscore_log2'].values,
# #     '2_Benign': df[df['category'] == '2_Benign']['zscore_log2'].values,
# #     '3_Conflicting': df[df['category'] == '3_Conflicting']['zscore_log2'].values,
# #     '4_VUS': df[df['category'] == '4_VUS']['zscore_log2'].values
# # }

# # Example with synthetic data (replace with your actual data):
# np.random.seed(42)
# data_dict = {
#     '1_Pathogenic': np.concatenate([
#         np.random.normal(-8, 2, 1000),  # Main peak
#         np.random.normal(-15, 1.5, 300)  # Extra peak
#     ]),
#     '2_Benign': np.random.normal(-5, 3, 1500),  # NULL distribution
#     '3_Conflicting': np.concatenate([
#         np.random.normal(-5, 3, 800),   # Similar to benign
#         np.random.normal(-12, 2, 400)   # Extra component
#     ]),
#     '4_VUS': np.concatenate([
#         np.random.normal(-6, 3.5, 1000),  # Main peak
#         np.random.normal(-2, 1, 200)      # Right tail
#     ])
# }

# # Initialize analyzer
# analyzer = DistributionAnalyzer(data_dict)

# # Fit GMM models for all categories
# print(\"Fitting Gaussian Mixture Models...\")
# analyzer.fit_all_categories()

# # Plot distributions with GMM fits
# print(\"\\nPlotting distributions with GMM overlays...\")
# analyzer.plot_distributions_with_gmm()

# # Generate comprehensive report
# analyzer.generate_summary_report()

# # Compare specific categories with threshold analysis
# print(\"\\n=== THRESHOLD ANALYSIS ===\")
# categories_to_compare = ['1_Pathogenic', '3_Conflicting', '4_VUS']

# for category in categories_to_compare:
#     print(f\"\\nAnalyzing {category} vs 2_Benign:\")
#     analyzer.plot_comparison_with_threshold('2_Benign', category)


# import numpy as np
# import pandas as pd
# import seaborn as sns
# import matplotlib.pyplot as plt
# from sklearn.metrics.pairwise import cosine_similarity
# from sklearn.cluster import AgglomerativeClustering
# from scipy.cluster.hierarchy import linkage, dendrogram, fcluster, optimal_leaf_ordering
# from scipy.spatial.distance import squareform, pdist
# from sklearn.preprocessing import StandardScaler
# from scipy.stats import rankdata
# import matplotlib.patches as mpatches
# from matplotlib.colors import ListedColormap
# import umap
# import warnings
# warnings.filterwarnings('ignore')


# def inverse_normal_transform(data):
#     """
#     Apply inverse normal transformation (rank-based normalization) to data.
#     """
#     from scipy.stats import norm
    
#     def transform_column(col):
#         valid_mask = ~pd.isna(col)
#         if valid_mask.sum() == 0:
#             return col
        
#         valid_data = col[valid_mask]
#         ranks = rankdata(valid_data, method='average')
#         percentiles = (ranks - 0.5) / len(ranks)
#         transformed_valid = norm.ppf(percentiles)
        
#         result = col.copy()
#         result[valid_mask] = transformed_valid
#         return result
    
#     if isinstance(data, pd.DataFrame):
#         return data.apply(transform_column, axis=0)
#     else:
#         result = np.zeros_like(data)
#         for i in range(data.shape[1]):
#             col = pd.Series(data[:, i])
#             result[:, i] = transform_column(col).values
#         return result


# def robust_cosine_clustering_samples_features(df_input, 
#                                             annotation_columns=None,
#                                             feature_columns=None,
#                                             n_clusters=None, 
#                                             method='average', 
#                                             use_inverse_normal=True,
#                                             min_clusters=2, 
#                                             max_clusters=20,
#                                             optimal_ordering=True):
#     """
#     Perform robust clustering using cosine similarities with samples x features format.
#     """
#     # Default annotation columns if not specified
#     if annotation_columns is None:
#         annotation_columns = ["clinvar_clnsig_clean", "Altered_95th_perc_GFP", "AUROC_Mean"]
    
#     # Identify feature columns (non-annotation columns)
#     if feature_columns is None:
#         feature_columns = [col for col in df_input.columns 
#                           if col not in annotation_columns]
    
#     # Extract feature data for clustering
#     feature_data = df_input[feature_columns].copy()
    
#     # Remove columns with all NaN or constant values
#     valid_feature_cols = []
#     for col in feature_columns:
#         col_data = feature_data[col].dropna()
#         if len(col_data) > 1 and col_data.nunique() > 1:
#             valid_feature_cols.append(col)
    
#     feature_data_clean = feature_data[valid_feature_cols]
#     print(f"Using {len(valid_feature_cols)} features for clustering (from {len(feature_columns)} total)")
    
#     # Data preprocessing
#     if use_inverse_normal:
#         print("Applying inverse normal transformation...")
#         data_transformed = inverse_normal_transform(feature_data_clean)
#     else:
#         print("Applying standard scaling...")
#         scaler = StandardScaler()
#         data_transformed = pd.DataFrame(
#             scaler.fit_transform(feature_data_clean), 
#             index=feature_data_clean.index, 
#             columns=feature_data_clean.columns
#         )
    
#     # Compute cosine distance between samples (rows)
#     data_values = data_transformed.values
#     cosine_dist = pdist(data_values, metric='cosine')
    
#     # Perform hierarchical clustering
#     linkage_matrix = linkage(cosine_dist, method=method)
    
#     # Optimal leaf ordering
#     if optimal_ordering:
#         opt_linkage = optimal_leaf_ordering(linkage_matrix, cosine_dist)
#     else:
#         opt_linkage = linkage_matrix
    
#     # Determine optimal number of clusters if not specified
#     if n_clusters is None:
#         n_clusters = find_optimal_clusters_pdist(
#             data_transformed, cosine_dist, linkage_matrix, 
#             min_clusters, max_clusters
#         )
    
#     # Get cluster labels
#     cluster_labels = fcluster(opt_linkage, n_clusters, criterion='maxclust')
    
#     # Create results dictionary
#     results = {
#         'cluster_labels': cluster_labels,
#         'n_clusters': n_clusters,
#         'linkage_matrix': opt_linkage,
#         'linkage_method': method,
#         'data_transformed': data_transformed,
#         'feature_columns': valid_feature_cols,
#         'annotation_columns': annotation_columns,
#         'original_data': df_input,
#         'cosine_distance': cosine_dist
#     }
#     return results


# def find_optimal_clusters_pdist(data, distance_vector, linkage_matrix, min_k, max_k):
#     """Find optimal number of clusters using multiple criteria with precomputed distances."""
#     from sklearn.metrics import silhouette_score, calinski_harabasz_score
    
#     scores = {'silhouette': [], 'calinski_harabasz': []}
#     k_range = range(min_k, min(max_k + 1, len(data)))
    
#     for k in k_range:
#         labels = fcluster(linkage_matrix, k, criterion='maxclust')
        
#         try:
#             distance_matrix = squareform(distance_vector)
#             sil_score = silhouette_score(distance_matrix, labels, metric='precomputed')
#             scores['silhouette'].append(sil_score)
            
#             ch_score = calinski_harabasz_score(data, labels)
#             scores['calinski_harabasz'].append(ch_score)
#         except:
#             scores['silhouette'].append(-1)
#             scores['calinski_harabasz'].append(0)
    
#     if scores['silhouette'] and max(scores['silhouette']) > -1:
#         optimal_k = k_range[np.argmax(scores['silhouette'])]
#     else:
#         optimal_k = min_k
#     return optimal_k


# def identify_annotation_types(df_input, annotation_columns):
#     """Automatically identify which annotations are categorical vs continuous."""
#     categorical_annotations = []
#     continuous_annotations = []
    
#     for col in annotation_columns:
#         if col in df_input.columns:
#             col_data = df_input[col].dropna()
#             if col_data.nunique() < 20 or col_data.dtype == 'object':
#                 categorical_annotations.append(col)
#             else:
#                 continuous_annotations.append(col)
#     return categorical_annotations, continuous_annotations


# def create_color_palettes(df_input, categorical_annotations, continuous_annotations, n_clusters):
#     """Create color palettes for annotations and clusters."""
#     palettes = {}
    
#     # Color-blind friendly palettes
#     categorical_colors = sns.color_palette("Set2", 10)
#     cluster_colors = sns.color_palette("tab10", n_clusters)
    
#     # Create palettes for categorical annotations
#     color_idx = 0
#     for ann_col in categorical_annotations:
#         if ann_col in df_input.columns:
#             unique_vals = sorted([v for v in df_input[ann_col].unique() if pd.notna(v)])
#             n_colors_needed = len(unique_vals)
            
#             if n_colors_needed <= len(categorical_colors):
#                 colors = categorical_colors[color_idx:color_idx + n_colors_needed]
#             else:
#                 colors = sns.color_palette("husl", n_colors_needed)
            
#             palettes[ann_col] = dict(zip(unique_vals, colors))
#             color_idx = (color_idx + n_colors_needed) % len(categorical_colors)
    
#     # Create palette for clusters
#     cluster_vals = list(range(1, n_clusters + 1))
#     palettes['Cluster'] = dict(zip(cluster_vals, cluster_colors))
    
#     # For continuous annotations, use colormaps
#     continuous_palettes = {}
#     for cont_ann in continuous_annotations:
#         if cont_ann in df_input.columns:
#             continuous_palettes[cont_ann] = 'viridis'
    
#     return palettes, continuous_palettes


# def create_row_colors_dataframe(df_input, clustering_results, palettes, continuous_palettes):
#     """Create a DataFrame with row colors for all annotations and clusters."""
#     categorical_annotations, continuous_annotations = identify_annotation_types(
#         df_input, clustering_results['annotation_columns']
#     )
    
#     row_colors_dict = {}
    
#     # Add categorical annotations
#     for ann_col in categorical_annotations:
#         if ann_col in df_input.columns:
#             # Map values to colors, fill NaN with gray
#             mapped_colors = df_input[ann_col].map(palettes[ann_col])
#             # mapped_colors = mapped_colors.fillna((0.7, 0.7, 0.7))  # Gray for NaN
#             row_colors_dict[ann_col] = mapped_colors
    
#     # Add cluster colors
#     cluster_series = pd.Series(clustering_results['cluster_labels'], index=df_input.index)
#     row_colors_dict['Cluster'] = cluster_series.map(palettes['Cluster'])
    
#     # Add continuous annotations
#     for cont_ann in continuous_annotations:
#         if cont_ann in df_input.columns:
#             values = pd.to_numeric(df_input[cont_ann], errors='coerce')
#             valid_mask = ~pd.isna(values)
            
#             if valid_mask.sum() > 0:
#                 min_val, max_val = values[valid_mask].min(), values[valid_mask].max()
                
#                 if max_val > min_val:
#                     normalized_values = (values - min_val) / (max_val - min_val)
#                 else:
#                     normalized_values = pd.Series(0.5, index=values.index)
                
#                 # Apply colormap
#                 cmap = plt.cm.get_cmap(continuous_palettes[cont_ann])
#                 colors = []
#                 for val in normalized_values:
#                     if pd.isna(val):
#                         colors.append((0.7, 0.7, 0.7))  # Gray for NaN
#                     else:
#                         colors.append(cmap(val))
                
#                 row_colors_dict[cont_ann] = colors
    
#     # Create DataFrame
#     row_colors_df = pd.DataFrame(row_colors_dict, index=df_input.index)
#     return row_colors_df


# def create_comprehensive_clustermap_samples_features(df_input, clustering_results, 
#                                                    figsize=(15, 12), cmap='vlag', 
#                                                    save_path=None, show_legends=True):
#     """
#     Create a comprehensive clustermap with samples as rows and features as columns.
#     Fixed version that properly handles row coloring.
#     """
    
#     # Get the clustering data (samples x features - NO transpose needed)
#     data = clustering_results['data_transformed']  # Keep as samples x features
#     linkage_matrix = clustering_results['linkage_matrix']
    
#     # Identify annotation types
#     categorical_annotations, continuous_annotations = identify_annotation_types(
#         df_input, clustering_results['annotation_columns']
#     )
    
#     # Create color palettes
#     palettes, continuous_palettes = create_color_palettes(
#         df_input, categorical_annotations, continuous_annotations,
#         clustering_results['n_clusters']
#     )
    
#     # Create row colors DataFrame (for samples)
#     row_colors_df = create_row_colors_dataframe(
#         df_input, clustering_results, palettes, continuous_palettes
#     )
    
#     print(f"Data shape: {data.shape}")
#     print(f"Row colors shape: {row_colors_df.shape}")
#     print(f"Annotation columns: {list(row_colors_df.columns)}")
    
#     # Create the clustermap - FIXED VERSION
#     g = sns.clustermap(
#         data,  # samples x features
#         row_linkage=linkage_matrix,  # Use precomputed linkage for samples
#         row_colors=row_colors_df,    # Colors for samples (rows)
#         col_cluster=True,            # Cluster features (columns)
#         row_cluster=False,           # Don't re-cluster samples, use precomputed
#         cmap=cmap,
#         figsize=figsize,
#         cbar_kws={'label': 'Transformed Feature Value'},
#         linewidths=0.01,
#         xticklabels=False,  # Hide feature labels
#         yticklabels=False   # Hide sample labels
#     )
    
#     # Customize the plot
#     g.ax_heatmap.set_ylabel('Samples (Clustered)', fontsize=12)
#     g.ax_heatmap.set_xlabel('Features', fontsize=12)
    
#     # Add legends if requested
#     if show_legends:
#         add_legends_to_clustermap(g, palettes, continuous_palettes, df_input, continuous_annotations)
    
#     plt.suptitle(f'Comprehensive Sample Clustering Analysis\n'
#                 f'Method: {clustering_results["linkage_method"]}, '
#                 f'Clusters: {clustering_results["n_clusters"]}, '
#                 f'Features: {len(clustering_results["feature_columns"])}', 
#                 y=0.98, fontsize=14, weight='bold')
    
#     if save_path:
#         plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
#     return g, palettes, continuous_palettes


# def add_legends_to_clustermap(g, palettes, continuous_palettes, df_input, continuous_annotations):
#     """Add legends to the clustermap."""
#     legend_x_start = 1.05
#     legend_y_start = 0.95
#     legend_spacing = 0.15
    
#     # Legend for categorical annotations and clusters
#     legend_count = 0
#     for ann_name, palette in palettes.items():
#         handles = [mpatches.Patch(color=color, label=label) 
#                   for label, color in palette.items()]
        
#         legend = plt.legend(handles=handles, title=ann_name,
#                           bbox_to_anchor=(legend_x_start, legend_y_start - legend_count * legend_spacing), 
#                           bbox_transform=g.ax_heatmap.transAxes,
#                           loc='upper left')
#         legend.set_title(ann_name, prop={'weight': 'bold'})
#         plt.gca().add_artist(legend)
#         legend_count += 1
    
#     # Add continuous annotation info
#     cont_info = []
#     for cont_ann in continuous_annotations:
#         if cont_ann in df_input.columns:
#             values = pd.to_numeric(df_input[cont_ann], errors='coerce')
#             valid_values = values.dropna()
#             if len(valid_values) > 0:
#                 cont_info.append(f"{cont_ann}: {valid_values.min():.2f}-{valid_values.max():.2f}")
    
#     if cont_info:
#         plt.figtext(0.02, 0.02, "Continuous annotations: " + ", ".join(cont_info), 
#                    fontsize=10, style='italic')


# def create_umap_visualization(df_input, clustering_results, 
#                             n_neighbors=15, min_dist=0.1, n_components=2,
#                             figsize=(12, 8), save_path=None):
#     """Create UMAP visualization of the samples colored by cluster assignments."""
    
#     # Get the transformed feature data
#     feature_data = clustering_results['data_transformed']
#     cluster_labels = clustering_results['cluster_labels']
    
#     # Identify annotation types
#     categorical_annotations, continuous_annotations = identify_annotation_types(
#         df_input, clustering_results['annotation_columns']
#     )
    
#     # Perform UMAP
#     print(f"Performing UMAP with n_neighbors={n_neighbors}, min_dist={min_dist}...")
#     reducer = umap.UMAP(
#         n_neighbors=n_neighbors,
#         min_dist=min_dist,
#         n_components=n_components,
#         random_state=42,
#         metric='cosine'
#     )
    
#     embedding = reducer.fit_transform(feature_data)
    
#     # Create UMAP DataFrame
#     umap_df = pd.DataFrame({
#         'UMAP1': embedding[:, 0],
#         'UMAP2': embedding[:, 1],
#         'Cluster': cluster_labels,
#         'Sample': df_input.index
#     })
    
#     # Add annotations to UMAP DataFrame
#     for ann_col in clustering_results['annotation_columns']:
#         if ann_col in df_input.columns:
#             umap_df[ann_col] = df_input[ann_col].values
    
#     # Create color palette for clusters
#     cluster_colors = sns.color_palette("tab10", clustering_results['n_clusters'])
#     cluster_palette = dict(zip(range(1, clustering_results['n_clusters'] + 1), cluster_colors))
    
#     # Create subplots
#     n_plots = 1 + len(categorical_annotations) + min(len(continuous_annotations), 2)
#     n_cols = min(3, n_plots)
#     n_rows = (n_plots + n_cols - 1) // n_cols
    
#     fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
#     if n_plots == 1:
#         axes = [axes]
#     elif n_rows == 1:
#         axes = [axes] if n_cols == 1 else list(axes)
#     else:
#         axes = axes.flatten()
    
#     plot_idx = 0
    
#     # Plot 1: Colored by clusters
#     scatter = axes[plot_idx].scatter(
#         umap_df['UMAP1'], umap_df['UMAP2'], 
#         c=[cluster_palette[c] for c in umap_df['Cluster']], 
#         alpha=0.7, s=50
#     )
#     axes[plot_idx].set_title(f'UMAP colored by Clusters (n={clustering_results["n_clusters"]})')
#     axes[plot_idx].set_xlabel('UMAP1')
#     axes[plot_idx].set_ylabel('UMAP2')
    
#     # Add cluster labels
#     for cluster_id in range(1, clustering_results['n_clusters'] + 1):
#         cluster_mask = umap_df['Cluster'] == cluster_id
#         if cluster_mask.sum() > 0:
#             center_x = umap_df.loc[cluster_mask, 'UMAP1'].mean()
#             center_y = umap_df.loc[cluster_mask, 'UMAP2'].mean()
#             axes[plot_idx].annotate(f'C{cluster_id}', (center_x, center_y), 
#                                   fontsize=12, fontweight='bold', 
#                                   ha='center', va='center',
#                                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
    
#     plot_idx += 1
    
#     # Plot categorical annotations
#     for ann_col in categorical_annotations:
#         if plot_idx >= len(axes):
#             break
#         if ann_col in umap_df.columns:
#             unique_vals = [v for v in umap_df[ann_col].unique() if pd.notna(v)]
#             if len(unique_vals) > 1:
#                 colors = sns.color_palette("Set2", len(unique_vals))
#                 color_map = dict(zip(unique_vals, colors))
                
#                 for val in unique_vals:
#                     mask = umap_df[ann_col] == val
#                     axes[plot_idx].scatter(
#                         umap_df.loc[mask, 'UMAP1'], umap_df.loc[mask, 'UMAP2'],
#                         c=[color_map[val]], label=str(val), alpha=0.7, s=50
#                     )
                
#                 axes[plot_idx].set_title(f'UMAP colored by {ann_col}')
#                 axes[plot_idx].set_xlabel('UMAP1')
#                 axes[plot_idx].set_ylabel('UMAP2')
#                 axes[plot_idx].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
                
#                 plot_idx += 1
    
#     # Plot continuous annotations (max 2)
#     for ann_col in continuous_annotations[:2]:
#         if plot_idx >= len(axes):
#             break
#         if ann_col in umap_df.columns:
#             values = pd.to_numeric(umap_df[ann_col], errors='coerce')
#             scatter = axes[plot_idx].scatter(
#                 umap_df['UMAP1'], umap_df['UMAP2'], 
#                 c=values, cmap='viridis', alpha=0.7, s=50
#             )
#             axes[plot_idx].set_title(f'UMAP colored by {ann_col}')
#             axes[plot_idx].set_xlabel('UMAP1')
#             axes[plot_idx].set_ylabel('UMAP2')
#             plt.colorbar(scatter, ax=axes[plot_idx])
            
#             plot_idx += 1
    
#     # Hide unused subplots
#     for i in range(plot_idx, len(axes)):
#         axes[i].set_visible(False)
    
#     plt.tight_layout()
    
#     if save_path:
#         plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
#     return fig, umap_df, reducer


# def analyze_cluster_characteristics_samples_features(df_input, clustering_results):
#     """Analyze and summarize characteristics of each cluster."""
    
#     categorical_annotations, continuous_annotations = identify_annotation_types(
#         df_input, clustering_results['annotation_columns']
#     )
    
#     cluster_summary = {}
    
#     for cluster_id in range(1, clustering_results['n_clusters'] + 1):
#         cluster_mask = clustering_results['cluster_labels'] == cluster_id
#         cluster_samples = df_input.index[cluster_mask]
        
#         cluster_summary[f'Cluster_{cluster_id}'] = {
#             'n_samples': len(cluster_samples),
#             'samples': list(cluster_samples)
#         }
        
#         # Categorical annotation summaries
#         for ann_col in categorical_annotations:
#             if ann_col in df_input.columns:
#                 cluster_values = df_input.loc[cluster_samples, ann_col]
#                 value_counts = cluster_values.value_counts()
#                 cluster_summary[f'Cluster_{cluster_id}'][ann_col] = value_counts.to_dict()
        
#         # Continuous annotation summaries
#         for cont_ann in continuous_annotations:
#             if cont_ann in df_input.columns:
#                 cluster_values = pd.to_numeric(df_input.loc[cluster_samples, cont_ann], errors='coerce')
#                 cluster_summary[f'Cluster_{cluster_id}'][cont_ann] = {
#                     'mean': cluster_values.mean(),
#                     'std': cluster_values.std(),
#                     'min': cluster_values.min(),
#                     'max': cluster_values.max()
#                 }
    
#     return cluster_summary


# def run_comprehensive_analysis_samples_features(df_input, 
#                                                annotation_columns=None,
#                                                feature_columns=None,
#                                                n_clusters=None,
#                                                method='average',
#                                                use_inverse_normal=True,
#                                                umap_params=None,
#                                                save_prefix=None):
#     """
#     Run the complete analysis pipeline for samples x features format.
#     """
    
#     print("Starting comprehensive clustering analysis (samples x features format)...")
    
#     # Set defaults
#     if annotation_columns is None:
#         annotation_columns = ["clinvar_clnsig_clean", "Altered_95th_perc_GFP", "AUROC"]
    
#     if umap_params is None:
#         umap_params = {'n_neighbors': 15, 'min_dist': 0.1}
    
#     # Perform clustering
#     print("Performing robust cosine similarity clustering...")
#     clustering_results = robust_cosine_clustering_samples_features(
#         df_input,
#         annotation_columns=annotation_columns,
#         feature_columns=feature_columns,
#         n_clusters=n_clusters,
#         method=method,
#         use_inverse_normal=use_inverse_normal
#     )
    
#     print(f"Optimal number of clusters: {clustering_results['n_clusters']}")
#     print(f"Using {len(clustering_results['feature_columns'])} features for clustering")
    
#     # Create comprehensive clustermap
#     print("Creating comprehensive clustermap...")
#     clustermap_save = f"{save_prefix}_clustermap.png" if save_prefix else None
#     g, palettes, continuous_palettes = create_comprehensive_clustermap_samples_features(
#         df_input, clustering_results, save_path=clustermap_save
#     )
    
#     # Create UMAP visualization
#     print("Creating UMAP visualization...")
#     umap_save = f"{save_prefix}_umap.png" if save_prefix else None
#     fig_umap, umap_df, reducer = create_umap_visualization(
#         df_input, clustering_results, save_path=umap_save, **umap_params
#     )
    
#     # Analyze cluster characteristics
#     print("Analyzing cluster characteristics...")
#     cluster_summary = analyze_cluster_characteristics_samples_features(df_input, clustering_results)
    
#     # Print summary
#     categorical_annotations, continuous_annotations = identify_annotation_types(
#         df_input, clustering_results['annotation_columns']
#     )
    
#     for cluster_name, summary in cluster_summary.items():
#         print(f"\n{cluster_name}: {summary['n_samples']} samples")
#         for ann_col in categorical_annotations:
#             if ann_col in summary:
#                 print(f"  {ann_col}: {summary[ann_col]}")
#         for cont_ann in continuous_annotations:
#             if cont_ann in summary:
#                 stats = summary[cont_ann]
#                 print(f"  {cont_ann}: mean={stats['mean']:.3f} ± {stats['std']:.3f}")
    
#     return {
#         'clustering_results': clustering_results,
#         'clustermap': g,
#         'palettes': palettes,
#         'continuous_palettes': continuous_palettes,
#         'cluster_summary': cluster_summary,
#         'umap_fig': fig_umap,
#         'umap_df': umap_df,
#         'umap_reducer': reducer
#     }