#!/home/shenrunx/software/anaconda3/envs/varchamp/bin python3

import pandas as pd
import numpy as np
from scipy.spatial.distance import mahalanobis
from scipy.linalg import pinv
import re


def standardize_location(location_str):
    """Standardize location names to lowercase and clean format"""
    if pd.isna(location_str) or location_str.strip() == '':
        return []
    # Split by semicolon and clean each location
    locations = []
    for loc in location_str.split(';'):
        # Remove GO IDs and parentheses content
        loc = re.sub(r'\s*\([^)]*\)', '', loc)
        # Clean whitespace and standardize
        loc = loc.strip().lower()
        if loc:
            locations.append(loc)
    return locations


def extract_unique_locations(df):
    """Extract all unique subcellular locations from the dataset"""
    all_locations = set()
    
    # Process main locations
    for loc_str in df['Metadata_hpa_main_loc'].dropna():
        all_locations.update(standardize_location(loc_str))
    
    # Process additional locations
    for loc_str in df['Metadata_hpa_additional_loc'].dropna():
        all_locations.update(standardize_location(loc_str))
    
    return sorted(list(all_locations))

    
def create_location_vectors(df, unique_locations):
    """Create binary location vectors for each allele"""
    location_matrix = np.zeros((len(df), len(unique_locations)))
    
    for idx, (_, row) in enumerate(df.iterrows()):
        # Get locations from both columns
        main_locs = standardize_location(row['Metadata_hpa_main_loc'])
        add_locs = standardize_location(row['Metadata_hpa_additional_loc'])
        all_locs = set(main_locs + add_locs)
        
        # Set binary indicators
        for loc in all_locs:
            if loc in unique_locations:
                loc_idx = unique_locations.index(loc)
                location_matrix[idx, loc_idx] = 1
    
    return location_matrix


def compute_mahalanobis_distances(features, target_features):
    """Compute Mahalanobis distances between target and reference features"""
    # Calculate covariance matrix
    cov_matrix = np.cov(features.T)
    
    # Use pseudo-inverse if matrix is singular
    try:
        inv_cov = np.linalg.inv(cov_matrix)
    except np.linalg.LinAlgError:
        inv_cov = pinv(cov_matrix)
    
    distances = []
    mean_features = np.mean(features, axis=0)
    
    for target in target_features:
        dist = mahalanobis(target, mean_features, inv_cov)
        distances.append(dist)
    
    return np.array(distances)


def annotate_unannotated_alleles(df, location_matrix, unique_locations, feature_cols, k=1):
    """Annotate alleles without location data using k-nearest neighbors"""
    # Identify annotated and unannotated alleles
    has_annotation = location_matrix.sum(axis=1) > 0
    annotated_idx = np.where(has_annotation)[0]
    unannotated_idx = np.where(~has_annotation)[0]
    
    if len(unannotated_idx) == 0:
        return location_matrix
    
    # Get features for annotated and unannotated alleles
    annotated_features = df.iloc[annotated_idx][feature_cols].values
    unannotated_features = df.iloc[unannotated_idx][feature_cols].values
    annotated_locations = location_matrix[annotated_idx]
    
    # For each unannotated allele, find k nearest neighbors
    updated_matrix = location_matrix.copy()
    
    for i, unanno_idx in enumerate(unannotated_idx):
        target_features = unannotated_features[i:i+1]
        
        # Compute distances to all annotated alleles
        distances = []
        for j, anno_features in enumerate(annotated_features):
            # Simple Euclidean distance (Mahalanobis can be computationally expensive)
            dist = np.linalg.norm(target_features[0] - anno_features)
            distances.append((dist, j))
        
        # Get k nearest neighbors
        distances.sort(key=lambda x: x[0])
        nearest_k = distances[:min(k, len(distances))]
        
        # Process k nearest neighbors
        if nearest_k:
            neighbor_locations = np.array([annotated_locations[j] for _, j in nearest_k])
            
            if k == 1:
                # Use single nearest neighbor directly
                updated_matrix[unanno_idx] = neighbor_locations[0]
            else:
                # Average the location vectors of k nearest neighbors
                avg_locations = np.mean(neighbor_locations, axis=0)
                # Convert to binary (>0.5 threshold)
                binary_locations = (avg_locations > 0.5).astype(int)
                updated_matrix[unanno_idx] = binary_locations
    
    return updated_matrix


def main():
    # Parameters
    K_NEIGHBORS = 1  # Number of nearest neighbors to use
    
    # Load data
    csv_path = "../../3_outputs/ref_var_allele_int_featsel_cp_feat_w_loc.csv"
    df = pd.read_csv(csv_path)
    
    print(f"Loaded data: {len(df)} alleles")
    
    # Get feature columns (non-Metadata columns)
    feature_cols = [col for col in df.columns if not col.startswith('Metadata_')]
    print(f"Feature columns: {len(feature_cols)}")
    
    # Extract unique locations
    unique_locations = extract_unique_locations(df)
    print(f"Unique subcellular locations found: {len(unique_locations)}")
    print("Locations:", unique_locations)
    
    # Create initial location vectors for annotated alleles
    initial_location_matrix = create_location_vectors(df, unique_locations)
    
    # Count annotated vs unannotated
    has_annotation = initial_location_matrix.sum(axis=1) > 0
    n_annotated = np.sum(has_annotation)
    n_unannotated = len(df) - n_annotated
    
    print(f"Annotated alleles: {n_annotated}")
    print(f"Unannotated alleles: {n_unannotated}")
    
    if n_unannotated > 0:
        print(f"Annotating unannotated alleles using {K_NEIGHBORS}-NN approach...")
        # Annotate unannotated alleles
        final_location_matrix = annotate_unannotated_alleles(
            df, initial_location_matrix, unique_locations, feature_cols, k=K_NEIGHBORS
        )
        
        # Count newly annotated
        newly_annotated = np.sum(final_location_matrix.sum(axis=1) > 0) - n_annotated
        print(f"Newly annotated alleles: {newly_annotated}")
    else:
        final_location_matrix = initial_location_matrix
        print("All alleles already have annotations")
    
    # Create output dataframe
    output_df = df[['Metadata_gene_allele']].copy()
    
    # Add location columns
    for i, loc in enumerate(unique_locations):
        col_name = f"loc_{loc.replace(' ', '_').replace('-', '_')}"
        output_df[col_name] = final_location_matrix[:, i]
    
    # Add original location annotations for reference
    output_df['original_main_loc'] = df['Metadata_hpa_main_loc']
    output_df['original_additional_loc'] = df['Metadata_hpa_additional_loc']
    
    # Save results
    output_path = "../../3_outputs/subcellular_localization_annotations.csv"
    output_df.to_csv(output_path, index=False)
    print(f"Results saved to: {output_path}")
    
    # Summary statistics
    print("\nSummary:")
    print(f"Total alleles processed: {len(df)}")
    print(f"Unique locations: {len(unique_locations)}")
    
    for i, loc in enumerate(unique_locations):
        count = np.sum(final_location_matrix[:, i])
        print(f"  {loc}: {count} alleles")
    
    return output_df, final_location_matrix, unique_locations


if __name__ == "__main__":
    output_df, location_matrix, locations = main()