"""
Enhanced Structural Features Integration Script
==============================================

This script demonstrates how to extract and combine all the new structural
features for mutation impact assessment. It integrates:
- B-factor gradients (structural flexibility transitions)
- Contact numbers (local packing density)
- Enhanced secondary structure (transitions, loop analysis)
- Distance to center (core vs surface positioning)
- Backbone dihedrals (conformational constraints)

Usage example for integrating with your existing VarChAMP analysis pipeline.

Author: Runxi Shen
Created: 2025-01-XX
"""

import os
import sys
import time
from pathlib import Path
from typing import List, Dict, Optional, Union
import numpy as np
import pandas as pd
import polars as pl

# Import all the new feature calculation modules
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import calculate_bfactor_gradients as bfactor_grad
import calculate_contact_numbers as contact_num
import calculate_secondary_structure as sec_struct
import calculate_distance_to_center as dist_center
import calculate_backbone_dihedrals as dihedrals


def extract_all_structural_features(pdb_files: List[str],
                                   output_dir: str,
                                   batch_size: int = 100,
                                   n_processes: int = None) -> Dict[str, List[str]]:
    """
    Extract all enhanced structural features from PDB files.
    
    Args:
        pdb_files: List of PDB file paths
        output_dir: Base output directory
        batch_size: Number of files per batch
        n_processes: Number of parallel processes
        
    Returns:
        Dictionary mapping feature types to parquet file lists
    """
    print(f"Processing {len(pdb_files)} PDB files...")
    print(f"Output directory: {output_dir}")
    
    # Create output subdirectories
    feature_dirs = {
        'bfactor_gradients': f"{output_dir}/bfactor_gradients",
        'contact_numbers': f"{output_dir}/contact_numbers", 
        'secondary_structure': f"{output_dir}/secondary_structure",
        'distance_to_center': f"{output_dir}/distance_to_center",
        'backbone_dihedrals': f"{output_dir}/backbone_dihedrals"
    }
    
    parquet_files = {}
    
    # 1. B-factor gradients (flexibility transitions)
    print("\\n1. Calculating B-factor gradients...")
    start_time = time.time()
    parquet_files['bfactor_gradients'] = bfactor_grad.process_and_save_batches(
        pdb_files, 
        feature_dirs['bfactor_gradients'],
        batch_size=batch_size,
        window_size=5,  # 5-residue window for gradient calculation
        n_processes=n_processes
    )
    print(f"   Completed in {time.time() - start_time:.2f} seconds")
    
    # 2. Contact numbers (local packing density)
    print("\\n2. Calculating contact numbers...")
    start_time = time.time()
    parquet_files['contact_numbers'] = contact_num.process_and_save_batches(
        pdb_files,
        feature_dirs['contact_numbers'], 
        batch_size=batch_size,
        cutoff_distances=[6.0, 8.0, 10.0],  # Different interaction distances
        atom_type='CA',
        n_processes=n_processes
    )
    print(f"   Completed in {time.time() - start_time:.2f} seconds")
    
    # 3. Enhanced secondary structure
    print("\\n3. Calculating enhanced secondary structure features...")
    start_time = time.time()
    parquet_files['secondary_structure'] = sec_struct.process_and_save_batches(
        pdb_files,
        feature_dirs['secondary_structure'],
        batch_size=batch_size,
        dssp_executable='mkdssp',
        n_processes=n_processes
    )
    print(f"   Completed in {time.time() - start_time:.2f} seconds")
    
    # 4. Distance to center (core vs surface)
    print("\\n4. Calculating distance to protein center...")
    start_time = time.time()
    parquet_files['distance_to_center'] = dist_center.process_and_save_batches(
        pdb_files,
        feature_dirs['distance_to_center'],
        batch_size=batch_size,
        center_types=['geometric', 'mass', 'ca_geometric'],
        n_processes=n_processes
    )
    print(f"   Completed in {time.time() - start_time:.2f} seconds")
    
    # 5. Backbone dihedral angles
    print("\\n5. Calculating backbone dihedral angles...")
    start_time = time.time()
    parquet_files['backbone_dihedrals'] = dihedrals.process_and_save_batches(
        pdb_files,
        feature_dirs['backbone_dihedrals'],
        batch_size=batch_size,
        include_omega=True,
        n_processes=n_processes
    )
    print(f"   Completed in {time.time() - start_time:.2f} seconds")
    
    print(f"\\nAll features extracted successfully!")
    return parquet_files


def consolidate_features_for_variants(parquet_files: Dict[str, List[str]],
                                    variants_df: pd.DataFrame,
                                    uniprot_col: str = 'Uniprot_acc',
                                    position_col: str = 'aapos',
                                    residue_col: str = 'aaref_3') -> pd.DataFrame:
    """
    Consolidate all structural features for a set of variants.
    
    Args:
        parquet_files: Dictionary of feature parquet files
        variants_df: DataFrame with variants (must have uniprot, position, residue)
        uniprot_col: Column name for UniProt accession
        position_col: Column name for amino acid position
        residue_col: Column name for reference residue
        
    Returns:
        DataFrame with all structural features added
    """
    print("Consolidating structural features for variants...")
    
    # Convert to Polars for efficient operations
    variants_pl = pl.from_pandas(variants_df)
    
    # For each feature type, create a lookup function
    feature_data = {}
    
    print("Loading feature data...")
    
    # 1. B-factor gradients
    bfactor_df = pl.concat([
        pl.read_parquet(f) for f in parquet_files['bfactor_gradients']
    ])
    feature_data['bfactor'] = bfactor_df
    
    # 2. Contact numbers  
    contact_df = pl.concat([
        pl.read_parquet(f) for f in parquet_files['contact_numbers']
    ])
    feature_data['contacts'] = contact_df
    
    # 3. Secondary structure
    ss_df = pl.concat([
        pl.read_parquet(f) for f in parquet_files['secondary_structure'] 
    ])
    feature_data['secondary'] = ss_df
    
    # 4. Distance to center
    distance_df = pl.concat([
        pl.read_parquet(f) for f in parquet_files['distance_to_center']
    ])
    feature_data['distance'] = distance_df
    
    # 5. Backbone dihedrals
    dihedral_df = pl.concat([
        pl.read_parquet(f) for f in parquet_files['backbone_dihedrals']
    ])
    feature_data['dihedrals'] = dihedral_df
    
    print("Joining features with variants...")
    
    # Create join keys - need to handle AF model naming
    variants_with_features = variants_pl.with_columns([
        pl.concat_str([
            pl.lit("AF-"),
            pl.col(uniprot_col),
            pl.lit("-F1-model_v4")
        ]).alias("af_pdb_id"),
        pl.col(position_col).cast(pl.Int64).alias("residue_pos")
    ])
    
    # Join each feature type
    for feature_name, feature_df in feature_data.items():
        print(f"  Joining {feature_name} features...")
        
        # Prepare feature df for joining
        feature_df = feature_df.with_columns([
            pl.col("residue_num").cast(pl.Int64)
        ])
        
        # Join on PDB ID and residue position
        variants_with_features = variants_with_features.join(
            feature_df,
            left_on=["af_pdb_id", "residue_pos"],
            right_on=["pdb_id", "residue_num"],
            how="left",
            suffix=f"_{feature_name}"
        )
    
    # Convert back to pandas
    result_df = variants_with_features.to_pandas()
    
    print(f"Feature consolidation complete. Added features for {len(result_df)} variants.")
    
    return result_df


def analyze_feature_correlations(consolidated_df: pd.DataFrame) -> pd.DataFrame:
    """
    Analyze correlations between different structural features.
    
    Args:
        consolidated_df: DataFrame with all features
        
    Returns:
        Correlation matrix DataFrame
    """
    print("Analyzing feature correlations...")
    
    # Select key structural features for correlation analysis
    feature_cols = [
        # B-factor gradients
        'absolute_gradient', 'relative_gradient', 'max_change',
        
        # Contact numbers
        'contacts_6A', 'contacts_8A', 'contacts_10A',
        'density_8A',
        
        # Distance metrics
        'ca_dist_geometric', 'normalized_dist_geometric',
        
        # Dihedral constraints
        'phi', 'psi', 'strain_score',
        
        # Secondary structure (will need encoding)
        'dist_to_structured', 'ss_run_length'
    ]
    
    # Filter to available columns
    available_cols = [col for col in feature_cols if col in consolidated_df.columns]
    
    if len(available_cols) < 2:
        print("Not enough numeric features for correlation analysis")
        return pd.DataFrame()
    
    correlation_matrix = consolidated_df[available_cols].corr()
    
    print("Feature correlation analysis complete.")
    return correlation_matrix


def create_feature_summary_stats(consolidated_df: pd.DataFrame) -> Dict[str, Dict]:
    """
    Create summary statistics for all structural features.
    
    Args:
        consolidated_df: DataFrame with all features
        
    Returns:
        Dictionary with summary statistics for each feature
    """
    print("Generating feature summary statistics...")
    
    # Group features by type
    feature_groups = {
        'flexibility': ['absolute_gradient', 'relative_gradient', 'max_change'],
        'packing': ['contacts_6A', 'contacts_8A', 'contacts_10A', 'density_8A'],
        'position': ['ca_dist_geometric', 'normalized_dist_geometric', 'radial_percentile_geometric'],
        'conformation': ['phi', 'psi', 'omega', 'strain_score'],
        'structure': ['dist_to_structured', 'ss_run_length', 'loop_length']
    }
    
    summary_stats = {}
    
    for group_name, feature_list in feature_groups.items():
        available_features = [f for f in feature_list if f in consolidated_df.columns]
        
        if available_features:
            group_stats = {}
            for feature in available_features:
                if consolidated_df[feature].dtype in ['float64', 'int64']:
                    group_stats[feature] = {
                        'count': consolidated_df[feature].count(),
                        'mean': consolidated_df[feature].mean(),
                        'std': consolidated_df[feature].std(),
                        'min': consolidated_df[feature].min(),
                        'max': consolidated_df[feature].max(),
                        'q25': consolidated_df[feature].quantile(0.25),
                        'q50': consolidated_df[feature].quantile(0.50), 
                        'q75': consolidated_df[feature].quantile(0.75)
                    }
            
            summary_stats[group_name] = group_stats
    
    return summary_stats


## Example usage and testing
def test_enhanced_features(sample_pdb_dir: str, output_dir: str, n_sample: int = 10):
    """
    Test the enhanced structural features on a sample of PDB files.
    
    Args:
        sample_pdb_dir: Directory containing sample PDB files
        output_dir: Output directory for results
        n_sample: Number of PDB files to process for testing
    """
    print("Testing enhanced structural features...")
    
    # Get sample PDB files
    pdb_files = list(Path(sample_pdb_dir).glob("*.pdb*"))[:n_sample]
    
    if not pdb_files:
        print(f"No PDB files found in {sample_pdb_dir}")
        return
    
    print(f"Testing with {len(pdb_files)} PDB files")
    
    # Extract features
    parquet_files = extract_all_structural_features(
        pdb_files, 
        output_dir,
        batch_size=5,
        n_processes=2
    )
    
    # Create sample variants for testing
    sample_variants = pd.DataFrame({
        'Uniprot_acc': ['X6R8D5', 'Q86YR6'],  # Sample UniProt IDs
        'aapos': [50, 100],
        'aaref_3': ['ALA', 'VAL']
    })
    
    # Test feature consolidation
    try:
        consolidated_df = consolidate_features_for_variants(
            parquet_files, 
            sample_variants
        )
        
        print(f"\\nConsolidated features shape: {consolidated_df.shape}")
        print("Available feature columns:")
        feature_cols = [col for col in consolidated_df.columns 
                       if any(feat in col for feat in ['gradient', 'contacts', 'dist_', 'phi', 'psi'])]
        for col in sorted(feature_cols):
            print(f"  {col}")
        
        # Test correlation analysis
        if len(feature_cols) > 1:
            correlations = analyze_feature_correlations(consolidated_df)
            print(f"\\nFeature correlation matrix shape: {correlations.shape}")
        
        # Test summary statistics
        summary_stats = create_feature_summary_stats(consolidated_df)
        print(f"\\nSummary statistics generated for {len(summary_stats)} feature groups")
        
        print("\\nTesting completed successfully!")
        
    except Exception as e:
        print(f"Error during testing: {e}")


if __name__ == "__main__":
    # Example usage
    sample_pdb_dir = "../../1_inputs/raw_data/UP000005640_9606_HUMAN_v4"
    output_dir = "../../3_outputs/processed_data/enhanced_structural_features"
    
    # Test with a small sample
    test_enhanced_features(sample_pdb_dir, output_dir, n_sample=5)