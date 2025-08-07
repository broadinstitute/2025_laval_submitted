"""
Distance to Protein Center Storage and Query System
=================================================

This module calculates various distance metrics from protein structure PDB files,
providing insights into:
- Global protein architecture
- Core vs. surface positioning
- Radial distribution of residues
- Geometric centrality measures

Distance to center metrics help understand mutation impact because:
- Core residues (close to geometric center) often have higher structural importance
- Surface residues may have different functional roles
- Radial position correlates with evolutionary conservation
- Distance gradients indicate protein compactness and organization

Key features:
- Multiple center definitions (geometric center, center of mass, domain centers)
- Both absolute distances and normalized metrics
- Parallel processing of PDB files using multiprocessing
- Efficient storage using Apache Parquet format
- Fast query capabilities with predicate pushdown
- Batch processing to handle large datasets with minimal memory overhead

Author: Runxi Shen
Created: 2025-01-XX
"""

import os
import gzip
import warnings
from pathlib import Path
from functools import partial
from multiprocessing import Pool, cpu_count
from typing import List, Dict, Tuple, Optional, Union
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from Bio.PDB import PDBParser
from io import StringIO


# Atomic masses for center of mass calculation
ATOMIC_MASSES = {
    'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'S': 32.065,
    'P': 30.974, 'SE': 78.96, 'FE': 55.845, 'ZN': 65.38, 'MG': 24.305,
    'CA': 40.078, 'MN': 54.938, 'CU': 63.546, 'NI': 58.693, 'CO': 58.933
}


def calculate_distance_to_center(pdb_file: Union[str, Path], 
                                center_types: List[str] = ['geometric', 'mass', 'ca_geometric'],
                                include_radial_stats: bool = True) -> Dict[Tuple[str, int], dict]:
    """
    Calculate distance to various protein centers for each residue.
    
    Args:
        pdb_file: Path to PDB file (can be gzipped)
        center_types: Types of centers to calculate ('geometric', 'mass', 'ca_geometric')
        include_radial_stats: Whether to include radial distribution statistics
        
    Returns:
        Dictionary mapping (chain_id, residue_number) to distance metrics
    """
    pdb_path = Path(pdb_file).resolve()
    
    try:
        # Handle gzipped files
        if pdb_path.suffix == '.gz':
            with gzip.open(pdb_path, 'rt') as f:
                pdb_content = f.read()
        else:
            with open(pdb_path, 'r') as f:
                pdb_content = f.read()
                
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("AF", StringIO(pdb_content))
        
        # Extract all atom coordinates and residue information
        all_coords = []
        all_masses = []
        ca_coords = []
        residue_data = []
        
        for residue in structure.get_residues():
            try:
                chain_id = residue.get_parent().id
                res_num = residue.id[1]
                res_name = residue.resname
                
                # Get CA coordinate
                ca_coord = residue["CA"].coord
                ca_coords.append(ca_coord)
                
                # Get all heavy atom coordinates and masses for this residue
                residue_coords = []
                residue_masses = []
                
                for atom in residue:
                    if atom.element != 'H':  # Skip hydrogens
                        coord = atom.coord
                        mass = ATOMIC_MASSES.get(atom.element, 12.0)  # Default to carbon
                        
                        all_coords.append(coord)
                        all_masses.append(mass)
                        residue_coords.append(coord)
                        residue_masses.append(mass)
                
                # Calculate residue center of mass
                if residue_masses:
                    residue_com = np.average(residue_coords, weights=residue_masses, axis=0)
                else:
                    residue_com = ca_coord
                
                residue_data.append({
                    'chain_id': chain_id,
                    'res_num': res_num,
                    'residue_name': res_name,
                    'ca_coord': ca_coord,
                    'residue_com': residue_com,
                    'residue_coords': residue_coords,
                    'num_heavy_atoms': len(residue_coords)
                })
                
            except KeyError:
                continue
        
        if not residue_data:
            return {}
        
        # Convert to numpy arrays
        all_coords = np.array(all_coords)
        all_masses = np.array(all_masses)
        ca_coords = np.array(ca_coords)
        
        # Calculate different protein centers
        centers = {}
        
        if 'geometric' in center_types:
            # Geometric center of all heavy atoms
            centers['geometric'] = np.mean(all_coords, axis=0)
        
        if 'mass' in center_types:
            # Center of mass of all heavy atoms
            centers['mass'] = np.average(all_coords, weights=all_masses, axis=0)
        
        if 'ca_geometric' in center_types:
            # Geometric center of CA atoms only
            centers['ca_geometric'] = np.mean(ca_coords, axis=0)
        
        # Calculate distances and metrics for each residue
        distances = {}
        
        for residue in residue_data:
            chain_id = residue['chain_id']
            res_num = residue['res_num']
            ca_coord = residue['ca_coord']
            residue_com = residue['residue_com']
            
            distance_metrics = {
                'residue_name': residue['residue_name'],
                'num_heavy_atoms': residue['num_heavy_atoms']
            }
            
            # Calculate distances to each center type
            for center_name, center_coord in centers.items():
                # Distance from CA atom
                ca_distance = np.linalg.norm(ca_coord - center_coord)
                distance_metrics[f'ca_dist_{center_name}'] = ca_distance
                
                # Distance from residue center of mass
                com_distance = np.linalg.norm(residue_com - center_coord)
                distance_metrics[f'com_dist_{center_name}'] = com_distance
            
            distances[(chain_id, res_num)] = distance_metrics
        
        # Add radial statistics if requested
        if include_radial_stats:
            # Calculate protein radius and other global metrics
            for center_name, center_coord in centers.items():
                ca_distances = [np.linalg.norm(r['ca_coord'] - center_coord) for r in residue_data]
                
                protein_radius = np.max(ca_distances)
                mean_radius = np.mean(ca_distances)
                radius_std = np.std(ca_distances)
                
                # Add normalized distances and radial position metrics
                for residue in residue_data:
                    chain_id = residue['chain_id']
                    res_num = residue['res_num']
                    
                    ca_dist = distances[(chain_id, res_num)][f'ca_dist_{center_name}']
                    
                    # Normalized distance (0 = center, 1 = surface)
                    normalized_dist = ca_dist / protein_radius if protein_radius > 0 else 0.0
                    
                    # Radial percentile (where this residue falls in radial distribution)
                    radial_percentile = np.sum(np.array(ca_distances) <= ca_dist) / len(ca_distances) * 100
                    
                    # Core/intermediate/surface classification
                    if normalized_dist <= 0.33:
                        radial_class = 'core'
                    elif normalized_dist <= 0.66:
                        radial_class = 'intermediate'
                    else:
                        radial_class = 'surface'
                    
                    # Update distance metrics
                    distances[(chain_id, res_num)].update({
                        f'normalized_dist_{center_name}': normalized_dist,
                        f'radial_percentile_{center_name}': radial_percentile,
                        f'radial_class_{center_name}': radial_class,
                        f'protein_radius_{center_name}': protein_radius,
                        f'mean_radius_{center_name}': mean_radius,
                        f'radius_std_{center_name}': radius_std
                    })
        
        return distances
        
    except Exception as e:
        raise RuntimeError(f"Error processing {pdb_path.name}: {e}")


def process_pdb_batch(pdb_files: List[Union[str, Path]], 
                      center_types: List[str] = ['geometric', 'mass', 'ca_geometric'],
                      include_radial_stats: bool = True,
                      n_processes: Optional[int] = None) -> Dict[str, dict]:
    """
    Process PDB files in parallel to calculate distance metrics.
    
    Args:
        pdb_files: List of PDB file paths
        center_types: Types of centers to calculate
        include_radial_stats: Whether to include radial statistics
        n_processes: Number of parallel processes
        
    Returns:
        Nested dictionary mapping filename to distance data
    """
    n_processes = n_processes or min(cpu_count(), len(pdb_files))
    
    process_func = partial(calculate_distance_to_center, 
                          center_types=center_types,
                          include_radial_stats=include_radial_stats)
    
    with Pool(n_processes) as pool:
        results = {}
        for filename, data in zip(pdb_files, pool.imap(process_func, pdb_files)):
            results[Path(filename).name] = data
        return results


def save_distances_batch(batch_results: Dict[str, dict], 
                        output_file: str,
                        compression: str = 'snappy') -> None:
    """
    Save distance results to Parquet format.
    
    Args:
        batch_results: Dictionary of distance data
        output_file: Output Parquet path
        compression: Compression algorithm
    """
    rows = []
    for pdb_file, residues in batch_results.items():
        pdb_id = Path(pdb_file).stem.split('.')[0]
        
        for (chain, res_num), data in residues.items():
            # Create base row
            row = {
                'pdb_id': pdb_id,
                'chain_id': chain,
                'residue_num': res_num,
                'residue_name': data['residue_name'],
                'num_heavy_atoms': data['num_heavy_atoms']
            }
            
            # Add all distance metrics
            for key, value in data.items():
                if key not in ['residue_name', 'num_heavy_atoms']:
                    row[key] = value
            
            rows.append(row)
    
    table = pa.Table.from_pandas(pd.DataFrame(rows))
    pq.write_table(table, output_file, compression=compression)


def process_and_save_batches(pdb_files: List[str], 
                           output_dir: str,
                           batch_size: int = 1000,
                           center_types: List[str] = ['geometric', 'mass', 'ca_geometric'],
                           n_processes: int = None) -> List[str]:
    """
    Process PDB files in memory-efficient batches.
    
    Args:
        pdb_files: List of input PDB paths
        output_dir: Output directory for Parquet files
        batch_size: Number of files per batch
        center_types: Types of centers to calculate
        n_processes: Parallel processes per batch
        
    Returns:
        List of created Parquet files
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    parquet_files = []
    
    for i in range(0, len(pdb_files), batch_size):
        batch = pdb_files[i:i+batch_size]
        batch_num = i // batch_size + 1
        
        print(f"Processing batch {batch_num} ({len(batch)} files)")
        results = process_pdb_batch(batch, 
                                   center_types=center_types,
                                   n_processes=n_processes)
        
        output_path = Path(output_dir) / f"distance_to_center_batch_{batch_num:04d}.parquet"
        save_distances_batch(results, str(output_path))
        
        parquet_files.append(str(output_path))
        print(f"Saved batch {batch_num} to {output_path}")
    
    return parquet_files


def query_distances(parquet_files: Union[str, List[str]],
                   pdb_ids: List[str] = None,
                   chains: List[str] = None,
                   radial_class: str = None,
                   min_distance: float = None,
                   max_distance: float = None,
                   distance_type: str = 'ca_dist_geometric') -> pd.DataFrame:
    """
    Query distance data with various filters.
    
    Args:
        parquet_files: Parquet file(s) to query
        pdb_ids: Filter by PDB IDs
        chains: Filter by chain IDs
        radial_class: Filter by radial class ('core', 'intermediate', 'surface')
        min_distance: Minimum distance
        max_distance: Maximum distance
        distance_type: Which distance metric to filter on
        
    Returns:
        Filtered DataFrame of distance data
    """
    if isinstance(parquet_files, str):
        parquet_files = [parquet_files]
    
    filters = []
    if pdb_ids:
        filters.append(('pdb_id', 'in', pdb_ids))
    if chains:
        filters.append(('chain_id', 'in', chains))
    if radial_class:
        filters.append(('radial_class_geometric', '=', radial_class))
    if min_distance is not None:
        filters.append((distance_type, '>=', min_distance))
    if max_distance is not None:
        filters.append((distance_type, '<=', max_distance))
    
    return pq.read_table(
        parquet_files,
        filters=filters or None
    ).to_pandas()


def get_core_residues(parquet_files: Union[str, List[str]],
                     percentile_threshold: float = 25,
                     center_type: str = 'geometric') -> pd.DataFrame:
    """
    Find core residues based on distance to protein center.
    
    Args:
        parquet_files: Parquet file(s) to query
        percentile_threshold: Percentile threshold for "core" residues
        center_type: Which center definition to use
        
    Returns:
        DataFrame of core residues
    """
    distance_col = f'ca_dist_{center_type}'
    
    # Get all data to calculate percentile
    all_data = query_distances(parquet_files)
    threshold = np.nanpercentile(all_data[distance_col], percentile_threshold)
    
    return query_distances(parquet_files, max_distance=threshold, distance_type=distance_col)


def get_surface_residues(parquet_files: Union[str, List[str]],
                        percentile_threshold: float = 75,
                        center_type: str = 'geometric') -> pd.DataFrame:
    """
    Find surface residues based on distance to protein center.
    
    Args:
        parquet_files: Parquet file(s) to query
        percentile_threshold: Percentile threshold for "surface" residues
        center_type: Which center definition to use
        
    Returns:
        DataFrame of surface residues
    """
    distance_col = f'ca_dist_{center_type}'
    
    all_data = query_distances(parquet_files)
    threshold = np.nanpercentile(all_data[distance_col], percentile_threshold)
    
    return query_distances(parquet_files, min_distance=threshold, distance_type=distance_col)


## Example usage
if __name__ == "__main__":
    # Configuration
    input_dir = "../../1_inputs/raw_data/UP000005640_9606_HUMAN_v4"
    output_dir = "../../3_outputs/processed_data/distance_to_center"
    
    # Get all compressed PDB files
    pdb_files = list(Path(input_dir).glob("*.pdb*"))
    
    # Process and save data
    parquet_files = process_and_save_batches(
        pdb_files, 
        output_dir,
        center_types=['geometric', 'mass', 'ca_geometric']
    )
    
    # Example queries
    # 1. Get distance data for a specific protein
    # df = query_distances(parquet_files, pdb_ids=["AF-P00533-F1"])
    
    # 2. Find core residues (closest 25% to center)
    # core = get_core_residues(parquet_files, percentile_threshold=25)
    
    # 3. Find surface residues (farthest 25% from center)
    # surface = get_surface_residues(parquet_files, percentile_threshold=75)
    
    # 4. Get residues in intermediate region
    # intermediate = query_distances(parquet_files, radial_class='intermediate')