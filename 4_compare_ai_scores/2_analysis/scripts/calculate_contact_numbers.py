"""
Residue Contact Number Storage and Query System
==============================================

This module calculates residue contact numbers and local packing density
from protein structure PDB files, providing insights into:
- Local structural constraints
- Packing density around each residue
- Burial vs. surface exposure (complementing RSA)
- Structural stability contributions

Contact number is the count of residues within a specified distance threshold,
indicating how tightly packed each residue is in its local environment.
Higher contact numbers typically indicate:
- More structural constraints
- Higher mutation sensitivity due to more interactions
- Buried/core residues vs. surface residues

Key features:
- Multiple distance cutoffs (6Å, 8Å, 10Å for different interaction types)
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


def calculate_contact_numbers(pdb_file: Union[str, Path], 
                             cutoff_distances: List[float] = [6.0, 8.0, 10.0],
                             atom_type: str = 'CA') -> Dict[Tuple[str, int], dict]:
    """
    Calculate residue contact numbers for different distance cutoffs.
    
    Args:
        pdb_file: Path to PDB file (can be gzipped)
        cutoff_distances: List of distance cutoffs in Angstroms
        atom_type: Which atom to use for distance calculations ('CA', 'CB', 'centroid')
        
    Returns:
        Dictionary mapping (chain_id, residue_number) to contact metrics
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
        
        # Extract residue coordinates and metadata
        residue_data = []
        for residue in structure.get_residues():
            try:
                chain_id = residue.get_parent().id
                res_num = residue.id[1]
                res_name = residue.resname
                
                # Get representative coordinate based on atom_type
                if atom_type == 'CA':
                    coord = residue["CA"].coord
                elif atom_type == 'CB':
                    # Use CB if available, otherwise CA (for Glycine)
                    if 'CB' in residue:
                        coord = residue["CB"].coord
                    else:
                        coord = residue["CA"].coord
                elif atom_type == 'centroid':
                    # Calculate centroid of all heavy atoms
                    coords = [atom.coord for atom in residue if atom.element != 'H']
                    coord = np.mean(coords, axis=0) if coords else residue["CA"].coord
                else:
                    coord = residue["CA"].coord
                
                residue_data.append({
                    'chain_id': chain_id,
                    'res_num': res_num,
                    'residue_name': res_name,
                    'coord': coord
                })
                
            except KeyError:
                continue
        
        # Calculate contact numbers
        contacts = {}
        
        for i, residue_i in enumerate(residue_data):
            chain_id = residue_i['chain_id']
            res_num = residue_i['res_num']
            coord_i = residue_i['coord']
            
            contact_metrics = {}
            
            # Calculate contacts for each distance cutoff
            for cutoff in cutoff_distances:
                contact_count = 0
                closest_distance = float('inf')
                distances = []
                
                for j, residue_j in enumerate(residue_data):
                    if i == j:  # Skip self
                        continue
                    
                    coord_j = residue_j['coord']
                    distance = np.linalg.norm(coord_i - coord_j)
                    distances.append(distance)
                    
                    if distance <= cutoff:
                        contact_count += 1
                    
                    if distance < closest_distance:
                        closest_distance = distance
                
                cutoff_key = f'contacts_{cutoff:.0f}A'
                contact_metrics[cutoff_key] = contact_count
            
            # Additional metrics
            if distances:
                contact_metrics.update({
                    'closest_distance': closest_distance,
                    'mean_distance': np.mean(distances),
                    'median_distance': np.median(distances),
                    'distances_std': np.std(distances),
                    'total_residues': len(residue_data),
                    'residue_name': residue_i['residue_name']
                })
            
            # Calculate packing density (contacts normalized by total residues)
            for cutoff in cutoff_distances:
                cutoff_key = f'contacts_{cutoff:.0f}A'
                density_key = f'density_{cutoff:.0f}A'
                if contact_metrics['total_residues'] > 1:
                    contact_metrics[density_key] = contact_metrics[cutoff_key] / (contact_metrics['total_residues'] - 1)
                else:
                    contact_metrics[density_key] = 0.0
            
            contacts[(chain_id, res_num)] = contact_metrics
            
        return contacts
        
    except Exception as e:
        raise RuntimeError(f"Error processing {pdb_path.name}: {e}")


def process_pdb_batch(pdb_files: List[Union[str, Path]], 
                      cutoff_distances: List[float] = [6.0, 8.0, 10.0],
                      atom_type: str = 'CA',
                      n_processes: Optional[int] = None) -> Dict[str, dict]:
    """
    Process PDB files in parallel to calculate contact numbers.
    
    Args:
        pdb_files: List of PDB file paths
        cutoff_distances: List of distance cutoffs
        atom_type: Which atom to use for calculations
        n_processes: Number of parallel processes
        
    Returns:
        Nested dictionary mapping filename to contact data
    """
    n_processes = n_processes or min(cpu_count(), len(pdb_files))
    
    process_func = partial(calculate_contact_numbers, 
                          cutoff_distances=cutoff_distances, 
                          atom_type=atom_type)
    
    with Pool(n_processes) as pool:
        results = {}
        for filename, data in zip(pdb_files, pool.imap(process_func, pdb_files)):
            results[Path(filename).name] = data
        return results


def save_contacts_batch(batch_results: Dict[str, dict], 
                       output_file: str,
                       compression: str = 'snappy') -> None:
    """
    Save contact number results to Parquet format.
    
    Args:
        batch_results: Dictionary of contact data
        output_file: Output Parquet path
        compression: Compression algorithm
    """
    rows = []
    for pdb_file, residues in batch_results.items():
        pdb_id = Path(pdb_file).stem.split('.')[0]
        
        for (chain, res_num), data in residues.items():
            row = {
                'pdb_id': pdb_id,
                'chain_id': chain,
                'residue_num': res_num,
                'residue_name': data['residue_name'],
                'contacts_6A': data.get('contacts_6A', 0),
                'contacts_8A': data.get('contacts_8A', 0),
                'contacts_10A': data.get('contacts_10A', 0),
                'density_6A': data.get('density_6A', 0.0),
                'density_8A': data.get('density_8A', 0.0),
                'density_10A': data.get('density_10A', 0.0),
                'closest_distance': data.get('closest_distance', np.nan),
                'mean_distance': data.get('mean_distance', np.nan),
                'median_distance': data.get('median_distance', np.nan),
                'distances_std': data.get('distances_std', np.nan),
                'total_residues': data.get('total_residues', 0)
            }
            rows.append(row)
    
    table = pa.Table.from_pandas(pd.DataFrame(rows))
    pq.write_table(table, output_file, compression=compression)


def process_and_save_batches(pdb_files: List[str], 
                           output_dir: str,
                           batch_size: int = 1000,
                           cutoff_distances: List[float] = [6.0, 8.0, 10.0],
                           atom_type: str = 'CA',
                           n_processes: int = None) -> List[str]:
    """
    Process PDB files in memory-efficient batches.
    
    Args:
        pdb_files: List of input PDB paths
        output_dir: Output directory for Parquet files
        batch_size: Number of files per batch
        cutoff_distances: Distance cutoffs for contact calculation
        atom_type: Which atom to use for calculations
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
                                   cutoff_distances=cutoff_distances, 
                                   atom_type=atom_type,
                                   n_processes=n_processes)
        
        output_path = Path(output_dir) / f"contact_numbers_batch_{batch_num:04d}.parquet"
        save_contacts_batch(results, str(output_path))
        
        parquet_files.append(str(output_path))
        print(f"Saved batch {batch_num} to {output_path}")
    
    return parquet_files


def query_contacts(parquet_files: Union[str, List[str]],
                  pdb_ids: List[str] = None,
                  chains: List[str] = None,
                  min_contacts: int = None,
                  max_contacts: int = None,
                  contact_type: str = 'contacts_8A') -> pd.DataFrame:
    """
    Query contact number data with various filters.
    
    Args:
        parquet_files: Parquet file(s) to query
        pdb_ids: Filter by PDB IDs
        chains: Filter by chain IDs
        min_contacts: Minimum contact number
        max_contacts: Maximum contact number
        contact_type: Which contact metric to filter on
        
    Returns:
        Filtered DataFrame of contact data
    """
    if isinstance(parquet_files, str):
        parquet_files = [parquet_files]
    
    filters = []
    if pdb_ids:
        filters.append(('pdb_id', 'in', pdb_ids))
    if chains:
        filters.append(('chain_id', 'in', chains))
    if min_contacts is not None:
        filters.append((contact_type, '>=', min_contacts))
    if max_contacts is not None:
        filters.append((contact_type, '<=', max_contacts))
    
    return pq.read_table(
        parquet_files,
        filters=filters or None
    ).to_pandas()


def get_buried_residues(parquet_files: Union[str, List[str]],
                       threshold_percentile: float = 75,
                       contact_type: str = 'contacts_8A') -> pd.DataFrame:
    """
    Find highly buried residues based on contact number.
    
    Args:
        parquet_files: Parquet file(s) to query
        threshold_percentile: Percentile threshold for "buried" residues
        contact_type: Which contact metric to use
        
    Returns:
        DataFrame of buried residues
    """
    # Get all data to calculate percentile
    all_data = query_contacts(parquet_files)
    threshold = np.nanpercentile(all_data[contact_type], threshold_percentile)
    
    # Filter for highly contacted (buried) residues
    return query_contacts(parquet_files, min_contacts=int(threshold), contact_type=contact_type)


def get_surface_residues(parquet_files: Union[str, List[str]],
                        threshold_percentile: float = 25,
                        contact_type: str = 'contacts_8A') -> pd.DataFrame:
    """
    Find surface residues based on low contact numbers.
    
    Args:
        parquet_files: Parquet file(s) to query
        threshold_percentile: Percentile threshold for "surface" residues
        contact_type: Which contact metric to use
        
    Returns:
        DataFrame of surface residues
    """
    all_data = query_contacts(parquet_files)
    threshold = np.nanpercentile(all_data[contact_type], threshold_percentile)
    
    return query_contacts(parquet_files, max_contacts=int(threshold), contact_type=contact_type)


## Example usage
if __name__ == "__main__":
    # Configuration
    input_dir = "../../1_inputs/raw_data/UP000005640_9606_HUMAN_v4"
    output_dir = "../../3_outputs/processed_data/contact_numbers"
    
    # Get all compressed PDB files
    pdb_files = list(Path(input_dir).glob("*.pdb*"))
    
    # Process and save data
    parquet_files = process_and_save_batches(
        pdb_files, 
        output_dir,
        cutoff_distances=[6.0, 8.0, 10.0],  # Different interaction distances
        atom_type='CA'  # Use CA atoms for distance calculations
    )
    
    # Example queries
    # 1. Get contact data for a specific protein
    # df = query_contacts(parquet_files, pdb_ids=["AF-P00533-F1"])
    
    # 2. Find highly buried residues (high contact numbers)
    # buried = get_buried_residues(parquet_files, threshold_percentile=90)
    
    # 3. Find surface residues (low contact numbers)
    # surface = get_surface_residues(parquet_files, threshold_percentile=10)
    
    # 4. Get residues with specific contact ranges
    # medium_contacts = query_contacts(parquet_files, min_contacts=5, max_contacts=15, contact_type='contacts_8A')