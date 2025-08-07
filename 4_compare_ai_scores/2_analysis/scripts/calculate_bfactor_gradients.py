"""
B-factor Gradients Storage and Query System
==========================================

This module calculates B-factor (pLDDT) gradients and local structural flexibility
metrics from AlphaFold-generated PDB files, providing insights into local structural
rigidity and mutation sensitivity.

B-factor gradients indicate how much the local confidence changes around each residue,
which can be indicative of:
- Structural flexibility/rigidity transitions
- Domain boundaries
- Functionally important regions with varying confidence
- Regions where mutations might have cascading structural effects

Key features:
- Parallel processing of PDB files using multiprocessing
- Multiple gradient metrics (absolute, relative, max neighborhood change)
- Efficient storage using Apache Parquet format
- Fast query capabilities with predicate pushdown
- Batch processing to handle large datasets with minimal memory overhead

Author: Runxi Shen
Created: 2025-01-XX
"""

import os
import gzip
import warnings
import tempfile
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


def calculate_bfactor_gradients(pdb_file: Union[str, Path], 
                               window_size: int = 3,
                               gradient_types: List[str] = ['absolute', 'relative', 'max_change']) -> Dict[Tuple[str, int], dict]:
    """
    Calculate B-factor (pLDDT) gradients for residues in a PDB file.
    
    Args:
        pdb_file: Path to PDB file (can be gzipped)
        window_size: Size of neighborhood window for gradient calculation
        gradient_types: Types of gradients to calculate
        
    Returns:
        Dictionary mapping (chain_id, residue_number) to gradient metrics
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
        
        # Extract residue data
        residue_data = []
        for residue in structure.get_residues():
            try:
                chain_id = residue.get_parent().id
                res_num = residue.id[1]
                bfactor = residue["CA"].bfactor  # pLDDT score
                
                residue_data.append({
                    'chain_id': chain_id,
                    'res_num': res_num,
                    'residue_name': residue.resname,
                    'bfactor': bfactor
                })
            except KeyError:
                continue
        
        # Sort by chain and residue number for proper gradient calculation
        residue_data.sort(key=lambda x: (x['chain_id'], x['res_num']))
        
        # Calculate gradients
        gradients = {}
        half_window = window_size // 2
        
        for i, residue in enumerate(residue_data):
            chain_id = residue['chain_id']
            res_num = residue['res_num']
            current_bfactor = residue['bfactor']
            
            # Get neighboring residues within the same chain
            neighbors = []
            for j in range(max(0, i - half_window), min(len(residue_data), i + half_window + 1)):
                if residue_data[j]['chain_id'] == chain_id:
                    neighbors.append(residue_data[j]['bfactor'])
            
            if len(neighbors) < 2:  # Need at least 2 points for gradient
                continue
                
            gradient_metrics = {}
            
            # Calculate different gradient types
            if 'absolute' in gradient_types:
                # Absolute gradient: standard deviation of neighborhood
                gradient_metrics['absolute_gradient'] = np.std(neighbors)
            
            if 'relative' in gradient_types:
                # Relative gradient: coefficient of variation
                mean_bfactor = np.mean(neighbors)
                if mean_bfactor > 0:
                    gradient_metrics['relative_gradient'] = np.std(neighbors) / mean_bfactor
                else:
                    gradient_metrics['relative_gradient'] = 0.0
            
            if 'max_change' in gradient_types:
                # Maximum change from current residue to any neighbor
                max_diff = max(abs(current_bfactor - neighbor) for neighbor in neighbors)
                gradient_metrics['max_change'] = max_diff
            
            # Additional metrics
            gradient_metrics.update({
                'neighborhood_mean': np.mean(neighbors),
                'neighborhood_min': np.min(neighbors),
                'neighborhood_max': np.max(neighbors),
                'residue_name': residue['residue_name'],
                'window_size': len(neighbors)
            })
            
            gradients[(chain_id, res_num)] = gradient_metrics
            
        return gradients
        
    except Exception as e:
        raise RuntimeError(f"Error processing {pdb_path.name}: {e}")


def process_pdb_batch(pdb_files: List[Union[str, Path]], 
                      window_size: int = 3,
                      gradient_types: List[str] = ['absolute', 'relative', 'max_change'],
                      n_processes: Optional[int] = None) -> Dict[str, dict]:
    """
    Process PDB files in parallel to calculate B-factor gradients.
    
    Args:
        pdb_files: List of PDB file paths
        window_size: Size of neighborhood window for gradient calculation
        gradient_types: Types of gradients to calculate
        n_processes: Number of parallel processes
        
    Returns:
        Nested dictionary mapping filename to gradient data
    """
    n_processes = n_processes or min(cpu_count(), len(pdb_files))
    
    # Create partial function with fixed parameters
    process_func = partial(calculate_bfactor_gradients, 
                          window_size=window_size, 
                          gradient_types=gradient_types)
    
    with Pool(n_processes) as pool:
        results = {}
        for filename, data in zip(pdb_files, pool.imap(process_func, pdb_files)):
            results[Path(filename).name] = data
        return results


def save_gradients_batch(batch_results: Dict[str, dict], 
                        output_file: str,
                        compression: str = 'snappy') -> None:
    """
    Save B-factor gradient results to Parquet format.
    
    Args:
        batch_results: Dictionary of gradient data
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
                'absolute_gradient': data.get('absolute_gradient', np.nan),
                'relative_gradient': data.get('relative_gradient', np.nan),
                'max_change': data.get('max_change', np.nan),
                'neighborhood_mean': data['neighborhood_mean'],
                'neighborhood_min': data['neighborhood_min'],
                'neighborhood_max': data['neighborhood_max'],
                'window_size': data['window_size']
            }
            rows.append(row)
    
    table = pa.Table.from_pandas(pd.DataFrame(rows))
    pq.write_table(table, output_file, compression=compression)


def process_and_save_batches(pdb_files: List[str], 
                           output_dir: str,
                           batch_size: int = 1000,
                           window_size: int = 3,
                           n_processes: int = None) -> List[str]:
    """
    Process PDB files in memory-efficient batches.
    
    Args:
        pdb_files: List of input PDB paths
        output_dir: Output directory for Parquet files
        batch_size: Number of files per batch
        window_size: Size of neighborhood window
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
        results = process_pdb_batch(batch, window_size=window_size, n_processes=n_processes)
        
        output_path = Path(output_dir) / f"bfactor_gradients_batch_{batch_num:04d}.parquet"
        save_gradients_batch(results, str(output_path))
        
        parquet_files.append(str(output_path))
        print(f"Saved batch {batch_num} to {output_path}")
    
    return parquet_files


def query_gradients(parquet_files: Union[str, List[str]],
                   pdb_ids: List[str] = None,
                   chains: List[str] = None,
                   min_gradient: float = None,
                   max_gradient: float = None,
                   gradient_type: str = 'absolute_gradient') -> pd.DataFrame:
    """
    Query B-factor gradient data with various filters.
    
    Args:
        parquet_files: Parquet file(s) to query
        pdb_ids: Filter by PDB IDs
        chains: Filter by chain IDs
        min_gradient: Minimum gradient value
        max_gradient: Maximum gradient value
        gradient_type: Which gradient metric to filter on
        
    Returns:
        Filtered DataFrame of gradient data
    """
    if isinstance(parquet_files, str):
        parquet_files = [parquet_files]
    
    filters = []
    if pdb_ids:
        filters.append(('pdb_id', 'in', pdb_ids))
    if chains:
        filters.append(('chain_id', 'in', chains))
    if min_gradient is not None:
        filters.append((gradient_type, '>=', min_gradient))
    if max_gradient is not None:
        filters.append((gradient_type, '<=', max_gradient))
    
    return pq.read_table(
        parquet_files,
        filters=filters or None
    ).to_pandas()


def get_high_gradient_residues(parquet_files: Union[str, List[str]],
                              threshold_percentile: float = 90,
                              gradient_type: str = 'absolute_gradient') -> pd.DataFrame:
    """
    Find residues with high gradient values (potential flexibility transition points).
    
    Args:
        parquet_files: Parquet file(s) to query
        threshold_percentile: Percentile threshold for "high" gradients
        gradient_type: Which gradient metric to use
        
    Returns:
        DataFrame of high-gradient residues
    """
    # First, get all data to calculate percentile
    all_data = query_gradients(parquet_files)
    threshold = np.nanpercentile(all_data[gradient_type], threshold_percentile)
    
    # Filter for high-gradient residues
    return query_gradients(parquet_files, min_gradient=threshold, gradient_type=gradient_type)


## Example usage
if __name__ == "__main__":
    # Configuration
    input_dir = "../../1_inputs/raw_data/UP000005640_9606_HUMAN_v4"
    output_dir = "../../3_outputs/processed_data/bfactor_gradients"
    
    # Get all compressed PDB files
    pdb_files = list(Path(input_dir).glob("*.pdb*"))
    
    # Process and save data
    parquet_files = process_and_save_batches(
        pdb_files, 
        output_dir,
        window_size=5  # Use 5-residue window for gradient calculation
    )
    
    # Example queries
    # 1. Get gradient data for a specific protein
    # df = query_gradients(parquet_files, pdb_ids=["AF-P00533-F1"])
    
    # 2. Find high-flexibility transition points
    # high_gradients = get_high_gradient_residues(parquet_files, threshold_percentile=95)
    
    # 3. Get residues with high relative gradients (normalized by mean confidence)
    # relative_gradients = query_gradients(parquet_files, min_gradient=0.5, gradient_type='relative_gradient')