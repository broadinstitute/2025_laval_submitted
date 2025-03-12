"""
pLDDT Storage and Query System
==============================

This module provides a solution for extracting, storing, and querying
pLDDT confidence scores from AlphaFold-generated PDB files.

Features:
- Parallel processing of compressed PDB files
- Column-oriented Parquet storage for fast access
- Query capabilities with predicate pushdown
- Batch processing optimized for large datasets

Author: Runxi Shen
Created: March 11, 2025
"""

import os
import gzip
import warnings
from pathlib import Path
from functools import partial
from multiprocessing import Pool, cpu_count
from typing import List, Dict, Tuple, Optional, Union
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from Bio.PDB import PDBParser
from io import StringIO


def extract_plddt_scores(pdb_file: Union[str, Path]) -> Dict[Tuple[str, int], float]:
    """
    Extract pLDDT scores from AlphaFold-generated PDB file.
    
    Args:
        pdb_file: Path to PDB file (.pdb or .pdb.gz)
        
    Returns:
        Dictionary mapping (chain_id, residue_number) to pLDDT scores
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
        
        plddt_scores = {}
        for residue in structure.get_residues():
            try:
                # Extract chain, residue number, and pLDDT (B-factor of CA atom)
                chain_id = residue.get_parent().id
                res_num = residue.id[1]
                plddt = residue["CA"].bfactor
                
                plddt_scores[(chain_id, res_num)] = {
                    'plddt': plddt,
                    'residue_name': residue.resname
                }
            except KeyError:
                continue  # Skip residues without CA atom
                
        return plddt_scores
        
    except Exception as e:
        raise RuntimeError(f"Error processing {pdb_path.name}: {e}")


def process_pdb_batch(pdb_files: List[Union[str, Path]], 
                      n_processes: Optional[int] = None) -> Dict[str, dict]:
    """
    Process PDB files in parallel to extract pLDDT scores.
    
    Args:
        pdb_files: List of PDB file paths
        n_processes: Number of parallel processes
        
    Returns:
        Nested dictionary mapping filename to pLDDT data
    """
    n_processes = n_processes or min(cpu_count(), len(pdb_files))
    
    with Pool(n_processes) as pool:
        results = {}
        for filename, data in zip(pdb_files, pool.imap(extract_plddt_scores, pdb_files)):
            results[Path(filename).name] = data
        return results


def save_plddt_batch(batch_results: Dict[str, dict], 
                    output_file: str,
                    compression: str = 'snappy') -> None:
    """
    Save pLDDT results to Parquet format.
    
    Args:
        batch_results: Dictionary of pLDDT data
        output_file: Output Parquet path
        compression: Compression algorithm
    """
    rows = []
    for pdb_file, residues in batch_results.items():
        pdb_id = Path(pdb_file).stem.split('.')[0]
        
        for (chain, res_num), data in residues.items():
            rows.append({
                'pdb_id': pdb_id,
                'chain_id': chain,
                'residue_num': res_num,
                'residue_name': data['residue_name'],
                'plddt': data['plddt']
            })
    
    table = pa.Table.from_pandas(pd.DataFrame(rows))
    pq.write_table(table, output_file, compression=compression)


def process_and_save_batches(pdb_files: List[str], 
                           output_dir: str,
                           batch_size: int = 1000,
                           n_processes: int = None) -> List[str]:
    """
    Process PDB files in memory-efficient batches.
    
    Args:
        pdb_files: List of input PDB paths
        output_dir: Output directory for Parquet files
        batch_size: Number of files per batch
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
        results = process_pdb_batch(batch, n_processes)
        
        output_path = Path(output_dir) / f"plddt_batch_{batch_num:04d}.parquet"
        save_plddt_batch(results, str(output_path))
        
        parquet_files.append(str(output_path))
        print(f"Saved batch {batch_num} to {output_path}")
    
    return parquet_files


def query_plddt(parquet_files: Union[str, List[str]],
               pdb_ids: List[str] = None,
               chains: List[str] = None,
               min_score: float = None,
               max_score: float = None) -> pd.DataFrame:
    """
    Query pLDDT data with various filters.
    
    Args:
        parquet_files: Parquet file(s) to query
        pdb_ids: Filter by PDB IDs
        chains: Filter by chain IDs
        min_score: Minimum pLDDT score
        max_score: Maximum pLDDT score
        
    Returns:
        Filtered DataFrame of pLDDT scores
    """
    if isinstance(parquet_files, str):
        parquet_files = [parquet_files]
    
    filters = []
    if pdb_ids:
        filters.append(('pdb_id', 'in', pdb_ids))
    if chains:
        filters.append(('chain_id', 'in', chains))
    if min_score is not None:
        filters.append(('plddt', '>=', min_score))
    if max_score is not None:
        filters.append(('plddt', '<=', max_score))
    
    return pq.read_table(
        parquet_files,
        filters=filters or None
    ).to_pandas()


## Example usage
if __name__ == "__main__":
    # Configuration
    input_dir = "../../1_inputs/raw_data/UP000005640_9606_HUMAN_v4"
    output_dir = "../../3_outputs/processed_data/plddt_scores"
    
    # Get all compressed PDB files
    pdb_files = list(Path(input_dir).glob("*.pdb*"))
    
    # Process and save data
    parquet_files = process_and_save_batches(pdb_files, output_dir)
    
    # Example queries
    # 1. Get all scores for a specific protein
    # df = query_plddt(parquet_files, pdb_ids=["AF-P00533-F1"])
    
    # 2. Find low-confidence regions across all proteins
    # low_conf = query_plddt(parquet_files, min_score=0, max_score=50)
    
    # 3. Get scores for specific chains
    # chain_b = query_plddt(parquet_files, chains=["B"])