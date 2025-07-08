"""
RSA Storage and Query System
============================

This module provides a comprehensive solution for calculating, storing, and querying
Relative Solvent Accessibility (RSA) data from protein structure PDB files.

The system is designed to efficiently process thousands of protein structures in parallel
and store the per-residue RSA values in an optimized column-oriented format (Parquet)
for fast querying and analysis.

Key features:
- Parallel processing of PDB files using multiprocessing
- Efficient storage using Apache Parquet format
- Fast query capabilities with predicate pushdown
- Batch processing to handle large datasets with minimal memory overhead
- Simple API for common access patterns

Author: Runxi Shen
Created: March 11, 2025
Last modified: March 11, 2025
"""

import os
import gzip
import warnings
import tempfile
import glob
from pathlib import Path
from functools import partial
from multiprocessing import Pool, cpu_count
from typing import List, Dict, Tuple, Optional, Union
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from Bio.PDB import PDBParser, DSSP, PDBIO

aa_single_to_three = {
    'A': 'Ala',  # Alanine
    'R': 'Arg',  # Arginine
    'N': 'Asn',  # Asparagine
    'D': 'Asp',  # Aspartic acid
    'C': 'Cys',  # Cysteine
    'Q': 'Gln',  # Glutamine
    'E': 'Glu',  # Glutamic acid
    'G': 'Gly',  # Glycine
    'H': 'His',  # Histidine
    'I': 'Ile',  # Isoleucine
    'L': 'Leu',  # Leucine
    'K': 'Lys',  # Lysine
    'M': 'Met',  # Methionine
    'F': 'Phe',  # Phenylalanine
    'P': 'Pro',  # Proline
    'S': 'Ser',  # Serine
    'T': 'Thr',  # Threonine
    'W': 'Trp',  # Tryptophan
    'Y': 'Tyr',  # Tyrosine
    'V': 'Val'   # Valine
}


def calculate_rsa(pdb_file: Union[str, Path], 
                 dssp_executable: str = 'mkdssp',
                 chain_id: Optional[str] = None) -> Dict[Tuple[str, int], float]:
    """
    Calculate relative solvent accessibility (RSA) for residues in a PDB file.
    
    Args:
        pdb_file (str or Path): Path to PDB file (can be gzipped)
        dssp_executable (str): Path to DSSP executable (default: 'mkdssp')
        chain_id (str, optional): If provided, only calculate RSA for this chain
        
    Returns:
        Dictionary mapping (chain_id, residue_number) to RSA values (0-100%)
    """
    # Create absolute path for the PDB file
    pdb_path = Path(pdb_file).resolve()
    
    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")
        
    # Determine if the file is gzipped
    is_gzipped = str(pdb_path).endswith('.gz')
    
    # Create a parser
    parser = PDBParser(QUIET=True)
    
    # Open and parse the PDB file
    try:
        if is_gzipped:
            with gzip.open(pdb_path, 'rt') as f:
                structure = parser.get_structure("protein", f)
        else:
            structure = parser.get_structure("protein", str(pdb_path))
    except Exception as e:
        raise ValueError(f"Error parsing PDB file: {e}")
    
    # Get the first model
    model = structure[0]
    
    # Create a temporary non-gzipped PDB file if needed
    temp_pdb_name = None
    try:
        if is_gzipped:
            with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as temp_pdb:
                temp_pdb_name = temp_pdb.name
                io = PDBIO()
                io.set_structure(structure)
                io.save(temp_pdb_name)
            pdb_file_for_dssp = temp_pdb_name
        else:
            pdb_file_for_dssp = str(pdb_path)
        
        # Suppress warnings during DSSP execution
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            dssp = DSSP(model, pdb_file_for_dssp, dssp=dssp_executable, file_type='PDB')
        
        # Extract RSA values
        rsa_values = {}
        for key in dssp.keys():
            dssp_chain_id = key[0]
            
            # Skip if we're only looking at a specific chain
            if chain_id is not None and dssp_chain_id != chain_id:
                continue
                
            res_id = key[1][1]  # Residue number
            dssp_data = dssp[key]
            res_name = dssp_data[1]
            rel_acc = dssp_data[3]  # Convert from fraction to percentage * 100 if percentage
            
            # Store the RSA value
            rsa_values[(dssp_chain_id, res_id, res_name)] = rel_acc
            
        return rsa_values
        
    except Exception as e:
        raise RuntimeError(f"Error running DSSP: {e}")
    
    finally:
        # Clean up the temporary file if it exists
        if temp_pdb_name and os.path.exists(temp_pdb_name):
            os.remove(temp_pdb_name)


def process_pdb_batch(pdb_files: List[Union[str, Path]], 
                      dssp_executable: str = 'mkdssp', 
                      chain_id: Optional[str] = None,
                      n_processes: Optional[int] = None) -> Dict[str, Dict[Tuple[str, int], float]]:
    """
    Process a batch of PDB files in parallel to calculate RSA values.
    
    Args:
        pdb_files (list): List of PDB file paths
        dssp_executable (str): Path to DSSP executable
        chain_id (str, optional): If provided, only calculate RSA for this chain
        n_processes (int, optional): Number of processes to use. Defaults to CPU count.
        
    Returns:
        Dictionary mapping PDB filename to RSA dictionaries
    """
    # Determine number of processes to use
    if n_processes is None:
        n_processes = min(cpu_count(), len(pdb_files))
    
    # Create a partial function with the fixed parameters
    process_func = partial(calculate_rsa, dssp_executable=dssp_executable, chain_id=chain_id)
    
    # Create a dictionary to store results
    results = {}
    
    # Process files in parallel
    with Pool(processes=n_processes) as pool:
        # Create a list of (filename, result) tuples
        for pdb_file, rsa_values in zip(pdb_files, pool.map(process_func, pdb_files)):
            results[os.path.basename(pdb_file)] = rsa_values
    
    return results


def save_rsa_batch(batch_results: Dict[str, Dict[Tuple[str, int], float]], 
                   output_file: str,
                #    include_residue_names: bool = True,
                   compression: str = 'snappy') -> None:
    """
    Save RSA results from multiple PDB files into a single Parquet file.
    
    Args:
        batch_results: Dictionary mapping PDB filenames to RSA dictionaries
        output_file: Path to output Parquet file
        include_residue_names: Whether to include residue names (requires additional PDB parsing)
        model_type: Type of model (for metadata)
        compression: Compression algorithm for Parquet ('snappy', 'gzip', 'brotli', etc.)
    """
    # Prepare data for DataFrame
    data = []
    
    for pdb_filename, rsa_values in batch_results.items():
        # Extract PDB ID from filename (remove path and extension)
        pdb_id = Path(pdb_filename).stem
        if pdb_id.endswith('.pdb'):
            pdb_id = pdb_id[:-4]
            
        # Add each residue's data
        for (chain_id, res_num, res_name), rsa in rsa_values.items():
            row = {
                'pdb_id': pdb_id,
                'chain_id': chain_id,
                'residue_num': res_num,
                'residue_name': aa_single_to_three[res_name].upper(),
                'rsa': float(rsa)
            }
            
            # # Add residue name if available and requested
            # if include_residue_names and hasattr(batch_results, 'residue_names'):
            #     res_name = batch_results.residue_names.get((pdb_id, chain_id, res_num), 'UNK')
            #     row['residue_name'] = res_name
            data.append(row)
    
    # Create DataFrame
    df = pd.DataFrame(data)
    
    # Save as Parquet with appropriate options
    table = pa.Table.from_pandas(df)
    pq.write_table(
        table, 
        output_file, 
        compression=compression
    )


def process_and_save_pdb_batches(pdb_files: List[str], 
                               output_dir: str,
                               batch_size: int = 1000,
                               n_processes: int = None) -> List[str]:
    """
    Process PDB files in batches and save results to Parquet files.
    
    Args:
        pdb_files: List of PDB file paths
        output_dir: Directory to save output Parquet files
        batch_size: Number of PDB files to process in each batch
        n_processes: Number of parallel processes to use
        
    Returns:
        List of paths to created Parquet files
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Default to using all available CPUs if not specified
    if n_processes is None:
        n_processes = cpu_count()
    
    parquet_files = []
    
    # Process files in batches
    for i in range(0, len(pdb_files), batch_size):
        batch = pdb_files[i:i+batch_size]
        batch_number = i // batch_size + 1
        
        print(f"Processing batch {batch_number}/{len(pdb_files)//batch_size + 1} ({len(batch)} files)")
        
        # Process the batch
        batch_results = process_pdb_batch(batch, n_processes=n_processes)
        
        # Save the results
        output_file = os.path.join(output_dir, f"rsa_batch_{batch_number:04d}.parquet")
        save_rsa_batch(batch_results, output_file)
        
        parquet_files.append(output_file)
        
        print(f"Saved batch {batch_number} to {output_file}")
        
    return parquet_files


def query_rsa_data(parquet_files: Union[str, List[str]], 
                  pdb_ids: Optional[List[str]] = None,
                  chains: Optional[List[str]] = None,
                  residue_range: Optional[Tuple[int, int]] = None) -> pd.DataFrame:
    """
    Query RSA data from Parquet files with efficient filtering.
    
    Args:
        parquet_files: Path to Parquet file or list of paths
        pdb_ids: Optional list of PDB IDs to filter by
        chains: Optional list of chain IDs to filter by
        residue_range: Optional tuple of (min_residue, max_residue) to filter by
        
    Returns:
        Filtered DataFrame with RSA data
    """
    # Handle single file or list of files
    if isinstance(parquet_files, str):
        parquet_files = [parquet_files]
    
    # Build filter expressions for PyArrow
    filters = []
    if pdb_ids:
        filters.append(('pdb_id', 'in', pdb_ids))
    if chains:
        filters.append(('chain_id', 'in', chains))
    if residue_range:
        min_res, max_res = residue_range
        filters.append(('residue_num', '>=', min_res))
        filters.append(('residue_num', '<=', max_res))
    
    # Convert filters to PyArrow format if any exist
    pa_filters = None
    if filters:
        pa_filters = filters if len(filters) == 1 else list(filters)
    
    # Read data with predicate pushdown (filtering at the storage level)
    df = pd.concat([
        pq.read_table(pf, filters=pa_filters).to_pandas()
        for pf in parquet_files
    ])
    
    return df


def get_structure_rsa(parquet_files: Union[str, List[str]], pdb_id: str) -> pd.DataFrame:
    """
    Get all RSA data for a specific structure.
    
    Args:
        parquet_files: Path to Parquet file or list of paths
        pdb_id: PDB ID to retrieve
        
    Returns:
        DataFrame with RSA data for the structure
    """
    return query_rsa_data(parquet_files, pdb_ids=[pdb_id])


def get_residue_rsa(parquet_files: Union[str, List[str]], 
                   pdb_id: str, 
                   chain_id: str, 
                   residue_num: int) -> float:
    """
    Get RSA value for a specific residue.
    
    Args:
        parquet_files: Path to Parquet file or list of paths
        pdb_id: PDB ID
        chain_id: Chain ID
        residue_num: Residue number
        
    Returns:
        RSA value or None if not found
    """
    df = query_rsa_data(
        parquet_files, 
        pdb_ids=[pdb_id], 
        chains=[chain_id], 
        residue_range=(residue_num, residue_num)
    )
    
    if len(df) == 0:
        return None
    return df.iloc[0]['rsa']


## Example usage
if __name__ == "__main__":
    # Process all PDB files in a directory
    pdb_dir = "../../1_inputs/raw_data/UP000005640_9606_HUMAN_v4"
    output_dir = "../../3_outputs/processed_data/rsa_scores"
    
    pdb_files = list(glob.glob(os.path.join(pdb_dir, "*.pdb*")))
    parquet_files = process_and_save_pdb_batches(pdb_files, output_dir)
    
    # Query examples
    # 1. Get RSA data for a specific PDB
    # df_6lzg = get_structure_rsa(parquet_files, "6lzg")
    # print(f"Retrieved {len(df_6lzg)} residues for 6LZG")
    
    # # 2. Get RSA for a specific residue
    # rsa = get_residue_rsa(parquet_files, "6lzg", "A", 42)
    # print(f"RSA for 6LZG chain A residue 42: {rsa}%")
    
    # # 3. Advanced query - all residues with RSA > 50%
    # all_data = query_rsa_data(parquet_files)
    # exposed_residues = all_data[all_data['rsa'] > 50.0]
    # print(f"Found {len(exposed_residues)} exposed residues (RSA > 50%)")