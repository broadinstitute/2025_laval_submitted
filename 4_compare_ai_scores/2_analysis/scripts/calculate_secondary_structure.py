"""
Enhanced Secondary Structure Storage and Query System
===================================================

This module provides comprehensive secondary structure analysis from protein
structure PDB files using DSSP, with enhanced features for mutation impact assessment:
- Standard secondary structure assignments (helix, sheet, loop)
- Secondary structure transitions and boundaries
- Loop classification and analysis
- Distance to stable structural elements
- Secondary structure stability metrics

Secondary structure context is crucial for understanding mutation impact because:
- Mutations in structured regions (helix/sheet) often have different effects than loops
- Transitions between secondary structures are often functionally important
- Loop regions have varying flexibility and functional roles
- Distance to stable elements indicates structural constraint levels

Key features:
- Parallel processing of PDB files using multiprocessing
- Enhanced DSSP analysis with additional structural metrics
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
from Bio.PDB import PDBParser, DSSP, PDBIO
from io import StringIO


# Secondary structure mapping
SS_MAPPING = {
    'H': 'alpha_helix',      # Alpha helix
    'B': 'beta_bridge',      # Beta bridge
    'E': 'beta_sheet',       # Extended strand (beta sheet)
    'G': '310_helix',        # 3-10 helix
    'I': 'pi_helix',         # Pi helix
    'T': 'turn',             # Hydrogen bonded turn
    'S': 'bend',             # Bend
    '-': 'loop'              # Loop or irregular
}

# Simplified structure types
SIMPLE_SS = {
    'H': 'helix',
    'G': 'helix',
    'I': 'helix',
    'E': 'sheet',
    'B': 'sheet',
    'T': 'loop',
    'S': 'loop',
    '-': 'loop'
}


def calculate_secondary_structure(pdb_file: Union[str, Path], 
                                 dssp_executable: str = 'mkdssp',
                                 analyze_transitions: bool = True) -> Dict[Tuple[str, int], dict]:
    """
    Calculate enhanced secondary structure features using DSSP.
    
    Args:
        pdb_file: Path to PDB file (can be gzipped)
        dssp_executable: Path to DSSP executable
        analyze_transitions: Whether to analyze secondary structure transitions
        
    Returns:
        Dictionary mapping (chain_id, residue_number) to structure features
    """
    pdb_path = Path(pdb_file).resolve()
    
    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")
        
    is_gzipped = str(pdb_path).endswith('.gz')
    parser = PDBParser(QUIET=True)
    
    try:
        # Parse structure
        if is_gzipped:
            with gzip.open(pdb_path, 'rt') as f:
                structure = parser.get_structure("protein", f)
        else:
            structure = parser.get_structure("protein", str(pdb_path))
    except Exception as e:
        raise ValueError(f"Error parsing PDB file: {e}")
    
    model = structure[0]
    
    # Create temporary PDB file for DSSP if needed
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
        
        # Run DSSP
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            dssp = DSSP(model, pdb_file_for_dssp, dssp=dssp_executable, file_type='PDB')
        
        # Extract basic DSSP data
        ss_data = {}
        residue_list = []
        
        for key in dssp.keys():
            chain_id = key[0]
            res_id = key[1][1]
            dssp_data = dssp[key]
            
            ss_code = dssp_data[2]  # Secondary structure
            phi = dssp_data[4]      # Phi angle
            psi = dssp_data[5]      # Psi angle
            
            residue_info = {
                'chain_id': chain_id,
                'res_num': res_id,
                'residue_name': dssp_data[1],
                'ss_code': ss_code,
                'ss_name': SS_MAPPING.get(ss_code, 'unknown'),
                'ss_simple': SIMPLE_SS.get(ss_code, 'unknown'),
                'phi': phi,
                'psi': psi,
                'accessibility': dssp_data[3],  # RSA from DSSP
            }
            
            ss_data[(chain_id, res_id)] = residue_info
            residue_list.append(residue_info)
        
        # Sort for transition analysis
        residue_list.sort(key=lambda x: (x['chain_id'], x['res_num']))
        
        # Analyze transitions and add enhanced features
        if analyze_transitions:
            ss_data = _analyze_ss_transitions(ss_data, residue_list)
        
        return ss_data
        
    except Exception as e:
        raise RuntimeError(f"Error running DSSP: {e}")
    
    finally:
        if temp_pdb_name and os.path.exists(temp_pdb_name):
            os.remove(temp_pdb_name)


def _analyze_ss_transitions(ss_data: Dict[Tuple[str, int], dict], 
                           residue_list: List[dict]) -> Dict[Tuple[str, int], dict]:
    """
    Analyze secondary structure transitions and add enhanced features.
    
    Args:
        ss_data: Basic secondary structure data
        residue_list: Sorted list of residues for transition analysis
        
    Returns:
        Enhanced secondary structure data with transition features
    """
    # Group by chain for proper transition analysis
    chains = {}
    for residue in residue_list:
        chain_id = residue['chain_id']
        if chain_id not in chains:
            chains[chain_id] = []
        chains[chain_id].append(residue)
    
    # Analyze each chain
    for chain_id, chain_residues in chains.items():
        for i, residue in enumerate(chain_residues):
            res_key = (residue['chain_id'], residue['res_num'])
            current_ss = residue['ss_simple']
            
            # Transition analysis
            is_transition = False
            transition_type = 'none'
            
            # Check neighbors
            prev_ss = chain_residues[i-1]['ss_simple'] if i > 0 else None
            next_ss = chain_residues[i+1]['ss_simple'] if i < len(chain_residues)-1 else None
            
            # Identify transition types
            if prev_ss and prev_ss != current_ss:
                is_transition = True
                transition_type = f"{prev_ss}_to_{current_ss}"
            elif next_ss and next_ss != current_ss:
                is_transition = True
                transition_type = f"{current_ss}_to_{next_ss}"
            
            # Loop analysis
            loop_length = 0
            loop_position = 'not_loop'
            if current_ss == 'loop':
                # Find loop boundaries
                loop_start = i
                loop_end = i
                
                # Find start of loop
                while loop_start > 0 and chain_residues[loop_start-1]['ss_simple'] == 'loop':
                    loop_start -= 1
                
                # Find end of loop
                while loop_end < len(chain_residues)-1 and chain_residues[loop_end+1]['ss_simple'] == 'loop':
                    loop_end += 1
                
                loop_length = loop_end - loop_start + 1
                loop_position_idx = i - loop_start
                
                if loop_position_idx == 0:
                    loop_position = 'loop_start'
                elif loop_position_idx == loop_length - 1:
                    loop_position = 'loop_end'
                else:
                    loop_position = 'loop_middle'
            
            # Distance to nearest structured element
            dist_to_helix = float('inf')
            dist_to_sheet = float('inf')
            
            for j, other_residue in enumerate(chain_residues):
                if i == j:
                    continue
                
                distance = abs(residue['res_num'] - other_residue['res_num'])
                
                if other_residue['ss_simple'] == 'helix' and distance < dist_to_helix:
                    dist_to_helix = distance
                elif other_residue['ss_simple'] == 'sheet' and distance < dist_to_sheet:
                    dist_to_sheet = distance
            
            # Secondary structure run length
            ss_run_length = 1
            run_position = 'single'
            
            if len([r for r in chain_residues if r['ss_simple'] == current_ss]) > 1:
                # Count consecutive residues with same SS
                run_start = i
                run_end = i
                
                while run_start > 0 and chain_residues[run_start-1]['ss_simple'] == current_ss:
                    run_start -= 1
                
                while run_end < len(chain_residues)-1 and chain_residues[run_end+1]['ss_simple'] == current_ss:
                    run_end += 1
                
                ss_run_length = run_end - run_start + 1
                run_pos_idx = i - run_start
                
                if run_pos_idx == 0:
                    run_position = 'ss_start'
                elif run_pos_idx == ss_run_length - 1:
                    run_position = 'ss_end'
                else:
                    run_position = 'ss_middle'
            
            # Update ss_data with enhanced features
            ss_data[res_key].update({
                'is_transition': is_transition,
                'transition_type': transition_type,
                'loop_length': loop_length,
                'loop_position': loop_position,
                'dist_to_helix': dist_to_helix if dist_to_helix != float('inf') else -1,
                'dist_to_sheet': dist_to_sheet if dist_to_sheet != float('inf') else -1,
                'dist_to_structured': min(dist_to_helix, dist_to_sheet) if min(dist_to_helix, dist_to_sheet) != float('inf') else -1,
                'ss_run_length': ss_run_length,
                'run_position': run_position
            })
    
    return ss_data


def process_pdb_batch(pdb_files: List[Union[str, Path]], 
                      dssp_executable: str = 'mkdssp',
                      analyze_transitions: bool = True,
                      n_processes: Optional[int] = None) -> Dict[str, dict]:
    """
    Process PDB files in parallel to calculate secondary structure features.
    
    Args:
        pdb_files: List of PDB file paths
        dssp_executable: Path to DSSP executable
        analyze_transitions: Whether to analyze transitions
        n_processes: Number of parallel processes
        
    Returns:
        Nested dictionary mapping filename to structure data
    """
    n_processes = n_processes or min(cpu_count(), len(pdb_files))
    
    process_func = partial(calculate_secondary_structure, 
                          dssp_executable=dssp_executable,
                          analyze_transitions=analyze_transitions)
    
    with Pool(n_processes) as pool:
        results = {}
        for filename, data in zip(pdb_files, pool.imap(process_func, pdb_files)):
            results[Path(filename).name] = data
        return results


def save_ss_batch(batch_results: Dict[str, dict], 
                  output_file: str,
                  compression: str = 'snappy') -> None:
    """
    Save secondary structure results to Parquet format.
    
    Args:
        batch_results: Dictionary of structure data
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
                'ss_code': data['ss_code'],
                'ss_name': data['ss_name'],
                'ss_simple': data['ss_simple'],
                'phi': data.get('phi', np.nan),
                'psi': data.get('psi', np.nan),
                'accessibility': data.get('accessibility', np.nan),
                'is_transition': data.get('is_transition', False),
                'transition_type': data.get('transition_type', 'none'),
                'loop_length': data.get('loop_length', 0),
                'loop_position': data.get('loop_position', 'not_loop'),
                'dist_to_helix': data.get('dist_to_helix', -1),
                'dist_to_sheet': data.get('dist_to_sheet', -1),
                'dist_to_structured': data.get('dist_to_structured', -1),
                'ss_run_length': data.get('ss_run_length', 1),
                'run_position': data.get('run_position', 'single')
            }
            rows.append(row)
    
    table = pa.Table.from_pandas(pd.DataFrame(rows))
    pq.write_table(table, output_file, compression=compression)


def process_and_save_batches(pdb_files: List[str], 
                           output_dir: str,
                           batch_size: int = 1000,
                           dssp_executable: str = 'mkdssp',
                           n_processes: int = None) -> List[str]:
    """
    Process PDB files in memory-efficient batches.
    
    Args:
        pdb_files: List of input PDB paths
        output_dir: Output directory for Parquet files
        batch_size: Number of files per batch
        dssp_executable: Path to DSSP executable
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
                                   dssp_executable=dssp_executable,
                                   n_processes=n_processes)
        
        output_path = Path(output_dir) / f"secondary_structure_batch_{batch_num:04d}.parquet"
        save_ss_batch(results, str(output_path))
        
        parquet_files.append(str(output_path))
        print(f"Saved batch {batch_num} to {output_path}")
    
    return parquet_files


def query_ss_data(parquet_files: Union[str, List[str]],
                  pdb_ids: List[str] = None,
                  chains: List[str] = None,
                  ss_types: List[str] = None,
                  transitions_only: bool = False) -> pd.DataFrame:
    """
    Query secondary structure data with various filters.
    
    Args:
        parquet_files: Parquet file(s) to query
        pdb_ids: Filter by PDB IDs
        chains: Filter by chain IDs
        ss_types: Filter by secondary structure types
        transitions_only: Only return transition residues
        
    Returns:
        Filtered DataFrame of secondary structure data
    """
    if isinstance(parquet_files, str):
        parquet_files = [parquet_files]
    
    filters = []
    if pdb_ids:
        filters.append(('pdb_id', 'in', pdb_ids))
    if chains:
        filters.append(('chain_id', 'in', chains))
    if ss_types:
        filters.append(('ss_simple', 'in', ss_types))
    if transitions_only:
        filters.append(('is_transition', '=', True))
    
    return pq.read_table(
        parquet_files,
        filters=filters or None
    ).to_pandas()


def get_loop_residues(parquet_files: Union[str, List[str]],
                     min_loop_length: int = 3,
                     max_loop_length: int = None) -> pd.DataFrame:
    """
    Get loop residues with specified length criteria.
    
    Args:
        parquet_files: Parquet file(s) to query
        min_loop_length: Minimum loop length
        max_loop_length: Maximum loop length (None for no limit)
        
    Returns:
        DataFrame of loop residues
    """
    filters = [
        ('ss_simple', '=', 'loop'),
        ('loop_length', '>=', min_loop_length)
    ]
    
    if max_loop_length:
        filters.append(('loop_length', '<=', max_loop_length))
    
    return pq.read_table(parquet_files, filters=filters).to_pandas()


## Example usage
if __name__ == "__main__":
    # Configuration
    input_dir = "../../1_inputs/raw_data/UP000005640_9606_HUMAN_v4"
    output_dir = "../../3_outputs/processed_data/secondary_structure"
    
    # Get all compressed PDB files
    pdb_files = list(Path(input_dir).glob("*.pdb*"))
    
    # Process and save data
    parquet_files = process_and_save_batches(pdb_files, output_dir)
    
    # Example queries
    # 1. Get all helix residues
    # helices = query_ss_data(parquet_files, ss_types=['helix'])
    
    # 2. Get transition residues only
    # transitions = query_ss_data(parquet_files, transitions_only=True)
    
    # 3. Get long loops (potentially flexible regions)
    # long_loops = get_loop_residues(parquet_files, min_loop_length=10)
    
    # 4. Get sheet residues far from helices
    # df = query_ss_data(parquet_files, ss_types=['sheet'])
    # far_from_helix = df[df['dist_to_helix'] > 10]