"""
Backbone Dihedral Angles Storage and Query System
===============================================

This module calculates backbone dihedral angles (phi, psi, omega) from protein
structure PDB files, providing insights into:
- Local conformational constraints
- Ramachandran plot regions and structural preferences
- Conformational flexibility and rigidity
- Mutation impact on backbone geometry

Backbone dihedral angles are crucial for understanding mutation impact because:
- Certain phi/psi combinations are energetically favorable (allowed regions)
- Mutations can force backbone into unfavorable conformations
- Proline and glycine have unique conformational constraints
- Dihedral angles correlate with secondary structure stability

Key features:
- Calculation of phi, psi, and omega dihedral angles
- Ramachandran region classification
- Conformational strain analysis
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
import math
from pathlib import Path
from functools import partial
from multiprocessing import Pool, cpu_count
from typing import List, Dict, Tuple, Optional, Union
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from Bio.PDB import PDBParser, calc_dihedral
from Bio.PDB.vectors import Vector
from io import StringIO


# Ramachandran regions (approximate boundaries in degrees)
RAMACHANDRAN_REGIONS = {
    'alpha_helix': {'phi': (-180, -30), 'psi': (-90, 50)},
    'beta_sheet': {'phi': (-180, -30), 'psi': (50, 180)},
    'left_helix': {'phi': (30, 180), 'psi': (-50, 90)},
    'extended': {'phi': (-90, 30), 'psi': (90, 180)},
    'bridge': {'phi': (-90, 30), 'psi': (-180, -90)}
}

def angle_in_range(angle: float, range_tuple: Tuple[float, float]) -> bool:
    """Check if angle is within a range, handling wraparound."""
    min_angle, max_angle = range_tuple
    if min_angle <= max_angle:
        return min_angle <= angle <= max_angle
    else:  # wraparound case
        return angle >= min_angle or angle <= max_angle


def classify_ramachandran_region(phi: float, psi: float) -> str:
    """
    Classify phi/psi angles into Ramachandran regions.
    
    Args:
        phi: Phi dihedral angle in degrees
        psi: Psi dihedral angle in degrees
        
    Returns:
        Region name or 'disallowed'
    """
    if math.isnan(phi) or math.isnan(psi):
        return 'unknown'
    
    # Check each region
    for region, bounds in RAMACHANDRAN_REGIONS.items():
        if (angle_in_range(phi, bounds['phi']) and 
            angle_in_range(psi, bounds['psi'])):
            return region
    
    return 'disallowed'


def calculate_backbone_dihedrals(pdb_file: Union[str, Path],
                                include_omega: bool = True,
                                analyze_strain: bool = True) -> Dict[Tuple[str, int], dict]:
    """
    Calculate backbone dihedral angles for residues in a PDB file.
    
    Args:
        pdb_file: Path to PDB file (can be gzipped)
        include_omega: Whether to calculate omega angles
        analyze_strain: Whether to analyze conformational strain
        
    Returns:
        Dictionary mapping (chain_id, residue_number) to dihedral data
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
        
        # Group residues by chain for proper dihedral calculation
        chains = {}
        for residue in structure.get_residues():
            chain_id = residue.get_parent().id
            if chain_id not in chains:
                chains[chain_id] = []
            chains[chain_id].append(residue)
        
        # Sort residues by residue number within each chain
        for chain_id in chains:
            chains[chain_id].sort(key=lambda r: r.id[1])
        
        dihedrals = {}
        
        # Calculate dihedrals for each chain
        for chain_id, residues in chains.items():
            for i, residue in enumerate(residues):
                try:
                    res_num = residue.id[1]
                    res_name = residue.resname
                    
                    dihedral_data = {
                        'residue_name': res_name,
                        'phi': np.nan,
                        'psi': np.nan,
                        'omega': np.nan,
                        'rama_region': 'unknown'
                    }
                    
                    # Calculate Phi angle (C[i-1] - N[i] - CA[i] - C[i])
                    if i > 0:  # Need previous residue
                        try:
                            prev_residue = residues[i-1]
                            phi = calc_dihedral(
                                prev_residue['C'].get_vector(),
                                residue['N'].get_vector(),
                                residue['CA'].get_vector(),
                                residue['C'].get_vector()
                            )
                            dihedral_data['phi'] = math.degrees(phi)
                        except KeyError:
                            pass  # Missing atoms
                    
                    # Calculate Psi angle (N[i] - CA[i] - C[i] - N[i+1])
                    if i < len(residues) - 1:  # Need next residue
                        try:
                            next_residue = residues[i+1]
                            psi = calc_dihedral(
                                residue['N'].get_vector(),
                                residue['CA'].get_vector(),
                                residue['C'].get_vector(),
                                next_residue['N'].get_vector()
                            )
                            dihedral_data['psi'] = math.degrees(psi)
                        except KeyError:
                            pass  # Missing atoms
                    
                    # Calculate Omega angle (CA[i-1] - C[i-1] - N[i] - CA[i])
                    if include_omega and i > 0:
                        try:
                            prev_residue = residues[i-1]
                            omega = calc_dihedral(
                                prev_residue['CA'].get_vector(),
                                prev_residue['C'].get_vector(),
                                residue['N'].get_vector(),
                                residue['CA'].get_vector()
                            )
                            dihedral_data['omega'] = math.degrees(omega)
                        except KeyError:
                            pass
                    
                    # Classify Ramachandran region
                    if not (math.isnan(dihedral_data['phi']) or math.isnan(dihedral_data['psi'])):
                        dihedral_data['rama_region'] = classify_ramachandran_region(
                            dihedral_data['phi'], dihedral_data['psi']
                        )
                    
                    # Analyze conformational strain
                    if analyze_strain:
                        strain_metrics = _analyze_conformational_strain(
                            dihedral_data, res_name
                        )
                        dihedral_data.update(strain_metrics)
                    
                    dihedrals[(chain_id, res_num)] = dihedral_data
                    
                except Exception as e:
                    # Skip problematic residues
                    continue
        
        return dihedrals
        
    except Exception as e:
        raise RuntimeError(f"Error processing {pdb_path.name}: {e}")


def _analyze_conformational_strain(dihedral_data: dict, res_name: str) -> dict:
    """
    Analyze conformational strain based on dihedral angles and residue type.
    
    Args:
        dihedral_data: Dictionary containing phi, psi, omega angles
        res_name: Three-letter residue name
        
    Returns:
        Dictionary with strain analysis metrics
    """
    phi = dihedral_data.get('phi', np.nan)
    psi = dihedral_data.get('psi', np.nan)
    omega = dihedral_data.get('omega', np.nan)
    rama_region = dihedral_data.get('rama_region', 'unknown')
    
    strain_data = {
        'is_strained': False,
        'strain_score': 0.0,
        'cis_peptide': False,
        'special_residue': False
    }
    
    if math.isnan(phi) or math.isnan(psi):
        return strain_data
    
    # Check for cis peptide bond (omega around 0°)
    if not math.isnan(omega):
        if abs(omega) < 30:  # Cis peptide
            strain_data['cis_peptide'] = True
            if res_name != 'PRO':  # Cis is more strained for non-proline
                strain_data['strain_score'] += 2.0
    
    # Special residue considerations
    if res_name == 'GLY':
        strain_data['special_residue'] = True
        # Glycine can access any phi/psi, so less strain
        strain_data['strain_score'] *= 0.5
    elif res_name == 'PRO':
        strain_data['special_residue'] = True
        # Proline is constrained to certain phi values
        if phi > -30:  # Proline should have phi around -60°
            strain_data['strain_score'] += 1.5
    
    # Ramachandran-based strain
    if rama_region == 'disallowed':
        strain_data['is_strained'] = True
        strain_data['strain_score'] += 3.0
    elif rama_region == 'bridge':  # Less favorable region
        strain_data['strain_score'] += 1.0
    
    # Additional strain metrics based on specific angle values
    # High energy regions in Ramachandran space
    if 0 < phi < 90 and -90 < psi < 0:  # Positive phi, negative psi region
        strain_data['strain_score'] += 1.5
    
    # Set strain flag based on total score
    strain_data['is_strained'] = strain_data['strain_score'] > 1.0
    
    return strain_data


def process_pdb_batch(pdb_files: List[Union[str, Path]], 
                      include_omega: bool = True,
                      analyze_strain: bool = True,
                      n_processes: Optional[int] = None) -> Dict[str, dict]:
    """
    Process PDB files in parallel to calculate dihedral angles.
    
    Args:
        pdb_files: List of PDB file paths
        include_omega: Whether to calculate omega angles
        analyze_strain: Whether to analyze conformational strain
        n_processes: Number of parallel processes
        
    Returns:
        Nested dictionary mapping filename to dihedral data
    """
    n_processes = n_processes or min(cpu_count(), len(pdb_files))
    
    process_func = partial(calculate_backbone_dihedrals, 
                          include_omega=include_omega,
                          analyze_strain=analyze_strain)
    
    with Pool(n_processes) as pool:
        results = {}
        for filename, data in zip(pdb_files, pool.imap(process_func, pdb_files)):
            results[Path(filename).name] = data
        return results


def save_dihedrals_batch(batch_results: Dict[str, dict], 
                        output_file: str,
                        compression: str = 'snappy') -> None:
    """
    Save dihedral angle results to Parquet format.
    
    Args:
        batch_results: Dictionary of dihedral data
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
                'phi': data['phi'],
                'psi': data['psi'],
                'omega': data['omega'],
                'rama_region': data['rama_region'],
                'is_strained': data.get('is_strained', False),
                'strain_score': data.get('strain_score', 0.0),
                'cis_peptide': data.get('cis_peptide', False),
                'special_residue': data.get('special_residue', False)
            }
            rows.append(row)
    
    table = pa.Table.from_pandas(pd.DataFrame(rows))
    pq.write_table(table, output_file, compression=compression)


def process_and_save_batches(pdb_files: List[str], 
                           output_dir: str,
                           batch_size: int = 1000,
                           include_omega: bool = True,
                           n_processes: int = None) -> List[str]:
    """
    Process PDB files in memory-efficient batches.
    
    Args:
        pdb_files: List of input PDB paths
        output_dir: Output directory for Parquet files
        batch_size: Number of files per batch
        include_omega: Whether to calculate omega angles
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
                                   include_omega=include_omega,
                                   n_processes=n_processes)
        
        output_path = Path(output_dir) / f"backbone_dihedrals_batch_{batch_num:04d}.parquet"
        save_dihedrals_batch(results, str(output_path))
        
        parquet_files.append(str(output_path))
        print(f"Saved batch {batch_num} to {output_path}")
    
    return parquet_files


def query_dihedrals(parquet_files: Union[str, List[str]],
                   pdb_ids: List[str] = None,
                   chains: List[str] = None,
                   rama_regions: List[str] = None,
                   strained_only: bool = False,
                   phi_range: Tuple[float, float] = None,
                   psi_range: Tuple[float, float] = None) -> pd.DataFrame:
    """
    Query dihedral angle data with various filters.
    
    Args:
        parquet_files: Parquet file(s) to query
        pdb_ids: Filter by PDB IDs
        chains: Filter by chain IDs
        rama_regions: Filter by Ramachandran regions
        strained_only: Only return conformationally strained residues
        phi_range: Filter by phi angle range (min, max)
        psi_range: Filter by psi angle range (min, max)
        
    Returns:
        Filtered DataFrame of dihedral data
    """
    if isinstance(parquet_files, str):
        parquet_files = [parquet_files]
    
    filters = []
    if pdb_ids:
        filters.append(('pdb_id', 'in', pdb_ids))
    if chains:
        filters.append(('chain_id', 'in', chains))
    if rama_regions:
        filters.append(('rama_region', 'in', rama_regions))
    if strained_only:
        filters.append(('is_strained', '=', True))
    if phi_range:
        filters.extend([('phi', '>=', phi_range[0]), ('phi', '<=', phi_range[1])])
    if psi_range:
        filters.extend([('psi', '>=', psi_range[0]), ('psi', '<=', psi_range[1])])
    
    return pq.read_table(
        parquet_files,
        filters=filters or None
    ).to_pandas()


def get_strained_residues(parquet_files: Union[str, List[str]],
                         min_strain_score: float = 1.0) -> pd.DataFrame:
    """
    Find conformationally strained residues.
    
    Args:
        parquet_files: Parquet file(s) to query
        min_strain_score: Minimum strain score threshold
        
    Returns:
        DataFrame of strained residues
    """
    filters = [('strain_score', '>=', min_strain_score)]
    return pq.read_table(parquet_files, filters=filters).to_pandas()


def get_ramachandran_outliers(parquet_files: Union[str, List[str]]) -> pd.DataFrame:
    """
    Find residues in disallowed Ramachandran regions.
    
    Args:
        parquet_files: Parquet file(s) to query
        
    Returns:
        DataFrame of Ramachandran outliers
    """
    return query_dihedrals(parquet_files, rama_regions=['disallowed'])


## Example usage
if __name__ == "__main__":
    # Configuration
    input_dir = "../../1_inputs/raw_data/UP000005640_9606_HUMAN_v4"
    output_dir = "../../3_outputs/processed_data/backbone_dihedrals"
    
    # Get all compressed PDB files
    pdb_files = list(Path(input_dir).glob("*.pdb*"))
    
    # Process and save data
    parquet_files = process_and_save_batches(
        pdb_files, 
        output_dir,
        include_omega=True  # Include omega angles
    )
    
    # Example queries
    # 1. Get dihedral data for a specific protein
    # df = query_dihedrals(parquet_files, pdb_ids=["AF-P00533-F1"])
    
    # 2. Find conformationally strained residues
    # strained = get_strained_residues(parquet_files, min_strain_score=2.0)
    
    # 3. Find Ramachandran outliers
    # outliers = get_ramachandran_outliers(parquet_files)
    
    # 4. Get alpha-helical residues
    # helices = query_dihedrals(parquet_files, rama_regions=['alpha_helix'])
    
    # 5. Find residues with specific phi/psi ranges
    # specific_angles = query_dihedrals(parquet_files, 
    #                                 phi_range=(-80, -40), 
    #                                 psi_range=(-60, -20))