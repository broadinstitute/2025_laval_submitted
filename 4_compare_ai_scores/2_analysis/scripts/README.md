# Enhanced Structural Features for Protein Mutation Impact Assessment

This directory contains a comprehensive suite of Python scripts for extracting enhanced structural features from AlphaFold protein structures. These features go beyond the traditional pLDDT and RSA scores to provide a multi-dimensional view of structural context that is crucial for understanding mutation impact.

## Overview

Traditional structural features like pLDDT (confidence) and RSA (solvent accessibility) provide valuable but limited information about mutation impact. This enhanced feature set captures:

- **Local structural flexibility and transitions**
- **Packing density and environmental constraints** 
- **Global protein architecture and positioning**
- **Conformational constraints and geometric stress**
- **Secondary structure context and boundaries**

## Feature Categories and Scripts

### 1. Structural Flexibility Analysis

#### `calculate_bfactor_gradients.py`

**Purpose**: Identifies regions of varying structural confidence and flexibility transitions that may be sensitive to mutations.

**Key Calculations**:

- **Absolute Gradient**: `np.std(neighborhood_plddt_values)`
  - Standard deviation of pLDDT scores in a local window (default: 5 residues)
  - High values indicate rapid confidence changes (flexibility transitions)
  - Units: pLDDT confidence units

- **Relative Gradient**: `np.std(neighborhood_plddt_values) / np.mean(neighborhood_plddt_values)`
  - Coefficient of variation, normalizes for local confidence level
  - Identifies regions where confidence variation is high relative to absolute confidence
  - Dimensionless ratio

- **Max Change**: `max(|current_plddt - neighbor_plddt|)`
  - Maximum pLDDT difference between current residue and any neighbor
  - Captures the most dramatic local confidence transitions
  - Units: pLDDT confidence units

**Output Features**:
- `absolute_gradient`: Local flexibility measure
- `relative_gradient`: Normalized flexibility measure  
- `max_change`: Maximum local confidence change
- `neighborhood_mean/min/max`: Local confidence statistics
- `window_size`: Number of neighbors used in calculation

**Biological Significance**:
- **High gradients** often indicate domain boundaries, loop regions, or functionally important transitions
- **Low gradients** suggest uniformly structured or unstructured regions
- **Mutations in high-gradient regions** may have cascading structural effects

---

### 2. Local Environment and Packing Analysis

#### `calculate_contact_numbers.py`

**Purpose**: Quantifies local packing density and structural constraints around each residue.

**Key Calculations**:

- **Contact Number**: Count of residues within distance cutoff
  ```python
  contacts_6A = sum(1 for neighbor in residues if distance(current, neighbor) <= 6.0)
  contacts_8A = sum(1 for neighbor in residues if distance(current, neighbor) <= 8.0)  
  contacts_10A = sum(1 for neighbor in residues if distance(current, neighbor) <= 10.0)
  ```
  - 6Å: Close contacts (mainly side-chain interactions)
  - 8Å: Medium-range contacts (backbone and extended side chains)
  - 10Å: Long-range contacts (domain interactions)

- **Packing Density**: `contact_number / (total_residues - 1)`
  - Normalized contact number accounting for protein size
  - Allows comparison across proteins of different sizes
  - Range: 0.0 (no contacts) to 1.0 (all residues in contact)

- **Distance Statistics**:
  - `closest_distance`: Distance to nearest neighbor
  - `mean_distance`: Average distance to all other residues
  - `median_distance`: Median distance to all other residues
  - `distances_std`: Standard deviation of distances

**Output Features**:
- `contacts_6A/8A/10A`: Contact counts at different distance thresholds
- `density_6A/8A/10A`: Normalized packing densities
- `closest_distance`: Nearest neighbor distance
- `mean/median_distance`: Distance distribution statistics
- `distances_std`: Distance variability measure

**Biological Significance**:
- **High contact numbers** indicate buried, constrained residues
- **Low contact numbers** suggest surface or flexible regions
- **Contact patterns** correlate with evolutionary conservation and mutation sensitivity
- **Different distance thresholds** capture different types of structural interactions

---

### 3. Secondary Structure Context Analysis

#### `calculate_secondary_structure.py`

**Purpose**: Provides detailed secondary structure context including transitions, loop analysis, and distance to stable structural elements.

**Key Calculations**:

- **Basic Secondary Structure**: Uses DSSP to assign standard classifications
  - H: α-helix, E: β-sheet, T: turn, S: bend, -: loop/coil
  - Simplified to: helix, sheet, loop

- **Transition Analysis**: Identifies boundaries between secondary structure elements
  ```python
  is_transition = (prev_ss != current_ss) or (current_ss != next_ss)
  transition_type = f"{prev_ss}_to_{current_ss}"
  ```

- **Loop Analysis**: For residues in loops
  - `loop_length`: Total length of the loop containing this residue
  - `loop_position`: Position within loop (start/middle/end)

- **Distance to Structured Elements**:
  ```python
  dist_to_helix = min(|current_pos - helix_pos| for helix_pos in helix_positions)
  dist_to_sheet = min(|current_pos - sheet_pos| for sheet_pos in sheet_positions)
  dist_to_structured = min(dist_to_helix, dist_to_sheet)
  ```

- **Secondary Structure Run Analysis**:
  - `ss_run_length`: Length of continuous secondary structure element
  - `run_position`: Position within the structural element (start/middle/end)

**Output Features**:
- `ss_code/ss_name/ss_simple`: Secondary structure assignments
- `phi/psi`: Backbone dihedral angles from DSSP
- `accessibility`: RSA values from DSSP
- `is_transition`: Boolean flag for structural transitions
- `transition_type`: Type of transition (helix_to_loop, etc.)
- `loop_length/loop_position`: Loop context information
- `dist_to_helix/sheet/structured`: Distance to stable elements
- `ss_run_length/run_position`: Position within structural elements

**Biological Significance**:
- **Transition regions** are often functionally important and mutation-sensitive
- **Loop regions** vary in flexibility and functional roles
- **Distance to structured elements** indicates degree of structural support
- **Position within elements** affects mutation tolerance (ends vs. middle)

---

### 4. Global Architecture and Positioning

#### `calculate_distance_to_center.py`

**Purpose**: Characterizes global positioning within the protein architecture to distinguish core vs. surface residues.

**Key Calculations**:

- **Protein Centers**: Multiple center definitions for robustness
  ```python
  geometric_center = np.mean(all_heavy_atom_coords, axis=0)
  center_of_mass = np.average(all_heavy_atom_coords, weights=atomic_masses, axis=0)
  ca_geometric_center = np.mean(ca_atom_coords, axis=0)
  ```

- **Distance Calculations**: From both CA atom and residue center-of-mass
  ```python
  ca_distance = ||ca_coord - center_coord||
  com_distance = ||residue_com - center_coord||
  ```

- **Normalized Distances**: Relative positioning within protein
  ```python
  normalized_distance = distance_to_center / protein_radius
  # 0.0 = at center, 1.0 = at surface
  ```

- **Radial Classification**: Categorical assignment
  ```python
  if normalized_dist <= 0.33: radial_class = 'core'
  elif normalized_dist <= 0.66: radial_class = 'intermediate'
  else: radial_class = 'surface'
  ```

- **Radial Percentile**: Statistical position in radial distribution
  ```python
  radial_percentile = (residues_closer_to_center / total_residues) * 100
  ```

**Output Features**:
- `ca_dist_geometric/mass/ca_geometric`: Distances from different centers
- `com_dist_geometric/mass/ca_geometric`: Center-of-mass distances
- `normalized_dist_*`: Relative radial positions (0-1)
- `radial_percentile_*`: Percentile in radial distribution
- `radial_class_*`: Categorical assignment (core/intermediate/surface)
- `protein_radius_*`: Overall protein size metrics

**Biological Significance**:
- **Core residues** are often evolutionarily conserved and mutation-intolerant
- **Surface residues** may have different functional roles (binding, stability)
- **Radial position** correlates with mutation impact and evolutionary constraint
- **Multiple center definitions** provide robustness against structural asymmetry

---

### 5. Conformational Constraints Analysis

#### `calculate_backbone_dihedrals.py`

**Purpose**: Analyzes backbone geometry and conformational constraints that affect mutation accommodation.

**Key Calculations**:

- **Dihedral Angles**: Standard protein geometry angles
  ```python
  phi = dihedral(C[i-1], N[i], CA[i], C[i])      # Backbone flexibility
  psi = dihedral(N[i], CA[i], C[i], N[i+1])     # Local conformation  
  omega = dihedral(CA[i-1], C[i-1], N[i], CA[i]) # Peptide bond geometry
  ```

- **Ramachandran Classification**: Assigns phi/psi combinations to allowed regions
  ```python
  regions = {
    'alpha_helix': phi ∈ (-180,-30), psi ∈ (-90,50),
    'beta_sheet': phi ∈ (-180,-30), psi ∈ (50,180),
    'left_helix': phi ∈ (30,180), psi ∈ (-50,90),
    'extended': phi ∈ (-90,30), psi ∈ (90,180),
    'disallowed': outside allowed regions
  }
  ```

- **Conformational Strain Analysis**: Identifies geometrically stressed regions
  ```python
  strain_score = 0
  if rama_region == 'disallowed': strain_score += 3.0
  if cis_peptide and residue != 'PRO': strain_score += 2.0
  if positive_phi_negative_psi: strain_score += 1.5
  is_strained = strain_score > 1.0
  ```

- **Special Residue Handling**:
  - **Glycine**: Reduced strain (can access any phi/psi)
  - **Proline**: Constrained phi angle (~-60°)
  - **Cis peptides**: Higher strain unless proline

**Output Features**:
- `phi/psi/omega`: Backbone dihedral angles (degrees)
- `rama_region`: Ramachandran plot classification
- `is_strained`: Boolean flag for conformational strain
- `strain_score`: Quantitative strain measure
- `cis_peptide`: Boolean flag for cis peptide bonds
- `special_residue`: Boolean flag for Gly/Pro

**Biological Significance**:
- **Allowed regions** indicate conformationally stable geometries
- **Disallowed regions** suggest high-energy conformations sensitive to mutation
- **Conformational strain** predicts difficulty of mutation accommodation
- **Special residues** have unique constraint profiles affecting mutation tolerance

---

## Integration and Usage

### `extract_enhanced_structural_features.py`

**Purpose**: Master script that orchestrates extraction of all enhanced features and provides integration utilities.

**Key Functions**:

- **`extract_all_structural_features()`**: 
  - Processes PDB files through all feature calculation pipelines
  - Handles batch processing and parallel execution
  - Creates organized output directory structure

- **`consolidate_features_for_variants()`**:
  - Joins all structural features with variant datasets
  - Handles AlphaFold naming conventions
  - Creates analysis-ready combined datasets

- **`analyze_feature_correlations()`**:
  - Computes correlation matrices between features
  - Identifies redundant or complementary features
  - Guides feature selection for downstream analysis

- **`create_feature_summary_stats()`**:
  - Generates comprehensive statistical summaries
  - Groups features by biological category
  - Provides data quality and distribution insights

---

## Biological Interpretation Framework

### Feature Combinations for Mutation Impact Assessment

**High-Impact Mutations (Likely Pathogenic)**:
- High contact numbers + low RSA (buried, constrained)
- High B-factor gradients (flexibility transitions)
- Core positioning + secondary structure transitions
- Conformational strain + disallowed Ramachandran regions

**Low-Impact Mutations (Likely Benign)**:
- Low contact numbers + high RSA (surface exposed)
- Low B-factor gradients (uniform confidence)
- Surface positioning + loop regions
- Allowed conformations + low strain scores

**Context-Dependent Mutations**:
- Intermediate contact numbers + intermediate positioning
- Loop regions with high structural support
- Transition regions with low strain
- Surface residues in functional domains

### Feature Scaling and Normalization

**Recommended Preprocessing**:
- **Distance features**: Normalize by protein size
- **Contact numbers**: Use density ratios for cross-protein comparison
- **Gradient features**: Consider protein-specific percentiles
- **Dihedral angles**: Use circular statistics where appropriate

### Integration with Existing Predictors

**Complementarity with AI Predictors**:
- **AlphaFold features** provide structural context for sequence-based predictions
- **ESM/ProtTrans embeddings** capture evolutionary patterns, structural features add geometric constraints
- **Combined models** can weight predictions by structural confidence and context

---

## Technical Implementation Details

### Dependencies

```python
# Core scientific computing
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

# Structural biology
from Bio.PDB import PDBParser, DSSP, calc_dihedral
from Bio.PDB.vectors import Vector

# Parallel processing
from multiprocessing import Pool, cpu_count
from functools import partial

# System utilities
import gzip
import tempfile
from pathlib import Path
```

### Performance Considerations

**Memory Management**:
- Batch processing prevents memory overflow
- Parquet format provides efficient storage and querying
- Lazy loading with Polars for large datasets

**Computational Efficiency**:
- Parallel processing across PDB files
- Vectorized numpy operations where possible
- Efficient distance calculations using spatial indexing

**Storage Format**:
- Columnar Parquet format for fast analytical queries
- Predicate pushdown for efficient filtering
- Compression reduces storage footprint

### Error Handling

**Robust Processing**:
- Graceful handling of missing atoms or residues
- Skip problematic structures with detailed logging
- Validation of input data integrity

**Quality Control**:
- Range checking for calculated values
- Consistency validation across features
- Missing data handling with appropriate defaults

---

## Example Usage Workflows

### 1. Extract Features for New Dataset

```python
import sys
sys.path.append('scripts/')
import extract_enhanced_structural_features as enhanced

# Process PDB files
pdb_files = list(Path("pdb_directory").glob("*.pdb.gz"))
parquet_files = enhanced.extract_all_structural_features(
    pdb_files, 
    "output_directory",
    batch_size=100,
    n_processes=8
)
```

### 2. Add Features to Variant Analysis

```python
# Load your variant dataset
variants_df = pd.read_csv("variants.csv")

# Consolidate all structural features
enhanced_variants = enhanced.consolidate_features_for_variants(
    parquet_files,
    variants_df,
    uniprot_col='UniProt_ID',
    position_col='Position', 
    residue_col='Reference_AA'
)

# Now enhanced_variants has all structural features joined
```

### 3. Feature Analysis and Selection

```python
# Analyze feature correlations
correlations = enhanced.analyze_feature_correlations(enhanced_variants)

# Get summary statistics
stats = enhanced.create_feature_summary_stats(enhanced_variants)

# Select non-redundant features
selected_features = correlations.abs() < 0.8  # Remove highly correlated
```

### 4. Integration with Machine Learning

```python
# Prepare feature matrix
structural_features = [
    'absolute_gradient', 'contacts_8A', 'ca_dist_geometric',
    'strain_score', 'dist_to_structured', 'normalized_dist_geometric'
]

X_structural = enhanced_variants[structural_features].fillna(0)
y = enhanced_variants['pathogenicity_label']

# Combine with sequence-based features
X_combined = pd.concat([X_sequence, X_structural], axis=1)
```

---

## Future Extensions

### Additional Features Under Development

1. **Network Centrality Measures**:
   - Betweenness centrality in residue interaction networks
   - Clustering coefficients for local connectivity
   - Shortest path analysis for allosteric networks

2. **Domain Architecture Features**:
   - Domain boundary detection and classification  
   - Inter-domain interaction analysis
   - Domain-specific evolutionary constraints

3. **Dynamic Properties**:
   - B-factor based flexibility predictions
   - Normal mode analysis for collective motions
   - Allosteric pathway identification

4. **Functional Site Integration**:
   - Distance to known functional sites
   - Binding pocket involvement
   - Catalytic site proximity analysis

### Performance Optimizations

1. **GPU Acceleration**: CUDA implementations for distance calculations
2. **Distributed Computing**: Spark/Dask integration for large-scale processing
3. **Incremental Updates**: Add new structures without reprocessing existing data

---

## References and Validation

### Experimental Validation

These features have been designed based on established principles from:
- Structural biology literature on mutation effects
- Protein folding and stability studies  
- Evolutionary constraint analysis
- Clinical mutation databases (ClinVar, HGMD)

### Benchmark Datasets

Recommended validation datasets:
- **ClinVar**: Clinically annotated variants
- **gnomAD**: Population frequency data
- **DMS experiments**: Deep mutational scanning results
- **Thermodynamic studies**: Stability measurements

### Citation

If you use these enhanced structural features in your research, please cite:

```
[Your VarChAMP publication]
Enhanced Structural Features for Protein Mutation Impact Assessment
https://github.com/[your-repo]/varchamp
```

---

## Support and Troubleshooting

### Common Issues

1. **DSSP Installation**: Ensure `mkdssp` is in PATH
2. **Memory Limitations**: Reduce batch_size for large proteins
3. **Missing PDB Files**: Check AlphaFold database organization
4. **Performance**: Adjust n_processes based on system capabilities

### Contact

For technical support or feature requests:
- Create issues in the GitHub repository
- Contact: [your contact information]

---

This comprehensive structural feature extraction suite provides a powerful foundation for understanding mutation impact through multi-scale structural analysis, significantly extending beyond traditional confidence and accessibility measures.