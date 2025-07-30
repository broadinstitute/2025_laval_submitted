# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This is the VarChAMP 2025 publication repository containing functional assay data analysis and manuscript generation code. The repository follows a structured organization where folders with prefixes 1_ through 6_ correspond to different manuscript sections.

### Repository Structure

Each analysis section follows a consistent pattern:
- **1_inputs**: Raw input data or annotations (e.g., ClinVar annotations, unnormalized assay data)
- **2_analysis**: Processing and analysis scripts (Jupyter notebooks, Python scripts, R code)
- **3_outputs**: Intermediate outputs that serve as inputs to other analyses
- **4_manuscript_figures**: High-resolution figures for main text and supplementary materials
- **5_manuscript_tables**: Formatted tables for manuscript inclusion

### Main Analysis Sections

1. **1_allele_collection/**: Allele collection and curation
2. **2_individual_assay_results/**: Individual functional assay analyses
   - `dual_ipa/`: Dual IPA assay analysis
   - `imaging/`: Cell painting imaging analysis with Snakemake pipeline
   - `ppi/`: Protein-protein interaction (Y2H) analysis
   - `threshold_calibration/`: Statistical threshold calibration using R
3. **3_integrated_assay_results/**: Cross-assay integration and comparison
4. **4_compare_ai_scores/**: AI predictor comparison and benchmarking
5. **5_mechanistic_case_studies/**: Deep-dive case studies
6. **6_ppi_analysis/**: Extended PPI analysis

## Development Environment Setup

### Primary Environment: `varchamp`
```bash
# Create and activate the main conda environment
conda env create -f 2_individual_assay_results/imaging/2_analysis/1_snakemake_pipeline/2025_varchamp_snakemake/env.yml
conda activate varchamp
```

### Specialized Environment: `ai4protein` 
For AI predictor analyses (section 4):
```bash
# Create environment for AI protein analyses
conda env create -f 4_compare_ai_scores/2_analysis/env.yml
conda activate ai4protein
```

### Cytoself Environment
For advanced imaging analysis:
```bash
cd 2_individual_assay_results/imaging/2_analysis/3_cytoself
./install_env.sh
conda activate cytoself
```

## Common Development Commands

### Jupyter Notebooks
Most analyses are implemented as Jupyter notebooks:
```bash
# Start Jupyter server for analysis
jupyter notebook
# or for lab interface
jupyter lab
```

### Snakemake Pipeline (Imaging Analysis)
```bash
cd 2_individual_assay_results/imaging/2_analysis/1_snakemake_pipeline/2025_varchamp_snakemake
# Run the complete imaging analysis pipeline
./run_snakemake_pipeline.sh
```

### R Statistical Analysis
```bash
# Run LLR threshold calibration
cd 2_individual_assay_results/threshold_calibration
./runLLR.sh  # Standard run
./runLLR_max.sh  # Maximum analysis
```

### Data Download Scripts
```bash
# Download external datasets
cd 4_compare_ai_scores/1_inputs && ./download_data.sh
cd 4_compare_ai_scores/3_outputs && ./download_processed_data.sh
cd 2_individual_assay_results/dual_ipa/1_inputs && ./download_raw_FACS.sh
```

## Key Architecture Components

### Analysis Pipeline Architecture
- **Input-Processing-Output Pattern**: Each section follows a consistent data flow from raw inputs through analysis scripts to processed outputs
- **Cross-section Dependencies**: Outputs from individual assays (section 2) feed into integrated analyses (section 3)
- **Manuscript Generation**: Figures and tables are generated directly from analysis outputs

### Technology Stack
- **Python**: Primary analysis language with pandas, numpy, matplotlib, seaborn
- **R**: Statistical modeling and LLR calculations
- **Jupyter**: Interactive analysis and visualization
- **Snakemake**: Workflow management for imaging pipeline
- **Conda**: Environment and dependency management
- **CUDA**: GPU acceleration for deep learning components

### Data Processing Patterns
- **Assay-specific Processing**: Each assay type has specialized processing pipelines
- **Statistical Calibration**: Threshold determination using log-likelihood ratios
- **Cross-assay Integration**: Standardized scoring across different assay modalities
- **AI Benchmarking**: Systematic comparison with computational predictors

## Git Workflow

Follow the standard fork-and-branch workflow:
```bash
# Clone your fork
git clone https://github.com/<YOUR-USERNAME>/2025_laval_submitted.git
cd 2025_laval_submitted
git checkout -b descriptive-branch-name

# Make changes, then commit and push
git add .
git commit -m "descriptive commit message"
git push  # Follow Git's suggestions for first push

# Create pull request through GitHub web interface
```

## Important Development Notes

- **Environment Isolation**: Use appropriate conda environments for different analysis sections
- **Data Dependencies**: Check that required input data exists before running analyses
- **Computational Requirements**: Imaging and AI analyses require significant computational resources
- **Reproducibility**: All analyses should produce consistent outputs across runs
- **Documentation**: Each analysis notebook should include clear explanations and parameter descriptions

## File Naming Conventions

- Analysis scripts are numbered sequentially (e.g., `1_process_data.ipynb`, `2_analyze_results.ipynb`)
- Output files include descriptive names indicating content and processing stage
- Environment files are named `env.yml` or include specific requirements (`env-reqs.txt`)
- Shell scripts use `.sh` extension and are often prefixed with `run_` or include specific actions (`download_`, `install_`)