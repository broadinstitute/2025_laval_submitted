#!/home/shenrunx/software/anaconda3/envs/ai4protein/bin/python3
"""
G2P Protein Features Mapping with Auto-Download and Detailed Failure Tracking
==============================================================================

This script maps G2P protein features to variants with:
- Automatic detection of existing G2P feature files
- Automatic download from G2P API for missing genes
- Comprehensive failure tracking for debugging
- Detailed summary reports

Usage:
    ./2_mapping_g2p_feats.py \
        --input_file "../../3_outputs/VarChAMP_Allele_Collection_1Perc_ClinVar_dbNSFP_pLDDT_RSA.parquet" \
        --output_file "../../3_outputs/VarChAMP_Allele_Collection_1Perc_ClinVar_dbNSFP_pLDDT_RSA_G2P.parquet" \
        --failure_report "../../3_outputs/g2p_mapping_failures.tsv"

Author: Runxi Shen (v2 with auto-download)
Created: 2025-01-07
"""

import os
import sys
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import time
from tqdm import tqdm
from typing import Dict, Tuple, List
import warnings
import re
import requests
import gzip
warnings.filterwarnings('ignore')

# Default settings
G2P_FEATURES_DIR = "../../3_outputs/processed_data/g2p_protein_features"


def three_letter_to_one_letter(three_letter: str) -> str:
    """Convert three-letter amino acid code to one-letter code."""
    aa_map = {
        'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
        'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
        'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
        'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
    }
    return aa_map.get(three_letter, '')


def one_letter_to_three_letter(one_letter: str) -> str:
    """Convert one-letter amino acid code to three-letter code."""
    aa_map = {
        'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
        'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
        'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
        'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'
    }
    return aa_map.get(one_letter, '')


def parse_aa_change(aa_change: str) -> Tuple[str, int, str]:
    """Parse amino acid change string like 'Leu359Pro' into (ref_aa, position, alt_aa)."""
    if pd.isna(aa_change) or not aa_change:
        return '', 0, ''

    match = re.match(r'^([A-Za-z]{3})(\d+)([A-Za-z]{3})$', str(aa_change))
    if match:
        ref_aa_3letter, pos, alt_aa_3letter = match.groups()
        ref_aa_1letter = three_letter_to_one_letter(ref_aa_3letter)
        return ref_aa_1letter, int(pos), alt_aa_3letter

    return '', 0, ''


def get_all_uniprot_swissprot_ids(protein_name: str) -> List[str]:
    """Query UniProt's REST search API to find ALL reviewed (Swiss-Prot) accessions."""
    url = "https://rest.uniprot.org/uniprotkb/search"
    query = f"gene:{protein_name} AND organism_id:9606 AND reviewed:true"
    params = {
        "query": query,
        "fields": "accession",
        "format": "json",
        "size": 500,
    }
    try:
        resp = requests.get(url, params=params, timeout=10)
        resp.raise_for_status()
        data = resp.json()
        results = data.get("results", [])
        return [result["primaryAccession"] for result in results] if results else []
    except Exception as e:
        print(f"  Error querying UniProt for {protein_name}: {e}")
        return []


def download_g2p_protein_features(gene_name: str, uniprot_id: str, output_dir: str) -> str:
    """Download protein features from G2P portal API."""
    url = f"https://g2p.broadinstitute.org/api/gene/{gene_name}/protein/{uniprot_id}/protein-features"
    headers = {'accept': 'application/csv'}

    try:
        response = requests.get(url, headers=headers, timeout=30)
        response.raise_for_status()

        if response.content:
            filename = f"G2P_{gene_name}_{uniprot_id}_protein_features.tsv"
            output_file = Path(output_dir) / filename

            # Handle gzipped content
            if response.headers.get('content-encoding') == 'gzip':
                try:
                    decompressed_content = gzip.decompress(response.content).decode('utf-8')
                    with open(output_file, 'w') as f:
                        f.write(decompressed_content)
                except gzip.BadGzipFile:
                    with open(output_file, 'wb') as f:
                        f.write(response.content)
            else:
                with open(output_file, 'wb') as f:
                    f.write(response.content)

            # Verify file has content
            if output_file.stat().st_size > 100:
                return str(output_file)

        return None

    except Exception as e:
        print(f"  Error downloading G2P features for {gene_name} ({uniprot_id}): {e}")
        return None


def find_or_download_g2p_features(gene_name: str, g2p_dir: str) -> str:
    """
    Find existing G2P features file for a gene, or download if not found.

    Returns:
        Path to G2P features file, or None if not available
    """
    g2p_path = Path(g2p_dir)

    # Check if file already exists
    existing_files = list(g2p_path.glob(f"G2P_{gene_name}_*_protein_features.tsv"))

    if existing_files:
        return str(existing_files[0])

    # File doesn't exist - try to download
    print(f"  G2P features not found for {gene_name}, attempting download...")

    # Get UniProt IDs
    uniprot_ids = get_all_uniprot_swissprot_ids(gene_name)

    if not uniprot_ids:
        print(f"  No UniProt IDs found for {gene_name}")
        return None

    print(f"  Found {len(uniprot_ids)} UniProt IDs: {uniprot_ids[:3]}...")

    # Try each UniProt ID
    for uniprot_id in uniprot_ids:
        print(f"  Trying {uniprot_id}...")
        protein_file = download_g2p_protein_features(gene_name, uniprot_id, g2p_dir)

        if protein_file:
            print(f"  ✓ Successfully downloaded G2P features for {gene_name}")
            time.sleep(0.5)  # Be respectful to API
            return protein_file

        time.sleep(0.5)

    print(f"  ✗ Could not download G2P features for {gene_name}")
    return None


def load_g2p_protein_features_with_autodownload(g2p_features_dir: str,
                                                 genes_needed: List[str]) -> Tuple[Dict[str, pd.DataFrame], List[str], Dict[str, str]]:
    """
    Load G2P protein features, automatically downloading missing genes.

    Returns:
        - Dictionary mapping gene symbol to protein features DataFrame
        - List of all feature column names (with g2p_ prefix)
        - Dictionary mapping gene symbol to UniProt ID used
    """
    g2p_data = {}
    gene_to_uniprot = {}
    all_feature_cols = set()
    features_dir = Path(g2p_features_dir)
    features_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nLoading/downloading G2P protein features...")

    # Track download statistics
    already_existed = 0
    newly_downloaded = 0
    failed_download = 0

    for gene in tqdm(genes_needed, desc="Processing genes"):
        # Find or download G2P features
        protein_file = find_or_download_g2p_features(gene, g2p_features_dir)

        if not protein_file:
            failed_download += 1
            continue

        # Extract UniProt ID from filename: G2P_{gene_name}_{uniprot_id}_protein_features.tsv
        filename = Path(protein_file).name
        match = re.match(rf'G2P_{re.escape(gene)}_([A-Z0-9]+)_protein_features\.tsv', filename)
        if match:
            gene_to_uniprot[gene] = match.group(1)

        # Check if this was newly downloaded
        if "downloaded" in protein_file or time.time() - Path(protein_file).stat().st_mtime < 3600:
            # File modified in last hour - likely newly downloaded
            if Path(protein_file).stat().st_mtime > time.time() - 60:
                newly_downloaded += 1
            else:
                already_existed += 1
        else:
            already_existed += 1

        try:
            # Load features
            df = pd.read_csv(protein_file, sep='\t')

            # Validate required columns
            if 'residueId' not in df.columns or 'AA' not in df.columns:
                continue

            # Set residueId as index for fast lookups
            df = df.set_index('residueId')

            # Store in dictionary
            g2p_data[gene] = df

            # Collect all feature column names (excluding AA since we use it for validation)
            feature_cols = [col for col in df.columns if col != 'AA']
            all_feature_cols.update(feature_cols)

        except Exception as e:
            print(f"  Error loading {protein_file}: {e}")

    # Create sorted list of feature columns with g2p_ prefix
    feature_columns = sorted([f'g2p_{col}' for col in all_feature_cols])

    print(f"\n📊 G2P Features Loading Summary:")
    print(f"  Already existed: {already_existed}")
    print(f"  Newly downloaded: {newly_downloaded}")
    print(f"  Failed to download: {failed_download}")
    print(f"  Total loaded: {len(g2p_data)} genes")
    print(f"  Feature columns: {len(feature_columns)}")

    return g2p_data, feature_columns, gene_to_uniprot


def map_g2p_features_with_tracking(variants_df: pd.DataFrame,
                                   g2p_data: Dict[str, pd.DataFrame],
                                   feature_columns: List[str],
                                   gene_to_uniprot: Dict[str, str]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Map G2P protein features to variants with comprehensive failure tracking.

    Returns:
        - Enhanced variants DataFrame
        - Failures DataFrame for debugging
    """
    print("\nParsing amino acid changes...")

    # Initialize tracking lists
    failure_records = []

    # Parse amino acid changes
    variants_df['ref_aa_1letter'] = ''
    variants_df['position'] = 0
    variants_df['alt_aa_3letter'] = ''

    for idx, row in variants_df.iterrows():
        ref_aa, position, alt_aa = parse_aa_change(row['aa_change'])
        variants_df.at[idx, 'ref_aa_1letter'] = ref_aa
        variants_df.at[idx, 'position'] = position
        variants_df.at[idx, 'alt_aa_3letter'] = alt_aa

    # Initialize all feature columns with empty strings
    print("Initializing feature columns...")
    for col in feature_columns:
        variants_df[col] = ''

    # Process by gene for efficiency
    print("\nMapping G2P features by gene...")
    genes_with_g2p = set(g2p_data.keys())
    genes_in_variants = set(variants_df['symbol'].dropna().unique())
    genes_to_process = genes_with_g2p & genes_in_variants

    print(f"Processing {len(genes_to_process)} genes with G2P data out of {len(genes_in_variants)} total genes")

    matched_count = 0
    total_processable = 0

    # Track failure reasons
    failure_stats = {
        'unparseable_aa_change': 0,
        'gene_missing_g2p': 0,
        'position_not_in_g2p': 0,
        'residue_mismatch': 0
    }

    for gene in tqdm(genes_to_process, desc="Processing genes"):
        # Get variants for this gene
        gene_mask = variants_df['symbol'] == gene
        gene_indices = variants_df.index[gene_mask]

        # Get G2P features for this gene
        g2p_df = g2p_data[gene]

        # Process each variant for this gene
        for idx in gene_indices:
            variant_id = variants_df.at[idx, 'ccsb_mutation_id'] if 'ccsb_mutation_id' in variants_df.columns else idx
            aa_change = variants_df.at[idx, 'aa_change']
            ref_aa = variants_df.at[idx, 'ref_aa_1letter']
            position = variants_df.at[idx, 'position']

            # Skip if parsing failed
            if not ref_aa or position == 0:
                failure_records.append({
                    'ccsb_mutation_id': variant_id,
                    'gene': gene,
                    'uniprot_id': gene_to_uniprot.get(gene, ''),
                    'aa_change': aa_change,
                    'G2P_ref_aa': '',
                    'failure_reason': 'unparseable_aa_change',
                    'details': f'Could not parse aa_change: {aa_change}'
                })
                failure_stats['unparseable_aa_change'] += 1
                continue

            total_processable += 1

            # Check if position exists in G2P data
            if position not in g2p_df.index:
                failure_records.append({
                    'ccsb_mutation_id': variant_id,
                    'gene': gene,
                    'uniprot_id': gene_to_uniprot.get(gene, ''),
                    'aa_change': aa_change,
                    'G2P_ref_aa': '',
                    'failure_reason': 'position_not_in_g2p',
                    'details': f'Position {position} not found in G2P data (protein length: {g2p_df.index.max()})'
                })
                failure_stats['position_not_in_g2p'] += 1
                continue

            # Get the residue row
            residue_row = g2p_df.loc[position]

            # Validate amino acid matches
            protein_aa = residue_row.get('AA', '')
            if protein_aa != ref_aa:
                # Convert G2P's one-letter AA to three-letter code
                g2p_ref_aa_3letter = one_letter_to_three_letter(protein_aa)
                failure_records.append({
                    'ccsb_mutation_id': variant_id,
                    'gene': gene,
                    'uniprot_id': gene_to_uniprot.get(gene, ''),
                    'aa_change': aa_change,
                    'G2P_ref_aa': g2p_ref_aa_3letter,
                    'failure_reason': 'residue_mismatch',
                    'details': f'Expected {ref_aa} at position {position}, but G2P has {protein_aa}'
                })
                failure_stats['residue_mismatch'] += 1
                continue

            # Match found - extract all features
            matched_count += 1

            # Map features to variant
            for orig_col in g2p_df.columns:
                if orig_col == 'AA':
                    continue

                g2p_col = f'g2p_{orig_col}'
                if g2p_col in feature_columns:
                    val = residue_row[orig_col]
                    variants_df.at[idx, g2p_col] = str(val) if pd.notna(val) else ''

    # Add failures for genes without G2P data
    genes_missing_g2p = genes_in_variants - genes_with_g2p
    for gene in genes_missing_g2p:
        gene_variants = variants_df[variants_df['symbol'] == gene]
        for idx, row in gene_variants.iterrows():
            variant_id = row['ccsb_mutation_id'] if 'ccsb_mutation_id' in variants_df.columns else idx
            failure_records.append({
                'ccsb_mutation_id': variant_id,
                'gene': gene,
                'uniprot_id': '',
                'aa_change': row['aa_change'],
                'G2P_ref_aa': '',
                'failure_reason': 'gene_missing_g2p',
                'details': 'No G2P protein features file found for this gene'
            })
            failure_stats['gene_missing_g2p'] += 1

    # Create failures DataFrame
    failures_df = pd.DataFrame(failure_records)

    print(f"\n✓ Mapping complete:")
    print(f"  Annotated: {matched_count}/{total_processable} parseable variants")
    print(f"  Match rate: {matched_count/total_processable:.2%}" if total_processable > 0 else "  No parseable variants")

    print(f"\n📊 Failure Statistics:")
    print(f"  Total failures: {len(failure_records)}")
    for reason, count in failure_stats.items():
        print(f"  - {reason}: {count}")

    # Drop temporary columns
    variants_df = variants_df.drop(columns=['ref_aa_1letter', 'position', 'alt_aa_3letter'])

    return variants_df, failures_df


def main():
    start_time = time.time()

    print("="*70)
    print("G2P PROTEIN FEATURES MAPPING (v2 with Auto-Download)")
    print("="*70)

    # Read input data
    print(f"\nReading input data from: {INPUT_FILE}")
    variants_df = pd.read_parquet(INPUT_FILE)
    print(f"Loaded {len(variants_df)} variants")

    # Get unique genes
    genes_needed = list(variants_df['symbol'].dropna().unique())
    print(f"Unique genes: {len(genes_needed)}")

    # Load/download G2P protein features
    g2p_data, feature_columns, gene_to_uniprot = load_g2p_protein_features_with_autodownload(
        G2P_FEATURES_DIR,
        genes_needed
    )

    if not g2p_data:
        print("ERROR: No G2P features loaded. Exiting.")
        sys.exit(1)

    # Map G2P features with failure tracking
    enhanced_df, failures_df = map_g2p_features_with_tracking(
        variants_df,
        g2p_data,
        feature_columns,
        gene_to_uniprot
    )

    # Count variants with features
    variants_with_features = 0
    if feature_columns:
        variants_with_features = len(enhanced_df[
            enhanced_df[feature_columns].apply(
                lambda row: any(str(val).strip() != '' for val in row), axis=1
            )
        ])

    # Save output
    print(f"\nWriting results to: {OUTPUT_FILE}")
    enhanced_df.to_parquet(OUTPUT_FILE, index=False)

    # Save failure report if requested
    if FAILURE_REPORT and not failures_df.empty:
        print(f"Writing failure report to: {FAILURE_REPORT}")
        failures_df.to_csv(FAILURE_REPORT, sep='\t', index=False)

    end_time = time.time()

    # Print summary
    print(f"\n{'='*70}")
    print("PROCESSING SUMMARY")
    print(f"{'='*70}")
    print(f"Total variants processed: {len(variants_df)}")
    print(f"Variants with G2P features: {variants_with_features}")
    print(f"Coverage rate: {variants_with_features/len(variants_df):.2%}")
    print(f"Failed variants: {len(failures_df)}")
    print(f"G2P feature columns added: {len(feature_columns)}")
    print(f"Processing time: {end_time - start_time:.2f} seconds")
    print(f"\nOutput files:")
    print(f"  Data: {OUTPUT_FILE}")
    if FAILURE_REPORT:
        print(f"  Failures: {FAILURE_REPORT}")
    print("="*70)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Map G2P protein features to variants with auto-download and failure tracking"
    )
    parser.add_argument("--input_file", type=str, required=True,
                        help="Input parquet file with variants")
    parser.add_argument("--output_file", type=str, required=True,
                        help="Output parquet file with G2P features")
    parser.add_argument("--g2p_dir", type=str, default=G2P_FEATURES_DIR,
                        help=f"Directory containing G2P protein features (default: {G2P_FEATURES_DIR})")
    parser.add_argument("--failure_report", type=str, default=None,
                        help="Output TSV file for detailed failure report (optional)")

    args = parser.parse_args()

    # Set global variables from args
    INPUT_FILE = args.input_file
    OUTPUT_FILE = args.output_file
    G2P_FEATURES_DIR = args.g2p_dir
    FAILURE_REPORT = args.failure_report

    main()
