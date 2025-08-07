"""
Comprehensive Protein Features Extraction Pipeline
==================================================

This script processes all genes/variants from the integrated assay results,
extracts G2P protein features with multi-UniProt fallback, validates amino acid matches,
and directly appends all protein feature data to each variant.

Combines all previous functionality into a single comprehensive script.

Author: Runxi Shen  
Created: 2025-01-XX
"""

import requests
import pandas as pd
import numpy as np
import re
import os
import gzip
from pathlib import Path
import time
import json
from typing import List, Dict, Tuple, Optional, Set
import warnings
warnings.filterwarnings('ignore')


def get_all_uniprot_swissprot_ids(protein_name: str) -> List[str]:
    """
    Query UniProt's REST search API to find ALL reviewed (Swiss-Prot) accessions
    for a given human gene/protein name.
    """
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
        print(f"Error querying UniProt for {protein_name}: {e}")
        return []


def extract_g2p_protein_features(gene_name: str, uniprot_id: str, output_dir: str) -> Optional[str]:
    """
    Extract protein features from G2P portal API and return the file path if successful.
    """
    url = f"https://g2p.broadinstitute.org/api/gene/{gene_name}/protein/{uniprot_id}/protein-features"
    headers = {'accept': 'application/csv'}
    
    try:
        response = requests.get(url, headers=headers, timeout=30)
        response.raise_for_status()
        
        if response.content:
            # Extract filename from content-disposition header
            content_disposition = response.headers.get('content-disposition', '')
            if 'filename=' in content_disposition and '"' in content_disposition:
                try:
                    filename = content_disposition.split('filename="')[1].split('"')[0]
                except IndexError:
                    filename = f"G2P_{gene_name}_{uniprot_id}_protein_features.tsv"
            else:
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
            
            return str(output_file)
        
    except Exception as e:
        print(f"Error extracting protein features for {gene_name} ({uniprot_id}): {e}")
    
    return None


def parse_aa_change(aa_change: str) -> Tuple[str, int, str]:
    """Parse amino acid change string like 'Leu359Pro' into (ref_aa, position, alt_aa)."""
    if pd.isna(aa_change) or not aa_change:
        return None, None, None
    
    match = re.match(r'^([A-Za-z]{3})(\d+)([A-Za-z]{3})$', aa_change)
    if match:
        ref_aa, pos, alt_aa = match.groups()
        return ref_aa, int(pos), alt_aa
    
    return None, None, None


def three_letter_to_one_letter(three_letter: str) -> str:
    """Convert three-letter amino acid code to one-letter code."""
    aa_map = {
        'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
        'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
        'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
        'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
    }
    return aa_map.get(three_letter, three_letter)


def validate_variants_against_protein_features(variants_df: pd.DataFrame, protein_features_file: str) -> Dict:
    """Validate that all variants' reference amino acids match the protein features file."""
    if not Path(protein_features_file).exists():
        return {"valid": False, "error": f"File not found: {protein_features_file}"}
    
    try:
        prot_df = pd.read_csv(protein_features_file, sep='\t')
        pos_to_aa = dict(zip(prot_df['residueId'], prot_df['AA']))
        
        mismatches = []
        valid_count = 0
        
        for _, variant_row in variants_df.iterrows():
            aa_change = variant_row['aa_change']
            ref_aa_3letter, position, alt_aa_3letter = parse_aa_change(aa_change)
            
            if ref_aa_3letter is None:
                continue
            
            ref_aa_1letter = three_letter_to_one_letter(ref_aa_3letter)
            
            if position in pos_to_aa:
                protein_aa = pos_to_aa[position]
                if protein_aa == ref_aa_1letter:
                    valid_count += 1
                else:
                    mismatches.append({
                        'variant': aa_change,
                        'position': position,
                        'expected_aa': ref_aa_1letter,
                        'protein_aa': protein_aa
                    })
            else:
                mismatches.append({
                    'variant': aa_change,
                    'position': position,
                    'expected_aa': ref_aa_1letter,
                    'protein_aa': 'POSITION_NOT_FOUND'
                })
        
        return {
            "valid": len(mismatches) == 0,
            "total_variants": len(variants_df),
            "valid_matches": valid_count,
            "mismatches": mismatches,
            "protein_length": len(pos_to_aa)
        }
        
    except Exception as e:
        return {"valid": False, "error": f"Error reading protein features: {e}"}


def extract_protein_features_for_variants(variants_df: pd.DataFrame, protein_features_file: str) -> pd.DataFrame:
    """
    Extract all protein features for each variant position and append to variants dataframe.
    """
    enhanced_df = variants_df.copy()
    
    if not Path(protein_features_file).exists():
        print(f"Warning: Protein features file not found: {protein_features_file}")
        return enhanced_df
    
    try:
        # Load protein features
        prot_df = pd.read_csv(protein_features_file, sep='\t')
        
        # Create mapping from position to all features
        protein_features_dict = {}
        for _, row in prot_df.iterrows():
            pos = int(row['residueId'])
            # Convert all feature values to strings to handle mixed types
            features = {col: str(val) if pd.notna(val) else '' for col, val in row.items() if col != 'residueId'}
            protein_features_dict[pos] = features
        
        # Get all possible feature columns
        all_feature_columns = set()
        for features in protein_features_dict.values():
            all_feature_columns.update(features.keys())
        all_feature_columns = sorted(all_feature_columns)
        
        # Initialize new columns with empty strings
        for col in all_feature_columns:
            enhanced_df[f'protein_{col}'] = ''
        
        # Add protein features for each variant
        for idx, row in enhanced_df.iterrows():
            aa_change = row['aa_change']
            ref_aa_3letter, position, alt_aa_3letter = parse_aa_change(aa_change)
            
            if position and position in protein_features_dict:
                features = protein_features_dict[position]
                for col in all_feature_columns:
                    enhanced_df.loc[idx, f'protein_{col}'] = features.get(col, '')
        
        print(f"  ✓ Added {len(all_feature_columns)} protein feature columns")
        return enhanced_df
        
    except Exception as e:
        print(f"  ✗ Error extracting protein features: {e}")
        return enhanced_df


def process_gene_with_fallback(gene_name: str, variants_df: pd.DataFrame, output_dir: str) -> Dict:
    """
    Process a gene with multiple UniProt ID fallback logic and extract detailed features.
    """
    print(f"\n{'='*70}")
    print(f"Processing gene: {gene_name} ({len(variants_df)} variants)")
    print(f"{'='*70}")
    
    # Get all possible UniProt IDs
    uniprot_ids = get_all_uniprot_swissprot_ids(gene_name)
    if not uniprot_ids:
        return {
            "gene": gene_name,
            "status": "failed",
            "reason": "no_uniprot_ids",
            "enhanced_variants": variants_df
        }
    
    print(f"Found {len(uniprot_ids)} UniProt IDs: {uniprot_ids}")
    
    # Try each UniProt ID until we find a perfect match
    for i, uniprot_id in enumerate(uniprot_ids):
        print(f"\n--- Trying UniProt ID {i+1}/{len(uniprot_ids)}: {uniprot_id} ---")
        
        # Extract protein features
        protein_file = extract_g2p_protein_features(gene_name, uniprot_id, output_dir)
        if not protein_file:
            continue
        
        # Validate variants
        validation_result = validate_variants_against_protein_features(variants_df, protein_file)
        
        print(f"  Validation: {'✓ Valid' if validation_result['valid'] else '✗ Invalid'}")
        print(f"  Matches: {validation_result.get('valid_matches', 0)}/{validation_result.get('total_variants', 0)}")
        
        if validation_result['valid']:
            # Perfect match found - extract detailed features
            print(f"  ✓ Perfect match! Extracting detailed protein features...")
            enhanced_variants = extract_protein_features_for_variants(variants_df, protein_file)
            
            return {
                "gene": gene_name,
                "status": "success",
                "uniprot_id": uniprot_id,
                "protein_file": protein_file,
                "validation_result": validation_result,
                "enhanced_variants": enhanced_variants
            }
        else:
            # Remove invalid file
            if Path(protein_file).exists():
                Path(protein_file).unlink()
        
        time.sleep(1)  # Be respectful to API
    
    return {
        "gene": gene_name,
        "status": "failed",
        "reason": "no_valid_uniprot_match",
        "enhanced_variants": variants_df
    }


def save_progress(results: Dict, progress_file: str):
    """Save processing progress to JSON file."""
    # Convert results to JSON-serializable format
    json_results = {}
    for gene, result in results.items():
        json_results[gene] = {
            "gene": result["gene"],
            "status": result["status"],
            "reason": result.get("reason", ""),
            "uniprot_id": result.get("uniprot_id", ""),
            "protein_file": result.get("protein_file", "")
        }
    
    with open(progress_file, 'w') as f:
        json.dump(json_results, f, indent=2)


def load_progress(progress_file: str) -> Dict:
    """Load processing progress from JSON file."""
    if Path(progress_file).exists():
        with open(progress_file, 'r') as f:
            return json.load(f)
    return {}


def main():
    """
    Main function to process all genes and extract comprehensive protein features.
    """
    print("Comprehensive Protein Features Extraction Pipeline")
    print("=" * 70)
    
    # Set up paths
    base_dir = Path(__file__).parent.parent.parent.parent
    input_csv = base_dir / "3_integrated_assay_results/3_outputs/all_gene_alleles_profiled_by_any_assay.csv"
    output_dir = base_dir / "4_compare_ai_scores/3_outputs/processed_data/g2p_protein_features"
    progress_file = base_dir / "4_compare_ai_scores/3_outputs/comprehensive_features_progress.json"
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Read input data
    print(f"\nReading input data from: {input_csv}")
    df = pd.read_csv(input_csv)
    unique_genes = df['symbol'].unique()
    
    print(f"Loaded {len(df)} variants for {len(unique_genes)} unique genes")
    
    # Load existing progress
    json_results = load_progress(progress_file)
    processed_genes = set(json_results.keys())
    print(f"Found {len(processed_genes)} previously processed genes")
    
    # Process remaining genes
    remaining_genes = [gene for gene in unique_genes if gene not in processed_genes]
    print(f"Processing {len(remaining_genes)} remaining genes")
    
    all_enhanced_variants = []
    full_results = {}
    
    for i, gene in enumerate(remaining_genes, 1):
        print(f"\n[{len(processed_genes) + i}/{len(unique_genes)}] Processing gene: {gene}")
        
        # Get variants for this gene
        gene_variants = df[df['symbol'] == gene].copy()
        
        # Process gene with fallback logic and feature extraction
        result = process_gene_with_fallback(gene, gene_variants, str(output_dir))
        full_results[gene] = result
        
        # Collect enhanced variants
        all_enhanced_variants.append(result['enhanced_variants'])
        
        # Print status
        if result['status'] == 'success':
            print(f"✓ Successfully processed {gene} with UniProt ID {result['uniprot_id']}")
        else:
            print(f"✗ Failed to process {gene}: {result['reason']}")
        
        # Save progress periodically
        if i % 10 == 0:
            save_progress(full_results, progress_file)
            print(f"  Progress saved ({len(processed_genes) + i}/{len(unique_genes)} genes)")
        
        time.sleep(1)  # Be respectful to APIs
    
    # Process previously successful genes (load their enhanced variants)
    print(f"\nProcessing {len(processed_genes)} previously successful genes...")
    for gene in processed_genes:
        gene_result = json_results[gene]
        gene_variants = df[df['symbol'] == gene].copy()
        
        if gene_result['status'] == 'success' and gene_result.get('protein_file'):
            protein_file = gene_result['protein_file']
            if Path(protein_file).exists():
                enhanced_variants = extract_protein_features_for_variants(gene_variants, protein_file)
                all_enhanced_variants.append(enhanced_variants)
                print(f"  ✓ Loaded features for {gene}")
            else:
                print(f"  ✗ Protein file missing for {gene}: {protein_file}")
                all_enhanced_variants.append(gene_variants)
        else:
            all_enhanced_variants.append(gene_variants)
    
    # Combine all enhanced variants
    print(f"\nCombining enhanced variants...")
    final_enhanced_df = pd.concat(all_enhanced_variants, ignore_index=True)
    
    # Save final results
    save_progress(full_results, progress_file)
    
    output_csv = base_dir / "4_compare_ai_scores/3_outputs/all_gene_alleles_with_comprehensive_protein_features.csv"
    final_enhanced_df.to_csv(output_csv, index=False)
    
    # Generate final statistics
    all_results = {**json_results, **full_results}
    successful_genes = [g for g, r in all_results.items() if r['status'] == 'success']
    failed_genes = [g for g, r in all_results.items() if r['status'] == 'failed']
    
    # Count variants with protein features
    protein_feature_cols = [col for col in final_enhanced_df.columns if col.startswith('protein_')]
    variants_with_features = 0
    if protein_feature_cols:
        # Count variants that have at least one non-empty protein feature
        variants_with_features = len(final_enhanced_df[
            final_enhanced_df[protein_feature_cols].apply(
                lambda row: any(str(val).strip() != '' for val in row), axis=1
            )
        ])
    
    print(f"\n{'='*70}")
    print("COMPREHENSIVE PROTEIN FEATURES SUMMARY")
    print(f"{'='*70}")
    print(f"Total genes: {len(unique_genes)}")
    print(f"Successfully processed: {len(successful_genes)}")
    print(f"Failed: {len(failed_genes)}")
    print(f"Success rate: {len(successful_genes)/len(unique_genes):.2%}")
    print(f"\nTotal variants: {len(final_enhanced_df)}")
    print(f"Variants with protein features: {variants_with_features}")
    print(f"Protein feature columns added: {len(protein_feature_cols)}")
    print(f"Coverage rate: {variants_with_features/len(final_enhanced_df):.2%}")
    
    print(f"\nOutput saved to: {output_csv}")
    print(f"Progress saved to: {progress_file}")
    print(f"\nSample successful genes: {successful_genes[:10]}")
    
    if failed_genes:
        print(f"\nFailed genes: {failed_genes[:10]}")


if __name__ == "__main__":
    main()