#!/usr/bin/env python3
"""
Notebook-style analysis of subcellular mislocalization patterns
"""

import pandas as pd
import numpy as np
from collections import defaultdict, Counter

def analyze_mislocalization_patterns(csv_path):
    """
    Analyze mislocalization patterns by comparing variant alleles to their reference genes
    """
    # Load data
    df = pd.read_csv(csv_path)
    
    # Get location columns
    loc_cols = [col for col in df.columns if col.startswith('loc_')]
    print(f"Found {len(loc_cols)} location categories")
    
    # Separate reference genes and variants
    df['gene_name'] = df['Metadata_gene_allele'].apply(lambda x: x.split('_')[0] if '_' in x else x)
    df['is_reference'] = ~df['Metadata_gene_allele'].str.contains('_')
    
    reference_genes = df[df['is_reference']].copy()
    variant_alleles = df[~df['is_reference']].copy()
    
    print(f"Reference genes: {len(reference_genes)}")
    print(f"Variant alleles: {len(variant_alleles)}")
    
    # Create reference gene lookup
    ref_patterns = {}
    for _, row in reference_genes.iterrows():
        gene_name = row['gene_name']
        pattern = tuple(row[loc_cols].values)
        ref_patterns[gene_name] = pattern
    
    print(f"Found reference patterns for {len(ref_patterns)} genes")
    
    # Analyze variants
    mislocalization_data = []
    
    for _, variant in variant_alleles.iterrows():
        gene_name = variant['gene_name']
        variant_name = variant['Metadata_gene_allele']
        variant_pattern = tuple(variant[loc_cols].values)
        
        if gene_name in ref_patterns:
            ref_pattern = ref_patterns[gene_name]
            
            # Check if patterns differ
            has_mislocalization = variant_pattern != ref_pattern
            
            mislocalization_data.append({
                'gene': gene_name,
                'variant': variant_name,
                'has_mislocalization': has_mislocalization,
                'reference_pattern': ref_pattern,
                'variant_pattern': variant_pattern
            })
    
    return mislocalization_data, loc_cols

def calculate_pattern_changes(mislocalization_data, loc_cols):
    """
    Calculate unique pattern changes and their frequencies
    """
    # Count alleles with changes
    total_variants = len(mislocalization_data)
    changed_variants = sum(1 for item in mislocalization_data if item['has_mislocalization'])
    
    print(f"Total variants analyzed: {total_variants}")
    print(f"Variants with mislocalization: {changed_variants}")
    print(f"Percentage with changes: {changed_variants/total_variants*100:.1f}%")
    
    # Analyze pattern changes for variants with mislocalization
    pattern_changes = []
    
    for item in mislocalization_data:
        if item['has_mislocalization']:
            ref_pattern = np.array(item['reference_pattern'])
            var_pattern = np.array(item['variant_pattern'])
            
            # Calculate change vector (gained - lost locations)
            change_vector = var_pattern - ref_pattern
            
            # Create change description
            gained_locs = [loc_cols[i].replace('loc_', '') for i, val in enumerate(change_vector) if val > 0]
            lost_locs = [loc_cols[i].replace('loc_', '') for i, val in enumerate(change_vector) if val < 0]
            
            change_desc = {
                'variant': item['variant'],
                'gene': item['gene'],
                'gained': gained_locs,
                'lost': lost_locs,
                'change_vector': tuple(change_vector),
                'n_gained': len(gained_locs),
                'n_lost': len(lost_locs)
            }
            
            pattern_changes.append(change_desc)
    
    return pattern_changes

def summarize_unique_patterns(pattern_changes):
    """
    Find unique mislocalization patterns and their frequencies
    """
    # Group by change vector (unique patterns)
    pattern_counts = Counter(item['change_vector'] for item in pattern_changes)
    
    print(f"\nFound {len(pattern_counts)} unique mislocalization patterns")
    print(f"Top 10 most frequent patterns:")
    
    for i, (pattern, count) in enumerate(pattern_counts.most_common(10)):
        print(f"{i+1:2d}. Pattern occurs {count:3d} times: {pattern}")
    
    # Group by semantic changes (gained/lost combinations)
    semantic_changes = defaultdict(list)
    
    for change in pattern_changes:
        gained_str = ";".join(sorted(change['gained'])) if change['gained'] else "none"
        lost_str = ";".join(sorted(change['lost'])) if change['lost'] else "none"
        change_key = f"Gained: {gained_str} | Lost: {lost_str}"
        semantic_changes[change_key].append(change['variant'])
    
    print(f"\nSemantic change patterns:")
    for i, (change_desc, variants) in enumerate(sorted(semantic_changes.items(), key=lambda x: len(x[1]), reverse=True)[:15]):
        print(f"{i+1:2d}. {change_desc}")
        print(f"    Count: {len(variants)}, Examples: {', '.join(variants[:3])}")
        print()
    
    return pattern_counts, semantic_changes

# Main analysis
if __name__ == "__main__":
    csv_path = "../../3_outputs/subcellular_localization_annotations.csv"
    
    print("=== Mislocalization Pattern Analysis ===\n")
    
    # Step 1: Load data and identify mislocalization
    mislocalization_data, loc_cols = analyze_mislocalization_patterns(csv_path)
    
    print("\n" + "="*50)
    
    # Step 2: Analyze pattern changes
    pattern_changes = calculate_pattern_changes(mislocalization_data, loc_cols)
    
    print("\n" + "="*50)
    
    # Step 3: Summarize unique patterns
    pattern_counts, semantic_changes = summarize_unique_patterns(pattern_changes)
    
    print(f"\nAnalysis complete! Found {len(pattern_changes)} variants with mislocalization changes.")