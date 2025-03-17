#!/home/shenrunx/software/anaconda3/envs/ai4protein/bin/ python3

import os
import sys
import polars as pl
import pandas as pd
import numpy as np
import concurrent.futures
import argparse
from tqdm import tqdm
import time

# Import custom modules
import extract_plddt_scores as get_plddt
import calculate_rsa_scores as get_rsa


# Define all paths and settings directly in the script
INPUT_FILE = f"../../3_outputs/processed_data/dbnsfp/dbNSFP5.0a_variant.clin_var_re-annot_dedup_variants.parquet"
OUTPUT_FILE = "../../3_outputs/processed_data/dbnsfp/dbNSFP5.0a_variant.clin_var_dedup_plddt_rsa.parquet"
PLDDT_DIR = "../../3_outputs/processed_data/plddt_scores"
RSA_DIR = "../../3_outputs/processed_data/rsa_scores"
CHUNK_SIZE = 10000
NUM_WORKERS = max(256, os.cpu_count())  # Use a reasonable default for parallelism


# Define amino acid conversion dictionary
aa_single_to_three = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
    'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'
}


def get_plddt_score(pdb_id, aapos, aaref, parquet_files_plddt):
    """Get pLDDT score for a specific residue."""
    df = get_plddt.query_plddt(parquet_files_plddt, pdb_ids=[f"AF-{pdb_id}-F1-model_v4"])
    if df is None or df.empty:
        return np.nan

    try:
        if int(aapos) in df.index and df.loc[int(aapos), "residue_name"] == aaref:
            return df.loc[int(aapos), "plddt"]
    except:
        pass
    return np.nan


def get_rsa_score(pdb_id, aapos, aaref, parquet_files_rsa):
    """Get RSA score for a specific residue."""
    df = get_rsa.get_structure_rsa(parquet_files_rsa, pdb_id=f"AF-{pdb_id}-F1-model_v4")
    if df is None or df.empty:
        return np.nan

    try:
        if int(aapos) in df.index and df.loc[int(aapos), "residue_name"] == aaref:
            return df.loc[int(aapos), "rsa"]
    except:
        pass
    return np.nan


def process_chunk(chunk_df, plddt_files, rsa_files):
    """Process a chunk of the dataframe to add pLDDT and RSA scores."""
    return (
        chunk_df.with_columns(
            pl.struct(pl.col(["Uniprot_acc", "aapos", "aaref_3"]))
              .map_elements(
                  lambda x: get_plddt_score(
                      pdb_id=x["Uniprot_acc"].split('-')[0] if x["Uniprot_acc"] else "",
                      aapos=int(x["aapos"]) if x["aapos"] and x["aapos"].isdigit() else 0,
                      aaref=x["aaref_3"].upper() if x["aaref_3"] else "",
                      parquet_files_plddt=plddt_files
                  ),
                  return_dtype=pl.Float64
              )
              .alias("plddt"),
            pl.struct(pl.col(["Uniprot_acc", "aapos", "aaref_3"]))
              .map_elements(
                  lambda x: get_rsa_score(
                      pdb_id=x["Uniprot_acc"].split('-')[0] if x["Uniprot_acc"] else "",
                      aapos=int(x["aapos"]) if x["aapos"] and x["aapos"].isdigit() else 0,
                      aaref=x["aaref_3"].upper() if x["aaref_3"] else "",
                      parquet_files_rsa=rsa_files
                  ),
                  return_dtype=pl.Float64
              )
              .alias("rsa")
        )
    )


# This function is now outside the main function so it can be pickled
def process_chunk_range(args):
    """Process a specific range of the dataframe.

    Args:
        args: Tuple containing (chunk_range, df_selected, plddt_files, rsa_files)
              where chunk_range is (start, end)
    """
    chunk_range, df_path, plddt_files, rsa_files = args
    start, end = chunk_range

    # Load the chunk using the provided path
    df_selected = pl.scan_parquet(df_path)
    chunk = df_selected.slice(start, end - start).collect()

    # Process the chunk
    return process_chunk(chunk, plddt_files, rsa_files)


def main():
    start_time = time.time()
    print(f"Starting processing with {NUM_WORKERS} workers and chunk size of {CHUNK_SIZE}")

    # Get parquet files for pLDDT and RSA scores
    parquet_files_plddt = [f"{PLDDT_DIR}/{par_f}" for par_f in os.listdir(PLDDT_DIR) if par_f.endswith('.parquet')]
    parquet_files_rsa = [f"{RSA_DIR}/{par_f}" for par_f in os.listdir(RSA_DIR) if par_f.endswith('.parquet')]

    print(f"Found {len(parquet_files_plddt)} pLDDT parquet files and {len(parquet_files_rsa)} RSA parquet files")

    # Read the ClinVar data
    print("Reading ClinVar data...")
    df_clinvar = pl.scan_parquet(INPUT_FILE)

    # Step 3: Count total rows to process
    total_rows = df_clinvar.select(pl.len()).collect().item()
    print(f"Processing {total_rows} rows in total")

    chunks_to_process = []
    # Create chunks
    for i in range(0, total_rows, CHUNK_SIZE):
        end = min(i + CHUNK_SIZE, total_rows)
        # Include the path to the parquet file in the arguments
        chunks_to_process.append(((i, end), INPUT_FILE, parquet_files_plddt, parquet_files_rsa))

    print(f"Divided data into {len(chunks_to_process)} chunks")

    # Process chunks in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=NUM_WORKERS) as executor:
        results = list(tqdm(
            executor.map(process_chunk_range, chunks_to_process),
            total=len(chunks_to_process),
            desc="Processing chunks"
        ))

    # Concatenate all processed chunks
    print("Combining results...")
    final_df = pl.concat(results)

    # Write to output file
    print(f"Writing results to {OUTPUT_FILE}")
    final_df.write_parquet(OUTPUT_FILE)

    end_time = time.time()
    print(f"Processing completed in {end_time - start_time:.2f} seconds")


if __name__ == "__main__":
    # Simple argparse just to allow basic customization without editing script
    parser = argparse.ArgumentParser(description="Process ClinVar variants with pLDDT and RSA scores in parallel")
    parser.add_argument("--workers", type=int, default=NUM_WORKERS,
                        help=f"Number of parallel workers (default: {NUM_WORKERS})")
    parser.add_argument("--chunk-size", type=int, default=CHUNK_SIZE,
                        help=f"Chunk size for processing (default: {CHUNK_SIZE})")
    args = parser.parse_args()

    # Override defaults with command line args if provided
    NUM_WORKERS = args.workers
    CHUNK_SIZE = args.chunk_size

    main()