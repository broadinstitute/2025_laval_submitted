import polars as pl
import re
from typing import Optional, Tuple


def process_accession_file(input_csv: str, output_parquet: str):
    """
    Process large accession CSV into optimized Parquet format with interval indexing
    [Uses Polars' streaming and chunk processing for memory efficiency]
    """
    # Read with schema optimization and streaming mode
    df = pl.scan_csv(
        input_csv,
        separator=',',
        has_header=False,
        new_columns=["uniprot_id","first_res_idx","last_res_idx","af_db_identifier","latest_version"]
    )

    # Add subfile index extraction (F1, F2 from AF-ID)
    processed = df.with_columns(
        pl.col("af_db_identifier")
        .str.extract(r"-F(\d+).*", 1)  # Capture all digits after F until end
        .alias("subfile_index")
    ).sort(["uniprot_id", "first_res_idx"])

    print(processed.collect())

    # Sort by uniprot_id and first_res_idx for interval search
    processed = processed.sort(["uniprot_id", "first_res_idx"])

    # Write optimized Parquet with partitioning
    processed.sink_parquet(
        output_parquet,
        row_group_size=10000,
        compression="zstd"
    )


process_accession_file("../../1_inputs/raw_data/alphafold/accession_ids.csv",
                       "../../3_outputs/processed_data/af_id_mapping/af_mapping.parquet")