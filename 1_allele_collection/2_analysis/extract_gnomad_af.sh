#!/bin/bash

# Set strict error handling
set -euo pipefail

# Prompt user for name variable
read -p "Enter your dataset name (used for directory and file naming): " name

# Define directories and files
VCF_DIR="vcf"
VARIANT_LIST="${name}_variants_list.txt"
OUTPUT_DIR="${name}_allele_frequencies"

# Create directories only if they don't exist
[[ ! -d "$VCF_DIR" ]] && mkdir -p "$VCF_DIR"
[[ ! -d "$OUTPUT_DIR" ]] && mkdir -p "$OUTPUT_DIR"

# Define list of chromosomes
CHROMOSOMES=({1..22} X Y)

# Check for missing VCF files
MISSING_VCF=()
for CHR in "${CHROMOSOMES[@]}"; do
    VCF_FILE="$VCF_DIR/gnomad.exomes.v4.1.sites.chr${CHR}.vcf.bgz"
    if [[ ! -f "$VCF_FILE" ]]; then
        MISSING_VCF+=("$CHR")
    fi
done

# Download missing VCF files
if [[ ${#MISSING_VCF[@]} -eq 0 ]]; then
    echo "All VCF files are present. Skipping download."
else
    echo "The following chromosome VCF files are missing and will be downloaded: ${MISSING_VCF[*]}"
    for CHR in "${MISSING_VCF[@]}"; do
        echo "⬇️ Downloading chr${CHR}..."
        aws s3 cp "s3://gnomad-public-us-east-1/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr${CHR}.vcf.bgz" "$VCF_DIR/"
    done
    echo "✅ Download of missing VCF files completed."
fi

# Process each chromosome
for CHR in "${CHROMOSOMES[@]}"; do
    VCF_FILE="$VCF_DIR/gnomad.exomes.v4.1.sites.chr${CHR}.vcf.bgz"
    CSI_FILE="${VCF_FILE}.csi"
    OUTPUT_FILE="$OUTPUT_DIR/chr${CHR}_allele_frequencies.txt"

    # Check if VCF file exists
    if [[ -f "$VCF_FILE" ]]; then
        echo "🧬 Processing chromosome $CHR..."

        # Index VCF file only if not already indexed
        if [[ ! -f "$CSI_FILE" ]]; then
            echo "🔍 Indexing $VCF_FILE..."
            bcftools index --threads 4 "$VCF_FILE"
        else
            echo "Index for $VCF_FILE already exists. Skipping indexing."
        fi

        # Extract allele frequencies for listed variants
        echo "📊 Extracting allele frequencies for chromosome $CHR..."
        bcftools view -R "$VARIANT_LIST" "$VCF_FILE" | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' > "$OUTPUT_FILE"

        echo "✅ Saved allele frequencies to $OUTPUT_FILE."
    else
        echo "⚠️ VCF file for chromosome $CHR not found, skipping..."
    fi
done

echo "🎉 All tasks completed."
