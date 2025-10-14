#!/bin/bash

set -e

echo -e "\033[34m SCRIPT FOR PAIRED-END METAGENOMICS \033[0m"

# ===================== #
# PIPELINE FOR SHOTGUN  #
# ===================== #

# == # INPUTS
INPUT_DIR="$1"
OUTPUT_DIR="$2"
THREADS="$3"

# == # INPUT VERIFICATION
if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$THREADS" ]; then
  echo -e "\033[31m❌ Usage: $0 <input_directory> <output_directory> <threads>\033[0m"
  exit 1
fi

set +e

# == # DATABASE VARIABLES (adjust paths as necessary)
KRAKEN_DB=~/dataAngelo/databases/k2_db/PlusPFP # JUST A EXAMPLE, PLEASE DOWNLOAD YOUR OWN KRAKEN2 DATABASE

# == # CREATION OF OUTPUT DIRECTORIES
echo -e "\033[36m-> Creating directory structure in $OUTPUT_DIR\033[0m"
mkdir -p "$OUTPUT_DIR"/{00_QC,01_PROFILING}
mkdir -p "$OUTPUT_DIR/00_QC"/FASTP
mkdir -p "$OUTPUT_DIR/01_PROFILING"/TAXONOMY


# == # PROCESSING EACH SAMPLE
for SAMPLE_DIR in "$INPUT_DIR"/*; do
  SAMPLE=$(basename "$SAMPLE_DIR")
  echo -e "\033[33m🚀 Iniciando processamento da amostra: $SAMPLE\033[0m"

  fastp \
    -i "$SAMPLE_DIR"/*_R1*.fastq.gz \
    -I "$SAMPLE_DIR"/*_R2*.fastq.gz \
    -f 1 -F 1 \
    -l 100 \
    -o "$OUTPUT_DIR/00_QC/FASTP/${SAMPLE}_R1_trimmed.fastq.gz" \
    -O "$OUTPUT_DIR/00_QC/FASTP/${SAMPLE}_R2_trimmed.fastq.gz" \
    -h "$OUTPUT_DIR/00_QC/FASTP/${SAMPLE}_fastp.html" \
    -j "$OUTPUT_DIR/00_QC/FASTP/${SAMPLE}_fastp.json" \
    -w "$THREADS"

  # == # TAXONOMIC AND FUNCTIONAL PROFILING (CONTIGS) ---
  echo -e "\033[35m KRAKEN2 & BRACKEN \033[0m"

  READ1="$OUTPUT_DIR/00_QC/FASTP/${SAMPLE}_R1_trimmed.fastq.gz"
  READ2="$OUTPUT_DIR/00_QC/FASTP/${SAMPLE}_R2_trimmed.fastq.gz"
  TAX_DIR="$OUTPUT_DIR/01_PROFILING/TAXONOMY"

  kraken2 \
    --db "$KRAKEN_DB" \
    --threads "$THREADS" \
    --paired "$READ1" "$READ2" \
    --report "$TAX_DIR/${SAMPLE}_kraken_report.txt" \
    --output "$TAX_DIR/${SAMPLE}_kraken_output.txt"

  for level in D P C O F G S; do
    echo "  → Running for level $level ..."
    bracken \
      -d "$KRAKEN_DB" \
      -i "$TAX_DIR/${SAMPLE}_kraken_report.txt" \
      -o "$TAX_DIR/${SAMPLE}_bracken_${level}.txt" \
      -r 150 \
      -l "$level" \
      -t 1
  done
     
  echo -e "\033[32m✅ Finished processing sample: $SAMPLE\033[0m"
done
