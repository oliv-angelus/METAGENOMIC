#!/bin/bash

set -e

echo -e "\033[34m SCRIPT FOR TAXONOMIC PROFILING (KRAKEN2 + BRACKEN) \033[0m"

# ===================== #
#       INPUTS          #
# ===================== #

INPUT_DIR="$1"
OUTPUT_DIR="$2"
THREADS="$3"

if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$THREADS" ]; then
  echo -e "\033[31m❌ Usage: $0 <input_directory_with_trimmed_reads> <output_directory> <threads>\033[0m"
  exit 1
fi

# ===================== #
#     DATABASE PATHS    #
# ===================== #

KRAKEN_DB=~/dataAngelo/databases/k2_db/PlusPFP # 

# ===================== #
#    OUTPUT STRUCTURE   #
# ===================== #

echo -e "\033[36m-> Creating taxonomy output structure in $OUTPUT_DIR\033[0m"
mkdir -p "$OUTPUT_DIR/TAXONOMY"
TAX_DIR="$OUTPUT_DIR/TAXONOMY"

# ===================== #
#   PROCESS EACH SAMPLE #
# ===================== #

for R1 in "$INPUT_DIR/"*_R1_trimmed.fastq.gz; do
  SAMPLE=$(basename "$R1" "_R1_trimmed.fastq.gz")
  R2="${INPUT_DIR}/${SAMPLE}_R2_trimmed.fastq.gz"

  if [ ! -f "$R2" ]; then
    echo -e "\033[33m⚠️ Paired file not found for $R1, skipping...\033[0m"
    continue
  fi

  echo -e "\033[33m🚀 Processing sample: $SAMPLE\033[0m"

  # --- [KRAKEN2] ---
  echo -e "\033[35m KRAKEN2 \033[0m"
  kraken2 \
    --db "$KRAKEN_DB" \
    --threads "$THREADS" \
    --paired "$R1" "$R2" \
    --use-names \
    --report "$TAX_DIR/${SAMPLE}_kraken_report.txt" \
    --output "$TAX_DIR/${SAMPLE}_kraken_output.txt"

  # --- [BRACKEN] ---
  echo -e "\033[35m BRACKEN \033[0m"
  for level in D P C O F G S; do
    echo "  → Bracken level: $level"
    bracken \
      -d "$KRAKEN_DB" \
      -i "$TAX_DIR/${SAMPLE}_kraken_report.txt" \
      -o "$TAX_DIR/${SAMPLE}_bracken_${level}.txt" \
      -r 150 \
      -l "$level" \
      -t 1
  done

  echo -e "\033[32m✅ Finished sample: $SAMPLE\033[0m"
done

echo -e "\033[32m✅ Taxonomic profiling completed! Results in: $TAX_DIR \033[0m"
