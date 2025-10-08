#!/bin/bash

set -e

echo -e "\033[34m SCRIPT FOR SINGLE-END NANOPORE METAGENOMICS \033[0m"

# ===================== #
# PIPELINE FOR SHOTGUN  #
# ===================== #

# == # INPUT VARIABLES
INPUT_DIR="$1"
OUTPUT_DIR="$2"
THREADS="$3"

# == # INPUT VERIFICATION
if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$THREADS" ]; then
  echo -e "\033[31m❌ Usage: $0 <input_directory> <output_directory> <threads>\033[0m"
  exit 1
fi

# == # DATABASE VARIABLES (adjust paths as necessary)
KRAKEN_DB=~/databases/kraken_database # EXAMPLE: Please download your own Kraken2 database
EGGNOG_DB=~/databases/eggnog_database # Path to the EggNOG database
CHECKM2_DB=~/databases/CheckM2_database # Path to the CheckM2 database

# == # CREATION OF OUTPUT DIRECTORIES
echo -e "\033[36m-> Creating directory structure in $OUTPUT_DIR\033[0m"
mkdir -p "$OUTPUT_DIR"/{00_QC,01_ASSEMBLY,02_MAPPING_AND_COVERAGE,03_PROFILING,04_BINNING,05_MAG_QC,06_MAG_ANNOTATION,REPORTS}
mkdir -p "$OUTPUT_DIR/00_QC"/{NANOPLOT,FILTLONG}
mkdir -p "$OUTPUT_DIR/01_ASSEMBLY"/{FLYE,QUAST}
mkdir -p "$OUTPUT_DIR/02_MAPPING_AND_COVERAGE"/{BAM,COVERM}
mkdir -p "$OUTPUT_DIR/03_PROFILING"/{TAXONOMY,CONTIG_ANNOTATION}
mkdir -p "$OUTPUT_DIR/03_PROFILING/CONTIG_ANNOTATION"/{PRODIGAL,EGGNOGMAPPER}
mkdir -p "$OUTPUT_DIR/03_PROFILING/CONTIG_ANNOTATION/PRODIGAL"/{PROTEINS,GENES_NT,GFF}
mkdir -p "$OUTPUT_DIR/04_BINNING"/{METABAT2}
mkdir -p "$OUTPUT_DIR/05_MAG_QC"/{CHECKM2}
mkdir -p "$OUTPUT_DIR/06_MAG_ANNOTATION"/{PRODIGAL,EGGNOGMAPPER}
mkdir -p "$OUTPUT_DIR/06_MAG_ANNOTATION/PRODIGAL"/{PROTEINS,GENES_NT,GFF}

# == # PROCESSING EACH SAMPLE
for RAW_FASTQ in "$INPUT_DIR"/*.fastq.gz; do
  SAMPLE=$(basename "$RAW_FASTQ" .fastq.gz)
  echo -e "\033[33m🚀 Starting processing for sample: $SAMPLE\033[0m"

  # --- [00] QUALITY CONTROL & FILTERING ---
  echo -e "\033[35m NANOPLOT \033[0m"
  NanoPlot \
    --fastq "$RAW_FASTQ" \
    -o "$OUTPUT_DIR/00_QC/NANOPLOT/$SAMPLE" \
    -t "$THREADS"

  echo -e "\033[35m FILTLONG \033[0m"
  FILTERED_FASTQ="$OUTPUT_DIR/00_QC/FILTLONG/${SAMPLE}_filtered.fastq.gz"
  filtlong --min_length 1000 "$RAW_FASTQ" | gzip > "$FILTERED_FASTQ"

  # --- [01] METAGENOME ASSEMBLY ---
  echo -e "\033[35m FLYE \033[0m"
  flye \
    --nano-raw "$FILTERED_FASTQ" \
    --out-dir "$OUTPUT_DIR/01_ASSEMBLY/FLYE/$SAMPLE" \
    --threads "$THREADS" \
    --meta
  
  echo -e "\033[35m QUAST \033[0m"
  CONTIGS_FA="$OUTPUT_DIR/01_ASSEMBLY/FLYE/$SAMPLE/assembly.fasta"
  quast \
    "$CONTIGS_FA" \
    -o "$OUTPUT_DIR/01_ASSEMBLY/QUAST/$SAMPLE" \
    -t "$THREADS"

  # --- [02] READ MAPPING AND COVERAGE ---
  BAM_DIR="$OUTPUT_DIR/02_MAPPING_AND_COVERAGE/BAM"
  SORTED_BAM="$BAM_DIR/${SAMPLE}.sorted.bam"
  DEPTH_FILE="$BAM_DIR/${SAMPLE}.depth.txt"
  
  echo -e "\033[35m MINIMAP2 \033[0m"
  minimap2 \
    -ax map-ont "$CONTIGS_FA" "$FILTERED_FASTQ" \
    -t "$THREADS" | samtools view -bS -@ "$THREADS" > "$BAM_DIR/${SAMPLE}.bam"
  
  echo -e "\033[35m SAMTOOLS \033[0m"
  samtools sort -@ "$THREADS" "$BAM_DIR/${SAMPLE}.bam" -o "$SORTED_BAM"
  samtools index "$SORTED_BAM"
  
  jgi_summarize_bam_contig_depths --outputDepth "$DEPTH_FILE" "$SORTED_BAM"

  echo -e "\033[35m COVERM \033[0m"
  coverm contig -b "$SORTED_BAM" > "$OUTPUT_DIR/02_MAPPING_AND_COVERAGE/COVERM/${SAMPLE}_coverage.tsv"

  # --- [03] TAXONOMIC AND FUNCTIONAL PROFILING (CONTIGS) ---
  echo -e "\033[35m KRAKEN2 & BRACKEN \033[0m"
  kraken2 --db "$KRAKEN_DB" --threads "$THREADS" --report "$OUTPUT_DIR/03_PROFILING/TAXONOMY/${SAMPLE}_kraken_report.txt" --output "$OUTPUT_DIR/03_PROFILING/TAXONOMY/${SAMPLE}_kraken_output.txt" "$CONTIGS_FA"
  bracken -d "$KRAKEN_DB" -i "$OUTPUT_DIR/03_PROFILING/TAXONOMY/${SAMPLE}_kraken_report.txt" -o "$OUTPUT_DIR/03_PROFILING/TAXONOMY/${SAMPLE}_bracken_species.txt" -r 150 -l G -t 1

  PRODIGAL_CONTIGS_DIR="$OUTPUT_DIR/03_PROFILING/CONTIG_ANNOTATION/PRODIGAL"
  CONTIGS_PROTEINS="$PRODIGAL_CONTIGS_DIR/PROTEINS/${SAMPLE}_contigs.faa"

  echo -e "\033[35m PRODIGAL - CONTIGS \033[0m"
  prodigal -i "$CONTIGS_FA" -p meta \
    -a "$CONTIGS_PROTEINS" \
    -d "$PRODIGAL_CONTIGS_DIR/GENES_NT/${SAMPLE}_contigs.fna" \
    -o "$PRODIGAL_CONTIGS_DIR/GFF/${SAMPLE}_contigs.gff" -f gff
  
  echo -e "\033[35m EGGNOG - CONTIGS \033[0m"
  emapper.py -i "$CONTIGS_PROTEINS" -m diamond --data_dir "$EGGNOG_DB" --cpu "$THREADS" -o "$OUTPUT_DIR/03_PROFILING/CONTIG_ANNOTATION/EGGNOGMAPPER/${SAMPLE}_contigs"

  # --- [04] METAGENOME BINNING ---
  echo -e "\033[35m METABAT2 \033[0m"
  METABAT2_OUT_DIR="$OUTPUT_DIR/04_BINNING/METABAT2/$SAMPLE"
  metabat2 -i "$CONTIGS_FA" -a "$DEPTH_FILE" -o "$METABAT2_OUT_DIR/bin" -t "$THREADS"

  # --- [05] MAG QUALITY CONTROL ---
  echo -e "\033[35m CHECKM2 \033[0m"
  checkm2 predict --threads "$THREADS" --input "$METABAT2_OUT_DIR" --output-directory "$OUTPUT_DIR/05_MAG_QC/CHECKM2/$SAMPLE" -x fa --force --database_path "$CHECKM2_DB"

  # --- [06] INDIVIDUAL ANNOTATION OF MAGs ---
  echo -e "\033[35m PRODIGAL & EGGNOG - BINS \033[0m"
  for MAG_FA in "$METABAT2_OUT_DIR"/bin.*.fa; do
    # Check if any bin files were created before proceeding
    [ -e "$MAG_FA" ] || continue
    MAG_NAME=$(basename "$MAG_FA" .fa)
    echo "  -> Annotating MAG: $MAG_NAME"

    # Prodigal for MAGs
    PRODIGAL_MAGS_DIR="$OUTPUT_DIR/06_MAG_ANNOTATION/PRODIGAL"
    MAG_PROTEINS="$PRODIGAL_MAGS_DIR/PROTEINS/${MAG_NAME}.faa"
    prodigal -i "$MAG_FA" -p meta \
        -a "$MAG_PROTEINS" \
        -d "$PRODIGAL_MAGS_DIR/GENES_NT/${MAG_NAME}.fna" \
        -o "$PRODIGAL_MAGS_DIR/GFF/${MAG_NAME}.gff" -f gff

    # EggNOG-mapper for MAGs
    emapper.py -i "$MAG_PROTEINS" -m diamond --data_dir "$EGGNOG_DB" --cpu "$THREADS" -o "$OUTPUT_DIR/06_MAG_ANNOTATION/EGGNOGMAPPER/${MAG_NAME}"
  done
     
  echo -e "\033[32m✅ Finished processing sample: $SAMPLE\033[0m"
done

# --- [07] FINAL REPORT ---
echo -e "\033[35m[07] STEP: MULTIQC\033[0m"
multiqc "$OUTPUT_DIR/00_QC" "$OUTPUT_DIR/01_ASSEMBLY" "$OUTPUT_DIR/05_MAG_QC/CHECKM2" -o "$OUTPUT_DIR/REPORTS"

echo -e "\033[32m✅ Pipeline finished successfully! Outputs are in: $OUTPUT_DIR \033[0m"