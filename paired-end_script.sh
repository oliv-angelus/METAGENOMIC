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
EGGNOG_DB=~/dataAngelo/databases/eggnog # PATH TO EGGNOG DATABASE
CHECKM2_DB=~/dataAngelo/databases/CheckM2_database/uniref100.KO.1.dmnd # PATH TO CHECKM2 DATABASE

# == # CREATION OF OUTPUT DIRECTORIES
echo -e "\033[36m-> Creating directory structure in $OUTPUT_DIR\033[0m"
mkdir -p "$OUTPUT_DIR"/{00_QC,01_ASSEMBLY,02_MAPPING_AND_COVERAGE,03_PROFILING,04_BINNING,05_MAG_QC,06_MAG_ANNOTATION,REPORTS}
mkdir -p "$OUTPUT_DIR/00_QC"/{FASTQC,FASTP}
mkdir -p "$OUTPUT_DIR/01_ASSEMBLY"/{MEGAHIT,QUAST}
mkdir -p "$OUTPUT_DIR/02_MAPPING_AND_COVERAGE"/{BAM,COVERM}
mkdir -p "$OUTPUT_DIR/03_PROFILING"/{TAXONOMY,CONTIG_ANNOTATION}
mkdir -p "$OUTPUT_DIR/03_PROFILING/CONTIG_ANNOTATION"/{PRODIGAL,EGGNOGMAPPER}
mkdir -p "$OUTPUT_DIR/03_PROFILING/CONTIG_ANNOTATION/PRODIGAL"/{PROTEINS,GENES_NT,GFF}
mkdir -p "$OUTPUT_DIR/04_BINNING"/METABAT2
mkdir -p "$OUTPUT_DIR/05_MAG_QC"/CHECKM2
mkdir -p "$OUTPUT_DIR/06_MAG_ANNOTATION"/{PRODIGAL,EGGNOGMAPPER}
mkdir -p "$OUTPUT_DIR/06_MAG_ANNOTATION/PRODIGAL"/{PROTEINS,GENES_NT,GFF}

# == # PROCESSING EACH SAMPLE
for SAMPLE_DIR in "$INPUT_DIR"/*; do
  SAMPLE=$(basename "$SAMPLE_DIR")
  echo -e "\033[33m🚀 Iniciando processamento da amostra: $SAMPLE\033[0m"

  # --- [00] QUALITY CONTROL ---
  echo -e "\033[35m FASTQC \033[0m"
  fastqc \
  "$SAMPLE_DIR"/*.fastq.gz \
  -t "$THREADS" \
  -o "$OUTPUT_DIR/00_QC/FASTQC"

  echo -e "\033[35m FASTP \033[0m"

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

  # --- [01] METAGENOME ASSEMBLY ---
  echo -e "\033[35m MEGAHIT \033[0m"
  megahit \
    -1 "$OUTPUT_DIR/00_QC/FASTP/${SAMPLE}_R1_trimmed.fastq.gz" \
    -2 "$OUTPUT_DIR/00_QC/FASTP/${SAMPLE}_R2_trimmed.fastq.gz" \
    -o "$OUTPUT_DIR/01_ASSEMBLY/MEGAHIT/$SAMPLE" \
    --min-contig-len 1500 \
    -t "$THREADS"
  
  echo -e "\033[35m QUAST \033[0m"

  CONTIGS_FA="$OUTPUT_DIR/01_ASSEMBLY/MEGAHIT/$SAMPLE/final.contigs.fa"
  quast \
  "$CONTIGS_FA" \
  -o "$OUTPUT_DIR/01_ASSEMBLY/QUAST/$SAMPLE" \
  -t "$THREADS"

  # --- [02] READS MAPPING AND COVERAGE ---
  BAM_DIR="$OUTPUT_DIR/02_MAPPING_AND_COVERAGE/BAM"
  SORTED_BAM="$BAM_DIR/${SAMPLE}.sorted.bam"
  DEPTH_FILE="$BAM_DIR/${SAMPLE}.depth.txt"
  
  echo -e "\033[35m MINIMAP \033[0m"

  minimap2 \
  -ax sr "$CONTIGS_FA" "$OUTPUT_DIR/00_QC/FASTP/${SAMPLE}_R1_trimmed.fastq.gz" "$OUTPUT_DIR/00_QC/FASTP/${SAMPLE}_R2_trimmed.fastq.gz" \
  -t "$THREADS" | samtools view -bS -@ "$THREADS" > "$BAM_DIR/${SAMPLE}.bam"
  
  echo -e "\033[35m SAMTOOLS \033[0m"

  samtools sort -@ "$THREADS" "$BAM_DIR/${SAMPLE}.bam" -o "$SORTED_BAM"
  samtools index "$SORTED_BAM"
  
  jgi_summarize_bam_contig_depths --outputDepth "$DEPTH_FILE" "$SORTED_BAM"

  echo -e "\033[35m COVERM \033[0m"

  coverm contig -b "$SORTED_BAM" > "$OUTPUT_DIR/02_MAPPING_AND_COVERAGE/COVERM/${SAMPLE}_coverage.tsv"

  # --- [03] TAXONOMIC AND FUNCTIONAL PROFILING (CONTIGS) ---
  echo -e "\033[35m KRAKEN2 & BRACKEN \033[0m"

  READ1="$OUTPUT_DIR/00_QC/FASTP/${SAMPLE}_R1_trimmed.fastq.gz"
  READ2="$OUTPUT_DIR/00_QC/FASTP/${SAMPLE}_R2_trimmed.fastq.gz"
  TAX_DIR="$OUTPUT_DIR/03_PROFILING/TAXONOMY"

  kraken2 \
    --db "$KRAKEN_DB" \
    --threads "$THREADS" \
    --paired "$READ1" "$READ2" \
    --report "$TAX_DIR/${SAMPLE}_kraken_report.txt" \
    --output "$TAX_DIR/${SAMPLE}_kraken_output.txt"

  for level in D P C O F G S; do
    echo "  → Rodando Bracken para nível $level ..."
    bracken \
      -d "$KRAKEN_DB" \
      -i "$TAX_DIR/${SAMPLE}_kraken_report.txt" \
      -o "$TAX_DIR/${SAMPLE}_bracken_${level}.txt" \
      -r 150 \
      -l "$level" \
      -t 1
  done

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
  metabat2 -i "$CONTIGS_FA" -a "$DEPTH_FILE" -o "$METABAT2_OUT_DIR/${SAMPLE}_bin" -t "$THREADS"

  # --- [05] QUALITY CONTROL ---
  echo -e "\033[35m CHECKM2 \033[0m"
    if compgen -G "$METABAT2_OUT_DIR/*.fa" > /dev/null; then
      checkm2 predict --threads "$THREADS" --input "$METABAT2_OUT_DIR" --output-directory "$OUTPUT_DIR/05_MAG_QC/CHECKM2/$SAMPLE" -x fa --force --database_path "$CHECKM2_DB"
    else
      echo -e "\033[33m⚠️ Nenhum bin encontrado para $SAMPLE, pulando CHECKM2.\033[0m"
    fi

  # --- [06] INDIVIDUAL ANNOTATION OF MAGs ---
  echo -e "\033[35m PRODIGAL & EGGNOG - BINS \033[0m"
  for MAG_FA in "$METABAT2_OUT_DIR"/*.fa; do
    # TEST IF FILE EXISTS
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

# --- [07] RELATÓRIO FINAL ---
echo -e "\033[35m[07] ETAPA: MULTIQC\033[0m"
multiqc "$OUTPUT_DIR/00_QC" "$OUTPUT_DIR/01_ASSEMBLY" "$OUTPUT_DIR/05_MAG_QC/CHECKM2" -o "$OUTPUT_DIR/REPORTS"

echo -e "\033[32m✅ Pipeline finished! Outputs in: $OUTPUT_DIR \033[0m"
