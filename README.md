# PEMS-Pipeline 🧬

### **P**aired-**E**nd **M**etagenomics **S**hotgun **Pipeline**

**PEMS-Pipeline** is an automated Bash workflow designed for the comprehensive analysis of paired-end shotgun metagenomic data. It streamlines the process from raw read quality control to the recovery and quality assessment of Metagenome-Assembled Genomes (MAGs).

## 📋 Workflow Overview

The pipeline executes the following sequential steps:

1.  **Quality Control:** Raw read processing using **FastQC** and **Fastp**.
2.  **Assembly:** De novo assembly using **MEGAHIT** and assessment with **QUAST**.
3.  **Mapping & Coverage:** Read mapping with **Minimap2**, sorting with **Samtools**, and coverage calculation with **CoverM**.
4.  **Profiling:** Taxonomic classification using **Kraken2/Bracken** and functional annotation with **Prodigal/EggNOG-mapper**.
5.  **Binning:** Genome binning using **MetaBAT2**.
6.  **MAG Quality:** Quality assessment of bins using **CheckM2**.

## 🛠️ Requirements & Installation

The only requirement to run PEMS-Pipeline is **Mamba** (Recommended for faster dependency resolution). All necessary bioinformatics tools are automatically installed via the provided `PEMS.yaml` file.

### 1. Installation

Clone the repository and create the environment using **mamba** with the `PEMS.yaml` file:

```bash
# Clone the repository
git clone [https://github.com/oliv-angelus/PEMS-Pipeline.git](https://github.com/oliv-angelus/PEMS-Pipeline.git)
cd PEMS-Pipeline

# Create the environment (using Mamba)
mamba env create -f PEMS.yaml

# Activate the environment
mamba activate PEMS
```

### 2. Databases setup

This project includes an automated script to download and configure the required databases (CheckM2 and EggNOG) and set up the directory structure. 

### ⚠️ Prerequisites
The script uses the `checkm2` tool to download its own specific database. Therefore, **you must activate the PEMS Mamba environment** where CheckM2 is installed before running the script.

### 🚀 Usage

*Option A: instaling at the default location (~/databases)*

```bash
chmod =x databases.sh
```

*Option B: instaling to a custom location provided as an argument*

```bash
./databases.sh /path/to/custom/database_dir
```

### 📝 Important Notes

*Disk Space:* The EggNOG database download is large. Ensure you have sufficient disk space before proceeding.

*Kraken2:* The script creates the directory structure for Kraken2 but DOES NOT download the database automatically (due to the large file size and variety of index options). Please follow the instructions printed at the end of the script execution to download the official index manually.


## 📂 Output Structure

After execution, **PEMS-Pipeline** organizes results into the following directory tree:

```text
results/
├── 00_QC/
│   ├── FASTQC/              # Quality reports before/after trimming
│   └── FASTP/               # Trimmed reads (*.fastq.gz) and HTML/JSON reports
├── 01_ASSEMBLY/
│   ├── MEGAHIT/             # Final contigs (final.contigs.fa)
│   └── QUAST/               # Assembly quality metrics
├── 02_MAPPING_AND_COVERAGE/
│   ├── BAM/                 # Sorted BAM files and depth calculations
│   └── COVERM/              # Coverage tables per sample
├── 03_PROFILING/
│   ├── TAXONOMY/            # Kraken2 reports and Bracken abundance estimations
│   └── CONTIG_ANNOTATION/   # Prodigal gene predictions (GFF/FAA/FNA) and EggNOG annotations
├── 04_BINNING/
│   └── METABAT2/            # Recovered bins (MAGs)
└── 05_MAG_QC/
    └── CHECKM2/             # Completeness and contamination reports for MAGs
```

### 📊 Analysis & Visualization

For downstream analysis, statistical testing, and generating plots from the data processed by this pipeline, please use the official companion application:

**OmniMeta**
🔗 https://github.com/oliv-angelus/OmniMeta.git

## 👤 Author

**Angelo Felipe Barbosa de Oliveira**

[![ORCID](https://img.shields.io/badge/ORCID-0000--0003--0831--447X-green.svg)](https://orcid.org/0000-0003-0831-447X)
[![Lattes](https://img.shields.io/badge/Lattes-Curriculum-blue.svg)](http://lattes.cnpq.br/5450775990055106)

## 📄 Citation

If you use **PEMS-Pipeline** in your research, please cite:

> **de Oliveira, A. F. B.** (2025). *PEMS-Pipeline: A Paired-End Metagenomics Shotgun Pipeline*. Available at: https://github.com/YOUR_USERNAME/PEMS-Pipeline

Alternatively, you can cite this repository using the `CITATION.cff` file provided.

## ⚖️ License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
