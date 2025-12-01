# METAGENOMIC

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![Language](https://img.shields.io/github/languages/top/YOUR-USERNAME/METAGENOMIC)
![Last Commit](https://img.shields.io/github/last-commit/YOUR-USERNAME/METAGENOMIC)

A repository containing a collection of scripts to automate and streamline shotgun metagenomics analysis workflows.

## 📖 About This Project

This project provides a suite of scripts designed to handle various stages of a metagenomic analysis pipeline, with tailored workflows for both Illumina (short-read) and Oxford Nanopore (long-read)(soon) sequencing data. The tools are built to manage the unique characteristics of each data type and can be adapted for hybrid assembly approaches. The pipeline covers all critical steps, from the initial quality control of raw reads to the final annotation and classification of high-quality Metagenome-Assembled Genomes (MAGs). The goal is to create a modular, reproducible, and easy-to-use toolkit for researchers.

## 🛠️ Core Technologies

The scripts and workflows in this repository are built using the following technologies:

* **Languages:**
    * `Shell (Bash)`
    * `R`

* **Environment Management:**
    * `Conda` and `Mamba` are used for robust and reproducible dependency management, ensuring that all tools run in a consistent environment.

## 🔬 Key Bioinformatic Tools

This repository leverages a range of standard and state-of-the-art bioinformatic tools, which are orchestrated by the scripts. The core tools utilized include:

#### **Quality Control**
* **FastQC**: Raw read quality assessment for illumina reads.
* **NanoPlot**: Raw read quality assessment for nanopore reads.
* **Fastp**: For trimming illumina reads.
* **Filtlong**: For trimming nanopore reads.
* **MultiQC**: Aggregates reports from multiple tools into a single HTML file.

#### **Assembly & Assessment**
* **MEGAHIT**: A fast and memory-efficient metagenome assembler for illumina reads.
* **Flye**: An assembler for long and noisy reads (useful for hybrid assembly), here we use for nanopore reads using -meta parameter.
* **Quast**: Evaluates the quality of genome assemblies.

#### **Metagenomic Binning**
* **MetaBAT2**: A tool for binning assembled contigs into MAGs based on sequence composition and coverage (i preffer this one).

#### **MAG Quality Assessment**
* **CheckM2**: The current standard for assessing the completeness and contamination of MAGs using a machine learning model.

#### **Gene Prediction & Functional Annotation**
* **Prodigal**: A fast and accurate gene prediction tool for microbial genomes.
* **eggNOG-mapper**: A tool for fast functional annotation of proteins.

#### **Taxonomic Classification**
* **Kraken2**: A k-mer-based tool for assigning taxonomic labels to sequences.
* **Bracken**: Estimates species/genus abundance from Kraken2 reports.

#### **Read Mapping & Coverage**
* **Minimap2**: A fast sequence mapper for DNA and mRNA sequences.
* **SAMtools**: A suite of tools for interacting with high-throughput sequencing data.
* **CoverM**: A tool for calculating coverage of contigs or genomes.

#### **Phylogenetic Analysis**
* **MAFFT**: A multiple sequence alignment program.
* **FastTree**: Infers approximately-maximum-likelihood phylogenetic trees from alignments.

## 🚀 Getting Started

This guide will walk you through setting up the computational environment and downloading the necessary databases to run the analyses in this repository.

1. Prerequisites
Before you begin, you must have Mamba (or Conda) installed on your system.

2. Clone the Repository
First, clone this repository to your local machine and navigate into the created directory.

```bash
git clone https://github.com/oliv-angelus/metagenomic
cd metagenomic
```
3. Create the Mamba Environment
The metagenomics.yaml file contains a list of all the required tools. Use the command below to create an isolated environment with everything installed:

```bash
mamba env create -f metagenomics.yaml
```
This may take a few minutes as Mamba will download and install dozens of bioinformatics packages. Once completed, activate the environment to start using it:

```bash
mamba activate metagenomics
```
4. Prepare the Databases
The databases.sh script automates the download and preparation of the CheckM2 and EggNOG databases.

Attention! Before running the script, make sure the metagenomics environment (created in the previous step) is active. The script requires the checkm2 command to function correctly.

a. Make the script executable (only needs to be done once):

```bash
chmod +x databases.sh
```

b. Run the script:

You can run the script in one of two ways:

* Option 1 (Default): Download the databases to the default directory (~/databases):

```bash
./databases.sh
```

* Option 2 (Custom): Specify a different directory (e.g., an external HDD):

```Bash
./databases.sh /path/to/your/directory
```

This process can be very time-consuming and may use tens of gigabytes of disk space and bandwidth, depending on your connection.

5. Download the Kraken2 Database (Manual Step)
As noted by the script, the Kraken2 database must be downloaded manually. This allows you to choose the version that is most appropriate for your analysis (e.g., Standard, PlusPF, etc.).

      1. Access the official Kraken2 database index: https://benlangmead.github.io/aws-indexes/k2.

      2. Choose and download the database of your preference.

      3. Unpack it and place the contents inside the kraken2_database directory that was created by the databases.sh script (${DB_PATH}/kraken2_database).

## 📄 License

This project is distributed under the MIT License. See the `LICENSE` file for more information.
