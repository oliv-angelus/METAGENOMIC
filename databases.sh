#!/bin/bash

# --- Definição das Cores ---
YELLOW='\033[1;33m'     # YELLOW
GREEN='\033[0;32m'      # GREEN
CYAN='\033[0;36m'       # CYAN
BOLD_CYAN='\033[1;36m'  # BOLD CYAN
RED='\033[1;31m'        # Bold Red
NC='\033[0m'            # NO COLOR (RESET)

# --- Standart Path for Databases ---
DB_PATH="${1:-$HOME/databases}"

# --- Critical Prerequisite Warning ---
echo -e "${RED}############################################################################${NC}"
echo -e "${RED}#                                                                          #${NC}"
echo -e "${RED}#   ATTENTION: PREREQUISITE REQUIRED                                       #${NC}"
echo -e "${RED}#                                                                          #${NC}"
echo -e "${RED}#   This script requires the 'checkm2' tool to be active in your           #${NC}"
echo -e "${RED}#   Conda/Mamba environment to download its database correctly.            #${NC}"
echo -e "${RED}#                                                                          #${NC}"
echo -e "${RED}#   => BEFORE PROCEEDING, PLEASE ACTIVATE THE CORRECT ENVIRONMENT.         #${NC}"
echo -e "${RED}#                                                                          #${NC}"
echo -e "${RED}#                                                                          #${NC}"
echo -e "${RED}############################################################################${NC}"
echo "" 

# --- INITIAL WARNING ---
echo -e "${CYAN}-------------------------------------------------------------------${NC}"
echo -e "${CYAN}NOTICE:${NC}"
echo -e "${CYAN}  This script will download and set up the databases in the default directory:${NC}"
echo -e "${BOLD_CYAN}    ${DB_PATH}${NC}"
echo ""
echo -e "${CYAN}  To change this location, please edit the path variable in the script.${NC}"
echo -e "${CYAN}-------------------------------------------------------------------${NC}"
echo ""


echo -e "${YELLOW}### 1. Criando diretórios para os bancos de dados... ###${NC}"

mkdir -p ${DB_PATH}/eggnog_database
mkdir -p ${DB_PATH}/kraken2_database

echo -e "${GREEN}Diretórios criados com sucesso.${NC}\n"

echo -e "${YELLOW}### 2. Baixando o banco de dados do CheckM2... ###${NC}"

checkm2 database --download --path ${DB_PATH}

echo -e "${GREEN}Download do CheckM2 concluído.${NC}\n"

echo -e "${YELLOW}### 3. Baixando os arquivos do banco de dados EggNOG... ###${NC}"

wget -P ${DB_PATH}/eggnog_database http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.db.gz
wget -P ${DB_PATH}/eggnog_database http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz
wget -P ${DB_PATH}/eggnog_database http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz
wget -P ${DB_PATH}/eggnog_database http://eggnog5.embl.de/download/emapperdb-5.0.2/mmseqs.tar.gz
wget -P ${DB_PATH}/eggnog_database http://eggnog5.embl.de/download/emapperdb-5.0.2/pfam.tar.gz
echo -e "${GREEN}Download do EggNOG concluído.${NC}\n"

echo -e "${YELLOW}### 4. Extraindo o banco de dados EggNOG... ###${NC}"
cd ${DB_PATH}/eggnog_database

tar -xzvf eggnog.taxa.tar.gz
tar -xzvf mmseqs.tar.gz
tar -xzvf pfam.tar.gz

gunzip eggnog.db.gz
gunzip eggnog_proteins.dmnd.gz

echo -e "${GREEN}Extração do EggNOG concluída.${NC}\n"

echo -e "${YELLOW}### 5. Limpando arquivos temporários (.tar.gz)... ###${NC}"

rm eggnog.taxa.tar.gz
rm mmseqs.tar.gz
rm pfam.tar.gz

echo -e "${GREEN}Arquivos temporários removidos.${NC}\n"

echo -e "${GREEN}############################################${NC}"
echo -e "${GREEN}### SETUP DOS BANCOS DE DADOS CONCLUÍDO! ###${NC}"
echo -e "${GREEN}############################################${NC}"

# --- WARNING K2 DATABASE ---
echo -e "${YELLOW}######################################################################${NC}"
echo -e "${YELLOW}#                                                                    #${NC}"
echo -e "${YELLOW}#   [ACTION REQUIRED] Kraken2 Database Setup                         #${NC}"
echo -e "${YELLOW}#                                                                    #${NC}"
echo -e "${YELLOW}#   The Kraken2 database was NOT downloaded automatically.           #${NC}"
echo -e "${YELLOW}#   You must choose and download a database that fits your needs.    #${NC}"
echo -e "${YELLOW}#                                                                    #${NC}"
echo -e "${YELLOW}#   Official Index URL: https://benlangmead.github.io/aws-indexes/k2   #${NC}"
echo -e "${YELLOW}#                                                                    #${NC}"
echo -e "${YELLOW}######################################################################${NC}"
