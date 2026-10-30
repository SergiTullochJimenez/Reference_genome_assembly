#!/bin/bash
# 01.1 Concatenate the raw read files of one species into one file per data type (and per read end).
# Expected layout: ${BASE_DIR}/${SPECIES}/data/{ONT,HIC,ILLUMINA}/ with the files as delivered by the sequencing centre.
set -euo pipefail

#### 0. Configuration ####
SPECIES="Species_name"
BASE_DIR="/path/to/project"          # contains one folder per species
THREADS=10

WD="${BASE_DIR}/${SPECIES}"
ONT="${WD}/data/ONT"
HIC="${WD}/data/HIC"
ILLUMINA="${WD}/data/ILLUMINA"

#### 1. ONT ####
zcat ${ONT}/*.gz | pigz -p ${THREADS} -c > ${ONT}/raw_ONT.fq.gz

#### 2. Hi-C and Illumina (R1 and R2 separately) ####
zcat ${HIC}/*_1.* | pigz -p ${THREADS} -c > ${HIC}/raw_HIC_1.fq.gz
zcat ${HIC}/*_2.* | pigz -p ${THREADS} -c > ${HIC}/raw_HIC_2.fq.gz
zcat ${ILLUMINA}/*_1.* | pigz -p ${THREADS} -c > ${ILLUMINA}/raw_ILLUMINA_1.fq.gz
zcat ${ILLUMINA}/*_2.* | pigz -p ${THREADS} -c > ${ILLUMINA}/raw_ILLUMINA_2.fq.gz
