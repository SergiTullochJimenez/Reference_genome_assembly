#!/bin/bash
# 01.2 Read filtering.
#   ONT: Filtlong v0.2.1 (min length 1 kb, min mean quality 80).
#   Hi-C: fastp v0.20 (Q30, polyG/polyX trimming, adapter detection, base correction), then the first 5 bp of every
#         read are removed for the Arima mapping pipeline (filtered_Arima_HIC_*).
#   Illumina WGS is filtered with the same fastp settings in analyses/psmc/01.Mapping.sh.
set -euo pipefail

#### 0. Configuration ####
SPECIES="Species_name"
BASE_DIR="/path/to/project"
FILTLONG="/path/to/Filtlong/bin/filtlong"
THREADS=15

WD="${BASE_DIR}/${SPECIES}"
ONT="${WD}/data/ONT"
HIC="${WD}/data/HIC"
outdir="${WD}/res/01.1.Trimming"
mkdir -p ${outdir}

#### 1. ONT ####
${FILTLONG} --min_length 1000 --min_mean_q 80 ${ONT}/raw_ONT.fq.gz | pigz -p ${THREADS} -c > ${outdir}/filtered_ONT.fq.gz

#### 2. Hi-C ####
fastp -i ${HIC}/raw_HIC_1.fq.gz -I ${HIC}/raw_HIC_2.fq.gz -o ${outdir}/filtered_HIC_1.fq.gz -O ${outdir}/filtered_HIC_2.fq.gz \
    --thread ${THREADS} --qualified_quality_phred 30 --correction --trim_poly_g --trim_poly_x --detect_adapter_for_pe

#### 3. Arima: remove the first 5 bp of each Hi-C read ####
zcat ${outdir}/filtered_HIC_1.fq.gz | awk '{ if(NR%2==0) {print substr($1,6) } else {print} }' | gzip > ${outdir}/filtered_Arima_HIC_1.fq.gz
zcat ${outdir}/filtered_HIC_2.fq.gz | awk '{ if(NR%2==0) {print substr($1,6) } else {print} }' | gzip > ${outdir}/filtered_Arima_HIC_2.fq.gz
