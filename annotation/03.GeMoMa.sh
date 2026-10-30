#!/bin/bash
# 03 Homology-based gene prediction with GeMoMa v1.9, projecting the annotations of three snake genomes
# (Crotalus adamanteus, Naja naja, Pantherophis guttatus) and using the RNA-seq alignments as intron evidence.
set -euo pipefail

#### 0. Configuration ####
SPECIES="Species_name"
BASE_DIR="/path/to/project"
REF_DIR="/path/to/reference_genomes"     # <ref>/genome.fna(.gz) + annotation.gff(.gz) for each reference
GEMOMA_ENV="GeMoMa"
GEMOMA_JAR="/path/to/GeMoMa-1.9.jar"
THREADS=40
MEM="240G"

WD="${BASE_DIR}/${SPECIES}"
OUT_DIR="${WD}/res/22.3.GeMoMa"
GENOME="${WD}/data/reference/${SPECIES}.masked.fa"
RNA_DIR="${WD}/res/22.2.Aligned_RNA_tissues"
mkdir -p ${OUT_DIR}

#### 1. One ERE.m= argument per tissue BAM ####
BAM_ARGS=$(ls ${RNA_DIR}/*.sorted.bam | sort | awk '{print "ERE.m="$1}' | tr '\n' ' ')

#### 2. GeMoMa pipeline ####
conda run -n ${GEMOMA_ENV} java -Xmx${MEM} -jar "${GEMOMA_JAR}" \
    CLI GeMoMaPipeline \
    threads=${THREADS} \
    outdir=${OUT_DIR} \
    GeMoMa.Score=ReAlign \
    AnnotationFinalizer.r=NO \
    o=false \
    t=${GENOME} \
    s=own i=crotalus_adamanteus g=${REF_DIR}/crotalus_adamanteus/genome.fna.gz a=${REF_DIR}/crotalus_adamanteus/annotation.gff.gz \
    s=own i=naja_naja g=${REF_DIR}/naja_naja/genome.fna a=${REF_DIR}/naja_naja/annotation.gff \
    s=own i=Pantherophis_guttatus g=${REF_DIR}/Pantherophis_guttatus/genome.fna a=${REF_DIR}/Pantherophis_guttatus/annotation.gff \
    r=MAPPED \
    ${BAM_ARGS}
