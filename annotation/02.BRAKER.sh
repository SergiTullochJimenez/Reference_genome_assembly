#!/bin/bash
# 02 Evidence-guided gene prediction with BRAKER3 v3.0.8 (container): RNA-seq alignments of all tissues merged into one
# BAM + a protein database (OrthoDB v12 Vertebrata proteins plus UniProt venom toxin sequences).
set -euo pipefail

#### 0. Configuration ####
SPECIES="Species_name"
BASE_DIR="/path/to/project"
BRAKER_SIF="/path/to/braker.sif"
AUGUSTUS_CONFIG="/path/to/augustus_config/config"   # writable copy of the AUGUSTUS config directory
PROTEIN_DB="/path/to/protein_evidence.fa"          # protein evidence (FASTA)
DB_DIR=$(dirname ${PROTEIN_DB})
THREADS=40

WD="${BASE_DIR}/${SPECIES}"
outdir="${WD}/res/21.Braker"
REF_GENOME="${WD}/data/reference/${SPECIES}.masked.fa"
BAM_DIR="${WD}/res/22.2.Aligned_RNA_tissues"
MERGED_BAM="${BAM_DIR}/merged_rnaseq.sorted.bam"
mkdir -p ${outdir}

#### 1. Merge the per-tissue BAMs ####
if [ ! -f "${MERGED_BAM}" ]; then
    samtools merge -@ ${THREADS} ${MERGED_BAM} ${BAM_DIR}/*.sorted.bam
    samtools index -@ ${THREADS} ${MERGED_BAM}
fi

#### 2. BRAKER3 ####
ulimit -n 65535
rm -rf ${outdir}/GeneMark-ETP                       # clean a previous broken run
apptainer exec \
    --bind ${WD}:${WD} \
    --bind ${AUGUSTUS_CONFIG}:/opt/Augustus/config \
    --bind ${DB_DIR}:${DB_DIR} \
    --env AUGUSTUS_CONFIG_PATH=/opt/Augustus/config \
    ${BRAKER_SIF} \
    braker.pl \
    --workingdir=${outdir} \
    --genome=${REF_GENOME} \
    --bam=${MERGED_BAM} \
    --prot_seq=${PROTEIN_DB} \
    --threads=${THREADS} \
    --species=${SPECIES}
