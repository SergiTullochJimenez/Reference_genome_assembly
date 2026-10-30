#!/bin/bash
# 01 RNA-seq evidence: per-tissue reads aligned with HISAT2 to the soft-masked genome (reads were not trimmed;
# HISAT2 soft-clips read ends), sorted and indexed with SAMtools.
# Expected layout: ${WD}/data/RNA/<tissue>/<lane files>_1.fq.gz / _2.fq.gz, one folder per tissue.
set -euo pipefail

#### 0. Configuration ####
SPECIES="Species_name"
BASE_DIR="/path/to/project"
HISAT2_ENV="HISAT2"                  # conda env with hisat2
THREADS=16
PARALLEL_JOBS=3                      # tissues aligned at the same time

WD="${BASE_DIR}/${SPECIES}"
GENOME="${WD}/data/reference/${SPECIES}.masked.fa"     # soft-masked genome (assembly/06.Masking.sh)
INDIR="${WD}/data/RNA"
REF_DIR="${WD}/res/20.0.RepeatMasking"
OUTDIR_CONCAT="${WD}/res/22.1.Concatenated_RNA_tissue"
OUTDIR_ALIGN="${WD}/res/22.2.Aligned_RNA_tissues"
mkdir -p ${REF_DIR} ${OUTDIR_CONCAT} ${OUTDIR_ALIGN}

#### 1. HISAT2 index ####
[ -f "${REF_DIR}/genome-idx.1.ht2" ] || conda run -n ${HISAT2_ENV} hisat2-build "${GENOME}" "${REF_DIR}/genome-idx"

#### 2. One R1 / R2 file per tissue (lanes concatenated) ####
mapfile -t TISSUES < <(find "$INDIR" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort)
echo "Found ${#TISSUES[@]} tissue folders: ${TISSUES[*]}"
for t in "${TISSUES[@]}"; do
    [ -s ${OUTDIR_CONCAT}/${t}_1.fastq.gz ] || cat ${INDIR}/${t}/*_1.f*q.gz > ${OUTDIR_CONCAT}/${t}_1.fastq.gz
    [ -s ${OUTDIR_CONCAT}/${t}_2.fastq.gz ] || cat ${INDIR}/${t}/*_2.f*q.gz > ${OUTDIR_CONCAT}/${t}_2.fastq.gz
done

#### 3. HISAT2 alignment, sort, index ####
for t in "${TISSUES[@]}"; do
(
    conda run -n ${HISAT2_ENV} hisat2 -p ${THREADS} -q -x ${REF_DIR}/genome-idx \
        -1 ${OUTDIR_CONCAT}/${t}_1.fastq.gz -2 ${OUTDIR_CONCAT}/${t}_2.fastq.gz \
        -S ${OUTDIR_ALIGN}/${t}.sam 1> ${OUTDIR_ALIGN}/${t}_hisat2.log 2> ${OUTDIR_ALIGN}/${t}_hisat2.err
    samtools view -@${THREADS} -bSh ${OUTDIR_ALIGN}/${t}.sam | samtools sort -@${THREADS} -o ${OUTDIR_ALIGN}/${t}.sorted.bam
    samtools index ${OUTDIR_ALIGN}/${t}.sorted.bam
    rm ${OUTDIR_ALIGN}/${t}.sam
) &
if [[ $(jobs -r -p | wc -l) -ge $PARALLEL_JOBS ]]; then wait -n; fi
done
wait
echo "All alignments finished."
