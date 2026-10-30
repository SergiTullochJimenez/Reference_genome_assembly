#!/bin/bash
# 05.3 Merge the InterProScan and DIAMOND results into the TSEBRA annotation with AGAT, producing the final
# functionally annotated GFF3/GTF, proteins and transcripts.
set -euo pipefail

#### 0. Configuration ####
SPECIES="Species_name"
BASE_DIR="/path/to/project"
AGAT_ENV="agat"
GFFREAD_ENV="tsebra"
BLAST_DB="/path/to/proteomes/colubroidea_all.fasta"     # same database as in 05.1

WD="${BASE_DIR}/${SPECIES}"
GTF="${WD}/res/23.TSEBRA/tsebra.gtf"
INTERPROSCAN_OUT="${WD}/res/24.Functional_annotation/2.Interproscan/tsebra_interproscan.tsv"
BLAST_OUT="${WD}/res/24.Functional_annotation/1.Blastp/tsebra_blastp_stitle.tsv"
REF_GENOME="${WD}/data/reference/${SPECIES}.masked.fa"
OUTDIR="${WD}/res/24.Functional_annotation/3.AGAT"
GFF_TMP="${OUTDIR}/tsebra.gff3"
OUT_GFF="${OUTDIR}/agat.gff3"
OUT_GTF="${OUTDIR}/agat.gtf"
OUT_LOG="${OUTDIR}/agat.log"
export TMPDIR="${WD}/tmp_agat"
mkdir -p "$TMPDIR" "$OUTDIR"
: > "$OUT_LOG"

#### 1. GTF -> GFF3 ####
conda run -n ${AGAT_ENV} agat_convert_sp_gxf2gxf.pl -g "$GTF" -o "$GFF_TMP" >> "$OUT_LOG" 2>&1

#### 2. Add InterProScan + DIAMOND annotations ####
conda run -n ${AGAT_ENV} agat_sp_manage_functional_annotation.pl \
    -f "$GFF_TMP" -i "$INTERPROSCAN_OUT" -b "$BLAST_OUT" -d "$BLAST_DB" --be 1e-5 -a -o "$OUT_GFF" >> "$OUT_LOG" 2>&1

#### 3. Final GTF (AGAT writes the annotated GFF3 inside the -o folder) ####
conda run -n ${AGAT_ENV} agat_convert_sp_gff2gtf.pl --gff "${OUT_GFF}/tsebra.gff3" -o "$OUT_GTF" >> "$OUT_LOG" 2>&1

#### 4. Proteins and transcripts ####
conda run -n ${GFFREAD_ENV} gffread "$OUT_GTF" -g "$REF_GENOME" \
    -y "${OUTDIR}/agat_proteins.fa" -w "${OUTDIR}/agat_transcripts.fa" >> "$OUT_LOG" 2>&1
