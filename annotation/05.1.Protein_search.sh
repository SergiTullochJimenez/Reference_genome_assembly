#!/bin/bash
# 05.1 Functional annotation, similarity search.
#   (a) Clean the TSEBRA proteins: strip a trailing stop character and drop proteins with an internal stop
#       (broken gene models).
#   (b) DIAMOND blastp against a custom Colubroidea protein database (e-value 1e-5, up to 5 targets).
set -euo pipefail

#### 0. Configuration ####
SPECIES="Species_name"
BASE_DIR="/path/to/project"
ENV="tsebra"                                           # conda env with diamond and seqkit
DB="/path/to/proteomes/colubroidea_all.fasta"          # custom Colubroidea protein database
EVALUE=1e-5
MAX_TARGET_SEQS=5
THREADS=32

WD="${BASE_DIR}/${SPECIES}"
RAW_PROTEINS="${WD}/res/23.TSEBRA/tsebra_proteins.fa"
CLEAN_PROTEINS="${WD}/res/23.TSEBRA/tsebra_proteins.clean.fa"
DROPPED_LOG="${WD}/res/23.TSEBRA/tsebra_proteins.dropped_internal_stop.txt"
OUT="${WD}/res/24.Functional_annotation/1.Blastp"
OUTFILE="${OUT}/tsebra_blastp_stitle.tsv"
mkdir -p "$OUT"

#### 1. Clean proteins ####
seqkit seq -w 0 "$RAW_PROTEINS" | awk -v dropped="$DROPPED_LOG" '
/^>/ {
    header = $0
    getline seq
    core = seq
    sub(/[.*]$/, "", core)          # strip a single trailing stop
    if (core ~ /[.*]/) {            # anything left = internal stop
        id = substr(header, 2); sub(/ .*/, "", id)
        print id >> dropped; n_dropped++; next
    }
    print header; print core; n_kept++
}
END { print "Kept: " n_kept > "/dev/stderr"; print "Dropped (internal stop): " n_dropped+0 > "/dev/stderr" }
' > "$CLEAN_PROTEINS"

#### 2. DIAMOND blastp ####
DB_BASE="${DB%.fasta}"
[ -f "${DB_BASE}.dmnd" ] || conda run -n ${ENV} diamond makedb --in "$DB" --db "$DB_BASE"
conda run -n ${ENV} diamond blastp \
    --query "$CLEAN_PROTEINS" --db "$DB_BASE" --out "$OUTFILE" \
    --evalue ${EVALUE} --threads ${THREADS} --max-target-seqs ${MAX_TARGET_SEQS} \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp stitle
