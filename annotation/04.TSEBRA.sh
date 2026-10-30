#!/bin/bash
# 04 Combine the BRAKER3 and GeMoMa gene sets with TSEBRA (configuration: tsebra_config.cfg), using RNA-seq intron
# hints (bam2hints) and the BRAKER hints; extract proteins and transcripts with gffread.
set -euo pipefail

#### 0. Configuration ####
SPECIES="Species_name"
BASE_DIR="/path/to/project"
TSEBRA_ENV="tsebra"                       # conda env with tsebra, gffread, bam2hints (AUGUSTUS tools)
TSEBRA_CFG="/path/to/this/repo/annotation/tsebra_config.cfg"

WD="${BASE_DIR}/${SPECIES}"
RNA_DIR="${WD}/res/22.2.Aligned_RNA_tissues"
gemoma="${WD}/res/22.3.GeMoMa"
OUT="${WD}/res/23.TSEBRA"
BRAKER_GTF="${WD}/res/21.Braker/braker.gtf"
BRAKER_HINTS="${WD}/res/21.Braker/hintsfile.gff"
GEMOMA_GTF="${gemoma}/gemoma.gtf"
RNA_HINTS="${gemoma}/rnaseq_hints.gff"
ref="${WD}/data/reference/${SPECIES}.masked.fa"
OUT_COMBINED="${OUT}/tsebra.gtf"
mkdir -p "$OUT"

#### 1. GeMoMa GFF -> GTF ####
conda run -n ${TSEBRA_ENV} gffread "${gemoma}/final_annotation.gff" -T -o "$GEMOMA_GTF"

#### 2. Intron hints from the per-tissue BAMs ####
: > "${RNA_HINTS}"
for bam in "${RNA_DIR}"/*.sorted.bam; do
    conda run -n ${TSEBRA_ENV} bam2hints --intronsonly --in="${bam}" --out=/dev/stdout | sed '/^$/d' >> "${RNA_HINTS}"
done

#### 3. TSEBRA ####
conda run -n ${TSEBRA_ENV} tsebra.py \
    -g "${BRAKER_GTF},${GEMOMA_GTF}" \
    -c "${TSEBRA_CFG}" \
    -e "${RNA_HINTS},${BRAKER_HINTS}" \
    -o "${OUT_COMBINED}"

#### 4. Proteins and transcripts ####
conda run -n ${TSEBRA_ENV} gffread "${OUT_COMBINED}" -g "$ref" -y "${OUT}/tsebra_proteins.fa" -w "${OUT}/tsebra_transcripts.fa"
