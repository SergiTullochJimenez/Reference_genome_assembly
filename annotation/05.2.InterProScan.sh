#!/bin/bash
# 05.2 Functional annotation, protein domains / families / GO terms with InterProScan 5.66-98.0.
set -euo pipefail

#### 0. Configuration ####
SPECIES="Species_name"
BASE_DIR="/path/to/project"
INTERPROSCAN="/path/to/interproscan-5.66-98.0/interproscan.sh"
THREADS=42
export _JAVA_OPTIONS="-Xmx32G"

WD="${BASE_DIR}/${SPECIES}"
PROTEINS="${WD}/res/23.TSEBRA/tsebra_proteins.clean.fa"      # from 05.1
OUT="${WD}/res/24.Functional_annotation/2.Interproscan"
mkdir -p "$OUT"

"$INTERPROSCAN" -i "$PROTEINS" -f TSV,XML,GFF3 -goterms -pa -dp -cpu ${THREADS} \
    -b "${OUT}/tsebra_interproscan" > "${OUT}/tsebra_interproscan.log" 2>&1
