#!/bin/bash
# 03.2 Hi-C scaffolding with YaHS v1.2 (Arima restriction motif GANTC, no contig error correction).
set -euo pipefail

#### 0. Configuration ####
SPECIES="Species_name"
BASE_DIR="/path/to/project"
YAHS="/path/to/yahs/yahs"

WD="${BASE_DIR}/${SPECIES}"
indir="${WD}/res/02.1.Assembly"
indir2="${WD}/res/03.1.Pre-Scaffolding"
outbase="${WD}/res/03.2.Scaffolding"

for version in consensus; do         # add hap1 hap2 to scaffold the haplotypes too
    ref="${indir}/assembly_${version}.fa"
    outdir="${outbase}/${version}/"
    mkdir -p "${outdir}"
    samtools faidx ${ref}
    ${YAHS} ${ref} ${indir2}/${version}/combined.dedup.bed -e GANTC --no-contig-ec -o ${outdir}
done
