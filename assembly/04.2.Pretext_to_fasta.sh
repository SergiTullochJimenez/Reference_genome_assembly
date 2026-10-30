#!/bin/bash
# 04.2 Apply the manual curation done in PretextView (exported AGP) to the scaffolds.
#   pretext_to_asm.py: from the Tree of Life agp-tpf-utils (https://github.com/sanger-tol/agp-tpf-utils)
#   rapid_join2.py:    helpers/ (joins the scaffolds following the TPF and writes the curated FASTA)
set -euo pipefail

#### 0. Configuration ####
SPECIES="Species_name"
BASE_DIR="/path/to/project"
AGP_TPF_UTILS="/path/to/agp-tpf-utils"           # cloned repository with its virtual environment in venv/
HELPERS="/path/to/this/repo/assembly/helpers"

WD="${BASE_DIR}/${SPECIES}"
indir="${WD}/res/03.2.Scaffolding/consensus"
indir2="${WD}/res/04.1.Pretext"
outbase="${WD}/res/04.2.Pretext_to_fasta"
mkdir -p ${outbase}

source ${AGP_TPF_UTILS}/venv/bin/activate
python ${AGP_TPF_UTILS}/src/tola/assembly/scripts/pretext_to_asm.py \
    -a ${indir}/_scaffolds_final.fa -p ${indir2}/consensus/out_consensus_fixed.pretext.agp_1 -o ${outbase}/pretextview_out.tpf
python ${HELPERS}/rapid_join2.py -f ${indir}/_scaffolds_final.fa -t ${outbase}/pretextview_out.tpf \
    -c ${outbase}/pretextview_out.chromosome.list.csv -o ${outbase}/curated.fasta
