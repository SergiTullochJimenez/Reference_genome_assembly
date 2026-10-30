#!/bin/bash
# 06.2 Annotation consistency with OMArk: OMAmer search against the LUCA database, then OMArk with the Squamata
# taxid (8509). SLURM array, one task per line of SPECIES_LIST.
#SBATCH --job-name=OMARK
#SBATCH --output=OMARK_%A_%a.log
#SBATCH --error=OMARK_%A_%a.err
#SBATCH --account=<account>
#SBATCH --partition=<partition>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=60G
#SBATCH --time=06:00:00
#SBATCH --array=1-<n_species>
set -euo pipefail

#### 0. Configuration ####
WD="/path/to/working_directory"
ANNOT_DIR="/path/to/Annotations"                 # <Species>_proteins.fa
SPECIES_LIST="/path/to/species.txt"
DB="/path/to/OMARK_DATABASE/LUCA.h5"
TAXID="8509"
# source activate omark                          # conda env with omamer and omark

species=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SPECIES_LIST} | tr -d '\r')
PROTEINS="${ANNOT_DIR}/${species}_proteins.fa"
OUT_DIR="${WD}/res/02.2.OMARK/${species}"
mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

#### 1. OMAmer ####
omamer search --db "${DB}" --query "${PROTEINS}" --out "${OUT_DIR}/${species}.omamer"

#### 2. OMArk ####
omark -f "${OUT_DIR}/${species}.omamer" -d "${DB}" -o "${OUT_DIR}" -t "${TAXID}"
