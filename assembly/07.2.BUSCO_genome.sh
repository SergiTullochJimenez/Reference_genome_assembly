#!/bin/bash
# 07.2 Gene-space completeness of each assembly: BUSCO v6.0.0 in genome mode with squamata_odb12 (container).
# SLURM array, one task per line of SPECIES_LIST.
#SBATCH --job-name=BUSCO_genome
#SBATCH --output=BUSCO_%A_%a.log
#SBATCH --error=BUSCO_%A_%a.err
#SBATCH --account=<account>
#SBATCH --partition=<partition>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=60G
#SBATCH --time=12:00:00
#SBATCH --array=1-<n_species>
set -euo pipefail

#### 0. Configuration ####
WD="/path/to/working_directory"
GENOMES_DIR="/path/to/Genomes"                   # <Species>/<Species>.masked.fa
SPECIES_LIST="/path/to/species.txt"
SIF_IMAGE="/path/to/busco_v6.0.0_cv1.sif"
DOWNLOAD_PATH="${WD}/data/BUSCO_downloads"
# module load apptainer                          # as appropriate on your system

unset BUSCO_CONFIG_FILE
species=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SPECIES_LIST} | tr -d '\r')
GENOME="${GENOMES_DIR}/${species}/${species}.masked.fa"
OUT_DIR="${WD}/res/01.BUSCO_Assessment/${species}"
mkdir -p "${OUT_DIR}" "${DOWNLOAD_PATH}"
cd "${OUT_DIR}"

apptainer exec --cleanenv "${SIF_IMAGE}" \
    busco -i "${GENOME}" -o "BUSCO_${species}" -m genome -l squamata_odb12 \
    -c "${SLURM_CPUS_PER_TASK}" --download_path "${DOWNLOAD_PATH}"
