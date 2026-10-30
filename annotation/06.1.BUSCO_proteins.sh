#!/bin/bash
# 06.1 Annotation completeness: BUSCO v6.0.0 in protein mode with squamata_odb12 on the predicted protein set
# (container). SLURM array, one task per line of SPECIES_LIST.
#SBATCH --job-name=BUSCO_proteins
#SBATCH --output=BUSCO_proteins_%A_%a.log
#SBATCH --error=BUSCO_proteins_%A_%a.err
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
SIF_IMAGE="/path/to/busco_v6.0.0_cv1.sif"
DOWNLOAD_PATH="${WD}/data/BUSCO_downloads"        # with lineages/squamata_odb12 already downloaded
# module load apptainer                          # as appropriate on your system

unset BUSCO_CONFIG_FILE
species=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SPECIES_LIST} | tr -d '\r')
PROTEINS="${ANNOT_DIR}/${species}_proteins.fa"
OUT_DIR="${WD}/res/02.BUSCO_Annotation/${species}"
mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

apptainer exec --cleanenv "${SIF_IMAGE}" \
    busco -i "${PROTEINS}" -o "BUSCO_${species}_proteins" -m proteins \
    -l ${DOWNLOAD_PATH}/lineages/squamata_odb12 -c "${SLURM_CPUS_PER_TASK}" -f \
    --download_path "${DOWNLOAD_PATH}" --offline
