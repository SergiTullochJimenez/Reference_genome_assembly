#!/bin/bash
#SBATCH --job-name=KimuraDivPartitions
#SBATCH --output=KimuraDivPartitions_%A_%a.log
#SBATCH --error=KimuraDivPartitions_%A_%a.err
#SBATCH --account=<account>
#SBATCH --partition=<partition>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --array=1-15%15
# Repeat landscapes per chromosome partition (autosomes, Z, W): for each species (SLURM array task), the RepeatMasker
# alignment file (.genome.fa.cat.gz, assembly/06.Masking.sh) is restricted to the scaffolds of each partition
# (build_partition_scaffold_list.py + filter_cat_by_scaffold.py), and calcDivergenceFromAlign.pl gives the
# CpG-adjusted Kimura divergence summary; the partition size (bp) is written alongside. Partitions with no scaffold
# (W in males, Z not identified) are skipped. Z/W scaffold names: analyses/macrosynteny/{z,w}_chromosome.txt.

#### 0. Configuration ####
SPECIES_LIST="/path/to/species.txt"
REPEATMASKER_DIR="/path/to/RepeatMasker"           # RepeatMasker 4.1.6 installation (provides util/calcDivergenceFromAlign.pl)
module purge
module load repeatmasker/4.1.6                     # as appropriate on your system

set -euo pipefail

CALC_DIV_SCRIPT="${REPEATMASKER_DIR}/util/calcDivergenceFromAlign.pl"
if [ ! -f "${CALC_DIV_SCRIPT}" ]; then
    echo "ERROR: calcDivergenceFromAlign.pl not found at ${CALC_DIV_SCRIPT}"
    echo "Set REPEATMASKER_DIR to your RepeatMasker installation."
    exit 1
fi

species=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SPECIES_LIST})

WD="/path/to/working_directory"
RM_DIR="${WD}/res/00.RepeatMasking/${species}"
OUT_DIR="${WD}/res/01.RepeatLandscape/${species}"

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

W_CHROM_FILE="${SCRIPT_DIR}/../macrosynteny/w_chromosome.txt"
Z_CHROM_FILE="${SCRIPT_DIR}/../macrosynteny/z_chromosome.txt"

echo "Processing task ${SLURM_ARRAY_TASK_ID}: ${species}..."

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

CAT_FILE="${RM_DIR}/${species}.genome.fa.cat.gz"
MASKED_FASTA="${RM_DIR}/${species}.genome.fa.masked"

if [ ! -f "${CAT_FILE}" ]; then
    echo "ERROR: ${CAT_FILE} not found — skipping ${species}"
    exit 1
fi
if [ ! -f "${MASKED_FASTA}" ]; then
    echo "ERROR: ${MASKED_FASTA} not found — skipping ${species}"
    exit 1
fi

FAI="${MASKED_FASTA}.fai"
if [ ! -s "${FAI}" ]; then
    samtools faidx "${MASKED_FASTA}"
fi

for partition in autosomes Z W; do
  echo "---- ${species} / ${partition} ----"

  SCAF_LIST="${OUT_DIR}/${species}.${partition}_scaffolds.txt"
  python3 "${SCRIPT_DIR}/build_partition_scaffold_list.py" \
      "${FAI}" "${W_CHROM_FILE}" "${Z_CHROM_FILE}" "${species}" "${partition}" "${SCAF_LIST}"

  N_SCAF=$(grep -vc '^#' "${SCAF_LIST}" || true)
  if [ "${N_SCAF}" -eq 0 ]; then
    echo "SKIP: ${species} / ${partition} -- empty scaffold list (see ${SCAF_LIST} for reason). No .divsum or genome-size fragment will be written for this combination."
    continue
  fi

  FILTERED_CAT="${OUT_DIR}/${species}.${partition}.cat.gz"
  python3 "${SCRIPT_DIR}/filter_cat_by_scaffold.py" \
      "${CAT_FILE}" "${SCAF_LIST}" "${FILTERED_CAT}"

  perl "${CALC_DIV_SCRIPT}" \
      -s "${species}.${partition}.divsum" \
      -a "${species}.${partition}.aligned.fa" \
      "${FILTERED_CAT}"

  if [ ! -s "${OUT_DIR}/${species}.${partition}.divsum" ]; then
      echo "ERROR: ${species}.${partition}.divsum was not created or is empty — calcDivergenceFromAlign.pl likely failed silently"
      exit 1
  fi
  echo "OK: ${species}.${partition}.divsum written."

  PARTITION_GENOME_SIZE=$(awk 'NR==FNR{if ($0 !~ /^#/ && NF) keep[$1]=1; next} ($1 in keep){sum+=$2} END{print sum+0}' \
      "${SCAF_LIST}" "${FAI}")
  echo "${species},${partition},${PARTITION_GENOME_SIZE}" > "${OUT_DIR}/${species}.${partition}.genome_size.csv"
  echo "${species} / ${partition}: genome size = ${PARTITION_GENOME_SIZE} bp"
done

echo "Finished all partitions for ${species}."
