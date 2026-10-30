#!/bin/bash
#SBATCH --job-name=MaskedFraction
#SBATCH --output=MaskedFraction_%A_%a.log
#SBATCH --error=MaskedFraction_%A_%a.err
#SBATCH --account=<account>
#SBATCH --partition=<partition>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=1-15%15
# Repeat content per chromosome partition (autosomes, Z, W): fraction of soft-masked bases of each partition in the
# masked assembly, per species (SLURM array task). Same partition definition as 01.KimuraDiv_partitions.sh.

#### 0. Configuration ####
SPECIES_LIST="/path/to/species.txt"
module purge
module load samtools                               # as appropriate on your system

set -euo pipefail

species=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SPECIES_LIST})

WD="/path/to/working_directory"
GENOME_DIR="${WD}/data/Genomes/${species}"
OUT_DIR="${WD}/res/01.RepeatLandscape/${species}"

MASKED_FASTA="${GENOME_DIR}/${species}.masked.fa"
FAI="${MASKED_FASTA}.fai"

if [ ! -f "${MASKED_FASTA}" ]; then
    echo "ERROR: ${MASKED_FASTA} not found — skipping ${species}"
    exit 1
fi
if [ ! -s "${FAI}" ]; then
    samtools faidx "${MASKED_FASTA}"
fi

echo "Processing task ${SLURM_ARRAY_TASK_ID}: ${species}..."

for partition in autosomes Z W; do
  echo "---- ${species} / ${partition} ----"

  SCAF_LIST="${OUT_DIR}/${species}.${partition}_scaffolds.txt"
  if [ ! -f "${SCAF_LIST}" ]; then
    echo "SKIP: ${species} / ${partition} -- no scaffold list at ${SCAF_LIST} " \
         "(run KimuraDiv_partitions.sh's build_partition_scaffold_list.py step first). "
    continue
  fi

  N_SCAF=$(grep -vc '^#' "${SCAF_LIST}" || true)
  if [ "${N_SCAF}" -eq 0 ]; then
    echo "SKIP: ${species} / ${partition} -- empty scaffold list (see ${SCAF_LIST} for reason)."
    continue
  fi

  REGION_FILE="${OUT_DIR}/${species}.${partition}_scaffolds.regions.tmp"
  grep -v '^#' "${SCAF_LIST}" | awk 'NF' > "${REGION_FILE}"

  MISSING=$(awk 'NR==FNR{fai[$1]=1; next} !( $1 in fai ){print $1}' "${FAI}" "${REGION_FILE}")
  if [ -n "${MISSING}" ]; then
    echo "ERROR: ${species} / ${partition} -- scaffold(s) in ${SCAF_LIST} not found in ${FAI}:"
    echo "${MISSING}"
    echo "This means the scaffold list was built against a DIFFERENT masked FASTA than" \
         "${MASKED_FASTA} (naming mismatch, e.g. a renamed/Old_ copy). Fix the path/name" \
         "before trusting this species/partition's number. Skipping for now."
    rm -f "${REGION_FILE}"
    continue
  fi

  PARTITION_FASTA="${OUT_DIR}/${species}.${partition}.subset.fa.tmp"
  samtools faidx "${MASKED_FASTA}" -r "${REGION_FILE}" > "${PARTITION_FASTA}"

  read -r TOTAL_BP MASKED_BP <<< "$(awk '
      /^>/ { next }
      {
          total += length($0)
          masked += gsub(/[acgt]/, "&")
      }
      END { print total+0, masked+0 }
  ' "${PARTITION_FASTA}")"

  rm -f "${REGION_FILE}" "${PARTITION_FASTA}"

  if [ "${TOTAL_BP}" -eq 0 ]; then
    echo "WARNING: ${species} / ${partition} -- 0 total bp extracted (scaffold name mismatch?). Skipping."
    continue
  fi

  PCT=$(python3 -c "print(f'{${MASKED_BP}/${TOTAL_BP}*100:.2f}')")

  OUT_CSV="${OUT_DIR}/${species}.${partition}.masked_fraction.csv"
  echo "${species},${partition},${TOTAL_BP},${MASKED_BP},${PCT}" > "${OUT_CSV}"
  echo "${species} / ${partition}: total=${TOTAL_BP} masked=${MASKED_BP} -> ${PCT}% -> ${OUT_CSV}"
done

echo "Finished all partitions for ${species}."
