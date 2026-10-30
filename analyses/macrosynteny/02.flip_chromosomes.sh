#!/bin/bash
#SBATCH --job-name=macrosynteny_flip
#SBATCH --account=<account>
#SBATCH --partition=<partition>
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=macrosynteny_flip_%j.out
#SBATCH --error=macrosynteny_flip_%j.err

# =============================================================================
# 02.flip_chromosomes.sh -- apply chromosome flips from flip_config.txt.
# Run after 01.build_inputs.sh, before 03.run_synteny.sh. No-op if
# flip_config.txt has no entries (nothing downstream changes).
# =============================================================================
set -euo pipefail
# Under sbatch, SLURM copies this script into a per-job spool directory
# (/var/spool/slurmd/jobNNNNNN/) and runs it FROM THERE -- so BASH_SOURCE-based
# directory detection resolves to that spool dir, not the real scripts folder,
# and sourcing sibling files like 00.config.sh fails ("No such file or
# directory") because only this one script was copied there, not its siblings.
# SLURM_SUBMIT_DIR (the directory `sbatch` was invoked FROM) is the reliable
# anchor under SLURM; fall back to BASH_SOURCE detection only when running
# this script directly (not via sbatch). Always `cd` into this scripts folder
# before running `sbatch <script>.sh` (relative path) so the two agree.
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${SCRIPT_DIR}/00.config.sh"

module purge
module load ${MODULE_MAMBA}
MRUN="mamba run -p ${JCVI_ENV_PREFIX}"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

if [[ ! -f "$FLIP_FILE" ]]; then
  log "No flip_config.txt found at ${FLIP_FILE} -- nothing to do."
  exit 0
fi

n_done=0
while IFS=$'\t' read -r SPECIES SCAFFOLDS; do
  [[ -z "${SPECIES:-}" || "$SPECIES" =~ ^# ]] && continue
  [[ -z "${SCAFFOLDS:-}" ]] && continue

  SP_OUT="${OUT_DIR}/${SPECIES}"
  BED="${SP_OUT}/${SPECIES}.longest.bed"
  FAI="${SP_OUT}/${SPECIES}.fa.fai"
  OUT_BED="${SP_OUT}/${SPECIES}.longest.flipped.bed"

  if [[ ! -f "$BED" || ! -f "$FAI" ]]; then
    echo "[${SPECIES}] SKIPPING flip -- run 01.build_inputs.sh first (missing $BED or $FAI)" >&2
    continue
  fi

  ${MRUN} python3 "${SCRIPT_DIR}/flip_bed.py" "$BED" "$FAI" "$SCAFFOLDS" "$OUT_BED"
  n_done=$((n_done + 1))
done < "$FLIP_FILE"

log "Flipped ${n_done} species (per flip_config.txt)."
