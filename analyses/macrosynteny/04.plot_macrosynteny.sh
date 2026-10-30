#!/bin/bash
#SBATCH --job-name=macrosynteny_plot
#SBATCH --account=<account>
#SBATCH --partition=<partition>
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=macrosynteny_plot_%j.out
#SBATCH --error=macrosynteny_plot_%j.err

# =============================================================================
# 04.plot_macrosynteny.sh -- final karyotype plot.
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

log() { echo "[$(date '+%H:%M:%S')] $*"; }

SEQIDS="${OUT_DIR}/seqids.txt"
LAYOUT="${OUT_DIR}/layout.txt"

for f in "$SEQIDS" "$LAYOUT"; do
  if [[ ! -s "$f" ]]; then
    echo "ERROR: missing or empty $f -- run 03.run_synteny.sh first." >&2
    exit 1
  fi
done

log "Plotting ${PLOT_OUT} ..."
# Read by karyotype.py's Z-alignment patch (Track.draw()/update_offsets()):
# per-species extra-gap deltas written by 03.run_synteny.sh, so every
# species' Z bar lines up at the same x position (mirrors the already-
# working, fixed-constant W right-alignment).
export KARYOTYPE_Z_GAP_FILE="${OUT_DIR}/z_gap_deltas.tsv"
mamba run -p "$JCVI_ENV_PREFIX" python -m jcvi.graphics.karyotype \
  --chrstyle=roundrect --basepair \
  --figsize="$PLOT_FIGSIZE" --dpi="$PLOT_DPI" --format=pdf \
  --outfile="$PLOT_OUT" \
  "$SEQIDS" "$LAYOUT"

log "Saved: ${PLOT_OUT}"
