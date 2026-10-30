#!/bin/bash
#SBATCH --job-name=macrosynteny_setup_env
#SBATCH --account=<account>
#SBATCH --partition=<partition>
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=macrosynteny_setup_env_%j.out
#SBATCH --error=macrosynteny_setup_env_%j.err
# 00.setup_env.sh -- one-time build of the conda environment used by every other script (python3, jcvi, diamond,
# samtools), at the prefix set in 00.config.sh. Safe to re-run: skips creation if the env already exists.
# Run sbatch from this folder: under SLURM the script is copied to a spool directory, so SLURM_SUBMIT_DIR is used
# to find 00.config.sh.
set -euo pipefail
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${SCRIPT_DIR}/00.config.sh"
module purge
module load ${MODULE_MAMBA}
log() { echo "[$(date '+%H:%M:%S')] $*"; }
mkdir -p "${CONDA_ENVS_DIRS}" "${CONDA_PKGS_DIRS}"

if [[ -d "${JCVI_ENV_PREFIX}" && -x "${JCVI_ENV_PREFIX}/bin/python3" ]]; then
  log "Env already exists at '${JCVI_ENV_PREFIX}' -- skipping creation."
else
  log "Creating env at '${JCVI_ENV_PREFIX}'..."
  mamba create -y -p "${JCVI_ENV_PREFIX}" -c conda-forge -c bioconda \
    python=3.10 jcvi diamond samtools numpy scipy matplotlib
fi
mamba run -p "${JCVI_ENV_PREFIX}" python3 -c "import jcvi; print('  jcvi import OK')"
mamba run -p "${JCVI_ENV_PREFIX}" diamond version
mamba run -p "${JCVI_ENV_PREFIX}" samtools --version | head -1
log "Done -- next: 01.build_inputs.sh"
