#!/usr/bin/env bash
# 00.config.sh -- shared configuration of the macrosynteny pipeline (MCscan / JCVI v1.6.6, DIAMOND v2.2.5).
# Sourced by the other scripts (`source 00.config.sh`), not run directly.
# Inputs per species:  GENOME_DIR/<Species>/<Species>.masked.fa
#                      ANNOT_DIR/<Species>.gtf (or .gff3)   ANNOT_DIR/<Species>_proteins.fa
#                      optional ANNOT_DIR/<Species>.seqid_map.tsv (annotation_id<TAB>genome_id) when the annotation
#                      uses different sequence names than the genome (e.g. RefSeq NC_* vs GenBank accessions)

#### Paths ####
GENOME_DIR="/path/to/Genomes"
ANNOT_DIR="/path/to/Annotations"
OUT_DIR="/path/to/output/Macrosynteny_results"
mkdir -p "$OUT_DIR"

SCRIPT_DIR_CONFIG="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
W_CHROM_FILE="${SCRIPT_DIR_CONFIG}/w_chromosome.txt"
Z_CHROM_FILE="${SCRIPT_DIR_CONFIG}/z_chromosome.txt"
FLIP_FILE="${SCRIPT_DIR_CONFIG}/flip_config.txt"
CHROM_COUNT_FILE="${SCRIPT_DIR_CONFIG}/chromosome_count.txt"   # species with an entry are plotted with their N
                                                               # largest scaffolds only, renamed 1..N / Z / W
MINSPAN_W=5

#### Species order (rows of the plot, top to bottom = ladderised nuclear phylogeny, outgroup excluded) ####
SPECIES_ORDER=(
  "Vipera_aspis"
  "Vipera_seoanei"
  "Vipera_latastei"
  "Malpolon_monspessulanus"
  "Natrix_maura"
  "Natrix_astreptophora"
  "Natrix_helvetica"
  "Hemorrhois_hippocrepis"
  "Hierophis_viridiflavus"
  "Macroprotodon_brevis"
  "Macroprotodon_mauritanicus"
  "Zamenis_scalaris"
  "Zamenis_longissimus"
  "Coronella_girondica"
  "Coronella_austriaca"
)

#### Parameters ####
THREADS="${SLURM_CPUS_PER_TASK:-32}"
DIAMOND_EVALUE=1e-5
DIAMOND_MAX_TARGETS=5
MINSPAN=50
N_PARALLEL=10

#### Conda environment (python3 + jcvi + diamond + samtools; built by 00.setup_env.sh) ####
JCVI_ENV="JCVI"
CONDA_ROOT="/path/to/conda"
unset CONDA_ENVS_PATH                      # conda refuses CONDA_ENVS_PATH and CONDA_ENVS_DIRS together
export CONDA_ENVS_DIRS="${CONDA_ROOT}/envs"
export CONDA_PKGS_DIRS="${CONDA_ROOT}/pkgs"
JCVI_ENV_PREFIX="${CONDA_ENVS_DIRS}/${JCVI_ENV}"
MODULE_MAMBA="miniconda3"                  # module providing mamba/conda on your system

#### SLURM ####
SLURM_ACCOUNT="<account>"

#### Plot ####
PLOT_FIGSIZE="15x12"
PLOT_DPI=500
PLOT_OUT="${OUT_DIR}/macrosynteny_Iberian_snakes.pdf"
