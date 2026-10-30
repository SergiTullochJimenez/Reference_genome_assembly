#!/bin/bash
# 06 Repeat annotation: species-specific library with RepeatModeler 2.0.5, soft-masking with RepeatMasker 4.1.6.
# SLURM array, one task per line of SPECIES_LIST.
#SBATCH --job-name=RepeatModeler
#SBATCH --output=RepeatModeler_%A_%a.log
#SBATCH --error=RepeatModeler_%A_%a.err
#SBATCH --account=<account>
#SBATCH --partition=<partition>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH --time=72:00:00
#SBATCH --array=1-<n_species>
set -euo pipefail

#### 0. Configuration ####
WD="/path/to/working_directory"
GENOMES_DIR="/path/to/Genomes"                   # <Species>/<Species>.genome.fa
SPECIES_LIST="/path/to/species.txt"              # one species per line
# module load repeatmodeler/2.0.5                # load RepeatModeler/RepeatMasker as appropriate on your system

species=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SPECIES_LIST})
OUT_DIR="${WD}/res/00.RepeatMasking/${species}"
THREADS=${SLURM_CPUS_PER_TASK}
GENOME="${GENOMES_DIR}/${species}/${species}.genome.fa"
MASKED_GENOME="${GENOMES_DIR}/${species}/${species}.masked.fa"
MASKED_BED="${GENOMES_DIR}/${species}/${species}.masked.bed"
DB="${species}_repeats"
mkdir -p ${OUT_DIR}
cd ${OUT_DIR}
rm -rf RM_* ${DB}*

#### 1. Repeat library ####
BuildDatabase -name ${DB} ${GENOME}
RepeatModeler -database ${DB} -threads ${THREADS}

#### 2. Soft-masking ####
CONS=$(ls RM_*/consensi.fa.classified | head -n 1)
RepeatMasker -pa ${THREADS} -lib ${CONS} -gff -dir ${OUT_DIR} -xsmall ${GENOME}

#### 3. Copy outputs next to the genome ####
cp "${OUT_DIR}/${species}.genome.fa.masked" "${MASKED_GENOME}"
awk 'BEGIN{OFS="\t"} !/^#/ {split($9,a,";"); print $1, $4-1, $5, a[1]}' "${OUT_DIR}/${species}.genome.fa.out.gff" > "${MASKED_BED}"
cp "${OUT_DIR}/${species}.genome.fa.out.gff" "${GENOMES_DIR}/${species}/${species}.masked.gff"
echo "Finished ${species}"
