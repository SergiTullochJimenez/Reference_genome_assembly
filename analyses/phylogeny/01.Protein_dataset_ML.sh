#!/bin/bash
#SBATCH --job-name=Nuclear_Phylogeny
#SBATCH --output=Nuclear_Phylogeny_%j.log
#SBATCH --error=Nuclear_Phylogeny_%j.err
#SBATCH --account=<account>
#SBATCH --partition=<partition>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=60G
#SBATCH --time=72:00:00
# Nuclear phylogeny, protein dataset: BUSCO genes complete and single-copy in all species (+ outgroup), aligned with
# MAFFT 7.525 (--auto), trimmed with trimAl 1.4.1 (-automated1), concatenated, and analysed with IQ-TREE 3
# (ModelFinder, 1,000 ultrafast bootstraps, rooted on Anolis sagrei).

set -euo pipefail

WD="/path/to/working_directory"
BUSCO_BASE="${WD}/res/01.BUSCO_Assessment"
LINEAGE="squamata_odb12"

SCRIPT_DIR="/path/to/this/repo/analyses/phylogeny"

SPECIES_LIST="/path/to/species_outgroup.txt"   # one name per line, outgroup included

OUT_DIR="${WD}/res/04.Nuclear_Phylogeny"
RAW_DIR="${OUT_DIR}/01.raw_loci"
ALN_DIR="${OUT_DIR}/02.mafft"
TRIM_DIR="${OUT_DIR}/03.trimal"
SUPERMATRIX_DIR="${OUT_DIR}/04.supermatrix"
TREE_DIR="${OUT_DIR}/05.iqtree"

mkdir -p "${RAW_DIR}" "${ALN_DIR}" "${TRIM_DIR}" "${SUPERMATRIX_DIR}" "${TREE_DIR}"

THREADS=${SLURM_CPUS_PER_TASK:-48}

if [[ ! -f "${SPECIES_LIST}" ]]; then
    echo "ERROR: species list not found: ${SPECIES_LIST}" >&2
    exit 1
fi
N_SPECIES_IN_LIST=$(grep -c . "${SPECIES_LIST}")
echo "Species list: ${SPECIES_LIST} (${N_SPECIES_IN_LIST} entries)"
cat "${SPECIES_LIST}"

module purge
module load python/3.14.0

echo "=== Step 1: extracting shared single-copy BUSCO loci ==="
python3 "${SCRIPT_DIR}/extract_shared_loci.py" \
    --busco-base "${BUSCO_BASE}" \
    --species-list "${SPECIES_LIST}" \
    --out-dir "${RAW_DIR}" \
    --lineage "${LINEAGE}"

N_LOCI=$(ls "${RAW_DIR}"/*.faa 2>/dev/null | wc -l)
echo "Found ${N_LOCI} shared loci"
if [[ "${N_LOCI}" -eq 0 ]]; then
    echo "ERROR: no shared loci extracted, aborting" >&2
    exit 1
fi

module purge
module load mafft/7.525
module load trimal/1.4.1

echo "=== Step 2: aligning + trimming ${N_LOCI} loci ==="

FAILED_LOG="${OUT_DIR}/failed_loci.txt"
: > "${FAILED_LOG}"

align_and_trim () {
    local raw=$1
    local locus
    locus=$(basename "${raw}" .faa)
    local aln="${ALN_DIR}/${locus}.aln.fasta"
    local trm="${TRIM_DIR}/${locus}.fasta"

    if ! mafft --auto --quiet "${raw}" > "${aln}" 2> "${ALN_DIR}/${locus}.mafft.log"; then
        echo "${locus}: mafft failed" >> "${FAILED_LOG}"
        return 0
    fi
    if ! trimal -in "${aln}" -out "${trm}" -automated1 2> "${TRIM_DIR}/${locus}.trimal.log"; then
        echo "${locus}: trimal failed" >> "${FAILED_LOG}"
        return 0
    fi

    [[ -s "${ALN_DIR}/${locus}.mafft.log" ]] || rm -f "${ALN_DIR}/${locus}.mafft.log"
    [[ -s "${TRIM_DIR}/${locus}.trimal.log" ]] || rm -f "${TRIM_DIR}/${locus}.trimal.log"
    return 0
}
export -f align_and_trim
export ALN_DIR TRIM_DIR FAILED_LOG

for raw in "${RAW_DIR}"/*.faa; do
    align_and_trim "${raw}" &
    if [[ $(jobs -r -p | wc -l) -ge ${THREADS} ]]; then
        wait -n
    fi
done
wait

if [[ -s "${FAILED_LOG}" ]]; then
    echo "WARNING: $(wc -l < "${FAILED_LOG}") loci failed alignment/trimming -- see ${FAILED_LOG}"
fi

N_TRIMMED=$(ls "${TRIM_DIR}"/*.fasta 2>/dev/null | wc -l)
echo "Trimmed alignments produced: ${N_TRIMMED} / ${N_LOCI}"

module purge
module load python/3.14.0

echo "=== Step 3: concatenating into supermatrix ==="
python3 "${SCRIPT_DIR}/concatenate_alignments.py" \
    --aln-dir "${TRIM_DIR}" \
    --species-list "${SPECIES_LIST}" \
    --out-fasta "${SUPERMATRIX_DIR}/supermatrix.faa" \
    --out-partition "${SUPERMATRIX_DIR}/partitions.txt"

module purge
module load iq-tree/3.1.3

echo "=== Step 4: IQ-TREE ==="

IQTREE_BIN=$(command -v iqtree3 || command -v iqtree2 || command -v iqtree || true)
if [[ -z "${IQTREE_BIN}" ]]; then
    echo "ERROR: no iqtree/iqtree2/iqtree3 binary found on PATH after module load" >&2
    exit 1
fi
echo "Using IQ-TREE binary: ${IQTREE_BIN}"

cd "${TREE_DIR}"
"${IQTREE_BIN}" \
    -s "${SUPERMATRIX_DIR}/supermatrix.faa" \
    -m MFP \
    -B 1000 \
    -T AUTO --threads-max "${THREADS}" \
    -o Anolis_sagrei \
    --prefix nuclear_busco_tree

echo "Done. Best tree: ${TREE_DIR}/nuclear_busco_tree.treefile"
echo "(partition file for a future partitioned re-run: ${SUPERMATRIX_DIR}/partitions.txt)"
