#!/bin/bash
#SBATCH --job-name=genomic_tree
#SBATCH --account=<account>
#SBATCH --partition=<partition>
#SBATCH --cpus-per-task=48
#SBATCH --mem=100G
#SBATCH --time=72:00:00
#SBATCH --output=genomic_tree_%j.out
#SBATCH --error=genomic_tree_%j.err
# Nuclear phylogeny, nucleotide dataset (independent of the protein dataset, no back-translation): the CDS of each
# shared BUSCO locus is extracted from each genome using the BUSCO/miniprot gene coordinates, aligned with MAFFT 7.525,
# trimmed with trimAl 1.4.1 (-automated1), concatenated, and analysed with IQ-TREE 3 (ModelFinder, 1,000 ultrafast
# bootstraps, rooted on Anolis sagrei).

set -uo pipefail   # no -e: per-locus align/trim failures are logged and skipped

SCRIPT_DIR="/path/to/this/repo/analyses/phylogeny"
WD="/path/to/working_directory"
SPECIES_LIST="/path/to/species_outgroup.txt"   # one name per line, outgroup included

BUSCO_BASE="${WD}/res/01.BUSCO_Assessment"
GENOME_DIR="${WD}/data/Genomes"
LINEAGE="squamata_odb12"

OUT_DIR="${WD}/res/04.Nuclear_Phylogeny/Genomic_tree"
RAW_DIR="${OUT_DIR}/01.loci_nuc_raw"
ALIGNED_DIR="${OUT_DIR}/02.loci_nuc_aligned"
TRIMMED_DIR="${OUT_DIR}/03.loci_nuc_trimmed"
SUPERMATRIX_DIR="${OUT_DIR}/04.supermatrix_nuc"
IQTREE_DIR="${OUT_DIR}/05.iqtree_nuc"
FAILED_LOG="${OUT_DIR}/failed_loci.txt"

mkdir -p "${RAW_DIR}" "${ALIGNED_DIR}" "${TRIMMED_DIR}" "${SUPERMATRIX_DIR}" "${IQTREE_DIR}"
: > "${FAILED_LOG}"

THREADS="${SLURM_CPUS_PER_TASK:-48}"

echo "=== Step 1: extracting nucleotide CDS for shared single-copy BUSCO loci ==="
module purge
module load python/3.14.0
module load samtools/1.19

python3 "${SCRIPT_DIR}/extract_nuc_loci.py" \
    --busco-base "${BUSCO_BASE}" \
    --species-list "${SPECIES_LIST}" \
    --lineage "${LINEAGE}" \
    --genome-dir "${GENOME_DIR}" \
    --out-dir "${RAW_DIR}"

N_LOCI=$(ls "${RAW_DIR}"/*.fna 2>/dev/null | wc -l)
echo "Extracted ${N_LOCI} raw nucleotide loci."

echo "=== Step 2: aligning + trimming ${N_LOCI} loci (nucleotide) ==="
module purge
module load mafft/7.525
module load trimal/1.4.1

align_and_trim() {
    local locus_path="$1"
    local locus
    locus="$(basename "${locus_path}" .fna)"
    local aligned="${ALIGNED_DIR}/${locus}.fna"
    local trimmed="${TRIMMED_DIR}/${locus}.fna"

    if ! mafft --auto --thread 1 "${locus_path}" > "${aligned}" 2>>"${OUT_DIR}/mafft.err"; then
        echo "ALIGN_FAIL: ${locus}" >> "${FAILED_LOG}"
        return 0
    fi
    if ! trimal -in "${aligned}" -out "${trimmed}" -automated1 2>>"${OUT_DIR}/trimal.err"; then
        echo "TRIM_FAIL: ${locus}" >> "${FAILED_LOG}"
        return 0
    fi
    if [ ! -s "${trimmed}" ]; then
        echo "EMPTY_AFTER_TRIM: ${locus}" >> "${FAILED_LOG}"
        return 0
    fi
    return 0
}
export -f align_and_trim
export ALIGNED_DIR TRIMMED_DIR FAILED_LOG OUT_DIR

ls "${RAW_DIR}"/*.fna | xargs -P "${THREADS}" -I{} bash -c 'align_and_trim "$@"' _ {}
wait

N_TRIMMED=$(ls "${TRIMMED_DIR}"/*.fna 2>/dev/null | wc -l)
echo "Trimmed alignments produced: ${N_TRIMMED} / ${N_LOCI}"

echo "=== Step 3: concatenating into supermatrix ==="
module purge
module load python/3.14.0

python3 "${SCRIPT_DIR}/concatenate_nuc_loci.py" \
    --aligned-dir "${TRIMMED_DIR}" \
    --species-list "${SPECIES_LIST}" \
    --out-fasta "${SUPERMATRIX_DIR}/supermatrix.fna" \
    --out-partition "${SUPERMATRIX_DIR}/partitions.txt"

echo "=== Step 4: IQ-TREE (nucleotide/codon model selection) ==="
module purge
module purge
module load iq-tree/3.1.3

echo "=== Step 4: IQ-TREE ==="

IQTREE_BIN=$(command -v iqtree3 || command -v iqtree2 || command -v iqtree || true)
if [[ -z "${IQTREE_BIN}" ]]; then
    echo "ERROR: no iqtree/iqtree2/iqtree3 binary found on PATH after module load" >&2
    exit 1
fi
echo "Using IQ-TREE binary: ${IQTREE_BIN}"

"${IQTREE_BIN}" \
    -s "${SUPERMATRIX_DIR}/supermatrix.fna" \
    -m MFP \
    -B 1000 \
    -T AUTO --threads-max "${THREADS}" \
    -o Anolis_sagrei \
    --prefix "${IQTREE_DIR}/nuclear_busco_tree_genomic"

echo "Done. Best tree: ${IQTREE_DIR}/nuclear_busco_tree_genomic.treefile"
echo "(partition file for a future partitioned re-run: ${SUPERMATRIX_DIR}/partitions.txt)"
echo "Failed loci (if any) logged to: ${FAILED_LOG}"
