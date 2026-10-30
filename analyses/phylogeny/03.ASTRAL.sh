#!/bin/bash
#SBATCH --job-name=astral_tree_both
#SBATCH --account=<account>
#SBATCH --partition=<partition>
#SBATCH --cpus-per-task=48
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --output=astral_tree_both_%j.out
#SBATCH --error=astral_tree_both_%j.err
# Coalescent species trees with ASTRAL 5.7.1 for both datasets. One ML gene tree per trimmed locus with IQ-TREE 3
# (fixed model: LG+G for proteins, GTR+G for nucleotides), then ASTRAL with full branch annotation (-t 2: quartet
# frequencies and local posterior probabilities of the three resolutions of every branch).
# Inputs: the trimmed loci of 01.Protein_dataset_ML.sh and 02.Nucleotide_dataset_ML.sh.

set -uo pipefail   # no -e: per-locus gene-tree failures are logged and skipped

WD="/path/to/working_directory"
ASTRAL_JAR="/path/to/ASTRAL-5.7.1/Astral/astral.5.7.1.jar"
THREADS="${SLURM_CPUS_PER_TASK:-48}"

build_gene_tree() {
    local locus_path="$1" ext="$2" gene_trees_dir="$3" failed_log="$4" astral_dir="$5" model="$6" st_flag="$7"
    local locus prefix
    locus="$(basename "${locus_path}" ".${ext}")"
    prefix="${gene_trees_dir}/${locus}"

    local st_args=()
    [ -n "${st_flag}" ] && st_args=(-st "${st_flag}")

    if ! iqtree3 -s "${locus_path}" "${st_args[@]}" -m "${model}" -T 1 --quiet --prefix "${prefix}" \
         > /dev/null 2>>"${astral_dir}/iqtree_genetrees.err"; then
        echo "GENE_TREE_FAIL: ${locus}" >> "${failed_log}"
        return 0
    fi
    if [ ! -s "${prefix}.treefile" ]; then
        echo "GENE_TREE_MISSING: ${locus}" >> "${failed_log}"
        return 0
    fi
    return 0
}
export -f build_gene_tree

run_astral_pipeline() {
    local label="$1" trimmed_dir="$2" ext="$3" astral_dir="$4" model="$5" st_flag="$6"
    local gene_trees_dir="${astral_dir}/gene_trees"
    local failed_log="${astral_dir}/failed_gene_trees.txt"
    local final_tree="${astral_dir}/astral_species_tree.tre"
    local tmp_tree="${astral_dir}/astral_species_tree.tre.tmp"

    echo ""
    echo "################################################################"
    echo "### ${label}"
    echo "################################################################"

    rm -rf "${astral_dir}"
    mkdir -p "${gene_trees_dir}"
    : > "${failed_log}"

    local n_loci
    n_loci=$(ls "${trimmed_dir}"/*."${ext}" 2>/dev/null | wc -l)
    if [ "${n_loci}" -eq 0 ]; then
        echo "ERROR [${label}]: no *.${ext} files found in ${trimmed_dir}. Aborting this part." >&2
        echo "Contents of ${trimmed_dir}:" >&2
        ls -la "${trimmed_dir}" >&2 || true
        return 1
    fi
    echo "[${label}] Found ${n_loci} loci in ${trimmed_dir} (*.${ext})."

    echo "=== [${label}] Step 1: one quick ML gene tree per locus (fixed ${model}, no ModelFinder) ==="
    module purge
    module load iq-tree/3.1.3

    export GENE_TREES_DIR="${gene_trees_dir}" FAILED_LOG="${failed_log}" ASTRAL_DIR="${astral_dir}"
    ls "${trimmed_dir}"/*."${ext}" | \
        xargs -P "${THREADS}" -I{} bash -c 'build_gene_tree "$1" "$2" "$3" "$4" "$5" "$6" "$7"' _ \
            {} "${ext}" "${gene_trees_dir}" "${failed_log}" "${astral_dir}" "${model}" "${st_flag}"
    wait

    local n_trees
    n_trees=$(ls "${gene_trees_dir}"/*.treefile 2>/dev/null | wc -l)
    echo "[${label}] Built ${n_trees} / ${n_loci} per-locus gene trees."
    if [ "${n_trees}" -eq 0 ]; then
        echo "ERROR [${label}]: zero gene trees were built -- refusing to run ASTRAL on empty input." >&2
        echo "Check ${failed_log} and ${astral_dir}/iqtree_genetrees.err" >&2
        return 1
    fi

    echo "=== [${label}] Step 2: collecting gene trees ==="
    cat "${gene_trees_dir}"/*.treefile > "${astral_dir}/gene_trees_all.tre"

    echo "=== [${label}] Step 3: running ASTRAL ==="
    module purge
    module load java

    rm -f "${tmp_tree}"
    if ! java -jar "${ASTRAL_JAR}" \
            -i "${astral_dir}/gene_trees_all.tre" \
            -o "${tmp_tree}" \
            -t 2 \
            2> "${astral_dir}/astral.log"; then
        echo "ERROR [${label}]: ASTRAL exited non-zero. See ${astral_dir}/astral.log" >&2
        return 1
    fi
    if [ ! -s "${tmp_tree}" ]; then
        echo "ERROR [${label}]: ASTRAL reported success but ${tmp_tree} is empty/missing." >&2
        return 1
    fi
    mv "${tmp_tree}" "${final_tree}"

    echo "[${label}] Done. Species tree: ${final_tree}"
    echo "[${label}] md5sum: $(md5sum "${final_tree}")"
    return 0
}

NUCLEAR_PHYLO_DIR="${WD}/res/04.Nuclear_Phylogeny"
TRIMMED_DIR_PROTEIN="${NUCLEAR_PHYLO_DIR}/03.trimal"
LOCUS_EXT_PROTEIN="fasta"
ASTRAL_DIR_PROTEIN="${NUCLEAR_PHYLO_DIR}/10.ASTRAL"

protein_status=0
run_astral_pipeline "PROTEIN dataset" \
    "${TRIMMED_DIR_PROTEIN}" "${LOCUS_EXT_PROTEIN}" "${ASTRAL_DIR_PROTEIN}" \
    "LG+G" "AA" || protein_status=1

GENOMIC_TREE_DIR="${NUCLEAR_PHYLO_DIR}/Genomic_tree"
TRIMMED_DIR_NUC="${GENOMIC_TREE_DIR}/03.loci_nuc_trimmed"
LOCUS_EXT_NUC="fna"
ASTRAL_DIR_NUC="${GENOMIC_TREE_DIR}/10.ASTRAL"

nuc_status=0
run_astral_pipeline "NUCLEOTIDE dataset" \
    "${TRIMMED_DIR_NUC}" "${LOCUS_EXT_NUC}" "${ASTRAL_DIR_NUC}" \
    "GTR+G" "" || nuc_status=1

echo ""
echo "################################################################"
echo "### Final check"
echo "################################################################"

if [ "${protein_status}" -ne 0 ] || [ "${nuc_status}" -ne 0 ]; then
    echo "At least one part FAILED (protein_status=${protein_status}, nuc_status=${nuc_status})." >&2
    exit 1
fi

TREE_PROTEIN="${ASTRAL_DIR_PROTEIN}/astral_species_tree.tre"
TREE_NUC="${ASTRAL_DIR_NUC}/astral_species_tree.tre"

if cmp -s "${TREE_PROTEIN}" "${TREE_NUC}"; then
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" >&2
    echo "ERROR: the protein and nucleotide ASTRAL trees are BYTE-IDENTICAL." >&2
    echo "That should be essentially impossible for two real independent runs" >&2
    echo "on different gene-tree sets. Do NOT use these outputs -- something is" >&2
    echo "still wrong (e.g. one input loci directory is empty/wrong and the" >&2
    echo "other's trees leaked in, or a filesystem/symlink issue)." >&2
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" >&2
    exit 1
fi

echo "OK: protein and nucleotide ASTRAL trees differ, as expected."
echo "  protein tree     : ${TREE_PROTEIN}"
echo "                     $(md5sum "${TREE_PROTEIN}")"
echo "  nucleotide tree  : ${TREE_NUC}"
echo "                     $(md5sum "${TREE_NUC}")"
echo ""
echo "NOTE: neither tree is rooted on the outgroup yet -- ASTRAL doesn't take"
echo "an outgroup flag. Root both .tre files on Anolis_sagrei afterwards, e.g. in R:"
echo "  library(ape)"
echo "  t <- read.tree('astral_species_tree.tre')"
echo "  t <- root(t, outgroup = 'Anolis_sagrei', resolve.root = TRUE)"
echo "  write.tree(t, 'astral_species_tree_rooted.tre')"
