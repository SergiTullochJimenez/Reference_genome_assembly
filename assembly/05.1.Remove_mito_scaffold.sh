#!/bin/bash
# 05.1 Identify the scaffold(s) carrying the mitochondrial genome in the curated assembly and remove them, so the
# nuclear assembly (Ref_NoMito.fa) can be decontaminated (05.2). The mitogenome itself is assembled in 08.Mitogenome.sh.
# Run in two passes: (1) MitoFinder + BLAST to find the mito scaffold(s); inspect mito_vs_ref.blast;
#                    (2) set MITO_SCAFFOLDS and rerun to write Ref_NoMito.fa.
set -euo pipefail

#### 0. Configuration ####
SPECIES="Species_name"
BASE_DIR="/path/to/project"
MITO_REFERENCE="/path/to/mito_reference.gb"      # GenBank mitogenome of a related species
MITO_SCAFFOLDS=""                                # e.g. "scaffold_42" (space-separated; can be more than one copy)
THREADS=24

WD="${BASE_DIR}/${SPECIES}"
ref="${WD}/res/04.2.Pretext_to_fasta/curated.fasta"
outbase="${WD}/res/05.1.Mitofinder"
mkdir -p ${outbase}
cd ${WD}

#### 1. Find mitochondrial genes in the assembly ####
if [[ ! -s ${outbase}/mito_vs_ref.blast ]]; then
    mitofinder -j comprobation -s ${ref} -r ${MITO_REFERENCE} -o 2 -p ${THREADS} -m 24
    mv comprobation/comprobation_MitoFinder_megahit_mitfi_Final_Results/* ${outbase}
    makeblastdb -in ${ref} -dbtype nucl -out ${outbase}/ref_db
    blastn -query ${outbase}/comprobation_final_genes_NT.fasta -db ${outbase}/ref_db -out ${outbase}/mito_vs_ref.blast \
        -evalue 1e-20 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
        -max_target_seqs 1 -num_threads ${THREADS}
fi
echo "Scaffolds hit by mitochondrial genes:"; cut -f2 ${outbase}/mito_vs_ref.blast | sort | uniq -c

#### 2. Remove the mito scaffold(s) ####
if [[ -z "${MITO_SCAFFOLDS}" ]]; then
    echo "Set MITO_SCAFFOLDS after checking ${outbase}/mito_vs_ref.blast, then rerun."; exit 0
fi
pattern=$(echo ${MITO_SCAFFOLDS} | tr ' ' '\n' | sed 's/.*/^&$/' | paste -sd'|' -)
seqkit grep -v -r -p "${pattern}" ${ref} > ${outbase}/Ref_NoMito.fa
seqkit grep -r -p "${pattern}" ${ref} > ${outbase}/Mitochondrial_${SPECIES}.fa
echo "Removed: ${MITO_SCAFFOLDS}"
