#!/bin/bash
# 05.2 Decontamination: BLASTn of the curated scaffolds (mitogenome removed) against NCBI nt (e-value < 1e-25),
# accession -> taxid mapping (helpers/fetch_taxids.py), BlobTools with ONT coverage, then keep the scaffolds
# whose phylum is in KEEP_CATEGORIES.
set -euo pipefail

#### 0. Configuration ####
SPECIES="Species_name"
BASE_DIR="/path/to/project"
HELPERS="/path/to/this/repo/assembly/helpers"   # fetch_taxids.py
BLAST_NT_DB="/path/to/blastdb/nt"
TAXDUMP_NAMES="/path/to/taxdump/names.dmp"
TAXDUMP_NODES="/path/to/taxdump/nodes.dmp"
THREADS=25

# Phylum categories to keep (column 6 of the BlobTools table). Inspect the blobplot before choosing.
KEEP_CATEGORIES=("Chordata")

WD="${BASE_DIR}/${SPECIES}"
outdir="${WD}/res/05.2.Decontamination"
ref="${WD}/res/05.1.Mitofinder/Ref_NoMito.fa"
ONT_FQ="${WD}/res/01.1.Trimming/filtered_ONT.fq.gz"
mkdir -p "$outdir"

#### 1. BLAST against nt ####
if [[ ! -s "${outdir}/output2.ncbi.blastn.out" ]]; then
    blastn -task megablast -db "$BLAST_NT_DB" -query "$ref" \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid" \
        -max_target_seqs 10 -max_hsps 1 -evalue 1e-25 -num_threads "$THREADS" \
        -out "${outdir}/output2.ncbi.blastn.out"
fi
[[ -s "${outdir}/cleaned_output2.blastn.out" ]] || sed -E 's/^([^\t]+)\t[^|]+\|[^|]+\|[^|]+\|([^|]+)\|/\1\t\2/' \
    "${outdir}/output2.ncbi.blastn.out" > "${outdir}/cleaned_output2.blastn.out"
[[ -s "${outdir}/sseqids.txt" ]] || cut -f2 "${outdir}/cleaned_output2.blastn.out" | sort | uniq > "${outdir}/sseqids.txt"

#### 2. Accession -> taxid ####
[[ -s "${outdir}/new_taxid_map.tsv" ]] || python "${HELPERS}/fetch_taxids.py" "${outdir}/sseqids.txt" "${outdir}/new_taxid_map.tsv"
[[ -s "${outdir}/new_taxid_map.tsv" ]] || { echo "ERROR: new_taxid_map.tsv is empty -- every scaffold would be 'undef'" >&2; exit 1; }

#### 3. Taxify BLAST hits ####
cd "$outdir"
[[ -s cleaned_output2.blastn.taxified.cleaned_output2.blastn.taxified.out ]] || \
    blobtools taxify -f cleaned_output2.blastn.out -m new_taxid_map.tsv -s 0 -t 1 -o cleaned_output2.blastn.taxified

#### 4. ONT coverage ####
if [[ ! -s "${outdir}/aln_ONT_sorted.bam" ]]; then
    minimap2 -ax map-ont -t "$THREADS" "$ref" "$ONT_FQ" | samtools sort -@5 -O BAM -o "${outdir}/aln_ONT_sorted.bam" -
fi
[[ -s "${outdir}/aln_ONT_sorted.bam.bai" ]] || samtools index "${outdir}/aln_ONT_sorted.bam"

#### 5. BlobDB ####
[[ -s blobdir.blobDB.json ]] || blobtools create -i "$ref" -b aln_ONT_sorted.bam \
    -t cleaned_output2.blastn.taxified.cleaned_output2.blastn.taxified.out -o blobdir \
    --names "$TAXDUMP_NAMES" --nodes "$TAXDUMP_NODES"
[[ -s blobdir.blobdir.blobDB.table.txt ]] || blobtools view -i blobdir.blobDB.json -o blobdir
ls blobdir*blobplot*.png >/dev/null 2>&1 || blobtools plot -i blobdir.blobDB.json -o blobdir

#### 6. Keep scaffolds by phylum (column 6 of the BlobTools table) ####
keep_pattern=$(IFS='|'; echo "${KEEP_CATEGORIES[*]}")
awk -F"\t" -v pat="$keep_pattern" 'BEGIN{split(pat, arr, "|"); for (i in arr) keep[arr[i]]=1} $6 in keep {print $1}' \
    blobdir.blobdir.blobDB.table.txt > keep_scaffolds.txt
echo "Keeping $(wc -l < keep_scaffolds.txt) scaffolds in categories: ${KEEP_CATEGORIES[*]}"
awk 'NR==FNR {keep[$1]; next} /^>/ {show = substr($0, 2) in keep} show' keep_scaffolds.txt "$ref" > Ref_NoMito_decontam.fa
