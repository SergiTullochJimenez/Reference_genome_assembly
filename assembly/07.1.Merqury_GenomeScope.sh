#!/bin/bash
# 07.1 Base accuracy and k-mer completeness (Merqury 1.3) and genome-size / heterozygosity estimates (GenomeScope2)
# from the fastp-filtered Illumina reads. k is chosen with Merqury's best_k.sh from the assembly size (k = 20 here).
# Use meryl 1.3 (the build bundled with Merqury 1.3): newer meryl builds produced wrong assembly k-mer sets
# (QV ~17) on these genomes.
set -uo pipefail

#### 0. Configuration ####
SPECIES=(Species_1 Species_2)
GENOME_DIR="/path/to/final_assemblies"            # <Species>.fa
READS_DIR="/path/to/fastp"                        # <Species>.1.fastq.gz / <Species>.2.fastq.gz
OUT_BASE="/path/to/output"
MERYL_BIN="/path/to/meryl-1.3/bin"
MERQURY_HOME="/path/to/merqury-1.3"
GENOMESCOPE_R="/path/to/genomescope.R"
THREADS=16
MERYL_MEMORY_GB=150

export PATH=${MERYL_BIN}:$PATH
MERQURY_OUT="${OUT_BASE}/Merqury"
GENOMESCOPE_OUT="${OUT_BASE}/GenomeScope"
mkdir -p "$MERQURY_OUT" "$GENOMESCOPE_OUT"

run_species() {
    local species="$1"
    local r1="${READS_DIR}/${species}.1.fastq.gz" r2="${READS_DIR}/${species}.2.fastq.gz"
    local genome="${GENOME_DIR}/${species}.fa"
    local sp_merqury="${MERQURY_OUT}/${species}" sp_genomescope="${GENOMESCOPE_OUT}/${species}"
    mkdir -p "$sp_merqury" "$sp_genomescope"
    [[ -s "${genome}.fai" ]] || samtools faidx "$genome"

    #### 1. k from the genome size ####
    local genome_size best_k
    genome_size=$(awk '{sum+=$2} END{print sum}' "${genome}.fai")
    if [[ -s "${sp_merqury}/best_k.txt" ]]; then
        best_k=$(tail -n1 "${sp_merqury}/best_k.txt")
    else
        best_k=$(printf '%.0f' "$("${MERQURY_HOME}/best_k.sh" "$genome_size" | tail -n1)")
        echo "$best_k" > "${sp_merqury}/best_k.txt"
    fi

    #### 2. Read k-mer database ####
    local meryl_db="${sp_merqury}/${species}.k${best_k}.meryl"
    [[ -d "$meryl_db" ]] || meryl k="$best_k" threads="$THREADS" memory="$MERYL_MEMORY_GB" count "$r1" "$r2" output "$meryl_db"

    #### 3. Merqury (QV, completeness) ####
    [[ -s "${sp_merqury}/${species}.qv" ]] || ( cd "$sp_merqury" && "${MERQURY_HOME}/merqury.sh" "$meryl_db" "$genome" "$species" )

    #### 4. GenomeScope2 ####
    local hist="${sp_genomescope}/${species}.k${best_k}.hist"
    [[ -s "$hist" ]] || meryl histogram "$meryl_db" > "$hist"
    [[ -s "${sp_genomescope}/summary.txt" ]] || "$GENOMESCOPE_R" -i "$hist" -o "$sp_genomescope" -k "$best_k" -p 2
    echo "[${species}] done"
}

for sp in "${SPECIES[@]}"; do
    run_species "$sp" || echo "[${sp}] failed, continuing" >&2
done
for sp in "${SPECIES[@]}"; do
    [[ -s "${MERQURY_OUT}/${sp}/${sp}.qv" ]] && printf "%-30s %s\n" "$sp" "$(cat "${MERQURY_OUT}/${sp}/${sp}.qv")"
done
