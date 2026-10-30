#!/bin/bash
# 04.1 Build the Hi-C contact map for manual curation in PretextView: Hi-C reads mapped to the YaHS scaffolds,
# PretextMap, plus four PretextGraph tracks (repeats, Hi-C coverage, telomeric repeats, gaps).
# Tools: bwa, samtools, bedtools, PretextMap/PretextGraph, tidk, RED (via helpers/repetitive_redmask.py).
set -euo pipefail

#### 0. Configuration ####
SPECIES="Species_name"
BASE_DIR="/path/to/project"
HELPERS="/path/to/this/repo/assembly/helpers"   # fasta_gaps.py, repetitive_redmask.py
PRETEXTMAP="/path/to/PretextMap"
PRETEXTGRAPH="/path/to/PretextGraph"
PYTHON="python"
RED_ENV="RED"                                    # conda env with RED for repetitive_redmask.py
THREADS=24

WD="${BASE_DIR}/${SPECIES}"
indir="${WD}/res/03.2.Scaffolding"
outbase="${WD}/res/04.1.Pretext"
fq1="${WD}/res/01.1.Trimming/filtered_Arima_HIC_1.fq.gz"
fq2="${WD}/res/01.1.Trimming/filtered_Arima_HIC_2.fq.gz"
mkdir -p "$outbase"

{
for version in consensus; do         # add hap1 hap2 to curate the haplotypes too
    ref="${indir}/${version}/_scaffolds_final.fa"
    outdir="${outbase}/${version}"
    mkdir -p "$outdir"
    echo "Processing version: ${version}"
    # keep only the final merged pretext file for this version
    trap 'rm -f "${outdir}/out.pretext" "${outdir}/out_1.pretext" "${outdir}/out_2.pretext" "${outdir}/out_3.pretext"' EXIT

    #### 1. Map the Hi-C reads to the scaffolds ####
    [[ -f "${ref}.bwt" ]] || bwa index "$ref"
    [[ -f "${ref}.fai" ]] || samtools faidx "$ref"
    if [[ ! -f "${outdir}/hic.bam" ]]; then
        bwa mem -t "$THREADS" "$ref" "$fq1" "$fq2" -R "@RG\tID:hic\tSM:hic" \
            | samtools sort -n -@ "$THREADS" -O bam -o "${outdir}/hic.bam" -
    else
        echo "hic.bam already exists, skipping mapping"
    fi

    #### 2. Contact map ####
    if [[ ! -f "${outdir}/out.pretext" ]]; then
        samtools view -h "${outdir}/hic.bam" \
            | "$PRETEXTMAP" -o "${outdir}/out.pretext" --sortby length --sortorder descend --mapq 10
    fi
    bedtools bamtobed -i "${outdir}/hic.bam" > "${outdir}/hic.bed"
    sort -k1,1 -k2,2n "${outdir}/hic.bed" > "${outdir}/hic_sorted.bed"

    #### 3. Gap track ####
    "$PYTHON" "${HELPERS}/fasta_gaps.py" --fasta "$ref" --gap "${outdir}/gaps.bed"
    [[ -s "${outdir}/gaps.bed" ]] || echo "WARNING: gaps.bed is empty" >&2
    chrom_sizes="${outdir}/chrom.sizes"
    sorted_gaps="${outdir}/sorted_gaps.bed"
    out_bedgraph="${outdir}/gaps_vs_background.bedgraph"
    cut -f1,2 "${ref}.fai" | sort -k1,1 > "$chrom_sizes"
    awk 'BEGIN{OFS="\t"} NR>1 {print $1, int($2), int($3)}' "${outdir}/gaps.bed" | sort -k1,1 -k2,2n > "$sorted_gaps"
    # background = genome minus gaps (score 1), gaps = score 100
    bedtools complement -i "$sorted_gaps" -g "$chrom_sizes" | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,1}' > "${outdir}/bg.bed"
    awk 'BEGIN{OFS="\t"}{print $1,$2,$3,100}' "$sorted_gaps" > "${outdir}/gaps_scored.bed"
    cat "${outdir}/bg.bed" "${outdir}/gaps_scored.bed" | sort -k1,1 -k2,2n > "$out_bedgraph"
    rm -f "${outdir}/bg.bed" "${outdir}/gaps_scored.bed"

    #### 4. Telomere track (most abundant repeat found by tidk explore) ####
    mkdir -p "${outdir}/explore"
    explore_tsv="${outdir}/explore/explore_out.tsv"
    (cd "${outdir}/explore" && tidk explore "$ref" --minimum 5 --maximum 12 --log) > "$explore_tsv"
    telomere_motif=""
    if [[ -s "$explore_tsv" ]]; then
        telomere_motif=$(awk -F'\t' 'NR>1 {s=0; for(i=2;i<=NF;i++) s+=$i; print $1, s}' "$explore_tsv" \
            | sort -k2,2 -nr | head -1 | cut -d' ' -f1)
    fi
    [[ -n "$telomere_motif" ]] || { echo "ERROR: no telomere motif found by tidk explore (${explore_tsv})" >&2; exit 1; }
    echo "Using telomere motif: ${telomere_motif}  (check it is the expected TTAGGG/CCCTAA in vertebrates)"
    tidk search "$ref" -d "${outdir}/search" -e tsv      --log -s "$telomere_motif" -o telomeres_out
    tidk search "$ref" -d "${outdir}/search" -e bedgraph --log -s "$telomere_motif" -o telomeres_out
    mkdir -p "${outdir}/plots/tidk-plot"
    tidk plot -t "${outdir}/search/telomeres_out_telomeric_repeat_windows.tsv" --height 120 -w 800 -o "${outdir}/plots/tidk-plot"

    #### 5. Hi-C coverage track ####
    bedtools genomecov -bga -i "${outdir}/hic_sorted.bed" -g "$chrom_sizes" > "${outdir}/coverage.bedgraph"

    #### 6. Repeat track (RED) ####
    mamba run -n "$RED_ENV" python "${HELPERS}/repetitive_redmask.py" -i "$ref" -o "${outdir}/${SPECIES}_rep_py2"
    awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$3-$2}' "${outdir}/${SPECIES}_rep_py2.repeats.bed" > "${outdir}/${SPECIES}_repeats.bedgraph"

    #### 7. Add all tracks to the pretext file ####
    "$PRETEXTGRAPH" -i "${outdir}/out.pretext"   -n repeats   -o "${outdir}/out_1.pretext" < "${outdir}/${SPECIES}_repeats.bedgraph"
    "$PRETEXTGRAPH" -i "${outdir}/out_1.pretext" -n coverage  -o "${outdir}/out_2.pretext" < "${outdir}/coverage.bedgraph"
    "$PRETEXTGRAPH" -i "${outdir}/out_2.pretext" -n telomeres -o "${outdir}/out_3.pretext" < "${outdir}/search/telomeres_out_telomeric_repeat_windows.bedgraph"
    "$PRETEXTGRAPH" -i "${outdir}/out_3.pretext" -n gaps      -o "${outdir}/out_${version}.pretext" < "$out_bedgraph"
done
echo "Finished ${SPECIES}"
} 2>&1 | tee "${outbase}/Pretext.log"
