#!/bin/bash
# =============================================================================
# hic_junction_check.sh -- are the fusion junctions real or scaffolding misjoins?
#
# For each junction (micro-chromosome fusions found by the breakpoint pipeline):
#  1. Hi-C "crossing" profile: read pairs with one mate in the 1 Mb left of a point
#     and the other in the 1 Mb right of it, for every 100 kb point along the
#     chromosome. A misjoin = a deep dip at the junction (contacts there fall to
#     between-chromosome levels); a real fusion = no deeper than other points.
#  2. ONT reads with alignments on BOTH sides of the junction (within 50 kb),
#     counting split/supplementary alignments: >0 = the join is physically spanned.
#  3. TTAGGG/CCCTAA count in the junction interval (interstitial telomeric repeats).
#
# Uses the Hi-C and ONT BAMs made by analyses/circos/01.Circos_tracks.sh:
#   bash 04.hic_junction_check.sh            (needs samtools; ~10-30 min)
# Output: ${OUT}/junction_summary.tsv and ${OUT}/hic_crossing_profiles.tsv
# =============================================================================
set -uo pipefail
BASE_ROOT="/path/to/project"                     # as in analyses/circos/01.Circos_tracks.sh
OUT="${PWD}/HiC_check"; mkdir -p "$OUT"
THREADS=8; MAPQ=30; BIN=100000; WIN=10          # crossing window = WIN bins (1 Mb) each side

# species  chromosome  junction_start  junction_end  (from the breakpoint pipeline's breakpoints.tsv)
# NOTE: the BAMs predate the C. austriaca split, so SUPER_SUPER_8 / SUPER_SUPER_11 are in the OLD
# coordinates (old SUPER_SUPER_8 = the known misjoin: it is the POSITIVE CONTROL -- expect a deep dip).
# No re-mapping is needed or done.
JUNCTIONS="
Coronella_austriaca     SUPER_SUPER_8   29222882  29264395
Coronella_austriaca     SUPER_SUPER_11  26728100  26769343
Coronella_austriaca     SUPER_SUPER_11  29376301  30378971
Hierophis_viridiflavus  scaffold_9      26114596  31441977
Hierophis_viridiflavus  scaffold_10     20553495  20966600
Hierophis_viridiflavus  scaffold_6      107695825 108415758
Natrix_astreptophora    SUPER_SUPER_10  13019212  13532961
Zamenis_scalaris        SUPER_SUPER_8   360068    1031051
"

echo -e "species\tchrom\tj_start\tj_end\thic_crossing_at_junction\thic_crossing_chrom_median\tratio\tpercentile_along_chrom\tont_reads_both_sides\ttelomeric_motifs_in_interval\tinterval_kb" > ${OUT}/junction_summary.tsv
echo -e "species\tchrom\tpos\tcrossing" > ${OUT}/hic_crossing_profiles.tsv

while read -r SP CHR JS JE <&3; do
  [[ -z "$SP" ]] && continue
  D="${BASE_ROOT}/${SP}/res/07.Circos/${SP}"
  HIC="${D}/bam/${SP}.hic.markdup.bam"; ONT="${D}/bam/${SP}.ont.bam"
  FA="${BASE_ROOT}/${SP}/res/07.Circos/${SP}.masked.fa"
  [[ -s "$HIC" ]] || { echo "missing $HIC"; continue; }
  echo "[$(date +%T)] $SP $CHR"
  # 1. Hi-C crossing profile (intra-chromosomal pairs, counted once)
  # R1 and R2 were mapped separately (bwa -5SP, no mate fields): pair them by read name
  samtools view -@ $THREADS -q $MAPQ -F 0xD04 "$HIC" "$CHR" | cut -f1,4 \
   | sort -k1,1 -S 2G -T "$OUT" \
   | awk '$1==n { a=(p<$2)?p:$2; b=(p<$2)?$2:p; print a"\t"b; n=""; next } { n=$1; p=$2 }' \
   | awk -v B=$BIN -v W=$WIN -v sp=$SP -v ch=$CHR '
      { i=int($1/B); j=int($2/B); if (j>i && j-i<=2*W) M[i","j]++; if (j>mx) mx=j }
      END {
        for (x=W; x<=mx-W; x++) {
          s=0; for (i=x-W; i<x; i++) for (j=x; j<x+W; j++) if ((i","j) in M) s+=M[i","j]
          print sp"\t"ch"\t"x*B"\t"s
        } }' > ${OUT}/tmp_profile.tsv
  cat ${OUT}/tmp_profile.tsv >> ${OUT}/hic_crossing_profiles.tsv
  mid=$(( (JS+JE)/2 ))
  atj=$(awk -v m=$mid '{ d=($3>m)?$3-m:m-$3; if (NR==1||d<bd){bd=d; a=$4} } END{print a+0}' ${OUT}/tmp_profile.tsv)
  med=$(cut -f4 ${OUT}/tmp_profile.tsv | sort -n | awk '{v[NR]=$1} END{print v[int((NR+1)/2)]+0}')
  pct=$(awk -v a=$atj '{ if ($4<a) b++ } END{ printf "%.1f", 100*b/NR }' ${OUT}/tmp_profile.tsv)
  ratio=$(awk -v a=$atj -v m=$med 'BEGIN{ if (m>0) printf "%.3f", a/m; else print "NA" }')
  # 2. ONT reads aligned within 50 kb on both sides of the junction
  ont=NA
  if [[ -s "$ONT" ]]; then
    L=$(( JS>50000 ? JS-50000 : 0 ))
    comm -12 <(samtools view -q 20 -F 0x104 "$ONT" "${CHR}:${L}-${JS}" | cut -f1 | sort -u) \
             <(samtools view -q 20 -F 0x104 "$ONT" "${CHR}:${JE}-$((JE+50000))" | cut -f1 | sort -u) > ${OUT}/tmp_ont.txt
    ont=$(wc -l < ${OUT}/tmp_ont.txt)
  fi
  # 3. telomeric motifs inside the junction interval
  tel=NA
  [[ -s "$FA" ]] && tel=$(samtools faidx "$FA" "${CHR}:${JS}-${JE}" | tail -n +2 | tr -d '\n' | grep -o -i -E "TTAGGG|CCCTAA" | wc -l)
  echo -e "${SP}\t${CHR}\t${JS}\t${JE}\t${atj}\t${med}\t${ratio}\t${pct}\t${ont}\t${tel}\t$(( (JE-JS)/1000 ))" >> ${OUT}/junction_summary.tsv
done 3< <(echo "$JUNCTIONS")
rm -f ${OUT}/tmp_profile.tsv ${OUT}/tmp_ont.txt
column -t ${OUT}/junction_summary.tsv 2>/dev/null || cat ${OUT}/junction_summary.tsv
echo "Results in ${OUT}/"
