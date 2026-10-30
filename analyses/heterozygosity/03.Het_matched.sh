#!/bin/bash
# Heterozygosity, step 3: additional filters matching Burriel-Carranza et al. (2024), applied to the filtered VCF of
# step 2: bcftools norm against the reference, bcftools filter --SnpGap 10, indels/MNPs removed (SNPs and invariant
# sites kept); heterozygosity recomputed per 100-kb window and genome-wide.
set -uo pipefail

WD="/path/to/project/PSMCs"
ROH_HET="${WD}/res/ROH_HET"                # from 21.roh_and_het.sh
REF_DIR="${WD}/data/genomes"               # <sp>.fa
SAMPLES="${WD}/data/samples_All.txt"
OUTBASE="${WD}/res/HET_MATCHED"
WINDOW=100000
N=4                                        # parallel species
KEEP_VCF=false                             # keep the re-filtered VCF

MODE=run
SPP_ARG=""
for a in "$@"; do
  case "$a" in
    --check) MODE=check ;;
    -*) echo "unknown option: $a"; exit 1 ;;
    *)  SPP_ARG="$SPP_ARG $a" ;;
  esac
done

mkdir -p "$OUTBASE"
if [[ -n "$SPP_ARG" ]]; then
  SPP="$SPP_ARG"
else
  SPP=$(awk '$2 != "x" && $3 != "x" && NF>=3 {print $1}' "$SAMPLES" | sort)
fi

echo "==================== PREFLIGHT ===================="
printf "%-28s %-6s %-5s %s\n" SPECIES vcf ref "old heterozygosity"
MISSING=0
for sp in $SPP; do
  v=$(ls "${ROH_HET}/${sp}"/${sp}_*.autosomes.filtered.vcf.gz 2>/dev/null | head -1)
  vs=--; [[ -n "$v" && -s "$v" ]] && vs=ok
  rs=--; [[ -s "${REF_DIR}/${sp}.fa" ]] && rs=ok
  old=$(awk -F'\t' '{print $4}' "${ROH_HET}/${sp}/results/${sp}.heterozygosity_genomewide.tsv" 2>/dev/null)
  [[ "$vs$rs" == okok ]] || MISSING=$((MISSING+1))
  printf "%-28s %-6s %-5s %s\n" "$sp" "$vs" "$rs" "${old:---}"
done
echo "---------------------------------------------------"
for t in bcftools tabix; do
  command -v "$t" >/dev/null && echo "$t: ok" || { echo "$t: MISSING"; MISSING=$((MISSING+1)); }
done
echo "==================================================="
[[ "$MISSING" -gt 0 ]] && echo "WARNING: $MISSING species/tool missing inputs (--)."
if [[ "$MODE" == check ]]; then echo "Preflight only. Nothing run."; exit 0; fi

for sp in $SPP; do
(
  REF="${REF_DIR}/${sp}.fa"
  FVCF=$(ls "${ROH_HET}/${sp}"/${sp}_*.autosomes.filtered.vcf.gz 2>/dev/null | head -1)
  [[ -n "$FVCF" && -s "$FVCF" && -s "$REF" ]] || { echo "[SKIP] $sp -- missing inputs"; exit 0; }

  OUT="${OUTBASE}/${sp}"; mkdir -p "$OUT"
  LOG="${OUT}/${sp}.log"
  echo "[START] $sp $(date +%F_%T)" | tee "$LOG"

  SAMPLE=$(bcftools query -l "$FVCF" | head -1)
  MVCF="${OUT}/${sp}_${SAMPLE}.matched.vcf.gz"

  if [[ ! -s "$MVCF" ]]; then
    echo "[INFO] norm + SnpGap 10 + indel removal" | tee -a "$LOG"
    bcftools norm -f "$REF" --check-ref w -Ou "$FVCF" 2>>"$LOG" \
      | bcftools filter --SnpGap 10 -Ou 2>>"$LOG" \
      | bcftools view -f PASS,. -V indels,mnps,other -Oz -o "$MVCF" 2>>"$LOG" || {
          echo "[ERROR] filtering failed for $sp" | tee -a "$LOG"; exit 1; }
    tabix -p vcf "$MVCF" 2>>"$LOG"
  fi

  n_before=$(bcftools index -n "$FVCF")
  n_after=$(bcftools index -n "$MVCF")
  echo "[INFO] records ${n_before} -> ${n_after}" | tee -a "$LOG"

  HETOUT="${OUT}/${sp}_${SAMPLE}_heterozygosity_${WINDOW}.matched.txt"
  bcftools query -f '%CHROM\t%POS\t[%GT]\n' "$MVCF" \
  | awk -v W="$WINDOW" 'BEGIN{OFS="\t"}
      { gt=$3; gsub(/\|/,"/",gt)
        if (gt=="0/0"||gt=="0/1"||gt=="1/0"||gt=="1/1") {
          w=int(($2-1)/W); k=$1 SUBSEP w
          cal[k]++
          if (gt=="0/1"||gt=="1/0") het[k]++
        } }
      END { for (k in cal) { split(k,a,SUBSEP)
              print a[1], a[2]*W + W/2, het[k]+0, cal[k], (het[k]+0)/cal[k] } }' \
  | sort -k1,1V -k2,2n \
  | cat <(printf "chrom\tmidpoint\theterozygous_sites\tcallable_sites\theterozygosity\n") - \
  > "$HETOUT"

  awk -F'\t' -v SP="$sp" -v NB="$n_before" -v NA="$n_after" 'NR>1{
        h+=$3; c+=$4;
        if ($4>0) { s+=$5; n++ }
      } END {
        if (c>0) printf "%s\t%d\t%d\t%d\t%d\t%.8f\t%.8f\t%d\n", SP, NB, NA, c, h, h/c, (n?s/n:0), n;
        else     printf "%s\t%d\t%d\tNA\tNA\tNA\tNA\t0\n", SP, NB, NA }' "$HETOUT" \
    > "${OUT}/${sp}.heterozygosity_matched.tsv"

  echo "[INFO] $(cat "${OUT}/${sp}.heterozygosity_matched.tsv")" | tee -a "$LOG"
  [[ "$KEEP_VCF" == "true" ]] || rm -f "$MVCF" "${MVCF}.tbi"
  echo "[DONE] $sp $(date +%F_%T)" | tee -a "$LOG"
) &
  if [[ $(jobs -r -p | wc -l) -ge "$N" ]]; then wait -n; fi
done
wait

SUM="${OUTBASE}/ALL_heterozygosity_matched.tsv"
printf "species\trecords_before\trecords_after\tcallable\thet_sites\thet_pooled\thet_window_mean\tn_windows\n" > "$SUM"
cat "${OUTBASE}"/*/*.heterozygosity_matched.tsv 2>/dev/null | sort -k6,6gr >> "$SUM"

CMP="${OUTBASE}/COMPARISON_old_vs_matched.tsv"
printf "species\told_pooled\tmatched_pooled\tpct_change\tmatched_window_mean\n" > "$CMP"
while IFS=$'\t' read -r sp nb na c h pooled wmean nw; do
  [[ "$sp" == "species" ]] && continue
  old=$(awk -F'\t' '{print $4}' "${ROH_HET}/${sp}/results/${sp}.heterozygosity_genomewide.tsv" 2>/dev/null)
  [[ -z "$old" || "$old" == "NA" ]] && old=NA
  pc=$(awk -v o="$old" -v n="$pooled" 'BEGIN{ if(o=="NA"||o+0==0) print "NA"; else printf "%+.2f", 100*(n-o)/o }')
  printf "%s\t%s\t%s\t%s\t%s\n" "$sp" "$old" "$pooled" "$pc" "$wmean" >> "$CMP"
done < "$SUM"

echo
echo "=========== OLD vs TALAVERA-MATCHED FILTERING ==========="
column -t "$CMP"
echo
echo "Wrote $SUM"
echo "Wrote $CMP"
echo
echo "Note: this matches their VARIANT FILTERING. Two differences remain by"
echo "design and should be stated in the methods rather than removed:"
echo "  * we exclude RepeatMasker intervals, they do not"
echo "  * we exclude sex-linked scaffolds via the autosome BED, they do not"
echo "  * our reads are mapped to each species' OWN assembly; their 24 samples"
echo "    are all mapped to the A. caudivolvulus reference"
