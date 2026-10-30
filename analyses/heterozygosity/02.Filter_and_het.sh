#!/bin/bash
# Heterozygosity, step 2: filter the all-sites VCFs and compute heterozygosity.
#   Exclude RepeatMasker intervals, keep sites with mean-depth/3 < DP < 2 x mean depth (per-species bounds from
#   samples_All.txt) and MQ > 20 (sites without MQ, i.e. invariant sites, are kept), concatenate the autosomes,
#   and count heterozygous / callable sites in non-overlapping 100-kb windows and genome-wide.
set -uo pipefail

WD="/path/to/project/PSMCs"
VCF_BASE="${WD}/res/GATK_ALLSITES"                  # from 20.allsites_vcf.sh
BED_DIR="${WD}/v2/data/autosome_beds"
GFF_DIR="${WD}/data/RepeatGFF"                      # <sp>.genome.fa.out.gff
SAMPLES="${WD}/data/samples_All.txt"
OUTBASE="${WD}/res/ROH_HET"
WINDOW=100000                                       # matches Talavera/Mochales-Riano
MIN_MQ=20
N=4                                                 # parallel species
KEEP_TMP=false                                      # keep per-contig filtered VCFs

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
printf "%-28s %-6s %-7s %-6s %s\n" SPECIES vcfs repeatGFF bed "DP range"
MISSING=0
for sp in $SPP; do
  nvcf=$(ls "${VCF_BASE}/${sp}"/${sp}_*.genotyped.allsites.vcf.gz 2>/dev/null | wc -l)
  g=--; [[ -s "${GFF_DIR}/${sp}.genome.fa.out.gff" ]] && g=ok
  b=--; [[ -s "${BED_DIR}/${sp}_autosomes.bed" ]] && b=ok
  read -r mn mx < <(awk -v s="$sp" '$1==s{print $2, $3}' "$SAMPLES")
  dp="--"; [[ -n "${mn:-}" && "$mn" != "x" ]] && dp="${mn}-${mx}"
  [[ "$nvcf" -gt 0 && "$b" == ok && "$dp" != "--" ]] || MISSING=$((MISSING+1))
  printf "%-28s %-6s %-7s %-6s %s\n" "$sp" "$nvcf" "$g" "$b" "$dp"
done
echo "==================================================="
echo "repeatGFF '--' means no RepeatMasker GFF found -- that species would run"
echo "WITHOUT repeat exclusion, which is NOT comparable to the others. Upload"
echo "the GFFs before the real run."
[[ "$MISSING" -gt 0 ]] && echo "WARNING: $MISSING species missing required inputs."
if [[ "$MODE" == check ]]; then echo "Preflight only. Nothing run."; exit 0; fi

for sp in $SPP; do
(
  echo "=================================================="
  echo "[START] $sp $(date +%F_%T)"

  VCFDIR="${VCF_BASE}/${sp}"
  BED="${BED_DIR}/${sp}_autosomes.bed"
  GFF="${GFF_DIR}/${sp}.genome.fa.out.gff"
  OUT="${OUTBASE}/${sp}"
  mkdir -p "${OUT}/tmp" "${OUT}/results"

  read -r MIN_DP MAX_DP < <(awk -v s="$sp" '$1==s{print $2, $3}' "$SAMPLES")
  if [[ -z "${MIN_DP:-}" || "$MIN_DP" == "x" ]]; then
    echo "[ERROR] no depth thresholds for $sp in $SAMPLES"; exit 1
  fi
  echo "[INFO] DP bounds ${MIN_DP}-${MAX_DP}, MQ>${MIN_MQ}"

  REPEAT_BED="${OUT}/${sp}_repeats.bed"
  if [[ -s "$GFF" ]]; then
    awk 'BEGIN{OFS="\t"} $1 !~ /^#/ {print $1, $4-1, $5}' "$GFF" \
      | sort -k1,1 -k2,2n > "$REPEAT_BED"
    echo "[INFO] repeat intervals: $(wc -l < "$REPEAT_BED")"
  else
    echo "[WARNING] no repeat GFF for $sp -- continuing WITHOUT repeat exclusion"
    : > "$REPEAT_BED"
  fi

  LIST="${OUT}/tmp/vcf_list.txt"; : > "$LIST"
  while read -r ch; do
    [[ -z "$ch" ]] && continue
    safe=$(printf '%s' "$ch" | tr -c 'A-Za-z0-9._-' '_')
    vcf="${VCFDIR}/${sp}_${safe}.genotyped.allsites.vcf.gz"
    [[ -s "$vcf" ]] || { echo "[WARN] missing $vcf"; continue; }
    [[ -s "${vcf}.tbi" ]] || tabix -p vcf "$vcf"
    tmpvcf="${OUT}/tmp/${safe}.filtered.vcf.gz"
    if [[ ! -s "$tmpvcf" ]]; then
      echo "[INFO] filtering $safe"
      if [[ -s "$REPEAT_BED" ]]; then
        bcftools view -T ^"$REPEAT_BED" \
          -i "INFO/DP>${MIN_DP} && INFO/DP<${MAX_DP} && (INFO/MQ=\".\" || INFO/MQ>${MIN_MQ})" \
          -Oz -o "$tmpvcf" "$vcf"
      else
        bcftools view \
          -i "INFO/DP>${MIN_DP} && INFO/DP<${MAX_DP} && (INFO/MQ=\".\" || INFO/MQ>${MIN_MQ})" \
          -Oz -o "$tmpvcf" "$vcf"
      fi
      tabix -p vcf "$tmpvcf"
    fi
    echo "$tmpvcf" >> "$LIST"
  done < <(cut -f1 "$BED" | sort -u)

  [[ -s "$LIST" ]] || { echo "[ERROR] no filtered VCFs for $sp"; exit 1; }

  SAMPLE=$(bcftools query -l "$(head -1 "$LIST")" | head -1)
  FVCF="${OUT}/${sp}_${SAMPLE}.autosomes.filtered.vcf.gz"
  if [[ ! -s "$FVCF" ]]; then
    echo "[INFO] concatenating $(wc -l < "$LIST") contigs"
    bcftools concat -f "$LIST" -Oz -o "$FVCF"
    tabix -p vcf "$FVCF"
  fi

  bcftools query -f '[%GT\n]' "$FVCF" | sort | uniq -c \
    > "${OUT}/results/${sp}_${SAMPLE}.genotype_counts.txt"
  echo "[INFO] genotype counts:"; sed 's/^/    /' "${OUT}/results/${sp}_${SAMPLE}.genotype_counts.txt"


  echo "[INFO] heterozygosity in ${WINDOW} bp windows"
  HETOUT="${OUT}/results/${sp}_${SAMPLE}_heterozygosity_${WINDOW}.txt"
  bcftools query -f '%CHROM\t%POS\t[%GT]\n' "$FVCF" \
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

  awk -F'\t' 'NR>1{h+=$3; c+=$4} END{
      if(c>0) printf "%s\t%d\t%d\t%.8f\n", SP, c, h, h/c
      else    printf "%s\tNA\tNA\tNA\n", SP }' SP="$sp" "$HETOUT" \
    > "${OUT}/results/${sp}.heterozygosity_genomewide.tsv"
  echo "[INFO] genome-wide: $(cat "${OUT}/results/${sp}.heterozygosity_genomewide.tsv")"

  [[ "$KEEP_TMP" == "true" ]] || rm -f "${OUT}"/tmp/*.filtered.vcf.gz*

  echo "[DONE] $sp $(date +%F_%T)"
) &
  if [[ $(jobs -r -p | wc -l) -ge "$N" ]]; then wait -n; fi
done
wait

echo
echo "=============== GENOME-WIDE HETEROZYGOSITY ==============="
printf "species\tcallable_sites\theterozygous_sites\theterozygosity\n" \
  > "${OUTBASE}/ALL_heterozygosity.tsv"
cat "${OUTBASE}"/*/results/*.heterozygosity_genomewide.tsv 2>/dev/null \
  | sort -k4,4gr >> "${OUTBASE}/ALL_heterozygosity.tsv"
column -t "${OUTBASE}/ALL_heterozygosity.tsv"
echo
echo "Wrote ${OUTBASE}/ALL_heterozygosity.tsv"
echo "Next: filter RG records by length >=100kb and quality >70, then classify"
echo "      short (100-500kb) / medium (500kb-1Mb) / long (>1Mb)."
