#!/bin/bash
# Heterozygosity, step 1: all-sites VCF per species and autosome with GATK v4.1.3.0 -- HaplotypeCaller
# (-ERC BP_RESOLUTION) on the filtered BAMs of analyses/psmc/01.Mapping.sh, then GenotypeGVCFs
# --include-non-variant-sites, so that confidently called homozygous positions enter the denominator.
set -uo pipefail

WD="/path/to/project/PSMCs"
GATK="/path/to/software/gatk4.1.3/gatk-4.1.3.0/gatk"
REF_DIR="${WD}/data/genomes"                        # <sp>.fa
BAM_DIR="${WD}/res/03.2.polished"                   # <sp>_mkdup_30.bam
BAM_SUFFIX="_mkdup_30.bam"
BED_DIR="${WD}/v2/data/autosome_beds"               # <sp>_autosomes.bed
SAMPLES="${WD}/data/samples_All.txt"                # sp min_dp max_dp mu g
OUTBASE="${WD}/res/GATK_ALLSITES"
N=10                                                # parallel contig jobs
CLEAN_GVCF=true                                     # delete gVCF after genotyping

REF_OVERRIDES="${WD}/data/reference_overrides.tsv"
resolve_ref() {
  local sp="$1" p=""
  if [[ -s "$REF_OVERRIDES" ]]; then
    p=$(awk -F'\t' -v s="$sp" '$1==s{print $2; exit}' "$REF_OVERRIDES")
  fi
  if [[ -n "$p" ]]; then printf '%s\n' "$p"; else printf '%s\n' "${REF_DIR}/${sp}.fa"; fi
}

MODE=run
SPP_ARG=""
for a in "$@"; do
  case "$a" in
    --check) MODE=check ;;
    -*) echo "unknown option: $a"; exit 1 ;;
    *)  SPP_ARG="$SPP_ARG $a" ;;
  esac
done

LOGBASE="${OUTBASE}/logs"; TMPBASE="${OUTBASE}/tmp"
mkdir -p "$OUTBASE" "$LOGBASE" "$TMPBASE"

if [[ -n "$SPP_ARG" ]]; then
  SPP="$SPP_ARG"
else
  SPP=$(awk '$2 != "x" && $3 != "x" && NF>=3 {print $1}' "$SAMPLES" | sort)
fi

echo "==================== PREFLIGHT ===================="
printf "%-28s %-4s %-4s %-4s %-4s %-5s %-6s %s\n" SPECIES fa bam bai bed dict contigs DP
MISSING=0
for sp in $SPP; do
  REF=$(resolve_ref "$sp"); BAM="${BAM_DIR}/${sp}${BAM_SUFFIX}"
  BED="${BED_DIR}/${sp}_autosomes.bed"; DICT="${REF%.fa}.dict"
  f=--; b=--; i=--; d=--; dc=--; nc=0; dp="--"
  [[ -s "$REF" ]] && f=ok
  [[ -s "$BAM" ]] && b=ok
  [[ -s "${BAM}.bai" ]] && i=ok
  [[ -s "$BED" ]] && { d=ok; nc=$(cut -f1 "$BED" | sort -u | wc -l); }
  [[ -s "$DICT" ]] && dc=ok
  read -r mn mx < <(awk -v s="$sp" '$1==s{print $2, $3}' "$SAMPLES")
  [[ -n "${mn:-}" && "$mn" != "x" ]] && dp="${mn}-${mx}"
  [[ "$f$b$d" == "okokok" ]] || MISSING=$((MISSING+1))
  printf "%-28s %-4s %-4s %-4s %-4s %-5s %-6s %s\n" "$sp" "$f" "$b" "$i" "$d" "$dc" "$nc" "$dp"
done
echo "---------------------------------------------------"
[[ -x "$GATK" ]] && echo "GATK: ok" || { echo "GATK: MISSING at $GATK"; MISSING=$((MISSING+1)); }
echo "Disk on $(df -h "$WD" | tail -1 | awk '{print $6": "$4" free of "$2}')"
echo "NOTE: BP_RESOLUTION gVCFs + all-sites VCFs are large."
echo "      Budget roughly 15-40 GB per species for final VCFs, plus gVCF"
echo "      intermediates transiently (removed when CLEAN_GVCF=true)."
echo "==================================================="
[[ "$MISSING" -gt 0 ]] && echo "WARNING: $MISSING species/tool missing inputs (--)."
if [[ "$MODE" == check ]]; then echo "Preflight only. Nothing run."; exit 0; fi

JOBLIST="${OUTBASE}/jobs.tsv"; : > "$JOBLIST"
MAPFILE="${OUTBASE}/contig_name_map.tsv"
[[ -s "$MAPFILE" ]] || printf "species\tsafe_name\treal_contig\n" > "$MAPFILE"

for sp in $SPP; do
  REF=$(resolve_ref "$sp"); BAM="${BAM_DIR}/${sp}${BAM_SUFFIX}"
  BED="${BED_DIR}/${sp}_autosomes.bed"
  for need in "$REF" "$BAM" "$BED"; do
    [[ -s "$need" ]] || { echo "[SKIP] $sp -- missing $need"; continue 2; }
  done

  OUTDIR="${OUTBASE}/${sp}"; TMPDIR="${TMPBASE}/${sp}"; ILDIR="${OUTDIR}/intervals"
  mkdir -p "$OUTDIR" "$TMPDIR" "$ILDIR"

  [[ -s "${BAM}.bai" ]] || { echo "[INFO] indexing BAM $sp"; samtools index "$BAM"; }
  [[ -s "${REF}.fai" ]] || { echo "[INFO] faidx $sp"; samtools faidx "$REF"; }
  DICT="${REF%.fa}.dict"
  [[ -s "$DICT" ]] || { echo "[INFO] dict $sp"; "$GATK" CreateSequenceDictionary -R "$REF" -O "$DICT"; }

  while read -r ch; do
    [[ -z "$ch" ]] && continue
    safe=$(printf '%s' "$ch" | tr -c 'A-Za-z0-9._-' '_')
    ilist="${ILDIR}/${safe}.list"
    printf "%s\n" "$ch" > "$ilist"
    grep -qP "^${sp}\t${safe}\t" "$MAPFILE" 2>/dev/null || \
      printf "%s\t%s\t%s\n" "$sp" "$safe" "$ch" >> "$MAPFILE"
    final="${OUTDIR}/${sp}_${safe}.genotyped.allsites.vcf.gz"
    if [[ -s "$final" && -s "${final}.tbi" ]]; then
      continue
    fi
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$sp" "$safe" "$ilist" "$REF" "$BAM" "$OUTDIR" "$TMPDIR" >> "$JOBLIST"
  done < <(cut -f1 "$BED" | sort -u)
done

echo "Jobs queued: $(wc -l < "$JOBLIST")"
[[ -s "$JOBLIST" ]] || { echo "Nothing to do -- all contigs already have final VCFs."; exit 0; }

while IFS=$'\t' read -r sp safe ilist ref bam outdir tmpdir; do
(
  log="${LOGBASE}/${sp}_${safe}.log"
  gvcf="${outdir}/${sp}_${safe}.bp_resolution.g.vcf.gz"
  final="${outdir}/${sp}_${safe}.genotyped.allsites.vcf.gz"

  echo "[START] $sp $safe $(date +%F_%T)" | tee "$log"

  if [[ -s "$gvcf" && -s "${gvcf}.tbi" ]]; then
    echo "[SKIP] gVCF exists" | tee -a "$log"
  else
    "$GATK" --java-options "-Djava.io.tmpdir=${tmpdir}" HaplotypeCaller \
      -R "$ref" -I "$bam" -O "$gvcf" \
      -ERC BP_RESOLUTION -L "$ilist" \
      --native-pair-hmm-threads 1 >> "$log" 2>&1 || {
        echo "[ERROR] HaplotypeCaller failed $sp $safe" | tee -a "$log"; exit 1; }
    [[ -s "${gvcf}.tbi" ]] || tabix -p vcf "$gvcf" >> "$log" 2>&1
  fi

  "$GATK" --java-options "-Djava.io.tmpdir=${tmpdir}" GenotypeGVCFs \
    -R "$ref" -V "$gvcf" --include-non-variant-sites \
    -L "$ilist" -O "$final" >> "$log" 2>&1 || {
      echo "[ERROR] GenotypeGVCFs failed $sp $safe" | tee -a "$log"; exit 1; }
  [[ -s "${final}.tbi" ]] || tabix -p vcf "$final" >> "$log" 2>&1

  if [[ "$CLEAN_GVCF" == "true" && -s "$final" && -s "${final}.tbi" ]]; then
    rm -f "$gvcf" "${gvcf}.tbi" "${gvcf}.idx"
  fi
  echo "[DONE] $sp $safe $(date +%F_%T)" | tee -a "$log"
) &
  if [[ $(jobs -r -p | wc -l) -ge "$N" ]]; then wait -n; fi
done < "$JOBLIST"
wait

echo "=================================================="
echo "[DONE] all-sites VCFs in ${OUTBASE}/<species>/"
echo "Contig name mapping: ${MAPFILE}"
