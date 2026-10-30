#!/bin/bash
#SBATCH --job-name=macrosynteny_build_inputs
#SBATCH --account=<account>
#SBATCH --partition=<partition>
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=macrosynteny_build_inputs_%j.out
#SBATCH --error=macrosynteny_build_inputs_%j.err

# =============================================================================
# 01.build_inputs.sh -- per-species prep, run once before 02.flip_chromosomes.sh.
# Run 00.setup_env.sh first (one-time), before this.
#
# For each species in SPECIES_ORDER (00.config.sh), from:
#   GENOME_DIR/<Species>/<Species>.masked.fa
#   ANNOT_DIR/<Species>.gtf  (or .gff3 / .gff -- auto-detected, see below)
#   ANNOT_DIR/<Species>_proteins.fa
# produces, in OUT_DIR/<Species>/:
#   <Species>.fa.fai   <Species>.bed   <Species>.longest.bed   <Species>.longest.faa
# Idempotent -- skips any step whose output already exists and is non-empty.
#
# Two robustness features:
#  1. ANNOTATION FORMAT: gtf_to_bed.py / longest_isoform.py already parse
#     GTF-style (`transcript_id "..."`) and GFF3-style (`ID=...;Parent=...`)
#     attributes automatically, per line -- the content, not the extension,
#     decides how it's parsed. This script just needs to FIND the file
#     regardless of which extension it was given, so it tries .gtf, then
#     .gff3, then .gff.
#  2. SEQID MISMATCH: some assemblies were submitted to both ENA/GenBank and
#     NCBI RefSeq and end up with two different chromosome-accession systems
#     -- your genome FASTA may say ">OZ187402.1 ... chromosome: 1" while a
#     RefSeq annotation for that same assembly uses "NC_138548.1" in its
#     seqid column for that same chromosome. If that happens, the BED
#     produced from the annotation won't match anything in the genome's
#     .fai, and every downstream step would silently see 0 overlap for that
#     species. This script checks for that after building the BED: if a
#     mapping file ANNOT_DIR/<Species>.seqid_map.tsv exists (two columns:
#     annotation_id<TAB>genome_id -- see remap_seqids.py's header for how to
#     build one from NCBI's sequence report), it's applied automatically. If
#     no mapping exists AND seqid overlap with the genome is poor, the
#     species is skipped with a clear error rather than silently proceeding
#     with broken (all-empty) synteny for it.
# =============================================================================
set -euo pipefail
# Under sbatch, SLURM copies this script into a per-job spool directory
# (/var/spool/slurmd/jobNNNNNN/) and runs it FROM THERE -- so BASH_SOURCE-based
# directory detection resolves to that spool dir, not the real scripts folder,
# and sourcing sibling files like 00.config.sh fails ("No such file or
# directory") because only this one script was copied there, not its siblings.
# SLURM_SUBMIT_DIR (the directory `sbatch` was invoked FROM) is the reliable
# anchor under SLURM; fall back to BASH_SOURCE detection only when running
# this script directly (not via sbatch). Always `cd` into this scripts folder
# before running `sbatch <script>.sh` (relative path) so the two agree.
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${SCRIPT_DIR}/00.config.sh"

module purge
module load ${MODULE_MAMBA}
MRUN="mamba run -p ${JCVI_ENV_PREFIX}"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

find_annotation() {
  # Echoes the first existing ANNOT_DIR/<species>.<ext> among gtf/gff3/gff.
  # Returns non-zero if none exist.
  local sp="$1"
  local ext
  for ext in gtf gff3 gff; do
    local cand="${ANNOT_DIR}/${sp}.${ext}"
    if [[ -f "$cand" ]]; then
      echo "$cand"
      return 0
    fi
  done
  return 1
}

for SPECIES in "${SPECIES_ORDER[@]}"; do
(
  SP_OUT="${OUT_DIR}/${SPECIES}"
  mkdir -p "$SP_OUT"

  MASKED="${GENOME_DIR}/${SPECIES}/${SPECIES}.masked.fa"
  PROTEINS="${ANNOT_DIR}/${SPECIES}_proteins.fa"

  if ! GTF=$(find_annotation "$SPECIES"); then
    echo "[${SPECIES}] SKIPPING -- no annotation file found (tried ${ANNOT_DIR}/${SPECIES}.{gtf,gff3,gff})" >&2
    exit 0
  fi

  for f in "$MASKED" "$PROTEINS"; do
    if [[ ! -f "$f" ]]; then
      echo "[${SPECIES}] SKIPPING -- missing input file: $f" >&2
      exit 0
    fi
  done

  log "[${SPECIES}] building inputs (annotation: $(basename "$GTF"))..."

  FAI="${SP_OUT}/${SPECIES}.fa.fai"
  if [[ ! -s "$FAI" ]]; then
    if [[ ! -s "${MASKED}.fai" ]]; then
      ${MRUN} samtools faidx "$MASKED"
    fi
    cp "${MASKED}.fai" "$FAI"
  fi

  BED_RAW="${SP_OUT}/${SPECIES}.raw.bed"
  if [[ ! -s "$BED_RAW" ]]; then
    ${MRUN} python3 "${SCRIPT_DIR}/gtf_to_bed.py" "$GTF" "$BED_RAW"
  fi

  BED="${SP_OUT}/${SPECIES}.bed"
  SEQID_MAP="${ANNOT_DIR}/${SPECIES}.seqid_map.tsv"
  if [[ -s "$SEQID_MAP" ]]; then
    if [[ ! -s "$BED" ]]; then
      ${MRUN} python3 "${SCRIPT_DIR}/remap_seqids.py" "$BED_RAW" "$SEQID_MAP" "$BED"
    fi
  else
    cp "$BED_RAW" "$BED"
  fi

  # ---- seqid consistency check: does the (possibly remapped) BED actually
  # overlap the genome's .fai? Catches the RefSeq-vs-ENA/GenBank accession
  # mismatch (or any other annotation/genome naming mismatch) before it
  # silently produces empty synteny for this species downstream.
  #
  # Weighted by GENE (BED row) count, not distinct-scaffold count: an
  # assembly can have a couple dozen real chromosomes plus hundreds of tiny
  # unplaced scaffolds with a gene or two each. A seqid_map.tsv covering
  # just the named chromosomes (the common case -- see remap_seqids.py) can
  # legitimately leave most of those hundreds of unplaced-scaffold IDs
  # unmapped while still covering the vast majority of actual genes. Scoring
  # by distinct-scaffold overlap would flag that as a failure; scoring by
  # gene-row overlap reflects what actually matters for the synteny call.
  if [[ -s "$BED" && -s "$FAI" ]]; then
    n_bed_rows=$(wc -l < "$BED")
    fai_ids="${SP_OUT}/.fai_ids.tmp"
    cut -f1 "$FAI" | sort -u > "$fai_ids"
    n_overlap_rows=$(cut -f1 "$BED" | sort | join -t $'\t' - "$fai_ids" | wc -l)
    rm -f "$fai_ids"
    if [[ "$n_bed_rows" -gt 0 ]]; then
      pct=$(( 100 * n_overlap_rows / n_bed_rows ))
      if [[ "$pct" -lt 50 ]]; then
        echo "[${SPECIES}] ERROR: only ${pct}% of annotated genes sit on a seqid " \
             "present in ${FAI} -- annotation and genome likely use different " \
             "chromosome accessions (e.g. RefSeq NC_* in the GTF vs ENA/GenBank " \
             "OZ*-style in the genome FASTA -- check for a " \
             "'#!genome-build-accession NCBI_Assembly:GCF_...' line near the top " \
             "of $GTF). Provide ${SEQID_MAP} (annotation_id<TAB>genome_id, see " \
             "remap_seqids.py) and rerun. Skipping ${SPECIES} for now." >&2
        exit 0
      fi
    fi
  fi

  LONGEST_BED="${SP_OUT}/${SPECIES}.longest.bed"
  LONGEST_FAA="${SP_OUT}/${SPECIES}.longest.faa"
  if [[ ! -s "$LONGEST_BED" || ! -s "$LONGEST_FAA" ]]; then
    ${MRUN} python3 "${SCRIPT_DIR}/longest_isoform.py" \
      "$GTF" "$PROTEINS" "$BED" "$LONGEST_FAA" "$LONGEST_BED"
  fi

  log "[${SPECIES}] done."
) &
  if [[ $(jobs -r -p | wc -l) -ge $N_PARALLEL ]]; then
    wait -n
  fi
done
wait

log "All species processed. Outputs under: ${OUT_DIR}/<Species>/"
