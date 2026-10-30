#!/bin/bash
#SBATCH --job-name=macrosynteny_synteny
#SBATCH --account=<account>
#SBATCH --partition=<partition>
#SBATCH --cpus-per-task=32
#SBATCH --mem=60G
#SBATCH --time=12:00:00
#SBATCH --output=macrosynteny_synteny_%j.out
#SBATCH --error=macrosynteny_synteny_%j.err

# =============================================================================
# 03.run_synteny.sh -- the core JCVI/MCscan synteny-calling step.
#
# Two kinds of comparisons:
#   1. STAIRCASE: every species vs its immediate neighbor in SPECIES_ORDER
#      (row i vs row i+1), built automatically from SPECIES_ORDER.
#   2. W-LINKS (CHAIN, not all-pairs): walk SPECIES_ORDER top to bottom;
#      each time a W-bearing species is found, link it ONLY to the NEXT
#      W-bearing species further down the order (extracting just their named
#      W scaffolds via extract_scaffold_subset.py and running diamond+jcvi on
#      that tiny subset), then keep walking from there. This mirrors the
#      staircase's own philosophy -- non-adjacent pairs are connected only
#      implicitly, through the chain, not with a direct edge -- so a species'
#      W is linked to at most 2 neighbors (previous and next W-bearer), NOT
#      to every other W-bearing species. (Earlier version of this script
#      linked EVERY pair of W-bearing species -- 28 edges for 8 females --
#      which was inconsistent with the staircase's own adjacent-only design
#      and, more importantly, unreadable in the plot. Fixed 2026-08-19.)
#
# SCAFFOLD ORDER WITHIN EACH ROW: seqids.txt now places each species' Z
# scaffold(s) (where known) and W scaffold(s) (where known) at the END of
# that species' row, after all autosomes -- instead of wherever they'd fall
# under plain sort -V. Added 2026-08-19.
#
# Produces, under OUT_DIR/:
#   pairs/<query>.<subject>.blast / .anchors / .anchors.simple      (staircase)
#   pairs/<query>.<subject>.W.blast / .anchors / .anchors.simple    (W-links)
#   seqids.txt   layoutY.txt   layoutEdges.txt   layout.txt
#
# Idempotent per pair (skips existing non-empty outputs).
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

log() { echo "[$(date '+%H:%M:%S')] $*" >&2; }

mkdir -p "${OUT_DIR}/pairs" "${OUT_DIR}/logs"

# Z-column alignment bookkeeping (see build_plot_seqids.py's docstring for
# the full explanation). z_alignment.tsv is appended to, once per species,
# by build_plot_seqids.py further down -- must start empty each run or a
# stale line from a previous run would corrupt the max() used below.
rm -f "${OUT_DIR}/z_alignment.tsv" "${OUT_DIR}/z_gap_deltas.tsv"

bed_for() {
  local sp="$1"
  local flipped="${OUT_DIR}/${sp}/${sp}.longest.flipped.bed"
  local plain="${OUT_DIR}/${sp}/${sp}.longest.bed"
  if [[ -s "$flipped" ]]; then echo "$flipped"; else echo "$plain"; fi
}
faa_for() { echo "${OUT_DIR}/$1/$1.longest.faa"; }

# run_pair SUBJECT QUERY SUBJECT_BED SUBJECT_FAA QUERY_BED QUERY_FAA LABEL
# Generic diamond+jcvi runner -- LABEL distinguishes staircase pairs from
# W-only pairs in the output filenames (so a full and a W-restricted
# comparison between the same two species never collide).
run_pair() {
  local subject="$1" query="$2" subject_bed="$3" subject_faa="$4" query_bed="$5" query_faa="$6" label="$7"
  local logf="${OUT_DIR}/logs/${query}.${subject}.${label}.log"
  # W-only comparisons are restricted to a single (often fragmented, low-
  # gene-count) chromosome, not a whole genome -- the global MINSPAN=50
  # tuned for full staircase comparisons is an unreasonable bar there and
  # silently zeroes out real W-W synteny for species whose W is split across
  # many small unlocalized contigs. Use a separate, much lower threshold for
  # label=="W" (MINSPAN_W, from 00.config.sh) instead.
  local minspan="$MINSPAN"
  # ${MINSPAN_W:-5} not $MINSPAN_W -- this script runs under `set -u`, so a
  # bare reference would abort the ENTIRE run with "unbound variable" if
  # 00.config.sh does not define it. Falls back to 5.
  [[ "$label" == "W" ]] && minspan="${MINSPAN_W:-5}"

  if [[ ! -f "$subject_faa" || ! -f "$query_faa" ]]; then
    log "[${query} vs ${subject} (${label})] SKIPPING -- missing protein FASTA"
    return 0
  fi
  if [[ ! -s "$subject_faa" || ! -s "$query_faa" ]]; then
    log "[${query} vs ${subject} (${label})] SKIPPING -- empty protein FASTA (e.g. W scaffold not found)"
    return 0
  fi

  local db="${OUT_DIR}/pairs/${subject}.${label}"
  local blast="${OUT_DIR}/pairs/${query}.${subject}.${label}.blast"
  local anchors="${OUT_DIR}/pairs/${query}.${subject}.${label}.anchors"
  local anchors_simple="${OUT_DIR}/pairs/${query}.${subject}.${label}.anchors.simple"

  {
    if [[ ! -s "${db}.dmnd" ]]; then
      ${MRUN} diamond makedb --in "$subject_faa" -d "$db"
    fi
    if [[ ! -s "$blast" ]]; then
      # jcvi.compara.synteny scan parses this file with jcvi's BlastLine,
      # which requires the FULL 12-column blast tabular format (qseqid
      # sseqid pident length mismatch gapopen qstart qend sstart send
      # evalue bitscore -- see jcvi/formats/pyblast.py's BlastLine.__init__,
      # which does int(args[4]) for mismatch, int(args[5]) for gapopen, and
      # int(args[6:10]) for the qstart/qend/sstart/send coordinates the
      # synteny scan itself needs). A 6-column outfmt (as an earlier version
      # of this script used) would hand BlastLine an e-value string where it
      # expects an integer mismatch count and crash on every single pair.
      ${MRUN} diamond blastp -q "$query_faa" -d "$db" \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        --evalue "$DIAMOND_EVALUE" --max-target-seqs "$DIAMOND_MAX_TARGETS" \
        --threads "$THREADS" > "$blast"
    fi
    if [[ ! -s "$anchors" ]]; then
      ${MRUN} python -m jcvi.compara.synteny scan \
        "$blast" "$anchors" --qbed="$query_bed" --sbed="$subject_bed" --no_strip_names
    fi
    if [[ ! -s "$anchors_simple" ]]; then
      ${MRUN} python -m jcvi.compara.synteny screen \
        "$anchors" "$anchors_simple" --qbed="$query_bed" --sbed="$subject_bed" \
        --simple --minspan="$minspan"
    fi
  } > "$logf" 2>&1 || { log "[${query} vs ${subject} (${label})] FAILED -- see $logf"; return 0; }

  if [[ ! -s "$anchors_simple" ]]; then
    log "[${query} vs ${subject} (${label})] WARNING: empty anchors.simple -- see $logf"
    return 0
  fi
  echo "$anchors_simple"
}

N_SPECIES=${#SPECIES_ORDER[@]}

# ---- 1. Staircase pairs (row i vs row i+1), full proteomes -----------------
declare -a EDGES   # each entry: "i,j,anchors_simple_path"
for ((i = 0; i < N_SPECIES - 1; i++)); do
  subject="${SPECIES_ORDER[$i]}"
  query="${SPECIES_ORDER[$((i + 1))]}"
  log "Staircase pair: ${query} vs ${subject} (rows $((i+1)),$i)"
  anchors_simple=$(run_pair "$subject" "$query" \
    "$(bed_for "$subject")" "$(faa_for "$subject")" \
    "$(bed_for "$query")" "$(faa_for "$query")" "full")
  if [[ -n "$anchors_simple" ]]; then
    EDGES+=("${i},$((i + 1)),${anchors_simple}")
  fi
done

# ---- 2. Load W and Z scaffold assignments (used for BOTH the W-chain below
#         AND the seqids.txt scaffold ordering in section 4) ----------------
declare -A W_SCAFFOLD Z_SCAFFOLD

read_sex_chrom_file() {
  local file="$1"
  local -n out_map="$2"   # nameref to the associative array to fill
  if [[ ! -f "$file" ]]; then return 0; fi
  while IFS=$'\t' read -r sp scaf; do
    if [[ -z "${sp:-}" || "$sp" =~ ^# ]]; then continue; fi
    if [[ -z "${scaf:-}" ]]; then continue; fi
    out_map["$sp"]="$scaf"
    fai="${OUT_DIR}/${sp}/${sp}.fa.fai"
    # scaf may be a comma-separated list (main sex chromosome + unlocalized
    # associated contigs) -- check each one individually against the .fai.
    if [[ -f "$fai" ]]; then
      IFS=',' read -r -a _scaf_list <<< "$scaf"
      for _s in "${_scaf_list[@]}"; do
        _s="${_s## }"; _s="${_s%% }"
        if ! grep -q -w "$_s" "$fai"; then
          log "WARNING: ${file} lists '${_s}' for ${sp}, but it's not in ${fai} -- typo, or not yet built?"
        fi
      done
    fi
  done < "$file"
}
read_sex_chrom_file "$W_CHROM_FILE" W_SCAFFOLD
read_sex_chrom_file "$Z_CHROM_FILE" Z_SCAFFOLD

# Species with a W, in SPECIES_ORDER order (not w_chromosome.txt's file order)
declare -a W_SPECIES
for sp in "${SPECIES_ORDER[@]}"; do
  [[ -n "${W_SCAFFOLD[$sp]:-}" ]] && W_SPECIES+=("$sp")
done

# ---- 3. W-chromosome CHAIN: link each W-bearing species only to the NEXT
#         W-bearing species further down SPECIES_ORDER (not all pairs) ------
n_w=${#W_SPECIES[@]}
if [[ $n_w -ge 2 ]]; then
  log "W-chromosome species (${n_w}, in phylogenetic order): ${W_SPECIES[*]} -- chaining consecutive pairs only"

  # Build the W-only bed/faa subset once per species (reused across the chain).
  declare -A W_BED W_FAA
  for sp in "${W_SPECIES[@]}"; do
    scaf="${W_SCAFFOLD[$sp]}"
    src_bed=$(bed_for "$sp")
    src_faa=$(faa_for "$sp")
    w_bed="${OUT_DIR}/${sp}/${sp}.W.bed"
    w_faa="${OUT_DIR}/${sp}/${sp}.W.faa"
    if [[ ! -s "$w_bed" || ! -s "$w_faa" ]]; then
      ${MRUN} python3 "${SCRIPT_DIR}/extract_scaffold_subset.py" "$src_bed" "$src_faa" "$scaf" "$w_bed" "$w_faa"
    fi
    W_BED["$sp"]="$w_bed"
    W_FAA["$sp"]="$w_faa"
  done

  for ((a = 0; a < n_w - 1; a++)); do
    sp_a="${W_SPECIES[$a]}"; sp_b="${W_SPECIES[$((a + 1))]}"
    idx_a=-1; idx_b=-1
    for ((k = 0; k < N_SPECIES; k++)); do
      [[ "${SPECIES_ORDER[$k]}" == "$sp_a" ]] && idx_a=$k
      [[ "${SPECIES_ORDER[$k]}" == "$sp_b" ]] && idx_b=$k
    done
    if [[ $idx_a -eq -1 || $idx_b -eq -1 ]]; then
      log "WARNING: ${sp_a}/${sp_b} from w_chromosome.txt not found in SPECIES_ORDER -- skipping"
      continue
    fi
    lo=$((idx_a < idx_b ? idx_a : idx_b))
    hi=$((idx_a < idx_b ? idx_b : idx_a))
    subject="${SPECIES_ORDER[$lo]}"; query="${SPECIES_ORDER[$hi]}"
    log "W-chain link: ${query} vs ${subject} (rows $hi,$lo) -- consecutive W-bearing pair, W-scaffold-only comparison"
    anchors_simple=$(run_pair "$subject" "$query" \
      "${W_BED[$subject]}" "${W_FAA[$subject]}" \
      "${W_BED[$query]}" "${W_FAA[$query]}" "W")
    if [[ -n "$anchors_simple" ]]; then
      EDGES+=("${lo},${hi},${anchors_simple}")
    else
      log "[${query} vs ${subject}] no W-W synteny block found (or comparison failed) -- no edge added"
    fi
  done
elif [[ $n_w -eq 1 ]]; then
  log "Only one species in w_chromosome.txt (${W_SPECIES[0]}) -- nothing to link yet."
else
  log "w_chromosome.txt has no entries yet -- no W-links added (list them in w_chromosome.txt)."
fi

# ---- 3. layoutEdges.txt ----------------------------------------------------
EDGES_FILE="${OUT_DIR}/layoutEdges.txt"
echo "# edges" > "$EDGES_FILE"
for e in "${EDGES[@]}"; do
  IFS=',' read -r ei ej epath <<< "$e"
  printf "e, %d, %d, %s\n" "$ei" "$ej" "$epath" >> "$EDGES_FILE"
done

# ---- 4. seqids.txt + layoutY.txt (based on FULL genome beds, not W-only) --
SEQIDS_FILE="${OUT_DIR}/seqids.txt"
LAYOUTY_FILE="${OUT_DIR}/layoutY.txt"
> "$SEQIDS_FILE"
echo "# y, xstart, xend, rotation, color, label, va, bed" > "$LAYOUTY_FILE"

y_start=0.90
y_step=$(LC_NUMERIC=C awk -v n="$N_SPECIES" 'BEGIN{ s=0.80/(n>1?n-1:1); if (s>0.07) s=0.07; printf "%.4f", s }')
# xstart widened from 0.10 -- the species label is drawn at
# (xstart - vpad/2) in karyotype.py, and at 0.10 it was crowding into the
# chromosome bars. 0.18 leaves clear room for the label text.
x_start=0.18
x_end=0.80
rot=0
color="#000000"
va="top"

# species<TAB>N line lookup from CHROM_COUNT_FILE (see 00.config.sh). A
# species with an entry gets a PLOT-ONLY bed (build_plot_seqids.py), renamed
# to clean "1".."N"/"Z"/"W" labels and restricted to its N largest
# scaffolds; a species without an entry keeps today's behaviour untouched.
declare -A CHROM_COUNT
if [[ -f "$CHROM_COUNT_FILE" ]]; then
  while IFS=$'\t' read -r sp n; do
    if [[ -z "${sp:-}" || "$sp" =~ ^# ]]; then continue; fi
    if [[ -z "${n:-}" ]]; then continue; fi
    CHROM_COUNT["$sp"]="$n"
  done < "$CHROM_COUNT_FILE"
fi

# Builds one species' seqids.txt row: all autosomes first (plain sort -V,
# same as before), then that species' Z scaffold(s) (if known), then its W
# scaffold(s) (if known) -- so Z and W always render at the end of the row
# instead of wherever they'd fall alphabetically/numerically among the
# autosomes. Scaffolds named in Z_SCAFFOLD/W_SCAFFOLD but not actually
# present in this species' BED are silently dropped (covers the common case
# of a species with no W/Z entry at all).
#
# NOTE on style: every conditional below uses explicit if/then rather than
# `[[ cond ]] && action`. That shorthand returns the exit status of `[[ ]]`
# itself whenever cond is false (the common case here -- most species have
# no Z, no W, or neither), and under `set -euo pipefail`, if that shorthand
# happens to be a function's last executed statement, the function's
# nonzero return silently aborts the ENTIRE calling loop the first time
# it's hit -- confirmed empirically while building this (a 3-species test
# loop using the shorthand died after species 1 with no other symptom).
# if/then does not have this problem.
ordered_seqids_line() {
  local species="$1" bed="$2"
  local all_scafs autosomes z_line w_line
  all_scafs=$(cut -f1 "$bed" | sort -V | uniq)

  local -A is_sex
  local scaf s
  for scaf in "${Z_SCAFFOLD[$species]:-}" "${W_SCAFFOLD[$species]:-}"; do
    if [[ -n "$scaf" ]]; then
      IFS=',' read -r -a _s <<< "$scaf"
      for s in "${_s[@]}"; do
        s="${s## }"; s="${s%% }"
        if [[ -n "$s" ]]; then
          is_sex["$s"]=1
        fi
      done
    fi
  done

  autosomes=$(echo "$all_scafs" | while read -r s; do
    if [[ -z "${is_sex[$s]:-}" ]]; then echo "$s"; fi
  done)

  z_line=""
  if [[ -n "${Z_SCAFFOLD[$species]:-}" ]]; then
    IFS=',' read -r -a _z <<< "${Z_SCAFFOLD[$species]}"
    for s in "${_z[@]}"; do
      s="${s## }"; s="${s%% }"
      if echo "$all_scafs" | grep -q -w -- "$s"; then
        z_line="${z_line:+$z_line,}${s}"
      fi
    done
  fi

  w_line=""
  if [[ -n "${W_SCAFFOLD[$species]:-}" ]]; then
    IFS=',' read -r -a _w <<< "${W_SCAFFOLD[$species]}"
    for s in "${_w[@]}"; do
      s="${s## }"; s="${s%% }"
      if echo "$all_scafs" | grep -q -w -- "$s"; then
        w_line="${w_line:+$w_line,}${s}"
      fi
    done
  fi

  {
    echo "$autosomes" | paste -sd, -
    if [[ -n "$z_line" ]]; then echo "$z_line"; fi
    if [[ -n "$w_line" ]]; then echo "$w_line"; fi
  } | paste -sd, - | sed 's/^,//; s/,$//; s/,,*/,/g'
}

i=0
for SPECIES in "${SPECIES_ORDER[@]}"; do
  bed=$(bed_for "$SPECIES")
  plot_bed="$bed"

  if [[ -s "$bed" ]]; then
    n_chrom="${CHROM_COUNT[$SPECIES]:-}"
    if [[ -n "$n_chrom" ]]; then
      fai="${OUT_DIR}/${SPECIES}/${SPECIES}.fa.fai"
      candidate_plot_bed="${OUT_DIR}/${SPECIES}/${SPECIES}.plot.bed"
      # Guarded explicitly: build_plot_seqids.py failing for ONE species (a
      # bug in it, a bad chromosome_count.txt entry, whatever) must NOT take
      # down the rest of this loop under set -e -- that silently leaves
      # every later species' seqids.txt row (and this whole script's final
      # layout.txt rebuild) stale/missing, which is much more confusing to
      # debug than one clearly-logged fallback. On failure, fall back to
      # the untouched ordered_seqids_line/original bed for this species only.
      seqids_line=""
      if seqids_line=$(${MRUN} python3 "${SCRIPT_DIR}/build_plot_seqids.py" \
          "$SPECIES" "$bed" "$fai" "$n_chrom" \
          "${Z_SCAFFOLD[$SPECIES]:-}" "${W_SCAFFOLD[$SPECIES]:-}" "$candidate_plot_bed" \
          2>>"${OUT_DIR}/logs/build_plot_seqids.log"); then
        echo "$seqids_line" >> "$SEQIDS_FILE"
        plot_bed="$candidate_plot_bed"
      else
        log "ERROR: build_plot_seqids.py FAILED for ${SPECIES} -- see ${OUT_DIR}/logs/build_plot_seqids.log. Falling back to its untouched, unfiltered bed for this run (plot_bed stays ${bed})."
        ordered_seqids_line "$SPECIES" "$bed" >> "$SEQIDS_FILE"
      fi
    else
      ordered_seqids_line "$SPECIES" "$bed" >> "$SEQIDS_FILE"
    fi
  else
    echo "" >> "$SEQIDS_FILE"
    log "WARNING: no BED for ${SPECIES} -- empty row in seqids.txt (run 01.build_inputs.sh)"
  fi

  y=$(LC_NUMERIC=C awk -v s="$y_start" -v d="$y_step" -v i="$i" 'BEGIN{printf "%.2f", s - i*d}')
  LC_NUMERIC=C printf " %.2f,     %.1f,    %.1f,       %d, %s, %s, %s, %s\n" \
    "$y" "$x_start" "$x_end" "$rot" "$color" "$SPECIES" "$va" "$plot_bed" >> "$LAYOUTY_FILE"
  i=$((i + 1))
done

cat "$LAYOUTY_FILE" "$EDGES_FILE" > "${OUT_DIR}/layout.txt"

# ---- 5. Z-column alignment: aggregate per-species natural Z-start
#         positions (appended by build_plot_seqids.py above, one line per
#         species with a resolved Z) into a per-species EXTRA-GAP delta,
#         so every species' Z bar can be shifted to a shared target x
#         position -- mirrors the already-working W right-alignment, which
#         needed no such computation since W is always the last bar in
#         every row. Consumed by karyotype.py's (separately patched)
#         Track.draw()/update_offsets() via the KARYOTYPE_Z_GAP_FILE env
#         var, set in 04.plot_macrosynteny.sh. -------------------------------
ALIGN_FILE="${OUT_DIR}/z_alignment.tsv"
DELTA_FILE="${OUT_DIR}/z_gap_deltas.tsv"
if [[ -s "$ALIGN_FILE" ]]; then
  awk -F'\t' '
    { z[$1] = $2; if ($2 > maxz) maxz = $2 }
    END {
      for (sp in z) printf "%s\t%.10f\n", sp, (maxz - z[sp])
    }
  ' "$ALIGN_FILE" > "$DELTA_FILE"
  log "Z-alignment: wrote per-species extra-gap deltas -> ${DELTA_FILE}"
  while IFS=$'\t' read -r sp delta; do
    log "  ${sp}: +${delta}"
  done < "$DELTA_FILE"
else
  log "NOTE: ${ALIGN_FILE} is empty/missing (no species had a resolved Z?) -- ${DELTA_FILE} not written; karyotype.py's Z-alignment patch will fall back to zero delta for every species."
  rm -f "$DELTA_FILE"
fi

log "Done. ${#EDGES[@]} edges (staircase + direct W-links). layout.txt + seqids.txt ready in ${OUT_DIR}/"
log "Next: 04.plot_macrosynteny.sh"
