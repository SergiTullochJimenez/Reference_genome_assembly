#!/usr/bin/env bash
# =============================================================================
# 01.Circos_tracks.sh -- data for the genome Circos plots (02.Circos_plot.R) and the Hi-C junction check
# Full pipeline: index → map → generate circos BED tracks
# Idempotent: each section checks its output before running.
# Re-run safely at any time — only missing or empty outputs are regenerated.
#
# USAGE:
#   bash 01.Circos_tracks.sh --species Species_name \
#       [--threads 16] [--binsize 1000000] [--minscaf 1000000]
#
# All file paths are derived from --species (see PATH LAYOUT below) --
# no other flags are needed for a normal run.
#
# PATH LAYOUT (fixed, one species per run):
#   WD      = /path/to/project/<species>
#   DIR     = ${WD}/res/07.Circos
#   masked  = ${DIR}/<species>.masked.fa           (softmasked assembly;
#                                                    used for BOTH indexing/
#                                                    mapping AND repeat density
#                                                    -- no separate unmasked
#                                                    fasta needed)
#   gtf     = ${DIR}/<species>.gtf
#   busco_asm = ${DIR}/<species>_full_table.tsv
#   busco_ann = ${DIR}/<species>_proteins_full_table.tsv
#   ONT   reads = ${WD}/res/01.1.Trimming/filtered_ONT.fq.gz
#   HiC   R1    = ${WD}/res/01.1.Trimming/filtered_Arima_HIC_1.fq.gz
#   HiC   R2    = ${WD}/res/01.1.Trimming/filtered_Arima_HIC_2.fq.gz
#   outdir      = ${DIR}/<species>/                (all BED/track outputs
#                                                    + bam/, logs/ land here)
#
# DEPENDENCIES:
#   samtools >= 1.17    bedtools >= 2.30    bwa >= 0.7 (classic, for HiC)
#   minimap2 >= 2.24    python3 (stdlib)    awk / sort
#
# OPTIONAL:
#   cooler   — HiC contact matrix export (section 11, commented out)
#
# OUTPUT BAMs (written to <outdir>/bam/):
#   <species>.ont.bam          — ONT mapped + sorted
#   <species>.hic.markdup.bam  — HiC mapped, chimera-filtered, markdup
# =============================================================================

set -euo pipefail

# ── Fixed root -- edit here if the mount point ever changes ───────────────────
BASE_ROOT="/path/to/project"

# ── Defaults ──────────────────────────────────────────────────────────────────
BINSIZE=1000000
MINSCAF=1000000
THREADS=20
SPECIES=""

# ── Argument parsing ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case $1 in
    --species)     SPECIES="$2";     shift 2 ;;
    --threads)     THREADS="$2";     shift 2 ;;
    --binsize)     BINSIZE="$2";     shift 2 ;;
    --minscaf)     MINSCAF="$2";     shift 2 ;;
    *) echo "ERROR: Unknown argument: $1"; exit 1 ;;
  esac
done

if [[ -z "$SPECIES" ]]; then
  echo "ERROR: --species is required, e.g.: bash $0 --species Macroprotodon_brevis"
  exit 1
fi

# ── Derive all paths from --species ────────────────────────────────────────────
SAMPLE="$SPECIES"
WD="${BASE_ROOT}/${SPECIES}"
DIR="${WD}/res/07.Circos"

MASKED="${DIR}/${SPECIES}.masked.fa"
FASTA="$MASKED"                      # same softmasked file used for both
GFF="${DIR}/${SPECIES}.gtf"          # GTF (not GFF3) -- see Section 8 parsing
BUSCO_ASM="${DIR}/${SPECIES}_full_table.tsv"
BUSCO_ANN="${DIR}/${SPECIES}_proteins_full_table.tsv"

READS_DIR="${WD}/res/01.1.Trimming"
ONT_READS="${READS_DIR}/filtered_ONT.fq.gz"
HIC_R1="${READS_DIR}/filtered_Arima_HIC_1.fq.gz"
HIC_R2="${READS_DIR}/filtered_Arima_HIC_2.fq.gz"

OUTDIR="${DIR}/${SPECIES}"

# ── Validate required inputs ──────────────────────────────────────────────────
for var in MASKED GFF BUSCO_ASM BUSCO_ANN ONT_READS HIC_R1 HIC_R2; do
  if [[ ! -f "${!var}" ]]; then
    echo "ERROR: File not found for ${var}: ${!var}"; exit 1
  fi
done

# ── Output directories ────────────────────────────────────────────────────────
mkdir -p "${OUTDIR}" "${OUTDIR}/bam" "${OUTDIR}/logs"
BAMDIR="${OUTDIR}/bam"
LOGDIR="${OUTDIR}/logs"

ONT_BAM="${BAMDIR}/${SAMPLE}.ont.bam"
HIC_BAM="${BAMDIR}/${SAMPLE}.hic.markdup.bam"

# =============================================================================
# HELPERS
# =============================================================================

log() { echo "[$(date '+%H:%M:%S')] $*"; }

# output_ok FILE [FILE2 ...]
# Returns 0 (true) if ALL listed files exist and are non-empty, else 1.
output_ok() {
  for f in "$@"; do
    [[ -f "$f" && -s "$f" ]] || return 1
  done
  return 0
}

# run_section LABEL OUTPUT_FILES... --- COMMAND
# Checks if all OUTPUT_FILES are present and non-empty.
# If yes: skips with a message. If no: runs COMMAND, then verifies outputs.
#
# Usage:
#   run_section "Section name" \
#       OUT:"file1.bed" OUT:"file2.bed" \
#       CMD: actual_command_here
#
# Because bash functions can't easily take mixed args, we use a simpler
# two-step pattern below: check_outputs + run_if_needed.

SECTION_OUTPUTS=()   # populated per section
SECTION_FAILED=0

check_section() {
  # Usage: check_section "Label" out1 [out2 ...]
  # Sets global SKIP_SECTION=1 if all outputs exist and are non-empty.
  local label="$1"; shift
  local missing=()
  for f in "$@"; do
    if [[ ! -f "$f" || ! -s "$f" ]]; then
      missing+=("$(basename "$f")")
    fi
  done
  if [[ ${#missing[@]} -eq 0 ]]; then
    log "  [SKIP] ${label} — outputs already present"
    SKIP_SECTION=1
  else
    log "  [RUN ] ${label} — missing/empty: ${missing[*]}"
    SKIP_SECTION=0
  fi
}

verify_outputs() {
  # Usage: verify_outputs "Label" out1 [out2 ...]
  local label="$1"; shift
  local failed=()
  for f in "$@"; do
    if [[ ! -f "$f" || ! -s "$f" ]]; then
      failed+=("$(basename "$f")")
    fi
  done
  if [[ ${#failed[@]} -gt 0 ]]; then
    log "  [WARN] ${label} produced empty/missing outputs: ${failed[*]}"
    SECTION_FAILED=1
  else
    local total_lines=0
    for f in "$@"; do
      n=$(wc -l < "$f" 2>/dev/null || echo 0)
      total_lines=$((total_lines + n))
    done
    log "  [OK  ] ${label} — ${total_lines} total lines across ${#@} file(s)"
  fi
}

# Temp file for scaffold name list (used in BUSCO filtering)
SCAF_LIST=$(mktemp)
trap "rm -f $SCAF_LIST" EXIT

# =============================================================================
# SECTION 0 — Scaffold sizes, filter, windows
# =============================================================================
log "========================================"
log "Section 0: Reference indexing & windows"
log "========================================"

# 0a. samtools faidx
check_section "samtools faidx" "${FASTA}.fai"
if [[ $SKIP_SECTION -eq 0 ]]; then
  samtools faidx "$FASTA"
  verify_outputs "samtools faidx" "${FASTA}.fai"
fi

# 0b. bwa index (classic, for HiC)
check_section "bwa index" "${FASTA}.bwt"
if [[ $SKIP_SECTION -eq 0 ]]; then
  bwa index "$FASTA" 2> "${LOGDIR}/bwa_index.log"
  verify_outputs "bwa index" "${FASTA}.bwt"
fi

# 0c. minimap2 index (cached .mmi for ONT)
MM2_IDX="${FASTA%.fasta}.ont.mmi"
check_section "minimap2 index" "$MM2_IDX"
if [[ $SKIP_SECTION -eq 0 ]]; then
  minimap2 -x map-ont -d "$MM2_IDX" "$FASTA" 2> "${LOGDIR}/mm2_index.log"
  verify_outputs "minimap2 index" "$MM2_IDX"
fi

# 0d. Scaffold size tables
check_section "chrom sizes" "${OUTDIR}/chrom_sizes.txt" "${OUTDIR}/genome.txt"
if [[ $SKIP_SECTION -eq 0 ]]; then
  awk '{print $1"\t"$2}' "${FASTA}.fai" > "${OUTDIR}/all_chrom_sizes.txt"
  awk -v minlen="$MINSCAF" '$2 >= minlen {print $1"\t"$2}' \
      "${OUTDIR}/all_chrom_sizes.txt" \
      | sort -k2,2rn > "${OUTDIR}/chrom_sizes.txt"
  cp "${OUTDIR}/chrom_sizes.txt" "${OUTDIR}/genome.txt"
  verify_outputs "chrom sizes" "${OUTDIR}/chrom_sizes.txt"
fi

NSCAF=$(wc -l < "${OUTDIR}/chrom_sizes.txt")
log "  Scaffolds >= $(( MINSCAF/1000000 )) Mb: ${NSCAF}"

# 0e. Windows BED
check_section "windows BED" "${OUTDIR}/windows.bed"
if [[ $SKIP_SECTION -eq 0 ]]; then
  bedtools makewindows \
      -g "${OUTDIR}/genome.txt" \
      -w "$BINSIZE" \
      > "${OUTDIR}/windows.bed"
  verify_outputs "windows BED" "${OUTDIR}/windows.bed"
fi

# Populate scaffold name list for BUSCO filtering (always refresh)
awk '{print $1}' "${OUTDIR}/chrom_sizes.txt" > "$SCAF_LIST"

# =============================================================================
# SECTION 1 — Map ONT reads
# =============================================================================
log "========================================"
log "Section 1: ONT mapping (minimap2)"
log "========================================"

check_section "ONT mapping" "$ONT_BAM"
if [[ $SKIP_SECTION -eq 0 ]]; then
  minimap2 \
      -ax map-ont \
      -t "$THREADS" \
      --secondary=no \
      -R "@RG\tID:ONT\tSM:${SAMPLE}\tPL:ONT" \
      "$MM2_IDX" \
      "$ONT_READS" \
      2> "${LOGDIR}/minimap2_ont.log" \
  | samtools sort \
      -@ "$THREADS" \
      -m 4G \
      --write-index \
      -o "$ONT_BAM" -
  verify_outputs "ONT mapping" "$ONT_BAM"
fi

samtools flagstat -@ "$THREADS" "$ONT_BAM" \
    > "${LOGDIR}/${SAMPLE}.ont.flagstat.txt" 2>/dev/null || true

# =============================================================================
# SECTION 2 — Map HiC reads (bwa mem classic, -5SP)
# =============================================================================
log "========================================"
log "Section 2: HiC mapping (bwa mem)"
log "========================================"

HIC_R1_RAW="${BAMDIR}/${SAMPLE}.hic.R1.bam"
HIC_R2_RAW="${BAMDIR}/${SAMPLE}.hic.R2.bam"
HIC_MERGED="${BAMDIR}/${SAMPLE}.hic.merged.bam"
HIC_NAMESORT="${BAMDIR}/${SAMPLE}.hic.namesorted.bam"

check_section "HiC mapping" "$HIC_BAM"
if [[ $SKIP_SECTION -eq 0 ]]; then

  log "  Mapping R1 ..."
  bwa mem -t "$THREADS" -5SP \
      -R "@RG\tID:HIC\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:hic" \
      "$FASTA" "$HIC_R1" \
      2> "${LOGDIR}/bwa_hic_R1.log" \
  | samtools view -@ "$THREADS" -b -F 4 -q 10 \
  | samtools sort -@ "$THREADS" -m 2G -n -o "$HIC_R1_RAW" -

  log "  Mapping R2 ..."
  bwa mem -t "$THREADS" -5SP \
      -R "@RG\tID:HIC\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:hic" \
      "$FASTA" "$HIC_R2" \
      2> "${LOGDIR}/bwa_hic_R2.log" \
  | samtools view -@ "$THREADS" -b -F 4 -q 10 \
  | samtools sort -@ "$THREADS" -m 2G -n -o "$HIC_R2_RAW" -

  log "  Merging + coordinate-sorting ..."
  samtools merge -@ "$THREADS" -n -f "$HIC_MERGED" "$HIC_R1_RAW" "$HIC_R2_RAW"
  rm -f "$HIC_R1_RAW" "$HIC_R2_RAW"
  samtools sort -@ "$THREADS" -m 4G -o "$HIC_NAMESORT" "$HIC_MERGED"
  rm -f "$HIC_MERGED"

  log "  Marking duplicates ..."
  samtools markdup \
      -@ "$THREADS" \
      --write-index \
      -f "${LOGDIR}/${SAMPLE}.hic.markdup.txt" \
      "$HIC_NAMESORT" "$HIC_BAM"
  rm -f "$HIC_NAMESORT"

  verify_outputs "HiC mapping" "$HIC_BAM"
fi

samtools flagstat -@ "$THREADS" "$HIC_BAM" \
    > "${LOGDIR}/${SAMPLE}.hic.flagstat.txt" 2>/dev/null || true

# =============================================================================
# SECTION 5 — GC content
# =============================================================================
log "========================================"
log "Section 5: GC content"
log "========================================"

check_section "GC content" "${OUTDIR}/gc_1mb.bed"
if [[ $SKIP_SECTION -eq 0 ]]; then
  bedtools nuc \
      -fi "$FASTA" \
      -bed "${OUTDIR}/windows.bed" \
  | awk 'BEGIN{OFS="\t"} NR>1 {print $1, $2, $3, $5}' \
  > "${OUTDIR}/gc_1mb.bed"
  verify_outputs "GC content" "${OUTDIR}/gc_1mb.bed"
fi

# =============================================================================
# SECTION 6 — Repeat density (softmasked FASTA)
# =============================================================================
log "========================================"
log "Section 6: Repeat density"
log "========================================"

check_section "Repeat density" "${OUTDIR}/repeat_1mb.bed"
if [[ $SKIP_SECTION -eq 0 ]]; then
  if [[ ! -s "$MASKED" ]]; then
    log "  [WARN] Masked FASTA is empty or missing: $MASKED"
    log "         Skipping repeat track — re-run with a valid --masked file"
    # Write an empty placeholder so the check knows we tried
    touch "${OUTDIR}/repeat_1mb.bed"
  else
    python3 - "$MASKED" "${OUTDIR}/windows.bed" "${OUTDIR}/repeat_1mb.bed" <<'PYEOF'
import sys, gzip
fasta_path, bed_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
seq = {}; cur = None
opener = gzip.open if fasta_path.endswith(".gz") else open
with opener(fasta_path, "rt") as fh:
    for line in fh:
        line = line.rstrip()
        if line.startswith(">"):
            cur = line[1:].split()[0]; seq[cur] = []
        elif cur:
            seq[cur].append(line)
seq = {k: "".join(v) for k, v in seq.items()}
with open(bed_path) as bed, open(out_path, "w") as out:
    for line in bed:
        p = line.rstrip().split("\t"); chrom, s, e = p[0], int(p[1]), int(p[2])
        if chrom not in seq: continue
        w = seq[chrom][s:e]
        if not w: continue
        out.write(f"{chrom}\t{s}\t{e}\t{sum(1 for c in w if c.islower())/len(w):.6f}\n")
PYEOF
    verify_outputs "Repeat density" "${OUTDIR}/repeat_1mb.bed"
  fi
fi

# =============================================================================
# SECTION 7 — Gene density
# =============================================================================
log "========================================"
log "Section 7: Gene density"
log "========================================"

check_section "Gene density" "${OUTDIR}/gene_density_1mb.bed"
if [[ $SKIP_SECTION -eq 0 ]]; then
  # Filter to this species' own scaffold whitelist (from chrom_sizes.txt) rather
  # than a hardcoded name pattern -- scaffold naming isn't consistent across all
  # 8 species (e.g. Hierophis_viridiflavus doesn't use "SUPER_SUPER*"), so a
  # regex like /^SUPER_SUPER/ would silently zero out this track for it.
  awk -v OFS='\t' '
      NR==FNR { keep[$1]=1; next }
      $3 == "gene" && ($1 in keep) { print $1, $4-1, $5 }
  ' "$SCAF_LIST" "$GFF" \
  | bedtools sort -g "${OUTDIR}/genome.txt" \
  | bedtools coverage \
      -a "${OUTDIR}/windows.bed" \
      -b - -counts -sorted \
      -g "${OUTDIR}/genome.txt" \
  | awk '{print $1"\t"$2"\t"$3"\t"$4}' \
  > "${OUTDIR}/gene_density_1mb.bed"
  verify_outputs "Gene density" "${OUTDIR}/gene_density_1mb.bed"
fi

# =============================================================================
# SECTION 8 — BUSCO hits (assembly + annotation)
# Uses grep -F -f for safe scaffold name matching (no regex fragility)
# =============================================================================
log "========================================"
log "Section 8: BUSCO hits"
log "========================================"

busco_to_bed() {
  # For genome-mode BUSCO: coordinates are directly in the table (col4=start, col5=end)
  local infile="$1" outfile="$2" label="$3"
  awk 'BEGIN{OFS="\t"}
       /^#/ { next }
       $2 == "Missing" { next }
       NF >= 6 {
         s = ($4 < $5) ? $4 : $5
         e = ($4 < $5) ? $5 : $4
         print $3, s, e, $2
       }' "$infile" \
  | grep -F -f "$SCAF_LIST" \
  | bedtools sort -g "${OUTDIR}/genome.txt" -i - \
  > "$outfile"

  local n; n=$(wc -l < "$outfile")
  if [[ "$n" -eq 0 ]]; then
    log "  [WARN] ${label}: 0 entries — check scaffold names match chrom_sizes.txt"
    log "         chrom_sizes example : '$(head -1 "$SCAF_LIST")'"
    log "         BUSCO table example : '$(awk '!/^#/ && $2!="Missing" && NF>=6 {print $3; exit}' "$infile")'"
  else
    log "  [OK  ] ${label}: ${n} entries"
  fi
}

busco_ann_to_bed() {
  local infile="$1" outfile="$2" gff="$3"

  python3 - "$infile" "$outfile" "$gff" "$SCAF_LIST" <<'PYEOF'
import sys, re

busco_path  = sys.argv[1]
out_path    = sys.argv[2]
gff_path    = sys.argv[3]
scaf_path   = sys.argv[4]

# ── Load scaffold whitelist ───────────────────────────────────────────────────
with open(scaf_path) as f:
    keep = set(l.split()[0] for l in f if l.strip())

# ── Parse annotation BUSCO table → transcript_id : status ────────────────────
# Columns: BuscoID  Status  TranscriptID  Score  Length  URL  Desc
# Strip surrounding pipes if present: |anno2.foo| -> anno2.foo
tx_status = {}
with open(busco_path) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        status = parts[1]
        if status == "Missing":
            continue
        tx_id = parts[2].strip("|")   # strip leading/trailing pipes
        # Keep only the most informative status per transcript
        # priority: Complete > Duplicated > Fragmented
        priority = {"Complete": 0, "Duplicated": 1, "Fragmented": 2}
        if tx_id not in tx_status or \
           priority.get(status, 9) < priority.get(tx_status[tx_id], 9):
            tx_status[tx_id] = status

print(f"  Annotation BUSCO: {len(tx_status)} transcript IDs to look up",
      flush=True)

# ── Parse GFF/GTF — extract transcript coords on kept scaffolds ──────────────
# Supports both attribute styles so this works regardless of which one the
# annotation file uses:
#   GFF3: ID=transcript123;Parent=gene123          -> id_re_gff3
#   GTF : gene_id "gene123"; transcript_id "tx123" -> id_re_gtf (preferred
#         when present, since "transcript_id" is what BUSCO's column-3
#         transcript IDs actually correspond to)
id_re_gff3 = re.compile(r'ID=([^;]+)')
id_re_gtf  = re.compile(r'transcript_id\s+"([^"]+)"')
results = []

with open(gff_path) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue
        feat_type = parts[2]
        if feat_type not in ("mRNA", "transcript"):
            continue
        scaffold = parts[0]
        if scaffold not in keep:
            continue
        attrs = parts[8]
        m = id_re_gtf.search(attrs)
        if m:
            tx_id = m.group(1)
        else:
            m = id_re_gff3.search(attrs)
            if not m:
                continue
            tx_id = m.group(1).strip("|")
        if tx_id not in tx_status:
            continue
        start = min(int(parts[3]), int(parts[4])) - 1   # BED is 0-based
        end   = max(int(parts[3]), int(parts[4]))
        results.append((scaffold, start, end, tx_status[tx_id]))

print(f"  Matched {len(results)} transcripts to SUPER_SUPER coordinates",
      flush=True)

if len(results) == 0:
    print("  WARNING: No matches found. Check that GFF transcript IDs match "
          "BUSCO table column 3.", flush=True)

# ── Sort by scaffold order then position ─────────────────────────────────────
# Use scaffold whitelist order (same as chrom_sizes.txt)
with open(scaf_path) as f:
    scaf_order = {l.split()[0]: i for i, l in enumerate(f) if l.strip()}

results.sort(key=lambda x: (scaf_order.get(x[0], 9999), x[1]))

with open(out_path, "w") as out:
    for scaffold, start, end, status in results:
        out.write(f"{scaffold}\t{start}\t{end}\t{status}\n")
PYEOF

  local n; n=$(wc -l < "$outfile")
  if [[ "$n" -eq 0 ]]; then
    log "  [WARN] Annotation BUSCO: 0 entries matched — check transcript ID format"
  else
    log "  [OK  ] Annotation BUSCO: ${n} entries"
  fi
}

check_section "BUSCO assembly" "${OUTDIR}/busco_assembly_hits.bed"
if [[ $SKIP_SECTION -eq 0 ]]; then
  busco_to_bed "$BUSCO_ASM" "${OUTDIR}/busco_assembly_hits.bed" "Assembly BUSCO"
fi

check_section "BUSCO annotation" "${OUTDIR}/busco_annotation_hits.bed"
if [[ $SKIP_SECTION -eq 0 ]]; then
  busco_ann_to_bed "$BUSCO_ANN" "${OUTDIR}/busco_annotation_hits.bed" "$GFF"
fi

# =============================================================================
# SECTION 9 — Coverage tracks
# =============================================================================
log "========================================"
log "Section 9: Coverage tracks"
log "========================================"

check_section "ONT coverage" "${OUTDIR}/ont_cov_1mb.bed"
if [[ $SKIP_SECTION -eq 0 ]]; then
  while read chrom length; do
    samtools depth -a \
        -r "${chrom}" \
        "$ONT_BAM" \
    | awk -v chrom="$chrom" -v binsize="$BINSIZE" '
        { bin = int($2 / binsize); depth[bin] += $3; count[bin]++ }
        END { for (b in depth) printf "%s\t%d\t%d\t%.7f\n",
              chrom, b*binsize, (b+1)*binsize, depth[b]/count[b] }'
  done < "${OUTDIR}/genome.txt" \
  | sort -k1,1 -k2,2n \
  > "${OUTDIR}/ont_cov_1mb.bed"
  verify_outputs "ONT coverage" "${OUTDIR}/ont_cov_1mb.bed"
fi

check_section "HiC coverage" "${OUTDIR}/hic_cov_1mb.bed"
if [[ $SKIP_SECTION -eq 0 ]]; then
  while read chrom length; do
    samtools depth -a \
        -r "${chrom}" \
        -G 1024 \
        "$HIC_BAM" \
    | awk -v chrom="$chrom" -v binsize="$BINSIZE" '
        { bin = int($2 / binsize); depth[bin] += $3; count[bin]++ }
        END { for (b in depth) printf "%s\t%d\t%d\t%.7f\n",
              chrom, b*binsize, (b+1)*binsize, depth[b]/count[b] }'
  done < "${OUTDIR}/genome.txt" \
  | sort -k1,1 -k2,2n \
  > "${OUTDIR}/hic_cov_1mb.bed"
  verify_outputs "HiC coverage" "${OUTDIR}/hic_cov_1mb.bed"
fi

# =============================================================================
# SECTION 10 — Gap / N density
# =============================================================================
log "========================================"
log "Section 10: Gap/N density"
log "========================================"

check_section "Gap/N density" "${OUTDIR}/gaps_1mb.bed"
if [[ $SKIP_SECTION -eq 0 ]]; then
  python3 - "$FASTA" "${OUTDIR}/windows.bed" "${OUTDIR}/gaps_1mb.bed" <<'PYEOF'
import sys, gzip
fasta_path, bed_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
seq = {}; cur = None
opener = gzip.open if fasta_path.endswith(".gz") else open
with opener(fasta_path, "rt") as fh:
    for line in fh:
        line = line.rstrip()
        if line.startswith(">"):
            cur = line[1:].split()[0]; seq[cur] = []
        elif cur:
            seq[cur].append(line)
seq = {k: "".join(v) for k, v in seq.items()}
with open(bed_path) as bed, open(out_path, "w") as out:
    for line in bed:
        p = line.rstrip().split("\t"); chrom, s, e = p[0], int(p[1]), int(p[2])
        if chrom not in seq: continue
        w = seq[chrom][s:e]
        if not w: continue
        out.write(f"{chrom}\t{s}\t{e}\t{sum(1 for c in w if c.upper()=='N')/len(w):.6f}\n")
PYEOF
  verify_outputs "Gap/N density" "${OUTDIR}/gaps_1mb.bed"
fi

# =============================================================================
# (OPTIONAL) SECTION 11 — HiC contact matrix
# Uncomment and set HIC_COOL to enable.
# =============================================================================
# HIC_COOL="path/to/contacts.mcool::resolutions/1000000"
# check_section "HiC matrix" "${OUTDIR}/hic_matrix.tsv"
# if [[ $SKIP_SECTION -eq 0 ]]; then
#   python3 - "$HIC_COOL" "${OUTDIR}/hic_matrix.tsv" "${OUTDIR}/chrom_sizes.txt" <<'PYEOF'
# import sys, cooler, numpy as np
# cool_path, out_tsv, chrom_file = sys.argv[1], sys.argv[2], sys.argv[3]
# with open(chrom_file) as f:
#     chroms = [l.split()[0] for l in f]
# c = cooler.Cooler(cool_path)
# mat = np.nan_to_num(c.matrix(balance=True).fetch(chroms), nan=0.0)
# np.savetxt(out_tsv, mat, delimiter="\t", fmt="%.6f")
# PYEOF
#   verify_outputs "HiC matrix" "${OUTDIR}/hic_matrix.tsv"
# fi

# =============================================================================
# FINAL SUMMARY
# =============================================================================
log "========================================"
log "Summary: output file status"
log "========================================"

ALL_OK=1

print_status() {
  local label="$1"; shift
  local ok=1
  for f in "$@"; do
    [[ -f "$f" && -s "$f" ]] || { ok=0; ALL_OK=0; }
  done
  if [[ $ok -eq 1 ]]; then
    local lines=0
    for f in "$@"; do lines=$(( lines + $(wc -l < "$f") )); done
    printf "  ✓  %-35s %d lines\n" "$label" "$lines"
  else
    for f in "$@"; do
      if [[ ! -f "$f" ]]; then
        printf "  ✗  %-35s MISSING\n" "$label ($(basename "$f"))"
      elif [[ ! -s "$f" ]]; then
        printf "  ✗  %-35s EMPTY\n"   "$label ($(basename "$f"))"
      fi
    done
  fi
}

print_status "chrom_sizes.txt"           "${OUTDIR}/chrom_sizes.txt"
print_status "windows.bed"               "${OUTDIR}/windows.bed"
print_status "gc_1mb.bed"                "${OUTDIR}/gc_1mb.bed"
print_status "repeat_1mb.bed"            "${OUTDIR}/repeat_1mb.bed"
print_status "gene_density_1mb.bed"      "${OUTDIR}/gene_density_1mb.bed"
print_status "busco_assembly_hits.bed"   "${OUTDIR}/busco_assembly_hits.bed"
print_status "busco_annotation_hits.bed" "${OUTDIR}/busco_annotation_hits.bed"
print_status "ont_cov_1mb.bed"           "${OUTDIR}/ont_cov_1mb.bed"
print_status "hic_cov_1mb.bed"           "${OUTDIR}/hic_cov_1mb.bed"
print_status "gaps_1mb.bed"              "${OUTDIR}/gaps_1mb.bed"
print_status "ONT BAM"                   "$ONT_BAM"
print_status "HiC BAM"                   "$HIC_BAM"

echo ""
if [[ $ALL_OK -eq 1 ]]; then
  log "All outputs present and non-empty. Ready to run: Rscript genome_circos.R"
else
  log "Some outputs are missing or empty (marked ✗ above)."
  log "Re-run this script — completed sections will be skipped automatically."
fi
if [[ $SECTION_FAILED -eq 1 ]]; then
  log "[WARN] One or more sections reported warnings — check the log above."
fi
