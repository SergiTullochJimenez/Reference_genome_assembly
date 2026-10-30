#!/bin/bash
# =============================================================================
# macro_micro_pipeline.sh  --  Steps 1 and 3 of the micro/macrochromosome plan
# (Step 0, the per-chromosome table with the macro/micro classification, is made by 00.chromosome_table.py)
#
#   Step 1  gene composition & structure, every annotated gene, 8 new genomes
#           (GC1/2/3, GC at fourfold sites, intron GC, ENC, CDS/intron lengths)
#   Step 2  OrthoFinder on the longest isoform of every gene (proteins translated
#           from the Step 1 CDS, so gene IDs are identical across steps):
#           core / lineage-specific / duplicated genes per chromosome
#   Step 3  dN, dS and omega (dN/dS) for OrthoFinder single-copy orthologs (default;
#           or genome-mode BUSCO genes, GENESET=busco) present in all 8,
#           located on the same chromosome class (micro/macro/Z) in all 8
#           -> codon alignment (MACSE, or MAFFT back-translated) -> codeml M0 on the fixed species tree
#   Collect tables, statistics and plots for both steps
#
# ONE COMMAND DOES EVERYTHING (after the one-time setup):
#   bash macro_micro_pipeline.sh setup     # once: builds the conda env (sbatch job)
#   bash macro_micro_pipeline.sh check     # preflight: lists what is found / missing
#   bash macro_micro_pipeline.sh submit    # submits all stages, chained with
#                                          # SLURM dependencies; nothing else to do
#   bash macro_micro_pipeline.sh status    # progress summary
#
# Resumable: every stage skips work whose output already exists.
# Final tables, statistics and plots: ${OUT_DIR}/results/
# =============================================================================
set -uo pipefail

# ============================== CONFIG ======================================
# Same genome and annotation inputs as analyses/macrosynteny
GENOME_DIR="/path/to/working_directory/data/Genomes"                # <Sp>/<Sp>.masked.fa
ANNOT_DIR="/path/to/Annotations"   # <Sp>.gtf
# Upload these two things here before 'check' (see README at the end of 'check'):
#   ${CHROM_TABLE} (per-chromosome table from 00.chromosome_table.py)
#   <Species>_full_table.tsv  x 8          (genome-mode BUSCO full tables, assembly/07.2.BUSCO_genome.sh)
INPUT_DIR="/path/to/working_directory/Macro_micro_inputs"
CHROM_TABLE="per_chromosome_summary_v2.csv"   # Step 0 table (00.2.Macro_micro_classes.py); passed to the SLURM jobs via env
OUT_DIR="/path/to/working_directory/Macro_micro_results"

SPECIES=(Coronella_austriaca Coronella_girondica Hierophis_viridiflavus
         Macroprotodon_brevis Macroprotodon_mauritanicus Natrix_astreptophora
         Zamenis_longissimus Zamenis_scalaris)

# Unrooted species tree for the 8, pruned from ML/nuclear_busco_tree_genomic.contree
# (Hemorrhois, the other Natrix, Malpolon, Vipera and Anolis removed). Short
# codes because PAML dislikes long names.
TREE="((HVIR,(MBRE,MMAU)),((ZSCA,ZLON),(CGIR,CAUS)),NAST);"

# SLURM
ACCOUNT="<account>"
PARTITION="<partition>"          # jobs of up to ~6 h
# Rooted version for OrthoFinder (-s); root between Natrix and the colubrids
TREE_ROOTED="(Natrix_astreptophora,((Hierophis_viridiflavus,(Macroprotodon_brevis,Macroprotodon_mauritanicus)),((Zamenis_scalaris,Zamenis_longissimus),(Coronella_girondica,Coronella_austriaca))));"
GENESET="orthofinder"  # gene set for Step 3: "orthofinder" (single-copy orthogroups) or "busco"

# OrthoFinder job (the heaviest one): adjust partition/time to your cluster limits
OF_PARTITION="${PARTITION}"
OF_CPUS=32
OF_MEM="128G"
OF_ANALYSIS_THREADS=2   # OrthoFinder -a: parallel gene-tree jobs. 8 ran out of RAM (2026-10-06)
OF_TIME="24:00:00"

N_CHUNKS=40           # step 3 array tasks (~250 genes each). 200 broke the per-user
                      # per-user job-submission limit of our cluster (QOSMaxSubmitJobPerUserLimit)
MAX_PARALLEL=40       # array tasks allowed to run at once
DNDS_PARTITION="<long_partition>"  # jobs > 6 h
DNDS_TIME="24:00:00"     # per array task (~270 genes: ~15 min with MAFFT, ~5 h with MACSE)
ALIGNER="mafft"          # codon aligner for step 3: "mafft" (protein alignment back-translated, ~3 s/gene)
                         # or "macse" (frameshift-aware, ~1 min/gene)

# Conda (same mechanism as 00.config.sh of the macrosynteny pipeline)
MODULE_MAMBA="miniconda3"   # module providing mamba/conda on your system
CONDA_ROOT="/path/to/conda"
unset CONDA_ENVS_PATH
export CONDA_ENVS_DIRS="${CONDA_ROOT}/envs"
export CONDA_PKGS_DIRS="${CONDA_ROOT}/pkgs"
ENV_PREFIX="${CONDA_ENVS_DIRS}/MICROMACRO"
OF_ENV_PREFIX="${CONDA_ENVS_DIRS}/ORTHOFINDER"   # separate env: avoids solver clashes
# =============================================================================

SCRIPTS="${OUT_DIR}/scripts"
LOGS="${OUT_DIR}/logs"
mkdir -p "$OUT_DIR" "$SCRIPTS" "$LOGS" "${OUT_DIR}/step1/proteins" "${OUT_DIR}/step2" "${OUT_DIR}/step3/genes" \
         "${OUT_DIR}/step3/chunks" "${OUT_DIR}/results"
SELF="$(readlink -f "$0")"
RUN="mamba run -p ${ENV_PREFIX}"
log() { echo "[$(date '+%F %T')] $*"; }

# ---------------------------------------------------------------------------
# Python helpers are embedded below and written to ${SCRIPTS} on every call,
# so this single file is the whole pipeline.
# ---------------------------------------------------------------------------
write_helpers() {

cat > "${SCRIPTS}/common.py" <<'PYEOF'
import gzip, re
CODES = {"Coronella_austriaca":"CAUS","Coronella_girondica":"CGIR","Hierophis_viridiflavus":"HVIR",
         "Macroprotodon_brevis":"MBRE","Macroprotodon_mauritanicus":"MMAU","Natrix_astreptophora":"NAST",
         "Zamenis_longissimus":"ZLON","Zamenis_scalaris":"ZSCA"}
COMP = str.maketrans("ACGTNacgtn","TGCANtgcan")
def revcomp(s): return s.translate(COMP)[::-1]

def read_fasta(path, keep=None):
    seqs, name, buf = {}, None, []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name is not None and (keep is None or name in keep):
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None and (keep is None or name in keep):
        seqs[name] = "".join(buf).upper()
    return seqs

GENCODE = {}
_b = "TCAG"; _aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i,a in enumerate(_b):
    for j,b in enumerate(_b):
        for k,c in enumerate(_b):
            GENCODE[a+b+c] = _aa[16*i+4*j+k]
FOURFOLD_PREFIX = {"CT","GT","TC","CC","AC","GC","CG","GG"}   # L V S P T A R G
PYEOF

# ----------------------------------------------------------------- STEP 1 ---
cat > "${SCRIPTS}/step1_composition.py" <<'PYEOF'
"""Per-gene composition for one species. Longest CDS transcript per gene.
Outputs <sp>_gene_composition.tsv and <sp>_longest_cds.fna (reused by step 3)."""
import sys, re, collections, math
sys.path.insert(0, sys.argv[5])
from common import *
sp, genome_path, gtf_path, outdir = sys.argv[1:5]

attr_tid = re.compile(r'transcript_id "([^"]+)"'); attr_gid = re.compile(r'gene_id "([^"]+)"')
attr_id  = re.compile(r'(?:^|;)\s*ID=([^;]+)'); attr_par = re.compile(r'(?:^|;)\s*Parent=([^;]+)')
cds = collections.defaultdict(list); t2g = {}; tinfo = {}
with open(gtf_path) as f:
    for line in f:
        if line.startswith("#"): continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9: continue
        seqid, ftype, s, e, strand, frame, attr = p[0], p[2], int(p[3]), int(p[4]), p[6], p[7], p[8]
        if ftype == "mRNA" or ftype == "transcript":
            m = attr_tid.search(attr); g = attr_gid.search(attr)
            if not m:
                m = attr_id.search(attr); g = attr_par.search(attr)
            if m and g: t2g[m.group(1)] = g.group(1)
            continue
        if ftype != "CDS": continue
        m = attr_tid.search(attr); g = attr_gid.search(attr)
        if not m:                                   # GFF3
            m = attr_par.search(attr); g = None
        if not m: continue
        tid = m.group(1).split(",")[0]
        if g: t2g[tid] = g.group(1)
        cds[tid].append((s, e, frame))
        tinfo[tid] = (seqid, strand)

# longest CDS per gene
best = {}
for tid, segs in cds.items():
    gid = t2g.get(tid, tid)
    L = sum(e - s + 1 for s, e, _ in segs)
    if gid not in best or L > best[gid][1]: best[gid] = (tid, L)

needed = {tinfo[t][0] for t, _ in best.values()}
genome = read_fasta(genome_path, keep=needed)

def gc(s):
    a = s.count("A") + s.count("T"); g = s.count("G") + s.count("C")
    return (g, a + g)

def enc(codons):
    fam = collections.defaultdict(collections.Counter)
    for c in codons:
        aa = GENCODE.get(c)
        if aa and aa not in "*MW": fam[aa][c] += 1
    deg = collections.Counter(GENCODE[c] for c in GENCODE if GENCODE[c] != "*")
    F = collections.defaultdict(list)
    for aa, cnt in fam.items():
        n = sum(cnt.values())
        if n < 2: continue
        F[deg[aa]].append((n * sum((v / n) ** 2 for v in cnt.values()) - 1) / (n - 1))
    Fm = {k: sum(v) / len(v) for k, v in F.items() if v and sum(v) > 0}
    if 2 not in Fm or 4 not in Fm: return float("nan")
    if 3 not in Fm: Fm[3] = (Fm[2] + Fm[4]) / 2
    if 6 not in Fm: Fm[6] = Fm[4]
    return min(61.0, 2 + 9 / Fm[2] + 1 / Fm[3] + 5 / Fm[4] + 3 / Fm[6])

out = open(f"{outdir}/{sp}_gene_composition.tsv", "w")
fa  = open(f"{outdir}/{sp}_longest_cds.fna", "w")
pa  = open(f"{outdir}/proteins/{sp}.faa", "w")
out.write("species\tgene_id\ttranscript_id\tchromosome\tstart\tend\tstrand\tcds_len\tn_cds_exons\t"
          "n_introns\tintron_len_total\tintron_len_mean\tgc_cds\tgc1\tgc2\tgc3\tgc4fold\tn_4fold\t"
          "gc_intron\tenc\tinternal_stop\n")
n_ok = 0
for gid, (tid, L) in best.items():
    seqid, strand = tinfo[tid]
    chrom = genome.get(seqid)
    if chrom is None: continue
    segs = sorted(cds[tid], key=lambda x: x[0])
    ordered = segs if strand == "+" else segs[::-1]
    pieces = [chrom[s - 1:e] for s, e, _ in segs]
    seq = "".join(pieces)
    if strand == "-": seq = revcomp(seq)
    fr = ordered[0][2]
    if fr in ("1", "2"): seq = seq[int(fr):]
    seq = seq[: len(seq) - len(seq) % 3]
    codons = [seq[i:i + 3] for i in range(0, len(seq), 3)]
    if codons and GENCODE.get(codons[-1]) == "*": codons = codons[:-1]
    internal_stop = int(any(GENCODE.get(c) == "*" for c in codons))
    pos = ["".join(c[k] for c in codons if len(c) == 3) for k in range(3)]
    g1, g2, g3 = (gc(x) for x in pos)
    four = [c[2] for c in codons if c[:2] in FOURFOLD_PREFIX and "N" not in c]
    g4 = gc("".join(four))
    introns = [(segs[i][1] + 1, segs[i + 1][0] - 1) for i in range(len(segs) - 1)
               if segs[i + 1][0] - segs[i][1] > 1]
    iseq = "".join(chrom[s - 1:e] for s, e in introns)
    gi = gc(iseq); gall = gc(seq)
    ilen = [e - s + 1 for s, e in introns]
    r = lambda t: f"{t[0]/t[1]:.5f}" if t[1] else "NA"
    out.write(f"{sp}\t{gid}\t{tid}\t{seqid}\t{segs[0][0]}\t{segs[-1][1]}\t{strand}\t{len(seq)}\t"
              f"{len(segs)}\t{len(introns)}\t{sum(ilen)}\t"
              f"{(sum(ilen)/len(ilen)) if ilen else 0:.1f}\t{r(gall)}\t{r(g1)}\t{r(g2)}\t{r(g3)}\t"
              f"{r(g4)}\t{len(four)}\t{r(gi)}\t{enc(codons):.3f}\t{internal_stop}\n")
    fa.write(f">{gid}\n{''.join(codons)}\n")
    prot = "".join(GENCODE.get(c, "X") for c in codons).replace("*", "X")
    if len(prot) >= 30: pa.write(f">{gid}\n{prot}\n")
    n_ok += 1
out.close(); fa.close(); pa.close()
print(f"{sp}: {len(best)} genes with CDS, {n_ok} written")
PYEOF

# ---------------------------------------------------------- STEP 3 prepare ---
cat > "${SCRIPTS}/step3_prepare.py" <<'PYEOF'
"""Build the Step 3 gene set and write one CDS multi-FASTA per locus.
geneset=orthofinder: OrthoFinder single-copy orthogroups (one gene in each of the 8).
geneset=busco      : genome-mode BUSCO hits mapped to the annotated gene they overlap.
Either way a locus is kept only if, in all 8 species, the gene has no internal stop,
a CDS >= 300 bp, and sits on the same chromosome class (micro, macro or Z);
class-switching loci are excluded."""
import sys, os, glob, collections, pandas as pd
sys.path.insert(0, sys.argv[1])
from common import *
scripts, input_dir, out_dir, n_chunks, geneset = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5]
species = list(CODES)
chrom = pd.read_csv(f"{input_dir}/{os.environ.get('CHROM_TABLE', 'per_chromosome_summary_v3.csv')}")
cls = {(r.species, r.chromosome): r["class"] for _, r in chrom.iterrows()}
size = {(r.species, r.chromosome): r.log10_Mb for _, r in chrom.iterrows()}
genes = {sp: pd.read_csv(f"{out_dir}/step1/{sp}_gene_composition.tsv", sep="\t").set_index("gene_id")
         for sp in species}
cds = {sp: read_fasta(f"{out_dir}/step1/{sp}_longest_cds.fna") for sp in species}

loci = {}          # locus -> {species: gene_id}
qc_rows = []
if geneset == "orthofinder":
    # -fg resumes write a second Results_* folder; use the one holding the orthogroups
    res = sorted(glob.glob(f"{out_dir}/step2/of/Results_*/Orthogroups/Orthogroups_SingleCopyOrthologues.txt"))[0]
    res = res.rsplit("/Orthogroups/", 1)[0]
    sco = {l.strip() for l in open(f"{res}/Orthogroups/Orthogroups_SingleCopyOrthologues.txt") if l.strip()}
    og = pd.read_csv(f"{res}/Orthogroups/Orthogroups.tsv", sep="\t", dtype=str).set_index("Orthogroup")
    for o in sorted(sco):
        loci[o] = {sp: og.at[o, sp].strip() for sp in species}
else:
    for sp in species:
        g = genes[sp].reset_index(); bychr = {c: x for c, x in g.groupby("chromosome")}
        for line in open(f"{input_dir}/{sp}_full_table.tsv"):
            if line.startswith("#"): continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 8 or p[1] != "Complete": continue
            bid, seqid = p[0], p[2]; s, e = sorted((int(p[3]), int(p[4]))); aa = int(p[7])
            x = bychr.get(seqid); status, gid = "no_gene", None
            if x is not None:
                ov = (x[["start", "end"]].clip(lower=s, upper=e).diff(axis=1)["end"] + 1)
                ov = ov.where((x.end >= s) & (x.start <= e), 0)
                if ov.max() > 0:
                    gid = x.at[ov.idxmax(), "gene_id"]
                    r = x.at[ov.idxmax(), "cds_len"] / (3 * aa)
                    status = "ok" if 0.7 <= r <= 1.3 else "length_mismatch"
            qc_rows.append(dict(species=sp, busco=bid, gene_id=gid, status=status))
            if status == "ok": loci.setdefault(bid, {})[sp] = gid
    loci = {k: v for k, v in loci.items() if len(v) == len(species)}
    pd.DataFrame(qc_rows).to_csv(f"{out_dir}/results/step3_busco_mapping_qc.tsv", sep="\t", index=False)

keep, drop = [], collections.Counter()
used = collections.Counter((sp, gid) for v in loci.values() for sp, gid in v.items())
for loc, v in sorted(loci.items()):
    if any(used[(sp, gid)] > 1 for sp, gid in v.items()): drop["gene_in_2_loci"] += 1; continue
    rows = [genes[sp].loc[gid] if gid in genes[sp].index else None for sp, gid in v.items()]
    if any(r is None for r in rows): drop["gene_not_in_step1"] += 1; continue
    if any(r.internal_stop == 1 for r in rows): drop["internal_stop"] += 1; continue
    if any(r.cds_len < 300 for r in rows): drop["cds_lt_300bp"] += 1; continue
    L = [r.cds_len for r in rows]
    if max(L) / min(L) > 1.5: drop["cds_length_ratio_gt_1.5"] += 1; continue
    classes = {cls.get((sp, r.chromosome), "unplaced") for sp, r in zip(v, rows)}
    if len(classes) != 1: drop["class_switching"] += 1; continue
    c = classes.pop()
    if c not in ("macro", "micro", "Z"): drop[f"class_{c}"] += 1; continue
    keep.append((loc, c, sum(size[(sp, r.chromosome)] for sp, r in zip(v, rows)) / len(rows), v))

gdir = f"{out_dir}/step3/genes"
with open(f"{out_dir}/step3/gene_list.tsv", "w") as f:
    f.write("locus\tclass\tmean_log10_Mb\t" + "\t".join(species) + "\n")
    for loc, c, ls, v in keep:
        f.write(f"{loc}\t{c}\t{ls:.4f}\t" + "\t".join(v[sp] for sp in species) + "\n")
        with open(f"{gdir}/{loc}.fna", "w") as fa:
            for sp in species: fa.write(f">{CODES[sp]}\n{cds[sp][v[sp]]}\n")
for i in range(n_chunks):
    with open(f"{out_dir}/step3/chunks/chunk_{i+1}.txt", "w") as f:
        f.write("".join(f"{k[0]}\n" for j, k in enumerate(keep) if j % n_chunks == i))
summary = (f"gene set: {geneset}\ncandidate loci (one gene in each of the 8): {len(loci)}\n"
           f"kept: {len(keep)}  [macro {sum(k[1]=='macro' for k in keep)}, "
           f"micro {sum(k[1]=='micro' for k in keep)}, Z {sum(k[1]=='Z' for k in keep)}]\n"
           f"dropped: {dict(drop)}\n")
open(f"{out_dir}/results/step3_gene_set_summary.txt", "w").write(summary); print(summary)
PYEOF

# ------------------------------------------------------ STEP 3 per chunk ---
cat > "${SCRIPTS}/step3_chunk.sh" <<'SHEOF'
#!/bin/bash
# args: chunk_file out_dir tree env_bin_dir aligner   (binaries called directly)
# The gene list is read on file descriptor 3 and every program gets </dev/null:
# (defensive: no program in the loop can consume the gene list)
set -uo pipefail
CHUNK="$1"; OUT="$2"; TREE="$3"; BIN="$4"; ALIGNER="${5:-macse}"
export PATH="${BIN}:${PATH}"
HERE="$(cd "$(dirname "$0")" && pwd)"
RES="${OUT}/step3/results_$(basename "$CHUNK" .txt).tsv"
[[ -s "$RES" ]] || echo -e "locus\tls\tlnL\tkappa\tomega\tdN_tree\tdS_tree\tstatus" > "$RES"
DONE=$(mktemp); cat "${OUT}"/step3/results_chunk_*.tsv 2>/dev/null | cut -f1 | grep -v '^locus$' > "$DONE"
while read -r B <&3; do
  [[ -z "$B" ]] && continue
  grep -qxF "$B" "$DONE" && continue      # already done in any results file (chunks may have been re-split)
  W="${OUT}/step3/work/${B}"; mkdir -p "$W"; cd "$W" || continue
  # 1. codon alignment
  rm -f nt.fasta aa.fasta nt_clean.fasta aa_clean.fasta aln.phy mlc
  if [[ "$ALIGNER" == "mafft" ]]; then     # translate -> MAFFT (protein) -> back-translate
    "${BIN}/python3" "${HERE}/codon_align.py" "${OUT}/step3/genes/${B}.fna" "${BIN}/mafft" nt_clean.fasta > align.log 2>&1 </dev/null
  else                                     # MACSE v2
    "${BIN}/macse" -prog alignSequences -seq "${OUT}/step3/genes/${B}.fna" \
        -out_NT nt.fasta -out_AA aa.fasta > align.log 2>&1 </dev/null
    "${BIN}/macse" -prog exportAlignment -align nt.fasta -codonForFinalStop --- \
        -codonForInternalStop NNN -codonForExternalFS --- -codonForInternalFS --- \
        -charForRemainingFS - -out_NT nt_clean.fasta -out_AA aa_clean.fasta >> align.log 2>&1 </dev/null
  fi
  if [[ ! -s nt_clean.fasta ]]; then echo -e "${B}\tNA\tNA\tNA\tNA\tNA\tNA\talign_failed" >> "$RES"; continue; fi
  # 2. PHYLIP + tree + control file
  "${BIN}/python3" - <<'PY'
seqs={};n=None
for l in open("nt_clean.fasta"):
    l=l.strip()
    if l.startswith(">"): n=l[1:]; seqs[n]=[]
    elif n: seqs[n].append(l.replace("!","-").upper())
seqs={k:"".join(v) for k,v in seqs.items()}
L=len(next(iter(seqs.values())))
with open("aln.phy","w") as f:
    f.write(f" {len(seqs)} {L}\n")
    for k,v in seqs.items(): f.write(f"{k}  {v}\n")
PY
  echo "$TREE" > tree.nwk
  cat > codeml.ctl <<CTL
seqfile = aln.phy
treefile = tree.nwk
outfile = mlc
noisy = 0
verbose = 0
runmode = 0
seqtype = 1
CodonFreq = 2
clock = 0
model = 0
NSsites = 0
icode = 0
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 0.4
cleandata = 1
getSE = 0
RateAncestor = 0
Small_Diff = .5e-6
method = 0
CTL
  "${BIN}/codeml" codeml.ctl > codeml.log 2>&1 </dev/null
  "${BIN}/python3" - "$B" >> "$RES" <<'PY'
import re,sys
b=sys.argv[1]
try: t=open("mlc").read()
except: print(f"{b}\tNA\tNA\tNA\tNA\tNA\tNA\tcodeml_failed"); sys.exit()
g=lambda p: (re.search(p,t).group(1) if re.search(p,t) else "NA")
ls=g(r"ls\s*=\s*(\d+)"); lnl=g(r"lnL\(.*?\):\s*(-?[\d.]+)")
kap=g(r"kappa \(ts/tv\)\s*=\s*([\d.]+)"); om=g(r"omega \(dN/dS\)\s*=\s*([\d.]+)")
dn=g(r"tree length for dN:\s*([\d.]+)"); ds=g(r"tree length for dS:\s*([\d.]+)")
st="ok" if "NA" not in (om,dn,ds) else "parse_failed"
print(f"{b}\t{ls}\t{lnl}\t{kap}\t{om}\t{dn}\t{ds}\t{st}")
PY
  rm -f "${W}/rst" "${W}/rst1" "${W}/rub" "${W}/2NG"* "${W}/4fold.nuc" "${W}/lnf"
done 3< "$CHUNK"
rm -f "$DONE"
SHEOF

cat > "${SCRIPTS}/codon_align.py" <<'PYEOF'
# translate CDS -> MAFFT protein alignment -> back-translate to codons
# (replaces MACSE: same principle -- gaps only between codons -- but ~100x faster)
import sys, subprocess, tempfile, os
fna, mafft, out = sys.argv[1:4]
code = {}
b = "TCAG"; aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
i = 0
for x in b:
    for y in b:
        for z in b: code[x+y+z] = aa[i]; i += 1
seqs = {}; n = None
for l in open(fna):
    l = l.strip()
    if l.startswith(">"): n = l[1:].split()[0]; seqs[n] = []
    elif n: seqs[n].append(l.upper())
codons, prot = {}, {}
for k, v in seqs.items():
    s = "".join(v); s = s[: len(s) // 3 * 3]
    cs = [s[i:i+3] for i in range(0, len(s), 3)]
    if cs and code.get(cs[-1]) == "*": cs = cs[:-1]          # terminal stop
    p = []
    for j, c in enumerate(cs):
        a = code.get(c, "X")
        if a == "*": a = "X"; cs[j] = "NNN"                   # internal stop / frameshift -> masked
        p.append(a)
    codons[k], prot[k] = cs, "".join(p)
if len(prot) < 4 or min(len(p) for p in prot.values()) < 50: sys.exit("too few / too short sequences")
with tempfile.NamedTemporaryFile("w", suffix=".faa", delete=False) as f:
    for k, p in prot.items(): f.write(f">{k}\n{p}\n")
    tmp = f.name
r = subprocess.run([mafft, "--auto", "--quiet", "--anysymbol", tmp], capture_output=True, text=True)
os.unlink(tmp)
aln = {}; n = None
for l in r.stdout.splitlines():
    if l.startswith(">"): n = l[1:].split()[0]; aln[n] = []
    elif n: aln[n].append(l.strip())
if set(aln) != set(prot): sys.exit("mafft failed: " + r.stderr[-500:])
with open(out, "w") as f:
    for k, a in aln.items():
        a = "".join(a).upper(); cs = iter(codons[k]); nt = []
        for ch in a: nt.append("---" if ch == "-" else next(cs))
        f.write(f">{k}\n{''.join(nt)}\n")
PYEOF

# ---------------------------------------------------------------- COLLECT ---
cat > "${SCRIPTS}/collect.py" <<'PYEOF'
"""Merge, test and plot Steps 1 and 3. Writes everything to OUT/results."""
import sys, os, glob, numpy as np, pandas as pd, warnings
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
out, input_dir = sys.argv[1], sys.argv[2]; R = f"{out}/results"
COL = {"Coronella_austriaca":"#8A00FF","Coronella_girondica":"#B19CD9","Hierophis_viridiflavus":"#E3D208",
       "Macroprotodon_brevis":"#FB8728","Macroprotodon_mauritanicus":"#FFA100","Natrix_astreptophora":"#94FF00",
       "Zamenis_longissimus":"#B23A48","Zamenis_scalaris":"#FF4880"}
chrom = pd.read_csv(f"{input_dir}/{os.environ.get('CHROM_TABLE', 'per_chromosome_summary_v3.csv')}")[["species","chromosome","class","log10_Mb"]]
lines = []

# ================================ STEP 1 ===================================
g = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob(f"{out}/step1/*_gene_composition.tsv"))])
g = g.merge(chrom, on=["species","chromosome"], how="left")
g["class"] = g["class"].fillna("unplaced")
g = g[(g.internal_stop == 0) & (g.cds_len >= 300)]
g.to_csv(f"{R}/step1_gene_composition_all.tsv", sep="\t", index=False)
vars1 = ["gc3","gc4fold","gc_intron","gc1","gc2","enc","cds_len","n_introns","intron_len_mean"]
auto = g[g["class"].isin(["macro","micro"])]
# chromosome means: the unit of replication (genes on a chromosome are not independent)
cm = (g[g["class"].isin(["macro","micro","Z"])]
      .groupby(["species","chromosome","class","log10_Mb"])[vars1].mean().reset_index())
cm["n_genes"] = g.groupby(["species","chromosome"]).size().reindex(
    list(zip(cm.species, cm.chromosome))).values
cm.to_csv(f"{R}/step1_chromosome_means.tsv", sep="\t", index=False)

rows = []
for v in vars1:
    a = cm[cm["class"].isin(["macro","micro"])]
    ma, mi = a[a["class"]=="macro"][v].dropna(), a[a["class"]=="micro"][v].dropna()
    u, p = stats.mannwhitneyu(ma, mi)
    rho, prho = stats.spearmanr(a.log10_Mb, a[v], nan_policy="omit")
    nsp = sum(1 for sp, s in a.groupby("species")
              if np.sign(s[s["class"]=="micro"][v].mean()-s[s["class"]=="macro"][v].mean())
              == np.sign(mi.mean()-ma.mean()))
    try:
        import statsmodels.formula.api as smf
        d = auto[["species","chromosome","log10_Mb",v]].dropna()
        m = smf.mixedlm(f"{v} ~ log10_Mb", d, groups="species",
                        vc_formula={"chrom":"0+C(chromosome)"}).fit(reml=True)
        slope, pslope = m.params["log10_Mb"], m.pvalues["log10_Mb"]
    except Exception:
        slope, pslope = np.nan, np.nan
    rows.append(dict(variable=v, macro_chrom_mean=ma.mean(), micro_chrom_mean=mi.mean(),
        mannwhitney_U_chrom=u, p_chrom=p, spearman_rho_vs_log10size=rho, p_spearman=prho,
        species_same_direction=f"{nsp}/8", mixed_slope_per_log10Mb=slope, p_mixed=pslope,
        n_macro_chrom=len(ma), n_micro_chrom=len(mi)))
st1 = pd.DataFrame(rows); st1.to_csv(f"{R}/step1_stats.tsv", sep="\t", index=False)
lines.append("STEP 1 (chromosome means, autosomes; mixed model on genes with species and chromosome random effects)")
lines.append(st1.round(5).to_string(index=False))

fig, ax = plt.subplots(1, 3, figsize=(12, 3.8))
for a_, v, lab in zip(ax, ["gc3","gc4fold","gc_intron"], ["GC3","GC at fourfold sites","Intron GC"]):
    for sp, s in cm[cm["class"].isin(["macro","micro"])].groupby("species"):
        a_.scatter(s.log10_Mb, s[v]*100, s=18, c=COL[sp], edgecolors="#333", linewidths=.3, label=sp.replace("_"," "))
    a_.axvline(np.log10(70), ls="--", c="#888", lw=.8); a_.set_xlabel("Chromosome length (log10 Mb)")
    a_.set_ylabel(f"{lab} (%)"); [a_.spines[k].set_visible(False) for k in ("top","right")]
ax[0].legend(fontsize=6, frameon=False, prop={"style":"italic","size":6})
fig.tight_layout(); fig.savefig(f"{R}/step1_GC3_GC4_intronGC_vs_size.png", dpi=250); plt.close(fig)

fig, ax = plt.subplots(1, 3, figsize=(12, 3.8))
for a_, v, lab in zip(ax, ["enc","cds_len","intron_len_mean"], ["ENC","CDS length (bp)","Mean intron length (bp)"]):
    for sp, s in cm[cm["class"].isin(["macro","micro"])].groupby("species"):
        a_.scatter(s.log10_Mb, s[v], s=18, c=COL[sp], edgecolors="#333", linewidths=.3)
    a_.axvline(np.log10(70), ls="--", c="#888", lw=.8); a_.set_xlabel("Chromosome length (log10 Mb)")
    a_.set_ylabel(lab); [a_.spines[k].set_visible(False) for k in ("top","right")]
fig.tight_layout(); fig.savefig(f"{R}/step1_ENC_genestructure_vs_size.png", dpi=250); plt.close(fig)

# ================================ STEP 2 ===================================
ogf = sorted(glob.glob(f"{out}/step2/of/Results_*/Orthogroups/Orthogroups.tsv"))
if ogf:
    res = ogf[0].rsplit("/Orthogroups/", 1)[0]; SP = list(COL)
    dupf = sorted(glob.glob(f"{out}/step2/of/Results_*/Gene_Duplication_Events/Duplications.tsv"))
    og = pd.read_csv(f"{res}/Orthogroups/Orthogroups.tsv", sep="\t", dtype=str)
    sco = {l.strip() for l in open(f"{res}/Orthogroups/Orthogroups_SingleCopyOrthologues.txt") if l.strip()}
    memb = {}                                   # (species, gene) -> orthogroup
    nsp = {}
    for _, r in og.iterrows():
        present = 0
        for sp in SP:
            if isinstance(r.get(sp), str) and r[sp].strip():
                present += 1
                for gname in r[sp].split(", "): memb[(sp, gname.strip())] = r.Orthogroup
        nsp[r.Orthogroup] = present
    dupg = set()
    dfile = dupf[0] if dupf else None
    if dfile is None:
        lines.append("STEP 2: gene trees not finished -> duplication fraction not available yet")
    if dfile and os.path.exists(dfile):
        dd = pd.read_csv(dfile, sep="\t", dtype=str)
        dd = dd[pd.to_numeric(dd["Support"], errors="coerce") >= 0.5]
        known = {(sp, gname) for sp, gname in memb}
        for col in ("Genes 1", "Genes 2"):
            for cell in dd[col].dropna():
                for tok in cell.split(", "):
                    tok = tok.strip()
                    for sp in SP:
                        if tok.startswith(sp + "_") and (sp, tok[len(sp)+1:]) in known:
                            dupg.add((sp, tok[len(sp)+1:])); break
                    else:
                        for sp in SP:
                            if (sp, tok) in known: dupg.add((sp, tok)); break
    g2 = g[["species","gene_id","chromosome","class","log10_Mb"]].copy()
    g2["orthogroup"] = [memb.get((s, x)) for s, x in zip(g2.species, g2.gene_id)]
    g2["og_n_species"] = g2.orthogroup.map(nsp).fillna(1)
    g2["core_all8"] = g2.og_n_species == 8
    g2["single_copy_all8"] = g2.orthogroup.isin(sco)
    g2["lineage_specific"] = g2.og_n_species <= 1
    g2["duplicated"] = [(s, x) in dupg for s, x in zip(g2.species, g2.gene_id)]
    g2.to_csv(f"{R}/step2_gene_orthology.tsv", sep="\t", index=False)
    v2 = ["core_all8","single_copy_all8","lineage_specific","duplicated"]
    c2 = (g2[g2["class"].isin(["macro","micro","Z"])]
          .groupby(["species","chromosome","class","log10_Mb"])[v2].mean().reset_index())
    c2.to_csv(f"{R}/step2_chromosome_fractions.tsv", sep="\t", index=False)
    rows = []
    a = c2[c2["class"].isin(["macro","micro"])]
    for v in v2:
        ma, mi = a[a["class"]=="macro"][v], a[a["class"]=="micro"][v]
        u, p = stats.mannwhitneyu(ma, mi); rho = stats.spearmanr(a.log10_Mb, a[v])
        rows.append(dict(variable=v, macro_chrom_mean=ma.mean(), micro_chrom_mean=mi.mean(),
                         Z_chrom_mean=c2[c2["class"]=="Z"][v].mean(), p_chrom=p,
                         spearman_rho_vs_log10size=rho.correlation, p_spearman=rho.pvalue))
    st2 = pd.DataFrame(rows); st2.to_csv(f"{R}/step2_stats.tsv", sep="\t", index=False)
    lines.append(f"\nSTEP 2 (OrthoFinder): {len(og)} orthogroups, {len(sco)} single-copy in all 8; "
                 f"fractions of genes per chromosome, autosomes")
    lines.append(st2.round(5).to_string(index=False))
    fig, ax = plt.subplots(1, 4, figsize=(15, 3.6))
    for a_, v, lab in zip(ax, v2, ["In orthogroups shared by all 8","Single-copy in all 8",
                                    "Lineage-specific","In a duplication event"]):
        for sp, s in a.groupby("species"):
            a_.scatter(s.log10_Mb, s[v]*100, s=18, c=COL[sp], edgecolors="#333", linewidths=.3)
        a_.axvline(np.log10(70), ls="--", c="#888", lw=.8); a_.set_xlabel("Chromosome length (log10 Mb)")
        a_.set_ylabel(f"{lab} (% genes)"); [a_.spines[k].set_visible(False) for k in ("top","right")]
    fig.tight_layout(); fig.savefig(f"{R}/step2_orthology_vs_size.png", dpi=250); plt.close(fig)
else:
    lines.append("\nSTEP 2: no OrthoFinder results found")

# ================================ STEP 3 ===================================
files = sorted(glob.glob(f"{out}/step3/results_chunk_*.tsv"))
if files:
    d = pd.concat([pd.read_csv(f, sep="\t") for f in files]).drop_duplicates("locus")
    gl = pd.read_csv(f"{out}/step3/gene_list.tsv", sep="\t")
    d = d.merge(gl[["locus","class","mean_log10_Mb"]], on="locus", how="left")
    d.to_csv(f"{R}/step3_dnds_per_gene_raw.tsv", sep="\t", index=False)
    ok = d[(d.status=="ok")].copy()
    for c in ["omega","dN_tree","dS_tree","ls"]: ok[c] = pd.to_numeric(ok[c], errors="coerce")
    n0 = len(ok)
    ok = ok[(ok.ls >= 100) & (ok.dS_tree > 0.001) & (ok.dS_tree < 2) & (ok.omega < 5)]
    ok.to_csv(f"{R}/step3_dnds_per_gene_filtered.tsv", sep="\t", index=False)
    lines.append(f"\nSTEP 3: {len(d)} genes run, {n0} ok, {len(ok)} after filters "
                 f"(filters: >=100 codons after cleandata; 0.001 < dS < 2; omega < 5)")
    rng = np.random.default_rng(1); rows = []
    for c, s in ok.groupby("class"):
        r = s.dN_tree.sum() / s.dS_tree.sum()
        bs = [s.sample(len(s), replace=True, random_state=int(rng.integers(1e9))) for _ in range(1000)]
        bs = np.array([b.dN_tree.sum()/b.dS_tree.sum() for b in bs])
        rows.append(dict(cls=c, n_genes=len(s), omega_ratio_of_sums=r, ci95_low=np.percentile(bs,2.5),
            ci95_high=np.percentile(bs,97.5), median_gene_omega=s.omega.median(),
            median_dN=s.dN_tree.median(), median_dS=s.dS_tree.median()))
    st3 = pd.DataFrame(rows)
    kw = stats.kruskal(*[s.omega for _, s in ok.groupby("class")])
    pair = {}
    for a1, a2 in [("macro","micro"),("macro","Z"),("micro","Z")]:
        x, y = ok[ok["class"]==a1], ok[ok["class"]==a2]
        if len(x) and len(y):
            pair[f"{a1}_vs_{a2}_omega_p"] = stats.mannwhitneyu(x.omega, y.omega).pvalue
            pair[f"{a1}_vs_{a2}_dS_p"] = stats.mannwhitneyu(x.dS_tree, y.dS_tree).pvalue
    au = ok[ok["class"].isin(["macro","micro"])]
    rho = stats.spearmanr(au.mean_log10_Mb, au.omega)
    st3.to_csv(f"{R}/step3_omega_by_class.tsv", sep="\t", index=False)
    pd.Series({"kruskal_omega_H": kw.statistic, "kruskal_omega_p": kw.pvalue,
               "spearman_omega_vs_log10size_rho": rho.correlation,
               "spearman_omega_vs_log10size_p": rho.pvalue, **pair}).to_csv(f"{R}/step3_tests.tsv", sep="\t")
    lines.append(st3.round(5).to_string(index=False))
    lines.append(f"Kruskal-Wallis omega across classes: H={kw.statistic:.2f}, p={kw.pvalue:.2e}")
    for k, v in pair.items(): lines.append(f"{k}: {v:.2e}")
    lines.append(f"Spearman omega vs log10 size (autosomes): rho={rho.correlation:.3f}, p={rho.pvalue:.2e}")

    fig, ax = plt.subplots(1, 3, figsize=(11, 3.6))
    order = [c for c in ["micro","macro","Z"] if c in set(ok["class"])]
    for a_, v, lab in zip(ax, ["omega","dN_tree","dS_tree"], ["omega (dN/dS)","dN (tree length)","dS (tree length)"]):
        a_.boxplot([ok[ok["class"]==c][v] for c in order], showfliers=False); a_.set_xticks(range(1, len(order)+1), order)
        a_.set_ylabel(lab); [a_.spines[k].set_visible(False) for k in ("top","right")]
    fig.tight_layout(); fig.savefig(f"{R}/step3_omega_dN_dS_by_class.png", dpi=250); plt.close(fig)
    fig, ax = plt.subplots(figsize=(4.5, 3.6))
    for c, col in [("micro","#555555"),("macro","#BBBBBB"),("Z","#1f5fbf")]:
        s = ok[ok["class"]==c]; ax.scatter(s.mean_log10_Mb, s.omega, s=4, c=col, alpha=.4, label=c)
    ax.set_yscale("log"); ax.set_xlabel("Mean chromosome length (log10 Mb)"); ax.set_ylabel("omega (log scale)")
    ax.legend(frameon=False, fontsize=7); fig.tight_layout()
    fig.savefig(f"{R}/step3_omega_vs_size.png", dpi=250); plt.close(fig)
else:
    lines.append("\nSTEP 3: no results yet")
open(f"{R}/SUMMARY.txt", "w").write("\n".join(lines) + "\n"); print("\n".join(lines))
PYEOF
}

write_sbatch() {
  # Step 1: one array task per species
  cat > "${SCRIPTS}/step1.sbatch" <<EOF
#!/bin/bash
#SBATCH --job-name=mm_step1
#SBATCH --account=<account>
#SBATCH --partition=<partition>
#SBATCH --array=0-$(( ${#SPECIES[@]} - 1 ))
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --output=${LOGS}/step1_%A_%a.out
module purge; module load ${MODULE_MAMBA}
unset CONDA_ENVS_PATH; export CONDA_ENVS_DIRS=${CONDA_ENVS_DIRS}; export CONDA_PKGS_DIRS=${CONDA_PKGS_DIRS}
SP=(${SPECIES[*]}); sp=\${SP[\$SLURM_ARRAY_TASK_ID]}
[[ -s ${OUT_DIR}/step1/\${sp}_gene_composition.tsv ]] && { echo "\$sp done"; exit 0; }
${RUN} python3 ${SCRIPTS}/step1_composition.py \$sp ${GENOME_DIR}/\$sp/\$sp.masked.fa \
  ${ANNOT_DIR}/\$sp.gtf ${OUT_DIR}/step1 ${SCRIPTS}
EOF
  echo "${TREE_ROOTED}" > "${SCRIPTS}/species_tree_rooted.nwk"
  cat > "${SCRIPTS}/step2.sbatch" <<EOF
#!/bin/bash
#SBATCH --job-name=mm_orthofinder
#SBATCH --account=<account>
#SBATCH --partition=<partition>
#SBATCH --cpus-per-task=${OF_CPUS}
#SBATCH --mem=${OF_MEM}
#SBATCH --time=${OF_TIME}
#SBATCH --output=${LOGS}/step2_orthofinder_%j.out
module purge; module load ${MODULE_MAMBA}
unset CONDA_ENVS_PATH; export CONDA_ENVS_DIRS=${CONDA_ENVS_DIRS}; export CONDA_PKGS_DIRS=${CONDA_PKGS_DIRS}
ulimit -n 20000 2>/dev/null || true
if ls -d ${OUT_DIR}/step2/of/Results_*/Gene_Duplication_Events >/dev/null 2>&1; then
  echo "OrthoFinder complete"; exit 0
elif [[ -s ${OUT_DIR}/step2/of/Results_mm/Orthogroups/Orthogroups.tsv ]]; then
  echo "Orthogroups exist -> resuming from them (gene trees, duplications only)"
  mamba run -p ${OF_ENV_PREFIX} orthofinder -fg ${OUT_DIR}/step2/of/Results_mm \
    -s ${SCRIPTS}/species_tree_rooted.nwk -t ${OF_CPUS} -a ${OF_ANALYSIS_THREADS} -n trees
else
  rm -rf ${OUT_DIR}/step2/of
  mamba run -p ${OF_ENV_PREFIX} orthofinder -f ${OUT_DIR}/step1/proteins \
    -s ${SCRIPTS}/species_tree_rooted.nwk -t ${OF_CPUS} -a ${OF_ANALYSIS_THREADS} -n mm -o ${OUT_DIR}/step2/of
fi
EOF
  cat > "${SCRIPTS}/step3_prepare.sbatch" <<EOF
#!/bin/bash
#SBATCH --job-name=mm_s3prep
#SBATCH --account=<account>
#SBATCH --partition=<partition>
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=${LOGS}/step3_prepare_%j.out
module purge; module load ${MODULE_MAMBA}
unset CONDA_ENVS_PATH; export CONDA_ENVS_DIRS=${CONDA_ENVS_DIRS}; export CONDA_PKGS_DIRS=${CONDA_PKGS_DIRS}
export CHROM_TABLE=${CHROM_TABLE}
${RUN} python3 ${SCRIPTS}/step3_prepare.py ${SCRIPTS} ${INPUT_DIR} ${OUT_DIR} ${N_CHUNKS} ${GENESET}
EOF
  cat > "${SCRIPTS}/step3_dnds.sbatch" <<EOF
#!/bin/bash
#SBATCH --job-name=mm_dnds
#SBATCH --account=<account>
#SBATCH --partition=<partition>
#SBATCH --array=1-${N_CHUNKS}%${MAX_PARALLEL}
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=${DNDS_TIME}
#SBATCH --output=${LOGS}/step3_dnds_%A_%a.out
module purge; module load ${MODULE_MAMBA}
unset CONDA_ENVS_PATH; export CONDA_ENVS_DIRS=${CONDA_ENVS_DIRS}; export CONDA_PKGS_DIRS=${CONDA_PKGS_DIRS}
bash ${SCRIPTS}/step3_chunk.sh ${OUT_DIR}/step3/chunks/chunk_\${SLURM_ARRAY_TASK_ID}.txt \
  ${OUT_DIR} "${TREE}" "${ENV_PREFIX}/bin" ${ALIGNER}
EOF
  cat > "${SCRIPTS}/collect.sbatch" <<EOF
#!/bin/bash
#SBATCH --job-name=mm_collect
#SBATCH --account=<account>
#SBATCH --partition=<partition>
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=${LOGS}/collect_%j.out
module purge; module load ${MODULE_MAMBA}
unset CONDA_ENVS_PATH; export CONDA_ENVS_DIRS=${CONDA_ENVS_DIRS}; export CONDA_PKGS_DIRS=${CONDA_PKGS_DIRS}
export CHROM_TABLE=${CHROM_TABLE}
${RUN} python3 ${SCRIPTS}/collect.py ${OUT_DIR} ${INPUT_DIR}
EOF
}

# =============================================================================
case "${1:-}" in
# -----------------------------------------------------------------------------
setup)
  cat > "${SCRIPTS}/setup.sbatch" <<EOF
#!/bin/bash
#SBATCH --job-name=mm_setup
#SBATCH --account=<account>
#SBATCH --partition=<partition>
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=${LOGS}/setup_%j.out
module purge; module load ${MODULE_MAMBA}
unset CONDA_ENVS_PATH; export CONDA_ENVS_DIRS=${CONDA_ENVS_DIRS}; export CONDA_PKGS_DIRS=${CONDA_PKGS_DIRS}
if [[ ! -x ${ENV_PREFIX}/bin/python3 ]]; then
mamba create -y -p ${ENV_PREFIX} -c conda-forge -c bioconda \
  python=3.11 pandas numpy scipy statsmodels matplotlib macse mafft paml
fi
mamba run -p ${ENV_PREFIX} python3 -c "import pandas, scipy, statsmodels; print('python OK')"
mamba run -p ${ENV_PREFIX} which macse mafft codeml
if [[ ! -x ${OF_ENV_PREFIX}/bin/orthofinder ]]; then
  mamba create -y -p ${OF_ENV_PREFIX} -c conda-forge -c bioconda orthofinder=2.5.5
fi
mamba run -p ${OF_ENV_PREFIX} orthofinder -h | head -3
EOF
  sbatch "${SCRIPTS}/setup.sbatch" ;;
# -----------------------------------------------------------------------------
check)
  write_helpers
  module purge >/dev/null 2>&1; module load ${MODULE_MAMBA} >/dev/null 2>&1
  printf "%-28s %-8s %-6s %-12s\n" SPECIES genome gtf busco_table
  for sp in "${SPECIES[@]}"; do
    g=--; a=--; b=--
    [[ -s "${GENOME_DIR}/${sp}/${sp}.masked.fa" ]] && g=ok
    [[ -s "${ANNOT_DIR}/${sp}.gtf" ]] && a=ok
    [[ -s "${INPUT_DIR}/${sp}_full_table.tsv" ]] && b=ok
    printf "%-28s %-8s %-6s %-12s\n" "$sp" "$g" "$a" "$b"
  done
  [[ -s "${INPUT_DIR}/${CHROM_TABLE}" ]] && echo "chromosome table: ok" || echo "chromosome table: MISSING"
  [[ -x "${ENV_PREFIX}/bin/python3" ]] && echo "conda env: ok" || echo "conda env: MISSING (run setup)"
  ls ${ENV_PREFIX}/bin/${ALIGNER} ${ENV_PREFIX}/bin/codeml 2>/dev/null || echo "${ALIGNER}/codeml: not found in env"
  [[ -x "${OF_ENV_PREFIX}/bin/orthofinder" ]] && echo "orthofinder env: ok" || echo "orthofinder env: MISSING (run setup)"
  echo "GTF chromosome names vs table (first species):"
  sp=${SPECIES[0]}; cut -f1 "${ANNOT_DIR}/${sp}.gtf" 2>/dev/null | grep -v "^#" | sort -u | head -5
  ;;
# -----------------------------------------------------------------------------
submit)
  write_helpers
  write_sbatch
  j1=$(sbatch --parsable "${SCRIPTS}/step1.sbatch")
  j2=$(sbatch --parsable --dependency=afterok:${j1} "${SCRIPTS}/step2.sbatch")
  # Step 3 needs only the orthogroups (written early by OrthoFinder), not the gene
  # trees, so it starts when OrthoFinder ends even if the tree stage fails.
  if [[ "$GENESET" == "orthofinder" ]]; then dep="afterany:${j2}"; else dep="afterok:${j1}"; fi
  j3=$(sbatch --parsable --dependency=${dep} "${SCRIPTS}/step3_prepare.sbatch")
  j4=$(sbatch --parsable --dependency=afterok:${j3} "${SCRIPTS}/step3_dnds.sbatch")
  j5=$(sbatch --parsable --dependency=afterany:${j4}:${j2} "${SCRIPTS}/collect.sbatch")
  log "submitted: step1=${j1}  orthofinder=${j2}  step3_prepare=${j3}  step3_dnds=${j4}  collect=${j5}"
  log "results will appear in ${OUT_DIR}/results/  (SUMMARY.txt first)" ;;
# -----------------------------------------------------------------------------
resume)      # after a failure: resubmit only what is missing, in the right order
  write_helpers; write_sbatch
  deps=""
  if ! ls -d ${OUT_DIR}/step2/of/Results_*/Gene_Duplication_Events >/dev/null 2>&1; then
    j2=$(sbatch --parsable "${SCRIPTS}/step2.sbatch"); deps="${j2}"
    log "OrthoFinder (resume from orthogroups): ${j2}"
  fi
  if [[ ! -s ${OUT_DIR}/step3/gene_list.tsv ]]; then
    j3=$(sbatch --parsable "${SCRIPTS}/step3_prepare.sbatch")
    j4=$(sbatch --parsable --dependency=afterok:${j3} "${SCRIPTS}/step3_dnds.sbatch") || { log "dN/dS submission failed"; exit 1; }
  else
    j4=$(sbatch --parsable "${SCRIPTS}/step3_dnds.sbatch") || { log "dN/dS submission failed"; exit 1; }
  fi
  log "step3 dN/dS: ${j4}"
  dep="afterany:${j4}"; [[ -n "$deps" ]] && dep="${dep}:${deps}"
  j5=$(sbatch --parsable --dependency=${dep} "${SCRIPTS}/collect.sbatch")
  log "collect: ${j5}  (results -> ${OUT_DIR}/results/)" ;;
# -----------------------------------------------------------------------------
step3)       # submit only step 3 (prepare -> dN/dS -> collect). Optional argument:
             # a job ID that collect should also wait for (e.g. the OrthoFinder job)
  write_helpers; write_sbatch
  module purge >/dev/null 2>&1; module load ${MODULE_MAMBA} >/dev/null 2>&1
  if [[ ! -x ${ENV_PREFIX}/bin/${ALIGNER} ]]; then
    log "installing ${ALIGNER} into ${ENV_PREFIX} ..."
    mamba install -y -p ${ENV_PREFIX} -c conda-forge -c bioconda ${ALIGNER} || { log "${ALIGNER} install failed"; exit 1; }
  fi
  # one aligner for every gene: results made with another aligner are set aside
  # (no marker file = the first run, which used MACSE)
  prev="$(cat ${OUT_DIR}/step3/ALIGNER 2>/dev/null || echo macse)"
  if [[ "$prev" != "$ALIGNER" ]] && ls ${OUT_DIR}/step3/results_chunk_*.tsv >/dev/null 2>&1; then
    old=${OUT_DIR}/step3/old_${prev}_$(date +%Y%m%d_%H%M); mkdir -p "$old"
    mv ${OUT_DIR}/step3/results_chunk_*.tsv "$old"/; log "previous (${prev}) results moved to $old"
  fi
  echo ${ALIGNER} > ${OUT_DIR}/step3/ALIGNER
  # the gene list must be split into exactly N_CHUNKS chunk files (the first run made
  # 200 chunks while the array ran only 40, so 80% of the genes were never submitted)
  if [[ -s ${OUT_DIR}/step3/gene_list.tsv ]]; then
    ng=$(( $(wc -l < ${OUT_DIR}/step3/gene_list.tsv) - 1 ))
    nc=$(ls ${OUT_DIR}/step3/chunks/chunk_*.txt 2>/dev/null | wc -l)
    inch=$(cat ${OUT_DIR}/step3/chunks/chunk_*.txt 2>/dev/null | grep -c .)
    if [[ $nc -ne $N_CHUNKS || $inch -ne $ng ]]; then
      rm -f ${OUT_DIR}/step3/chunks/chunk_*.txt
      tail -n +2 ${OUT_DIR}/step3/gene_list.tsv | cut -f1 \
        | awk -v n=${N_CHUNKS} -v d=${OUT_DIR}/step3/chunks '{ print > (d "/chunk_" ((NR-1) % n + 1) ".txt") }'
      log "re-split ${ng} genes into ${N_CHUNKS} chunks (were ${nc} chunk files holding ${inch} genes)"
    fi
  fi
  if [[ ! -s ${OUT_DIR}/step3/gene_list.tsv ]]; then
    j3=$(sbatch --parsable "${SCRIPTS}/step3_prepare.sbatch") || exit 1
    j4=$(sbatch --parsable --dependency=afterok:${j3} "${SCRIPTS}/step3_dnds.sbatch") || exit 1
  else
    j4=$(sbatch --parsable "${SCRIPTS}/step3_dnds.sbatch") || exit 1
  fi
  dep="afterany:${j4}"; [[ -n "${2:-}" ]] && dep="${dep}:${2}"
  j5=$(sbatch --parsable --dependency=${dep} "${SCRIPTS}/collect.sbatch") || exit 1
  log "step3_prepare=${j3:-skipped}  dnds=${j4}  collect=${j5}" ;;
# -----------------------------------------------------------------------------
collect)     # re-run only the statistics/plots, e.g. after a partial step 3
  write_helpers; sbatch "${SCRIPTS}/collect.sbatch" ;;
# -----------------------------------------------------------------------------
status)
  echo "step1 tables : $(ls ${OUT_DIR}/step1/*_gene_composition.tsv 2>/dev/null | wc -l) / ${#SPECIES[@]}"
  echo "orthofinder  : $(ls ${OUT_DIR}/step2/of/Results_*/Orthogroups/Orthogroups.tsv >/dev/null 2>&1 && echo done || echo 'not yet')"
  echo "step3 genes  : $(($(cat ${OUT_DIR}/step3/gene_list.tsv 2>/dev/null | wc -l) - 1)) prepared"
  echo "step3 done   : $(cat ${OUT_DIR}/step3/results_chunk_*.tsv 2>/dev/null | grep -vc '^locus')  (ok: $(cat ${OUT_DIR}/step3/results_chunk_*.tsv 2>/dev/null | grep -c 'ok$'))"
  squeue -u "$USER" -o "%.12i %.12j %.8T %.10M" 2>/dev/null | grep -E "mm_|JOBID" ;;
*)
  sed -n '2,24p' "$SELF" ;;
esac
