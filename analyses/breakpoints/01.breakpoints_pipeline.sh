#!/bin/bash
# =============================================================================
# breakpoints_pipeline.sh -- chromosome breakpoints across the 15 Iberian snakes
#
# Same tools and parameters as the macrosynteny analysis (DIAMOND blastp,
# e-value 1e-5, 5 targets -> JCVI synteny scan), reusing its longest-isoform
# inputs, but every species is compared against two fixed REFERENCES instead of
# its neighbour in the figure. Each species' chromosomes then become runs of
# reference-chromosome blocks ("painting"); a candidate breakpoint is where two
# consecutive blocks come from different reference chromosomes.
#
#   Stage A  paint     DIAMOND + JCVI scan, 15 species x 2 references
#   Stage B  features  100 kb windows from each soft-masked genome: GC, repeat
#                      fraction, assembly gaps, telomeric repeats, gene count
#   Stage C  analyze   blocks, breakpoints, shared regions placed on the tree,
#                      end-to-end joins, matched-control environment tests, plots
#
#   bash breakpoints_pipeline.sh check     # preflight
#   bash breakpoints_pipeline.sh submit    # A and B in parallel, then C
#   bash breakpoints_pipeline.sh analyze   # re-run only C (e.g. new parameters)
#   bash breakpoints_pipeline.sh status
#
# Resumable. Final tables and plots: ${OUT_DIR}/results/
# =============================================================================
set -uo pipefail

# ============================== CONFIG ======================================
MS_OUT="/path/to/working_directory/Macrosynteny_results"   # <Sp>/<Sp>.longest.bed/.faa/.fa.fai
GENOME_DIR="/path/to/working_directory/data/Genomes"        # <Sp>/<Sp>.masked.fa (soft-masked)
INPUT_DIR="/path/to/working_directory/Breakpoints_inputs"   # per-species input tables (see check)
OUT_DIR="/path/to/working_directory/Breakpoints_results"

SPECIES=(Vipera_aspis Vipera_seoanei Vipera_latastei Malpolon_monspessulanus
         Natrix_maura Natrix_astreptophora Natrix_helvetica Hemorrhois_hippocrepis
         Hierophis_viridiflavus Macroprotodon_brevis Macroprotodon_mauritanicus
         Zamenis_scalaris Zamenis_longissimus Coronella_girondica Coronella_austriaca)
# Two references: the most contiguous colubrid (male, 2n=36) and a viper with the
# standard 2n=36 karyotype. Every other species is painted against both.
REFS=(Zamenis_longissimus Vipera_seoanei)

# Analysis parameters (stage C only; change and run 'analyze' again)
MIN_RUN=15          # genes a block needs to count (filters paralogue noise)
MIN_CHROM_MB=5      # scaffolds smaller than this are ignored
MATCH_TOL_MB=3      # breaks from different species closer than this = same region
TERMINAL_MB=3       # a break within this distance of a reference end = end-to-end join
WINDOW=100000       # feature window
FLANK_MB=2          # environment profiles span +/- this around each breakpoint

DIAMOND_EVALUE=1e-5
DIAMOND_MAX_TARGETS=5
THREADS=32
ACCOUNT="<account>"
PARTITION="<partition>"

MODULE_MAMBA="miniconda3"   # module providing mamba/conda on your system
CONDA_ROOT="/path/to/conda"
unset CONDA_ENVS_PATH
export CONDA_ENVS_DIRS="${CONDA_ROOT}/envs"
export CONDA_PKGS_DIRS="${CONDA_ROOT}/pkgs"
JCVI_ENV="${CONDA_ENVS_DIRS}/JCVI"           # from the macrosynteny pipeline (diamond + jcvi)
PY_ENV="${CONDA_ENVS_DIRS}/MICROMACRO"       # from the macro/micro pipeline (pandas, matplotlib)
# =============================================================================

SCRIPTS="${OUT_DIR}/scripts"; LOGS="${OUT_DIR}/logs"
mkdir -p "$SCRIPTS" "$LOGS" "${OUT_DIR}/pairs" "${OUT_DIR}/features" "${OUT_DIR}/results"
log() { echo "[$(date '+%F %T')] $*"; }
ENVLINES="module purge; module load ${MODULE_MAMBA}
unset CONDA_ENVS_PATH; export CONDA_ENVS_DIRS=${CONDA_ENVS_DIRS}; export CONDA_PKGS_DIRS=${CONDA_PKGS_DIRS}"

write_helpers() {
cat > "${SCRIPTS}/features.py" <<'PYEOF'
"""Per-species 100 kb window features straight from the soft-masked genome:
GC, repeat fraction (lower-case bases), assembly gap fraction (N), telomeric
repeat density (TTAGGG/CCCTAA per kb), gene count. Also writes gap BED."""
import sys, re, collections
sp, genome, bed, out, W = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5])
genes = collections.defaultdict(list)
for l in open(bed):
    p = l.split("\t")
    if len(p) >= 3: genes[p[0]].append((int(p[1]) + int(p[2])) // 2)
tel = re.compile(r"(?=(TTAGGG|CCCTAA))")
fo = open(f"{out}/{sp}.windows.tsv", "w")
fo.write("species\tchrom\tstart\tend\tgc\trepeat\tgap\ttelomeric_per_kb\tgenes\tsoftmask_seen\n")
fg = open(f"{out}/{sp}.gaps.bed", "w")
def flush(name, seq):
    if not seq: return
    s = "".join(seq); L = len(s)
    for m in re.finditer(r"[Nn]{10,}", s): fg.write(f"{name}\t{m.start()}\t{m.end()}\n")
    gpos = sorted(genes.get(name, [])); gi = 0
    for a in range(0, L, W):
        w = s[a:a + W]; n = w.count("N") + w.count("n"); acgt = len(w) - n
        lower = sum(w.count(c) for c in "acgt")
        up = w.upper(); gc = up.count("G") + up.count("C")
        t = len(tel.findall(up))
        while gi < len(gpos) and gpos[gi] < a: gi += 1
        gj = gi
        while gj < len(gpos) and gpos[gj] < a + W: gj += 1
        fo.write(f"{sp}\t{name}\t{a}\t{a+len(w)}\t{(gc/acgt) if acgt else float('nan'):.4f}\t"
                 f"{(lower/acgt) if acgt else float('nan'):.4f}\t{n/len(w):.4f}\t"
                 f"{t/(len(w)/1000):.3f}\t{gj-gi}\t{int(lower>0)}\n")
name, seq = None, []
with open(genome) as f:
    for l in f:
        if l.startswith(">"):
            flush(name, seq); name, seq = l[1:].split()[0], []
        else: seq.append(l.strip())
flush(name, seq); fo.close(); fg.close()
print(sp, "done")
PYEOF
cat > "${SCRIPTS}/analyze.py" <<'PYEOF'
"""Breakpoint calling, event placement and plots from gene-anchor paintings.

For every species Q and reference R:
  anchors (Q gene, R gene) -> each Q chromosome becomes a run of R-chromosome
  blocks -> a BREAKPOINT is where two consecutive blocks come from different
  R chromosomes. Each breakpoint is also projected onto R: the R positions of
  the last gene before and first gene after it. A projection near an R
  chromosome end is 'terminal' (end-to-end join, e.g. a micro fusion); otherwise
  it is an 'internal' break of that R chromosome.
  Internal breaks from different species that fall within MATCH_TOL_MB on the
  same R chromosome are clustered -> shared breakpoint regions. The set of
  species sharing a region is compared with the clades of the species tree:
  one clade = one event on that clade's stem; several clades = convergent
  candidates (independent breaks of the same R region).
Outputs tables + figures in OUT/results."""
import sys, os, glob, json, collections, itertools
import numpy as np, pandas as pd
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
mpl.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Helvetica", "TeX Gyre Heros", "Liberation Sans", "DejaVu Sans"],
                     "font.size": 7, "pdf.fonttype": 42})

cfg = json.loads(sys.argv[1])
OUT, MS, INP = cfg["out"], cfg["ms_out"], cfg["inputs"]
REFS, SPECIES = cfg["refs"], cfg["species"]
MIN_RUN, MIN_CHROM = cfg["min_run"], cfg["min_chrom_mb"] * 1e6
TOL, TERM = cfg["match_tol_mb"] * 1e6, cfg["terminal_mb"] * 1e6
FLANK = cfg.get("flank_mb", 2) * 1e6
R = f"{OUT}/results"; os.makedirs(R, exist_ok=True)

def readlist(f):
    d = {}
    if os.path.exists(f):
        for l in open(f):
            if l.startswith("#") or not l.strip(): continue
            p = l.split(); d[p[0]] = p[1].split(",") if len(p) > 1 else []
    return d
ZL, WL = readlist(f"{INP}/z_chromosome.txt"), readlist(f"{INP}/w_chromosome.txt")
COL = {}
for l in open(f"{INP}/Species_colours.txt"):
    p = l.split()
    if len(p) >= 2: COL[p[0]] = p[1].strip('"')
TREE = open(f"{INP}/species_tree.nwk").read().strip()

# ---------------------------------------------------------------- tree clades
def clades(nwk):
    s = nwk.rstrip(";"); stack, out, i, cur = [], [], 0, ""
    leaves_stack = []
    while i < len(s):
        c = s[i]
        if c == "(":
            leaves_stack.append([])
        elif c in ",)":
            if cur.strip():
                name = cur.split(":")[0].strip()
                if name and leaves_stack: leaves_stack[-1].append(name)
            cur = ""
            if c == ")":
                grp = leaves_stack.pop(); out.append(frozenset(grp))
                if leaves_stack: leaves_stack[-1].extend(grp)
                j = i + 1
                while j < len(s) and s[j] not in ",)": j += 1
                i = j - 1
        else:
            cur += c
        i += 1
    leaves = max(out, key=len)
    return [frozenset([x]) for x in leaves] + out, leaves
CLADES, ALL = clades(TREE)

def place(shared, universe):
    """Smallest set of clades (restricted to 'universe' = species analysed vs this
    reference) whose union is exactly 'shared'."""
    # Full clades only: the reference belongs to its clade but by definition does not
    # show a break against itself, so a clade containing it cannot explain the pattern.
    shared = frozenset(shared)
    cands = sorted({c for c in CLADES if c <= shared}, key=len, reverse=True)
    picked, left = [], set(shared)
    for c in cands:
        if c <= left:
            picked.append(c); left -= c
        if not left: break
    return picked

def clade_name(c):
    c = sorted(c)
    if len(c) == 1: return c[0]
    gen = sorted({x.split("_")[0] for x in c})
    return (gen[0] if len(gen) == 1 else "+".join(gen)) + f" ({len(c)} spp)"

# --------------------------------------------------------- genomes / chromosomes
def chrom_info(sp):
    bed = pd.read_csv(f"{MS}/{sp}/{sp}.longest.bed", sep="\t", header=None,
                      usecols=[0, 1, 2, 3], names=["chrom", "start", "end", "gene"]).drop_duplicates("gene")
    bed["mid"] = (bed.start + bed.end) // 2
    fai = f"{MS}/{sp}/{sp}.fa.fai"
    if os.path.exists(fai):
        L = pd.read_csv(fai, sep="\t", header=None, usecols=[0, 1], index_col=0).iloc[:, 0]
    else:
        L = bed.groupby("chrom").end.max()
    gb = set(bed.chrom)
    keep = [c for c in L.sort_values(ascending=False).index if L[c] >= MIN_CHROM and c in gb]
    z = [c for c in ZL.get(sp, []) if c in keep]; w = [c for c in WL.get(sp, []) if c in keep]
    auto = [c for c in keep if c not in z and c not in w]
    lab = {c: str(i + 1) for i, c in enumerate(auto)}
    for i, c in enumerate(z): lab[c] = "Z" if i == 0 else f"Z{i+1}"
    for i, c in enumerate(w): lab[c] = "W" if i == 0 else f"W{i+1}"
    order = auto + z            # W excluded: degenerate, its Z-gametolog genes paint as false breaks
    return bed[bed.chrom.isin(order)].set_index("gene"), order, lab, L

INFO = {sp: chrom_info(sp) for sp in SPECIES}
def label(sp, c): return INFO[sp][2].get(c, c)
def chlen(sp, c): return int(INFO[sp][3][c])

def anchors(q, r):
    f = f"{OUT}/pairs/{q}.{r}.anchors"
    rows = [l.split()[:2] for l in open(f) if not l.startswith("#") and l.strip()]
    return pd.DataFrame(rows, columns=["gq", "gr"])

# ------------------------------------------------------------------ painting
def paint(q, r):
    bq, oq, _, _ = INFO[q]; br, orr, _, _ = INFO[r]
    a = anchors(q, r)
    a = a[a.gq.isin(bq.index) & a.gr.isin(br.index)].drop_duplicates("gq")
    a["qchrom"] = bq.loc[a.gq, "chrom"].values; a["qpos"] = bq.loc[a.gq, "mid"].values
    a["rchrom"] = br.loc[a.gr, "chrom"].values; a["rpos"] = br.loc[a.gr, "mid"].values
    return a

def blocks_of(g):
    """g: anchors on one query chromosome sorted by qpos -> list of blocks."""
    labs = g.rchrom.tolist(); runs = []
    for k, lab in enumerate(labs):
        if runs and runs[-1][0] == lab: runs[-1][2] = k
        else: runs.append([lab, k, k])
    runs = [x for x in runs if x[2] - x[1] + 1 >= MIN_RUN]
    merged = []
    for x in runs:
        if merged and merged[-1][0] == x[0]: merged[-1][2] = x[2]
        else: merged.append(list(x))
    out = []
    for lab, i, j in merged:
        sub = g.iloc[i:j + 1]
        rp = sub.rpos.values
        orient = "+" if np.corrcoef(sub.qpos, rp)[0, 1] >= 0 else "-" if len(sub) > 2 else "?"
        out.append(dict(rchrom=lab, n_genes=j - i + 1, q_start=int(sub.qpos.iloc[0]), q_end=int(sub.qpos.iloc[-1]),
                        r_min=int(rp.min()), r_max=int(rp.max()), r_first=int(rp[0]), r_last=int(rp[-1]),
                        orient=orient))
    return out

# --------------------------------------------------------------- features
FEAT = {}
for sp in SPECIES:
    f = f"{OUT}/features/{sp}.windows.tsv"
    if os.path.exists(f): FEAT[sp] = pd.read_csv(f, sep="\t")
GAPS = {}
for sp in SPECIES:
    f = f"{OUT}/features/{sp}.gaps.bed"
    if os.path.exists(f) and os.path.getsize(f):
        GAPS[sp] = pd.read_csv(f, sep="\t", header=None, names=["chrom", "start", "end"])

def gap_in(sp, c, s, e):
    g = GAPS.get(sp)
    if g is None: return np.nan
    g = g[(g.chrom == c) & (g.end >= s) & (g.start <= e)]
    return int(len(g))

# ---------------------------------------------------------- main loop
blocks_rows, bp_rows, paint_cache = [], [], {}
for r in REFS:
    for q in SPECIES:
        if q == r: continue
        if not os.path.exists(f"{OUT}/pairs/{q}.{r}.anchors"): continue
        a = paint(q, r); paint_cache[(q, r)] = a
        for qc in INFO[q][1]:
            g = a[a.qchrom == qc].sort_values("qpos")
            bl = blocks_of(g)
            for k, b in enumerate(bl):
                blocks_rows.append(dict(reference=r, species=q, q_chrom=qc, q_label=label(q, qc), block=k + 1,
                                        r_chrom=b["rchrom"], r_label=label(r, b["rchrom"]), **b))
            for b1, b2 in zip(bl, bl[1:]):
                s, e = b1["q_end"], b2["q_start"]
                ends = []
                for side, b, pos in (("left", b1, b1["r_last"]), ("right", b2, b2["r_first"])):
                    Lr = chlen(r, b["rchrom"])
                    term = "start" if pos <= TERM else "end" if Lr - pos <= TERM else "internal"
                    ends.append((b["rchrom"], pos, term))
                Lq = chlen(q, qc)
                bp_rows.append(dict(reference=r, species=q, q_chrom=qc, q_label=label(q, qc),
                    q_start=s, q_end=e, interval_kb=round((e - s) / 1e3, 1),
                    dist_to_q_end_mb=round(min(s, Lq - e) / 1e6, 2), q_chrom_mb=round(Lq / 1e6, 1),
                    q_class="micro" if Lq < 70e6 else "macro",
                    left_r_chrom=ends[0][0], left_r_label=label(r, ends[0][0]), left_r_pos=ends[0][1], left_type=ends[0][2],
                    right_r_chrom=ends[1][0], right_r_label=label(r, ends[1][0]), right_r_pos=ends[1][1], right_type=ends[1][2],
                    genes_left=b1["n_genes"], genes_right=b2["n_genes"],
                    gaps_in_interval=gap_in(q, qc, s, e),
                    small_scaffold_involved=int(Lq < 10e6 or chlen(r, ends[0][0]) < 10e6 or chlen(r, ends[1][0]) < 10e6)))
BL = pd.DataFrame(blocks_rows); BP = pd.DataFrame(bp_rows)
BL.to_csv(f"{R}/blocks.tsv", sep="\t", index=False)

# ------------------------------------------------- splits (fissions seen from the reference side)
# A fission in species Q leaves each daughter chromosome as one clean block, so it
# produces no junction INSIDE any Q chromosome. It shows up instead as one reference
# chromosome whose blocks are spread over two Q chromosomes. Walk each reference
# chromosome in reference order: consecutive blocks on different Q chromosomes = split.
# v1.2: walk the ANCHORS of each reference chromosome in reference order (not the
# blocks: a block's r_min/r_max span is stretched by stray anchors, which put v1.1
# split positions at meaningless chromosome midpoints).
def runs_of(labels):
    runs = []
    for k, lab in enumerate(labels):
        if runs and runs[-1][0] == lab: runs[-1][2] = k
        else: runs.append([lab, k, k])
    runs = [x for x in runs if x[2] - x[1] + 1 >= MIN_RUN]
    merged = []
    for x in runs:
        if merged and merged[-1][0] == x[0]: merged[-1][2] = x[2]
        else: merged.append(list(x))
    return merged
split_rows = []
for (q, r), a in paint_cache.items():
    for rc, g in a.groupby("rchrom"):
        g = g.sort_values("rpos"); rp = g.rpos.values
        m = runs_of(g.qchrom.tolist())
        for (l1, i1, j1), (l2, i2, j2) in zip(m, m[1:]):
            le, rs = int(rp[j1]), int(rp[i2]); pos = (le + rs) // 2; Lr = chlen(r, rc)
            if pos <= TERM or Lr - pos <= TERM: continue
            split_rows.append(dict(reference=r, species=q, r_chrom=rc, r_label=label(r, rc),
                r_left_end=le, r_right_start=rs, r_pos=pos, interval_kb=round((rs - le) / 1e3, 1),
                q_chrom_left=label(q, l1), q_chrom_right=label(q, l2),
                genes_left=j1 - i1 + 1, genes_right=j2 - i2 + 1,
                small_scaffold_involved=int(chlen(q, l1) < 10e6 or chlen(q, l2) < 10e6)))
SP_ = pd.DataFrame(split_rows)
if not len(SP_): SP_ = pd.DataFrame(columns=["reference","species","r_chrom","r_label","r_left_end","r_right_start","r_pos","interval_kb","q_chrom_left","q_chrom_right","genes_left","genes_right","small_scaffold_involved"])
SP_.to_csv(f"{R}/splits.tsv", sep="\t", index=False)

def interpret(sp, r):
    """single event / event in the reference's own clade / convergent candidate."""
    shared = frozenset(sp); pl = place(sp, frozenset(s for s in SPECIES if s != r))
    comp = frozenset(ALL) - shared
    if len(pl) == 1: return pl, f"single event: {clade_name(pl[0])}"
    if comp in set(CLADES):
        return pl, (f"single event in the reference lineage: {clade_name(comp)}" if r in comp
                    else f"single event: {clade_name(comp)}")
    return pl, "convergent candidate"

# ------------------------------------------------- shared breakpoint regions
# Each query breakpoint gives up to two reference projections. Internal ones are
# breaks of a reference chromosome; terminal ones are end-to-end joins.
proj = []
for i, b in BP[BP.small_scaffold_involved == 0].iterrows():
    for side in ("left", "right"):
        proj.append(dict(bp_index=i, reference=b.reference, species=b.species, side=side,
                         r_chrom=b[f"{side}_r_chrom"], r_label=b[f"{side}_r_label"],
                         r_pos=b[f"{side}_r_pos"], r_lo=b[f"{side}_r_pos"], r_hi=b[f"{side}_r_pos"],
                         kind=b[f"{side}_type"]))
for i, x in SP_[SP_.small_scaffold_involved == 0].iterrows():
    proj.append(dict(bp_index=f"split{i}", reference=x.reference, species=x.species, side="split",
                     r_chrom=x.r_chrom, r_label=x.r_label, r_pos=x.r_pos,
                     r_lo=x.r_left_end, r_hi=x.r_right_start, kind="internal"))
PJ = pd.DataFrame(proj)
if not len(PJ): PJ = pd.DataFrame(columns=["bp_index","reference","species","side","r_chrom","r_label","r_pos","r_lo","r_hi","kind"])
regions = []
for (r, rc), g in PJ[PJ.kind == "internal"].groupby(["reference", "r_chrom"]):
    # same region = uncertainty intervals overlap or lie within TOL of each other
    g = g.sort_values("r_lo"); cl, cur = [], [g.iloc[0]]; hi = g.iloc[0].r_hi
    for _, x in g.iloc[1:].iterrows():
        if x.r_lo - hi <= TOL: cur.append(x); hi = max(hi, x.r_hi)
        else: cl.append(cur); cur = [x]; hi = x.r_hi
    cl.append(cur)
    for c in cl:
        c = pd.DataFrame(c); sp = sorted(set(c.species))
        pl, interp = interpret(sp, r)
        regions.append(dict(reference=r, r_chrom=rc, r_label=label(r, rc),
            r_start=int(c.r_lo.min()), r_end=int(c.r_hi.max()), n_species=len(sp),
            species=",".join(sp), seen_as=",".join(sorted(set(c.side.replace({"left": "junction", "right": "junction"})))),
            placement=" | ".join(clade_name(x) for x in pl),
            n_independent_events_min=1 if not interp.startswith("convergent") else len(pl),
            interpretation=interp))
RG = pd.DataFrame(regions)
if len(RG):
    RG = RG.sort_values(["reference", "r_chrom", "r_start"]); RG["region_id"] = [f"R{i+1:03d}" for i in range(len(RG))]
RG.to_csv(f"{R}/shared_breakpoint_regions.tsv", sep="\t", index=False)
# terminal joins: which reference chromosome ends are glued together, by species
term = []
for i, b in BP[(BP.left_type != "internal") & (BP.right_type != "internal") & (BP.small_scaffold_involved == 0)].iterrows():
    key = tuple(sorted([f"{b.left_r_label}:{b.left_type}", f"{b.right_r_label}:{b.right_type}"]))
    term.append(dict(reference=b.reference, join=" + ".join(key), species=b.species, q_chrom=b.q_label))
TJ = pd.DataFrame(term)
if len(TJ):
    TJ = (TJ.groupby(["reference", "join"]).agg(species=("species", lambda x: ",".join(sorted(set(x)))),
                                                n_species=("species", "nunique")).reset_index())
    TJ["interpretation"] = [interpret(s.split(","), r)[1] for s, r in zip(TJ.species, TJ.reference)]
TJ.to_csv(f"{R}/terminal_joins.tsv", sep="\t", index=False)
BP.to_csv(f"{R}/breakpoints.tsv", sep="\t", index=False)
open(f"{R}/breakpoints_query_coords.bed", "w").write(
    "".join(f"{b.q_chrom}\t{int(b.q_start)}\t{int(b.q_end)}\t{b.species}|{b.reference}|{i}\n" for i, b in BP.iterrows()))

# ----------------------------------------------------------------- FIGURES
def short(sp): g, s = sp.split("_"); return f"{g[0]}. {s}"
def ref_palette(r):
    base = plt.get_cmap("tab20").colors + plt.get_cmap("tab20b").colors + plt.get_cmap("tab20c").colors
    return {c: base[i % len(base)] for i, c in enumerate(INFO[r][1])}

# 1. dot plots: small multiples per reference
for r in REFS:
    qs = [q for q in SPECIES if (q, r) in paint_cache]
    n = len(qs); nc = 4; nr = int(np.ceil(n / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(nc * 3.2, nr * 3.2)); axs = np.array(axs).ravel()
    _, orr, _, Lr = INFO[r]; offr = dict(zip(orr, np.cumsum([0] + [Lr[c] for c in orr[:-1]])))
    for ax, q in zip(axs, qs):
        _, oq, _, Lq = INFO[q]; offq = dict(zip(oq, np.cumsum([0] + [Lq[c] for c in oq[:-1]])))
        a = paint_cache[(q, r)]
        ax.scatter((a.rchrom.map(offr) + a.rpos) / 1e6, (a.qchrom.map(offq) + a.qpos) / 1e6,
                   s=.15, c="#222222", lw=0, rasterized=True)
        for c in orr: ax.axvline(offr[c] / 1e6, c="#dddddd", lw=.3)
        for c in oq: ax.axhline(offq[c] / 1e6, c="#dddddd", lw=.3)
        ax.set_xlim(0, sum(Lr[c] for c in orr) / 1e6); ax.set_ylim(0, sum(Lq[c] for c in oq) / 1e6)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_title(short(q), style="italic", fontsize=7, color=COL.get(q, "k"))
        sub = BP[(BP.species == q) & (BP.reference == r)]
        for _, b in sub.iterrows():
            ax.plot((offr[b.left_r_chrom] + b.left_r_pos) / 1e6, (offq[b.q_chrom] + (b.q_start + b.q_end) / 2) / 1e6,
                    marker="o", ms=2.5, mfc="none", mec="#d62728", mew=.6)
    for ax in axs[n:]: ax.axis("off")
    fig.suptitle(f"Gene-anchor dot plots against {short(r)} (x); red circles = candidate breakpoints",
                 fontsize=8, style="italic")
    fig.tight_layout(); fig.savefig(f"{R}/fig1_dotplots_vs_{r}.png", dpi=250); fig.savefig(f"{R}/fig1_dotplots_vs_{r}.pdf")
    plt.close(fig)

# 2. painted karyotypes: overview (one row per species) + one page per species
for r in REFS:
    pal = ref_palette(r); xmax = 0
    qs = [q for q in SPECIES if (q, r) in paint_cache]
    fig, ax = plt.subplots(figsize=(11, 0.45 * len(qs) + 1.2))
    for i, q in enumerate(qs):
        y = len(qs) - i; x0 = 0
        for qc in INFO[q][1]:
            L = chlen(q, qc) / 1e6
            ax.add_patch(plt.Rectangle((x0, y - .3), L, .6, fc="#f0f0f0", ec="#bbbbbb", lw=.3))
            sub = BL[(BL.reference == r) & (BL.species == q) & (BL.q_chrom == qc)]
            for _, b in sub.iterrows():
                ax.add_patch(plt.Rectangle((x0 + b.q_start / 1e6, y - .3), (b.q_end - b.q_start) / 1e6, .6,
                                           fc=pal[b.r_chrom], lw=0))
            x0 += L + 8
        ax.text(-14, y, short(q), ha="right", va="center", style="italic", fontsize=7)
        ax.plot(-8, y, marker="o", ms=5, c=COL.get(q, "#555"), mec="#333", mew=.3)
        xmax = max(xmax, x0)
    ax.set_xlim(-140, xmax); ax.set_ylim(.4, len(qs) + .8); ax.axis("off")
    hs = [plt.Rectangle((0, 0), 1, 1, fc=pal[c]) for c in INFO[r][1]]
    ax.legend(hs, [label(r, c) for c in INFO[r][1]], ncol=12, fontsize=5.5, frameon=False,
              loc="upper center", bbox_to_anchor=(.5, 1.06), title=f"{short(r)} chromosome", title_fontsize=6)
    fig.tight_layout(); fig.savefig(f"{R}/fig2_painted_karyotypes_vs_{r}.png", dpi=250)
    fig.savefig(f"{R}/fig2_painted_karyotypes_vs_{r}.pdf"); plt.close(fig)

    with PdfPages(f"{R}/fig2b_painting_per_species_vs_{r}.pdf") as pdf:
        for q in qs:
            oq = INFO[q][1]; a = paint_cache[(q, r)]
            fig, ax = plt.subplots(figsize=(9, 0.32 * len(oq) + 1))
            for i, qc in enumerate(oq):
                y = len(oq) - i; g = a[a.qchrom == qc]
                ax.add_patch(plt.Rectangle((0, y - .3), chlen(q, qc) / 1e6, .6, fc="#f4f4f4", ec="#999999", lw=.3))
                ax.scatter(g.qpos / 1e6, np.full(len(g), y), marker="|", s=25, lw=.5,
                           c=[pal.get(x, "#cccccc") for x in g.rchrom])
                ax.text(-2, y, label(q, qc), ha="right", va="center", fontsize=6)
                for _, b in BP[(BP.species == q) & (BP.reference == r) & (BP.q_chrom == qc)].iterrows():
                    ax.plot((b.q_start + b.q_end) / 2e6, y - .45, marker="^", ms=3.5, c="#d62728", mec="none")
                for _, b in BL[(BL.reference == r) & (BL.species == q) & (BL.q_chrom == qc)].iterrows():
                    ax.text((b.q_start + b.q_end) / 2e6, y + .33, b.r_label, ha="center", fontsize=4.5, color="#333333")
            ax.set_ylim(.3, len(oq) + .9); ax.set_yticks([]); ax.set_xlabel("Position (Mb)")
            for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
            ax.set_title(f"{short(q)} painted by {short(r)} chromosomes (labels = {short(r)} chromosome; "
                         f"red = candidate breakpoint)", fontsize=7, loc="left", style="italic")
            fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

# 3. breakpoint map in reference coordinates (recurrence)
for r in REFS:
    orr = INFO[r][1]; Lr = INFO[r][3]
    sub = PJ[(PJ.reference == r) & (PJ.kind == "internal")]
    fig, ax = plt.subplots(figsize=(10, 0.42 * len(orr) + 1.2))
    for i, c in enumerate(orr):
        y = len(orr) - i
        ax.add_patch(plt.Rectangle((0, y - .15), Lr[c] / 1e6, .3, fc="#e8e8e8", ec="#999999", lw=.3))
        ax.text(-3, y, label(r, c), ha="right", va="center", fontsize=6)
        g = sub[sub.r_chrom == c]
        for j, (pos, gg) in enumerate(g.groupby((g.r_pos // TOL).astype(int))):
            for k, (_, x) in enumerate(gg.drop_duplicates("species").iterrows()):
                ax.plot(x.r_pos / 1e6, y + .28 + k * .17, marker="o", ms=3.2, c=COL.get(x.species, "#555"),
                        mec="#333333", mew=.25)
        tj = PJ[(PJ.reference == r) & (PJ.kind != "internal") & (PJ.r_chrom == c)]
        for k, (_, x) in enumerate(tj.drop_duplicates(["species", "kind"]).iterrows()):
            xe = 0 if x.kind == "start" else Lr[c] / 1e6
            ax.plot(xe, y - .3 - k * .15, marker=">" if x.kind == "start" else "<", ms=3.2,
                    c=COL.get(x.species, "#555"), mec="#333333", mew=.25)
    ax.set_ylim(.2, len(orr) + 1.4); ax.set_yticks([]); ax.set_xlabel(f"Position on {short(r)} chromosome (Mb)")
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    used = sorted(set(PJ[PJ.reference == r].species))
    hs = [mpl.lines.Line2D([], [], marker="o", ls="", mfc=COL.get(s, "#555"), mec="#333", mew=.25, ms=4) for s in used]
    ax.legend(hs, [short(s) for s in used], fontsize=5.5, frameon=False, ncol=5, loc="upper center",
              bbox_to_anchor=(.5, 1.08), prop={"style": "italic", "size": 5.5})
    ax.set_title(f"Breaks projected onto {short(r)}: dots = internal breaks (stacked = shared region); triangles at ends = end-to-end joins",
                 fontsize=7.5, loc="left", y=1.08)
    fig.tight_layout(); fig.savefig(f"{R}/fig3_breakpoint_map_vs_{r}.png", dpi=250)
    fig.savefig(f"{R}/fig3_breakpoint_map_vs_{r}.pdf"); plt.close(fig)

# 4. genomic environment around breakpoints vs matched controls
envrows = []
if FEAT:
    rng = np.random.default_rng(1)
    feats = ["gc", "repeat", "genes", "telomeric_per_kb", "gap"]
    W = int(FEAT[next(iter(FEAT))].eval("end-start").max())
    offs = np.arange(-FLANK, FLANK + 1, W)
    prof = collections.defaultdict(list); ctrl = collections.defaultdict(list)
    seen = set()
    for _, b in BP[BP.small_scaffold_involved == 0].iterrows():
        key = (b.species, b.q_chrom, int((b.q_start + b.q_end) // 2 // 1e6))
        if key in seen or b.species not in FEAT: continue     # same break seen from both references
        seen.add(key)
        fw = FEAT[b.species]; fc = fw[fw.chrom == b.q_chrom].set_index("start")
        if fc.empty: continue
        centre = (b.q_start + b.q_end) // 2
        def take(tab, c):
            v = []
            for o in offs:
                st = int((c + o) // W * W)
                v.append(tab.loc[st, feats].values if st in tab.index else [np.nan] * len(feats))
            return np.array(v, dtype=float)
        prof[b.q_class].append(take(fc, centre))
        # matched controls: same species, same chromosome class, same distance-to-end bin
        same = fw[fw.chrom.isin([c for c in INFO[b.species][1]
                                 if (chlen(b.species, c) < 70e6) == (b.q_class == "micro")])]
        d_end = min(centre, chlen(b.species, b.q_chrom) - centre)
        for _ in range(50):
            row = same.sample(1, random_state=int(rng.integers(1e9))).iloc[0]
            Lc = chlen(b.species, row.chrom)
            if Lc < 2 * d_end + 2 * FLANK: continue
            pos = d_end if rng.random() < .5 else Lc - d_end
            ctrl[b.q_class].append(take(fw[fw.chrom == row.chrom].set_index("start"), pos))
        envrows.append(dict(species=b.species, q_chrom=b.q_label, centre=centre, q_class=b.q_class,
                            **{f: fc.loc[int(centre // W * W), f] if int(centre // W * W) in fc.index else np.nan
                               for f in feats}))
    ENV = pd.DataFrame(envrows); ENV.to_csv(f"{R}/breakpoint_window_features.tsv", sep="\t", index=False)
    labs = {"gc": "GC fraction", "repeat": "Repeat fraction (soft-masked)", "genes": "Genes per window",
            "telomeric_per_kb": "TTAGGG/CCCTAA per kb", "gap": "Assembly gap fraction"}
    fig, axs = plt.subplots(2, len(feats), figsize=(3 * len(feats), 5.4), sharex=True)
    stats_rows = []
    for row, cls in enumerate(["macro", "micro"]):
        for j, f in enumerate(feats):
            ax = axs[row, j]
            for data, col, lab in ((prof[cls], "#d62728", "breakpoints"), (ctrl[cls], "#777777", "matched controls")):
                if not data: continue
                m = np.array(data)[:, :, j]
                mu = np.nanmean(m, 0)
                bs = [np.nanmean(m[rng.integers(0, len(m), len(m))], 0) for _ in range(300)]
                lo, hi = np.nanpercentile(bs, 2.5, 0), np.nanpercentile(bs, 97.5, 0)
                ax.plot(offs / 1e6, mu, c=col, lw=1, label=f"{lab} (n={len(m)})")
                ax.fill_between(offs / 1e6, lo, hi, color=col, alpha=.2, lw=0)
            if prof[cls] and ctrl[cls]:
                ci = len(offs) // 2
                a_ = np.array(prof[cls])[:, ci, j]; b_ = np.array(ctrl[cls])[:, ci, j]
                a_, b_ = a_[~np.isnan(a_)], b_[~np.isnan(b_)]
                if len(a_) > 1 and len(b_) > 1:
                    obs = a_.mean() - b_.mean(); pool = np.r_[a_, b_]
                    perm = [rng.permutation(pool)[:len(a_)].mean() - b_.mean() for _ in range(2000)]
                    p = (np.sum(np.abs(perm) >= abs(obs)) + 1) / 2001
                    stats_rows.append(dict(chrom_class=cls, feature=f, n_breakpoints=len(a_), bp_mean=a_.mean(),
                                           control_mean=b_.mean(), difference=obs, permutation_p=p))
            ax.axvline(0, c="#333", lw=.5, ls=":")
            ax.set_title(f"{labs[f]} ({cls})", fontsize=7)
            if row == 1: ax.set_xlabel("Distance from breakpoint (Mb)")
            for s in ("top", "right"): ax.spines[s].set_visible(False)
        axs[row, 0].legend(fontsize=5.5, frameon=False)
    fig.tight_layout(); fig.savefig(f"{R}/fig4_breakpoint_environment.png", dpi=250)
    fig.savefig(f"{R}/fig4_breakpoint_environment.pdf"); plt.close(fig)
    pd.DataFrame(stats_rows).to_csv(f"{R}/breakpoint_environment_tests.tsv", sep="\t", index=False)

# ----------------------------------------------------------------- SUMMARY
lines = [f"Breakpoint pipeline summary  (MIN_RUN={MIN_RUN} genes, chromosomes >= {cfg['min_chrom_mb']} Mb, "
         f"shared-region tolerance {cfg['match_tol_mb']} Mb, terminal = within {cfg['terminal_mb']} Mb of an end)", ""]
for r in REFS:
    b = BP[BP.reference == r]
    lines.append(f"Reference {r}: {len(b)} candidate breakpoints in {b.species.nunique()} species "
                 f"({int(b.small_scaffold_involved.sum())} involve scaffolds < 10 Mb -> likely assembly fragments; "
                 f"{int((b.gaps_in_interval > 0).sum()) if 'gaps_in_interval' in b else 'NA'} contain an assembly gap)")
    sp_ = SP_[SP_.reference == r]
    qs_ = [q for q in SPECIES if q != r and (q, r) in paint_cache]
    ok_b, ok_s = b[b.small_scaffold_involved == 0], sp_[sp_.small_scaffold_involved == 0]
    tab = pd.DataFrame({
        "autosomes_ge5Mb": [sum(1 for c in INFO[q][1] if not label(q, c).startswith("Z")) for q in qs_],
        "junctions_inside_chromosomes": [int((ok_b.species == q).sum()) for q in qs_],
        "reference_chromosome_splits": [int((ok_s.species == q).sum()) for q in qs_],
        "fragment_events_excluded": [int(((b.species == q) & (b.small_scaffold_involved == 1)).sum()
                                         + ((sp_.species == q) & (sp_.small_scaffold_involved == 1)).sum()) for q in qs_]},
        index=pd.Index(qs_, name="species"))
    lines.append(tab.to_string()); lines.append("")
if len(RG):
    lines.append("Shared internal breakpoint regions (reference coordinates):")
    lines.append(RG[["region_id", "reference", "r_label", "r_start", "r_end", "n_species", "seen_as",
                     "interpretation"]].to_string(index=False))
if len(TJ):
    lines.append("\nEnd-to-end joins (e.g. microchromosome fusions):"); lines.append(TJ.to_string(index=False))
open(f"{R}/SUMMARY.txt", "w").write("\n".join(lines) + "\n"); print("\n".join(lines))
PYEOF
}

write_sbatch() {
  cat > "${SCRIPTS}/A_paint.sbatch" <<EOF
#!/bin/bash
#SBATCH --job-name=bp_paint
#SBATCH --account=<account>
#SBATCH --partition=<partition>
#SBATCH --cpus-per-task=${THREADS}
#SBATCH --mem=60G
#SBATCH --time=24:00:00
#SBATCH --output=${LOGS}/A_paint_%j.out
${ENVLINES}
for r in ${REFS[*]}; do
  db=${OUT_DIR}/pairs/\${r}.dmnd
  [[ -s \$db ]] || mamba run -p ${JCVI_ENV} diamond makedb --in ${MS_OUT}/\${r}/\${r}.longest.faa -d \${db%.dmnd}
  for q in ${SPECIES[*]}; do
    [[ \$q == \$r ]] && continue
    blast=${OUT_DIR}/pairs/\${q}.\${r}.blast; anch=${OUT_DIR}/pairs/\${q}.\${r}.anchors
    [[ -s \$anch ]] && { echo "\$q vs \$r done"; continue; }
    if [[ ! -s \$blast ]]; then
      mamba run -p ${JCVI_ENV} diamond blastp -q ${MS_OUT}/\${q}/\${q}.longest.faa -d \${db%.dmnd} \\
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \\
        --evalue ${DIAMOND_EVALUE} --max-target-seqs ${DIAMOND_MAX_TARGETS} --threads ${THREADS} > \$blast
    fi
    mamba run -p ${JCVI_ENV} python -m jcvi.compara.synteny scan \$blast \$anch \\
      --qbed=${MS_OUT}/\${q}/\${q}.longest.bed --sbed=${MS_OUT}/\${r}/\${r}.longest.bed --no_strip_names
    echo "\$q vs \$r: \$(grep -vc '^#' \$anch) anchor pairs"
  done
done
EOF
  cat > "${SCRIPTS}/B_features.sbatch" <<EOF
#!/bin/bash
#SBATCH --job-name=bp_features
#SBATCH --account=<account>
#SBATCH --partition=<partition>
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=${LOGS}/B_features_%j.out
${ENVLINES}
for sp in ${SPECIES[*]}; do
  [[ -s ${OUT_DIR}/features/\${sp}.windows.tsv ]] && { echo "\$sp done"; continue; }
  mamba run -p ${PY_ENV} python3 ${SCRIPTS}/features.py \$sp ${GENOME_DIR}/\${sp}/\${sp}.masked.fa \\
    ${MS_OUT}/\${sp}/\${sp}.longest.bed ${OUT_DIR}/features ${WINDOW}
done
EOF
  CFG=$(printf '{"out":"%s","ms_out":"%s","inputs":"%s","refs":[%s],"species":[%s],"min_run":%s,"min_chrom_mb":%s,"match_tol_mb":%s,"terminal_mb":%s,"flank_mb":%s}' \
    "$OUT_DIR" "$MS_OUT" "$INPUT_DIR" "$(printf '"%s",' "${REFS[@]}" | sed 's/,$//')" \
    "$(printf '"%s",' "${SPECIES[@]}" | sed 's/,$//')" "$MIN_RUN" "$MIN_CHROM_MB" "$MATCH_TOL_MB" "$TERMINAL_MB" "$FLANK_MB")
  cat > "${SCRIPTS}/C_analyze.sbatch" <<EOF
#!/bin/bash
#SBATCH --job-name=bp_analyze
#SBATCH --account=<account>
#SBATCH --partition=<partition>
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=${LOGS}/C_analyze_%j.out
${ENVLINES}
mamba run -p ${PY_ENV} python3 ${SCRIPTS}/analyze.py '${CFG}'
EOF
}

case "${1:-}" in
check)
  write_helpers
  printf "%-28s %-5s %-5s %-5s %-7s %s\n" SPECIES bed faa fai genome softmasked
  for sp in "${SPECIES[@]}"; do
    b=--; f=--; i=--; g=--; sm=--
    [[ -s ${MS_OUT}/${sp}/${sp}.longest.bed ]] && b=ok
    [[ -s ${MS_OUT}/${sp}/${sp}.longest.faa ]] && f=ok
    [[ -s ${MS_OUT}/${sp}/${sp}.fa.fai ]] && i=ok
    if [[ -s ${GENOME_DIR}/${sp}/${sp}.masked.fa ]]; then
      g=ok
      lc=$(head -c 3000000 ${GENOME_DIR}/${sp}/${sp}.masked.fa | grep -v '>' | tr -cd 'acgt' | wc -c)
      [[ $lc -gt 1000 ]] && sm=yes || sm="NO (hard-masked? repeat feature unusable)"
    fi
    printf "%-28s %-5s %-5s %-5s %-7s %s\n" "$sp" "$b" "$f" "$i" "$g" "$sm"
  done
  for x in z_chromosome.txt w_chromosome.txt Species_colours.txt species_tree.nwk; do
    [[ -s ${INPUT_DIR}/$x ]] && echo "input $x: ok" || echo "input $x: MISSING (place it in ${INPUT_DIR})"
  done
  [[ -x ${JCVI_ENV}/bin/diamond ]] && echo "JCVI env: ok" || echo "JCVI env: diamond not found"
  [[ -x ${PY_ENV}/bin/python3 ]] && echo "python env: ok" || echo "python env: MISSING" ;;
submit)
  write_helpers; write_sbatch
  ja=$(sbatch --parsable "${SCRIPTS}/A_paint.sbatch") || exit 1
  jb=$(sbatch --parsable "${SCRIPTS}/B_features.sbatch") || exit 1
  jc=$(sbatch --parsable --dependency=afterok:${ja},afterany:${jb} "${SCRIPTS}/C_analyze.sbatch") || exit 1
  log "paint=${ja}  features=${jb}  analyze=${jc}   results -> ${OUT_DIR}/results/" ;;
analyze)
  write_helpers; write_sbatch; sbatch "${SCRIPTS}/C_analyze.sbatch" ;;
status)
  echo "anchor files : $(ls ${OUT_DIR}/pairs/*.anchors 2>/dev/null | wc -l) / $(( ${#SPECIES[@]} * ${#REFS[@]} - ${#REFS[@]} ))"
  echo "feature files: $(ls ${OUT_DIR}/features/*.windows.tsv 2>/dev/null | wc -l) / ${#SPECIES[@]}"
  echo "results      : $(ls ${OUT_DIR}/results 2>/dev/null | wc -l) files"
  squeue -u "$USER" -o "%.12i %.14j %.8T %.10M" 2>/dev/null | grep -E "bp_|JOBID" ;;
*) sed -n '2,25p' "$0" ;;
esac
