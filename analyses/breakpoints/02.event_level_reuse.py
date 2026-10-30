"""Event-level breakpoint-reuse test.

UNIT = independent evolutionary event, not a species row. A break seen in the same place in the three
Natrix is ONE event (placed on the Natrix stem by shared_breakpoint_regions.tsv: 'Natrix (3 spp)'),
so it is counted once. Events are built from the rows of breakpoints.tsv (internal junction projections)
and splits.tsv on the primary reference (Zamenis_longissimus), macro autosomes 1-7 only.

STAT   = number of PAIRS of events from different genera (independent lineages) whose reference
         intervals come within TOL of each other (default 0.5 Mb, as before).
NULL   = every event re-placed uniformly at random on the reference macro autosomes (chromosome chosen
         proportionally to its length), same width and same number of events per genus; 5,000 draws.
Run: python3 02.event_level_reuse.py   (needs pandas, numpy; reads ./results/*.tsv from 01.breakpoints_pipeline.sh)"""
import sys, re, itertools, numpy as np, pandas as pd
REF = "Zamenis_longissimus"; R = "results"; TOLS = [0.5e6, 1e6, 3e6]; NPERM = 5000
rng = np.random.default_rng(7)

bp = pd.read_csv(f"{R}/breakpoints.tsv", sep="\t"); sp = pd.read_csv(f"{R}/splits.tsv", sep="\t")
rg = pd.read_csv(f"{R}/shared_breakpoint_regions.tsv", sep="\t"); rg = rg[rg.reference == REF]
bp = bp[(bp.reference == REF) & (bp.small_scaffold_involved == 0)]; sp = sp[(sp.reference == REF) & (sp.small_scaffold_involved == 0)]

# projections of breaks onto the reference (internal only), same definition as the pipeline
rows = []
for _, b in bp.iterrows():
    for side in ("left", "right"):
        if b[f"{side}_type"] == "internal":
            p = b[f"{side}_r_pos"]; rows.append((b.species, b[f"{side}_r_chrom"], b[f"{side}_r_label"], p, p))
for _, x in sp.iterrows():
    rows.append((x.species, x.r_chrom, x.r_label, x.r_left_end, x.r_right_start))
P = pd.DataFrame(rows, columns=["species", "chrom", "label", "lo", "hi"]).drop_duplicates()
MACRO = [str(i) for i in range(1, 8)]
P = P[P.label.astype(str).isin(MACRO)]
genus = lambda s: s.split("_")[0]

# assign each projection to an event = (region, clade) using the pipeline's own placement
def species_of_entry(entry, region_species):
    m = re.match(r"^(.*) \((\d+) spp\)$", entry)
    if not m: return [entry]
    gens = m.group(1).split("+"); return [s for s in region_species if genus(s) in gens]
events = []
for _, r in rg.iterrows():
    if str(r.r_label) not in MACRO: continue
    rs = r.species.split(","); claimed = set()
    for entry in [e.strip() for e in r.placement.split("|")]:
        members = [s for s in species_of_entry(entry, rs) if s not in claimed]; claimed |= set(members)
        sub = P[(P.chrom == r.r_chrom) & (P.species.isin(members)) & (P.hi >= r.r_start) & (P.lo <= r.r_end)]
        if len(sub):
            events.append(dict(region=r.region_id, chrom=r.r_chrom, label=str(r.r_label), clade=entry, genus=genus(members[0]),
                               lo=sub.lo.min(), hi=sub.hi.max(), n_species=len(members)))
E = pd.DataFrame(events)
E["width"] = E.hi - E.lo
L = (pd.read_csv("../Circos_plot/Zamenis_longissimus/chrom_sizes.txt", sep="\t", header=None, names=["c", "len"]).set_index("c").len)
lens = {c: int(L[c]) for c in E.chrom.unique() if c in L.index}
E = E[E.chrom.isin(lens)].reset_index(drop=True)
E.to_csv(f"{R}/reuse_events_macro.tsv", sep="\t", index=False)
print(f"{len(E)} independent events on {len(lens)} reference macro autosomes "
      f"(from {len(P)} species-level projections); by genus: {E.genus.value_counts().to_dict()}")

def pairs(chrom, lo, hi, gen, tol):
    n = 0
    for i, j in itertools.combinations(range(len(chrom)), 2):
        if gen[i] != gen[j] and chrom[i] == chrom[j] and lo[i] - tol <= hi[j] and lo[j] - tol <= hi[i]: n += 1
    return n
def clusters(chrom, lo, hi, gen, tol, mid=False):
    """number of hotspots = connected groups (>=2 events from different genera) -- a hotspot hit by 3
    lineages counts once, not as 3 pairs. mid=True compares event CENTRES (robust to wide merged intervals)."""
    n = len(chrom); par = list(range(n))
    def f(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    c = (np.asarray(lo) + np.asarray(hi)) / 2
    for i, j in itertools.combinations(range(n), 2):
        if gen[i] == gen[j] or chrom[i] != chrom[j]: continue
        ok = abs(c[i] - c[j]) <= tol if mid else (lo[i] - tol <= hi[j] and lo[j] - tol <= hi[i])
        if ok: par[f(i)] = f(j)
    grp = {}
    for i in range(n): grp.setdefault(f(i), []).append(i)
    return sum(1 for g in grp.values() if len(g) >= 2 and len({gen[k] for k in g}) >= 2)
def stat(E, tol): return pairs(E.chrom.values, E.lo.values, E.hi.values, E.genus.values, tol)

cl = np.array(list(lens)); w = np.array([lens[c] for c in cl], float); w /= w.sum()
out = []
for tol in TOLS:
    obs = stat(E, tol); null = []
    for _ in range(NPERM):
        c = rng.choice(cl, len(E), p=w)
        lo = np.array([rng.uniform(0, max(lens[x] - wd, 1)) for x, wd in zip(c, E.width)])
        null.append(pairs(c, lo, lo + E.width.values, E.genus.values, tol))
    null = np.array(null)
    out.append(dict(stat="pairs, interval overlap", tol_Mb=tol / 1e6, observed=obs, null_mean=null.mean(), null_95=np.percentile(null, 95),
                    p=(np.sum(null >= obs) + 1) / (NPERM + 1)))
    for mid in (False, True):
        o = clusters(E.chrom.values, E.lo.values, E.hi.values, E.genus.values, tol, mid); nl = []
        for _ in range(NPERM):
            c = rng.choice(cl, len(E), p=w)
            lo = np.array([rng.uniform(0, max(lens[x] - wd, 1)) for x, wd in zip(c, E.width)])
            nl.append(clusters(c, lo, lo + E.width.values, E.genus.values, tol, mid))
        nl = np.array(nl)
        out.append(dict(stat="hotspots, " + ("centre distance" if mid else "interval overlap"), tol_Mb=tol / 1e6, observed=o,
                        null_mean=nl.mean(), null_95=np.percentile(nl, 95), p=(np.sum(nl >= o) + 1) / (NPERM + 1)))
print(pd.DataFrame(out).round(4).to_string(index=False))
pd.DataFrame(out).to_csv(f"{R}/reuse_test_event_level.tsv", sep="\t", index=False)

# which pairs (for the paper's examples)
print("\noverlapping event pairs at 0.5 Mb:")
for i, j in itertools.combinations(range(len(E)), 2):
    a, b = E.loc[i], E.loc[j]
    if a.genus != b.genus and a.chrom == b.chrom and a.lo - 0.5e6 <= b.hi and b.lo - 0.5e6 <= a.hi:
        print(f"  Zl{a.label}: {a.clade} [{a.lo/1e6:.2f}-{a.hi/1e6:.2f} Mb]  x  {b.clade} [{b.lo/1e6:.2f}-{b.hi/1e6:.2f} Mb]")
