"""Event-level breakpoint-environment test, Z. longissimus reference only.
Inputs (from 01.breakpoints_pipeline.sh): results/breakpoints.tsv, results/shared_breakpoint_regions.tsv and
results/terminal_joins.tsv, keeping only the rows with reference == Zamenis_longissimus, e.g.
   awk -F"\t" 'NR==1 || $1=="Zamenis_longissimus"' <pipeline_results>/breakpoints.tsv > results/breakpoints.tsv
(same for the other two files), and features/<species>.windows.tsv (stage B of the pipeline).
Only junction breakpoints (breakpoints.tsv) are tested; fissions (splits.tsv) are not.
Breakpoint-environment test with the EVENT as the unit (fixes the pseudoreplication of the pipeline's test,
which counted a break shared by three Natrix as three independent breakpoints).

Per species-level breakpoint (junction inside a chromosome, as breakpoints.tsv; the same junction seen from
both references counted once):
   value   = feature (GC, repeat fraction, genes, telomeric repeats, gap fraction) in the 100 kb window at the breakpoint
   control = 50 windows from the SAME species, SAME chromosome class (micro <70 Mb / macro) and the SAME distance to the
             chromosome end (same construction as breakpoints_pipeline.sh)
   z       = value - mean(control)            pct = fraction of its own controls that are <= value
Breakpoints are then grouped into independent EVENTS (same event = same reference region AND same clade placement in
shared_breakpoint_regions.tsv, e.g. 'Natrix (3 spp)'; for end-to-end fusions: same join and same genus) and the
z / pct of the members are AVERAGED, so every event counts once.
TEST   : mean(pct) against 0.5 (uniform under the null) by simulation; mean z with a sign-permutation test.
Run: python3 03.event_level_environment.py"""
import re, glob, os, itertools, collections, numpy as np, pandas as pd
R = "results"; F = "features"; REFS = ["Zamenis_longissimus", "NONE"]; PRIMARY = REFS[0]
W = 100000; FLANK = 2e6; NCTRL = 50; NSIM = 20000; MICRO_MB = 70
FEATS = ["gc", "repeat", "genes", "telomeric_per_kb", "gap"]
rng = np.random.default_rng(11)
genus = lambda s: s.split("_")[0]

bp = pd.read_csv(f"{R}/breakpoints.tsv", sep="\t"); bp = bp[bp.small_scaffold_involved == 0].copy()
rg = pd.read_csv(f"{R}/shared_breakpoint_regions.tsv", sep="\t")
tj = pd.read_csv(f"{R}/terminal_joins.tsv", sep="\t")
FEAT = {os.path.basename(f).split(".")[0]: pd.read_csv(f, sep="\t") for f in glob.glob(f"{F}/*.windows.tsv")}
CL = {sp: f.groupby("chrom").end.max() for sp, f in FEAT.items()}

# ---------------------------------------------------------------- event keys
def species_of_entry(entry, members):
    m = re.match(r"^(.*) \((\d+) spp\)$", entry)
    return [entry] if not m else [s for s in members if genus(s) in m.group(1).split("+")]
REG = []  # (ref, chrom, lo, hi, {species: clade})
for _, r in rg.iterrows():
    members = r.species.split(","); claimed = set(); cmap = {}
    for entry in [e.strip() for e in r.placement.split("|")]:
        sel = [s for s in species_of_entry(entry, members) if s not in claimed]; claimed |= set(sel)
        for s in sel: cmap[s] = entry
    REG.append((r.reference, r.r_chrom, r.r_start, r.r_end, cmap, r.region_id))
def event_key(b):
    """key from the reference projection of this breakpoint (prefer the primary reference)."""
    for side in ("left", "right"):
        if b[f"{side}_type"] == "internal":
            for ref, ch, lo, hi, cmap, rid in REG:
                if ref == b.reference and ch == b[f"{side}_r_chrom"] and lo - 1e6 <= b[f"{side}_r_pos"] <= hi + 1e6 and b.species in cmap:
                    return ("region", rid, cmap[b.species])
    ends = sorted([f"{b.left_r_label}:{b.left_type}", f"{b.right_r_label}:{b.right_type}"])
    return ("join", b.reference, " + ".join(ends), genus(b.species))
bp["centre"] = ((bp.q_start + bp.q_end) // 2).astype(int)
bp["junc"] = bp.species + "|" + bp.q_chrom + "|" + (bp.centre // 1_000_000).astype(str)
bp["ref_rank"] = bp.reference.map({PRIMARY: 0, REFS[1]: 1})
bp = bp.sort_values("ref_rank")
keys = {}
for _, b in bp.iterrows():
    keys.setdefault(b.junc, event_key(b))           # primary-reference key wins
bp = bp.drop_duplicates("junc").copy(); bp["event"] = bp.junc.map(keys)
bp["event_id"] = bp.event.astype(str)

# ------------------------------------------------------------ window features
def win(fc, pos):
    st = int(pos // W * W); return fc.loc[st, FEATS].values.astype(float) if st in fc.index else np.full(len(FEATS), np.nan)
rows = []
for _, b in bp.iterrows():
    fw = FEAT[b.species]; fc = fw[fw.chrom == b.q_chrom].set_index("start")
    if fc.empty: continue
    Lq = float(CL[b.species][b.q_chrom]); cls = "micro" if Lq < MICRO_MB * 1e6 else "macro"
    v = win(fc, b.centre); d_end = min(b.centre, Lq - b.centre)
    chroms = [c for c in CL[b.species].index if (CL[b.species][c] < MICRO_MB * 1e6) == (cls == "micro")]
    same = fw[fw.chrom.isin(chroms)]; ctrl = []
    for _ in range(NCTRL * 4):
        if len(ctrl) >= NCTRL: break
        row = same.sample(1, random_state=int(rng.integers(1e9))).iloc[0]; Lc = float(CL[b.species][row.chrom])
        if Lc < 2 * d_end + 2 * FLANK: continue
        pos = d_end if rng.random() < .5 else Lc - d_end
        ctrl.append(win(fw[fw.chrom == row.chrom].set_index("start"), pos))
    ctrl = np.array(ctrl, float)
    if len(ctrl) < 10: continue
    rec = dict(species=b.species, q_chrom=b.q_chrom, centre=b.centre, cls=cls, event=b.event_id)
    for k, f in enumerate(FEATS):
        c = ctrl[:, k]; c = c[~np.isnan(c)]
        rec[f"z_{f}"] = v[k] - c.mean() if len(c) and not np.isnan(v[k]) else np.nan
        rec[f"pct_{f}"] = ((c < v[k]).sum() + 0.5 * (c == v[k]).sum()) / len(c) if len(c) and not np.isnan(v[k]) else np.nan
    rows.append(rec)
B = pd.DataFrame(rows); B.to_csv(f"{R}/breakpoint_level_z.tsv", sep="\t", index=False)
agg = {c: "mean" for c in B.columns if c.startswith(("z_", "pct_"))}
E = B.groupby(["cls", "event"]).agg(**{k: (k, "mean") for k in agg}, n_breakpoints=("species", "size"),
                                     species=("species", lambda x: ",".join(sorted(set(x))))).reset_index()
E.to_csv(f"{R}/event_level_environment_events.tsv", sep="\t", index=False)

def test(df, f):
    pct = df[f"pct_{f}"].dropna().values; z = df[f"z_{f}"].dropna().values
    if len(pct) < 3: return None
    sim = rng.random((NSIM, len(pct))).mean(1)
    p_pct = (np.sum(np.abs(sim - .5) >= abs(pct.mean() - .5)) + 1) / (NSIM + 1)
    sg = rng.choice([-1, 1], (NSIM, len(z))); p_z = (np.sum(np.abs((sg * z).mean(1)) >= abs(z.mean())) + 1) / (NSIM + 1)
    return dict(n=len(pct), mean_pct=pct.mean(), mean_diff=z.mean(), p_percentile=p_pct, p_signflip=p_z)
out = []
for unit, df in (("breakpoint (species rows)", B), ("EVENT", E)):
    for cls in ("macro", "micro"):
        for f in FEATS:
            t = test(df[df.cls == cls], f)
            if t: out.append(dict(unit=unit, cls=cls, feature=f, **t))
T = pd.DataFrame(out); T.to_csv(f"{R}/event_level_environment_tests.tsv", sep="\t", index=False)
pd.set_option("display.width", 200)
print(T.round(4).to_string(index=False))
print("\nevents:", E.groupby("cls").size().to_dict(), " breakpoints:", B.groupby("cls").size().to_dict())
