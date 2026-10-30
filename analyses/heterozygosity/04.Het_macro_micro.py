"""Heterozygosity of macro- vs micro-chromosomes, from the 100-kb window tables of 03.Het_matched.sh (windows/*.txt) and the
reference .fai files (fai/<Species>.fa.fai). Macro = autosomes >= 70 Mb (same threshold as the 8 new genomes), micro = autosomes
5-70 Mb; windows with >= 10 kb callable; sex chromosomes are not in the autosome BED.
UNIT of inference = chromosome (windows within a chromosome are not independent); species = 14 (primary) / 15 (sensitivity), sign/Wilcoxon across species.
Vipera_aspis (assembly with 10 macro-sized scaffolds) is excluded from the primary test and reported as sensitivity.
The paper reports the all-15-species result (printed as 'ALL' and 'SENSITIVITY (all 15, incl. V. aspis)').
Run: python3 04.Het_macro_micro.py  (pandas, numpy, scipy)"""
import glob, os, numpy as np, pandas as pd
from scipy import stats
MINCAL = 10000; MICRO_MIN = 5e6; MACRO_MIN = 70e6; rng = np.random.default_rng(3)
rows_c, rows_s, rows_w = [], [], []
for f in sorted(glob.glob("windows/*.txt")):
    b = os.path.basename(f).split("_heterozygosity")[0]; sp = b[:len(b)//2]
    fai = pd.read_csv(f"fai/{sp}.fa.fai", sep="\t", header=None, usecols=[0, 1], names=["c", "L"]).set_index("c").L
    w = pd.read_csv(f, sep="\t"); w = w[w.callable_sites >= MINCAL].copy()
    cs = w.groupby("chrom").agg(het=("heterozygous_sites", "sum"), cal=("callable_sites", "sum"), n=("chrom", "size"))
    cs["L"] = fai[cs.index].values; cs = cs[cs.L >= MICRO_MIN].sort_values("L", ascending=False)
    cs["cls"] = np.where(cs.L >= MACRO_MIN, "macro", "micro")
    cs["het_rate"] = cs.het / cs.cal; cs["callable_frac"] = cs.cal / (cs.n * 1e5); cs["species"] = sp
    cs = cs.reset_index(); rows_c.append(cs)
    w = w[w.chrom.isin(cs.chrom)].merge(cs[["chrom", "cls"]], on="chrom"); w["species"] = sp; rows_w.append(w)
    M, m = cs[cs.cls == "macro"], cs[cs.cls == "micro"]
    pm, pm2 = M.het.sum() / M.cal.sum(), m.het.sum() / m.cal.sum()
    wm = w[w.cls == "macro"].heterozygosity; wi = w[w.cls == "micro"].heterozygosity
    u_c = stats.mannwhitneyu(m.het_rate, M.het_rate, alternative="two-sided").pvalue
    u_w = stats.mannwhitneyu(wi, wm, alternative="two-sided").pvalue
    rho, prho = stats.spearmanr(np.log10(cs.L), cs.het_rate)
    rows_s.append(dict(species=sp, n_macro=len(M), n_micro=len(m), het_macro=pm, het_micro=pm2, ratio=pm2 / pm,
                       p_chrom_MWU=u_c, p_window_MWU=u_w, callable_macro=M.callable_frac.mean(), callable_micro=m.callable_frac.mean(),
                       rho_size=rho, p_rho=prho))
C = pd.concat(rows_c); S = pd.DataFrame(rows_s); W = pd.concat(rows_w)
C.to_csv("het_by_chromosome_70Mb.tsv", sep="\t", index=False); S.to_csv("het_macro_vs_micro_per_species_70Mb.tsv", sep="\t", index=False)
pd.set_option("display.width", 220); print(S.round(5).to_string(index=False))
P = S[S.species != "Vipera_aspis"]
for lab, d in (("PRIMARY (no V. aspis)", P), ("ALL", S)):
    r = np.log(d.ratio); n_up = (d.ratio > 1).sum()
    print(f"\n{lab}: micro>macro in {n_up}/{len(d)}; sign test p={stats.binomtest(n_up, len(d), .5).pvalue:.4g}; "
          f"Wilcoxon on log ratio p={stats.wilcoxon(r).pvalue:.4g}; median ratio {np.exp(np.median(r)):.3f}; mean ratio {np.exp(r.mean()):.3f}")
# does it survive accounting for callable fraction (repeat-masking proxy)?  chromosome-level, species fixed effects
def model(Cc, lab):
    Cc = Cc.copy(); Cc["ln_rel"] = np.log(Cc.het_rate / Cc.groupby("species").het_rate.transform("mean"))
    Cc["cf"] = Cc.callable_frac - Cc.groupby("species").callable_frac.transform("mean"); Cc["m"] = (Cc.cls == "micro").astype(float)
    X = np.column_stack([np.ones(len(Cc)), Cc.m, Cc.cf]); y = Cc.ln_rel.values
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    spp = Cc.species.values; us = np.unique(spp); bs = []; r = np.random.default_rng(3)
    for _ in range(2000):  # cluster bootstrap by species
        pick = r.choice(us, len(us)); idx = np.concatenate([np.where(spp == s)[0] for s in pick])
        bs.append(np.linalg.lstsq(X[idx], y[idx], rcond=None)[0])
    bs = np.array(bs)
    print(f"\n{lab}: chromosome-level model ln(het/species mean) ~ micro + callable_fraction (species-demeaned), {len(us)} spp, {len(Cc)} chromosomes:")
    print(f"  micro effect = {beta[1]:+.3f} (x{np.exp(beta[1]):.2f}); 95% species-bootstrap CI {np.percentile(bs[:,1],2.5):+.3f}..{np.percentile(bs[:,1],97.5):+.3f} (x{np.exp(np.percentile(bs[:,1],2.5)):.2f}-{np.exp(np.percentile(bs[:,1],97.5)):.2f})")
    print(f"  callable_fraction slope = {beta[2]:+.3f}")
    return Cc
Cc = model(C[C.species != "Vipera_aspis"], "PRIMARY (no V. aspis)")
model(C, "SENSITIVITY (all 15, incl. V. aspis)")
Cm = C[C.species != "Vipera_aspis"]
print("\ncallable fraction macro vs micro (mean):", Cm.groupby("cls").callable_frac.mean().round(3).to_dict())
# within micro and within macro: size trend
for c in ("macro", "micro"):
    d = Cc[Cc.cls == c]; r, p = stats.spearmanr(np.log10(d.L), d.ln_rel); print(f"within {c}: Spearman het~size rho={r:+.3f} p={p:.3g} (n={len(d)})")
