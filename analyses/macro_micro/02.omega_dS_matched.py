"""Is the lower omega of micro-chromosome genes an artefact of their ~2x higher dS?
Per-locus omega/dN/dS (tree-wide, codeml M0) from results/step3_dnds_per_gene_filtered.tsv (macro / micro / Z autosome classes).
Tests: (1) dS-stratified comparison; (2) 1:1 nearest-dS matching (caliper 5%) macro vs micro, matched pairs Wilcoxon;
(3) restriction to low dS (<0.5, <0.3) ; (4) regression log(omega) ~ micro + log(dS) + log(codons); (5) dN comparison at matched dS.
Run from the pipeline output folder: python3 02.omega_dS_matched.py  (pandas numpy scipy)"""
import numpy as np, pandas as pd
from scipy import stats
d = pd.read_csv("results/step3_dnds_per_gene_filtered.tsv", sep="\t"); d = d[d.status == "ok"].dropna(subset=["omega", "dS_tree", "dN_tree"])
d = d[d["class"].isin(["macro", "micro"])].copy(); d["m"] = (d["class"] == "micro").astype(int); d["lds"] = np.log10(d.dS_tree); d["lw"] = np.log10(d.omega.clip(lower=1e-4))
print("genes:", d.groupby("class").size().to_dict()); print("median dS", d.groupby("class").dS_tree.median().round(3).to_dict(), " median omega", d.groupby("class").omega.median().round(3).to_dict())
mw = lambda a, b: stats.mannwhitneyu(a, b).pvalue
print("\n(0) raw: omega micro vs macro p=%.2g" % mw(d[d.m == 1].omega, d[d.m == 0].omega))
# (1) strata of dS
print("\n(1) dS strata (quantile bins of dS over all genes)")
d["bin"] = pd.qcut(d.dS_tree, 6)
rows = []
for b, g in d.groupby("bin", observed=True):
    a, c = g[g.m == 1], g[g.m == 0]
    rows.append((str(b), len(c), len(a), c.omega.median(), a.omega.median(), mw(a.omega, c.omega)))
S = pd.DataFrame(rows, columns=["dS_bin", "n_macro", "n_micro", "omega_macro", "omega_micro", "p"]); print(S.round(4).to_string(index=False))
# Cochran-Mantel-Haenszel-like: stratified permutation of class labels within bins
rng = np.random.default_rng(1); obs = d[d.m == 1].omega.median() - d[d.m == 0].omega.median()
def stat(lbl): return np.mean([ (g.omega[lbl[g.index] == 1].median() - g.omega[lbl[g.index] == 0].median()) for _, g in d.groupby("bin", observed=True)])
lab = d.m.copy(); o = stat(lab); null = []
for _ in range(2000):
    l2 = lab.copy()
    for _, g in d.groupby("bin", observed=True): l2.loc[g.index] = rng.permutation(g.m.values)
    null.append(stat(l2))
null = np.array(null); print("stratified: mean within-bin median difference (micro-macro) = %.4f ; permutation p = %.4f" % (o, (np.sum(np.abs(null) >= abs(o)) + 1) / 2001))
# (2) nearest-dS matching
mi = d[d.m == 1].sort_values("dS_tree"); ma = d[d.m == 0].sort_values("dS_tree").copy(); pairs = []; used = set(); arr = ma.dS_tree.values; idx = ma.index.values
for i, r in mi.sample(frac=1, random_state=3).iterrows():
    j = np.searchsorted(arr, r.dS_tree); best = None
    for k in (j - 3, j - 2, j - 1, j, j + 1, j + 2):
        if 0 <= k < len(arr) and idx[k] not in used and abs(arr[k] - r.dS_tree) <= 0.05 * r.dS_tree:
            if best is None or abs(arr[k] - r.dS_tree) < abs(arr[best] - r.dS_tree): best = k
    if best is not None: used.add(idx[best]); pairs.append((i, idx[best]))
P = pd.DataFrame(pairs, columns=["mi", "ma"]); a = d.loc[P.mi]; c = d.loc[P.ma]
print("\n(2) matched on dS (+-5%%): %d pairs of %d micro genes; dS medians %.3f vs %.3f" % (len(P), len(mi), a.dS_tree.median(), c.dS_tree.median()))
for v in ("omega", "dN_tree"):
    print("   %s median micro %.4f vs macro %.4f ; Wilcoxon paired p=%.2g ; micro<macro in %.1f%% of pairs" % (v, a[v].median(), c[v].median(), stats.wilcoxon(a[v].values, c[v].values).pvalue, 100 * (a[v].values < c[v].values).mean()))
rs = (a.omega.values.sum()) / 1; 
print("   ratio of sums omega (dN/dS): micro %.4f macro %.4f" % (a.dN_tree.sum() / a.dS_tree.sum(), c.dN_tree.sum() / c.dS_tree.sum()))
# (3) low dS only
print("\n(3) restricting dS")
for cut in (1.0, 0.5, 0.3):
    g = d[d.dS_tree < cut]; x, y = g[g.m == 1], g[g.m == 0]
    print("   dS<%.1f: n %d/%d  omega micro %.4f vs macro %.4f  p=%.2g ; median dS %.3f vs %.3f" % (cut, len(y), len(x), x.omega.median(), y.omega.median(), mw(x.omega, y.omega), x.dS_tree.median(), y.dS_tree.median()))
# (4) regression
X = np.column_stack([np.ones(len(d)), d.m, d.lds, np.log10(d.ls)]); b, *_ = np.linalg.lstsq(X, d.lw.values, rcond=None)
res = d.lw.values - X @ b; sig2 = res @ res / (len(d) - 4); cov = sig2 * np.linalg.inv(X.T @ X); se = np.sqrt(np.diag(cov))
print("\n(4) OLS log10(omega) ~ micro + log10(dS) + log10(codons):")
for n_, bb, ss in zip(["intercept", "micro", "log10 dS", "log10 codons"], b, se): print("   %-13s %+.4f (se %.4f)  t=%.1f" % (n_, bb, ss, bb / ss))
print("   -> micro effect at equal dS and length: omega x%.3f" % 10 ** b[1])
print("   (note: omega and dS share dS in the denominator, so omega~dS has a negative bias; the matched design (2) is the cleaner test)")
