"""00.2 Step 0 of the macro/micro comparison: chromosome classes. Autosomes >= 70 Mb = macro, < 70 Mb = micro (70 Mb lies in
the size gap shared by the 8 genomes); Z and W are classed by name, and the Z of the three male genomes (named as ordinary
scaffolds) is listed in MALE_Z. Also writes the per-species counts and pooled / per-species Mann-Whitney tests.
Run: python3 00.2.Macro_micro_classes.py <per_chromosome_summary.csv from 00.1> <output_dir>
The output per_chromosome_summary_v2.csv is the CHROM_TABLE of 01.macro_micro_pipeline.sh."""
import sys, numpy as np, pandas as pd
from scipy import stats
IN, OUT = sys.argv[1], sys.argv[2]
THR = 70e6
MALE_Z = {("Hierophis_viridiflavus","scaffold_4"),("Macroprotodon_brevis","SUPER_SUPER_4"),
          ("Zamenis_longissimus","SUPER_SUPER_4")}
d = pd.read_csv(IN)
def cls(r):
    if (r.species, r.chromosome) in MALE_Z or r.chromosome.endswith("_Z"): return "Z"
    if r.chromosome.endswith("_W"): return "W"
    return "macro" if r.length >= THR else "micro"
d["class"] = d.apply(cls, axis=1)
d["log10_Mb"] = np.log10(d.length/1e6)
d.to_csv(f"{OUT}/per_chromosome_summary_v2.csv", index=False)
a = d[d["class"].isin(["macro","micro"])]
M = {"gc":"GC content","gene_density":"Gene density (genes / Mb window)","repeat":"Repeat content (fraction)"}
cnt=[]
for sp,g in a.groupby("species"):
    ma=g[g["class"]=="macro"].length; mi=g[g["class"]=="micro"].length
    cnt.append(dict(species=sp,n_macro=len(ma),n_micro=len(mi),macro_min_Mb=ma.min()/1e6,micro_max_Mb=mi.max()/1e6,gap_Mb=(ma.min()-mi.max())/1e6))
cnt=pd.DataFrame(cnt); cnt.to_csv(f"{OUT}/macro_micro_counts_v2.csv",index=False); print(cnt.to_string())
pool=[]
for k,l in M.items():
    x=a[a["class"]=="macro"][k].dropna(); y=a[a["class"]=="micro"][k].dropna()
    u,p=stats.mannwhitneyu(x,y,alternative="two-sided")
    pool.append(dict(metric=l,macro_n=len(x),micro_n=len(y),macro_mean=x.mean(),macro_median=x.median(),micro_mean=y.mean(),micro_median=y.median(),mannwhitney_U=u,p_value=p))
pool=pd.DataFrame(pool); pool.to_csv(f"{OUT}/macro_vs_micro_stats_pooled_v2.csv",index=False); print(pool.to_string())
per=[]
for sp,g in a.groupby("species"):
    for k,l in M.items():
        x=g[g["class"]=="macro"][k].dropna(); y=g[g["class"]=="micro"][k].dropna()
        u,p=stats.mannwhitneyu(x,y,alternative="two-sided")
        per.append(dict(species=sp,metric=l,macro_mean=x.mean(),micro_mean=y.mean(),
            direction="micro > macro" if y.mean()>x.mean() else "macro > micro",mannwhitney_U=u,p_value=p))
per=pd.DataFrame(per); per.to_csv(f"{OUT}/macro_vs_micro_stats_per_species_v2.csv",index=False); print(per.to_string())
