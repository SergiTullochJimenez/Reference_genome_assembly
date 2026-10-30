"""00.1 Step 0 of the macro/micro comparison: per-chromosome means (length-weighted over the 1-Mb windows) of GC content,
gene density and repeat content, from the Circos tracks of the 8 new genomes (analyses/circos/01.Circos_tracks.sh).
Writes <CIRCOS_BASE>/Chrom_comparisson/per_chromosome_summary.csv; the class column is overwritten by 00.2.
Run: CIRCOS_BASE=/path/to/Circos_plot python3 00.1.Per_chromosome_means.py"""
import os
import numpy as np
import pandas as pd
from scipy import stats

BASE = os.path.expanduser(os.environ.get("CIRCOS_BASE"))   # folder with one sub-folder of 1-Mb tracks per species
OUTDIR = os.path.join(BASE, "Chrom_comparisson")
os.makedirs(OUTDIR, exist_ok=True)

SPECIES = ["Coronella_austriaca","Coronella_girondica","Hierophis_viridiflavus",
           "Macroprotodon_brevis","Macroprotodon_mauritanicus","Natrix_astreptophora",
           "Zamenis_longissimus","Zamenis_scalaris"]

MACRO_THRESHOLD = 50_000_000

# H. viridiflavus has 16 chromosomes; scaffold_17 / scaffold_18 are unplaced scaffolds (its scaffolds are not named
# SUPER_SUPER_*, so the name filter in load_chrom_sizes does not remove them).
EXCLUDE = {
    "Hierophis_viridiflavus": {"scaffold_17", "scaffold_18"},
}

TRACKS = {
    "gc": "gc_1mb.bed",
    "gene_density": "gene_density_1mb.bed",
    "repeat": "repeat_1mb.bed",
    "hic_cov": "hic_cov_1mb.bed",
    "ont_cov": "ont_cov_1mb.bed",
}
metric_labels = {
    "gc": "GC content",
    "gene_density": "Gene density (genes / Mb window)",
    "repeat": "Repeat content (fraction)",
    "hic_cov": "Hi-C coverage",
    "ont_cov": "ONT coverage",
}

# only GC, gene density and repeat content are written to the table
PLOT_METRICS = ["gc", "gene_density", "repeat"]
STATS_METRIC_LABELS = {k: metric_labels[k] for k in PLOT_METRICS}

def load_chrom_sizes(path):
    sizes = {}
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                sizes[parts[0]] = int(parts[1])
    ss = {k: v for k, v in sizes.items() if k.startswith("SUPER_SUPER")}
    return ss if ss else sizes

def classify(name, length):
    if name.endswith("_Z"):
        return "Z"
    if name.endswith("_W"):
        return "W"
    return "macro" if length >= MACRO_THRESHOLD else "micro"

def load_track_means(path, chrom_names):
    sums = {c: 0.0 for c in chrom_names}
    weights = {c: 0.0 for c in chrom_names}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom, start, end, val = parts[0], parts[1], parts[2], parts[3]
            if chrom not in sums:
                continue
            try:
                start = int(start); end = int(end); val = float(val)
            except ValueError:
                continue
            w = max(end - start, 1)
            sums[chrom] += val * w
            weights[chrom] += w
    return {c: (sums[c] / weights[c] if weights[c] > 0 else np.nan) for c in chrom_names}

rows = []
for sp in SPECIES:
    sp_dir = os.path.join(BASE, sp)
    sizes_path = os.path.join(sp_dir, "chrom_sizes.txt")
    if not os.path.exists(sizes_path):
        print(f"WARNING: missing chrom_sizes.txt for {sp}")
        continue
    sizes = load_chrom_sizes(sizes_path)
    excl = EXCLUDE.get(sp, set())
    sizes = {k: v for k, v in sizes.items() if k not in excl}
    chrom_names = list(sizes.keys())

    track_means = {}
    missing_tracks = []
    for key, fname in TRACKS.items():
        fpath = os.path.join(sp_dir, fname)
        if not os.path.exists(fpath):
            missing_tracks.append(fname)
            track_means[key] = {c: np.nan for c in chrom_names}
            continue
        track_means[key] = load_track_means(fpath, chrom_names)
    if missing_tracks:
        print(f"WARNING: {sp} missing tracks: {missing_tracks}")

    for c in chrom_names:
        length = sizes[c]
        cls = classify(c, length)
        row = {"species": sp, "chromosome": c, "length": length, "class": cls}
        for key in PLOT_METRICS:
            row[key] = track_means[key][c]
        rows.append(row)

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUTDIR, "per_chromosome_summary.csv"), index=False)
print(f"Wrote per_chromosome_summary.csv ({len(df)} chromosomes across {df['species'].nunique()} species)")
