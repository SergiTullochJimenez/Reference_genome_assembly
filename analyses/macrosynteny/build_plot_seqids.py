#!/usr/bin/env python3
"""
build_plot_seqids.py -- for one species with a known chromosome count, build
a PLOTTING-ONLY bed (gene rows unchanged, seqid column renamed) and print the
seqids.txt row to use for it.

Why this exists: the real .longest.bed uses raw ENA/GenBank accessions or
similar for the seqid column (e.g. "ENA_OZ203524_OZ203524.1"), and includes
every scaffold in the assembly, not just true chromosomes. For the final
karyotype plot we only want to draw the top-N largest scaffolds (the real
chromosomes), labeled cleanly as "1","2",...,"Z","W" instead of the raw
accession. This does NOT touch gene names/order/anything used by diamond or
jcvi's synteny scan -- it only renames the seqid column for the copy of the
bed that the FINAL PLOT's layout.txt points to. Genes on excluded (non-
chromosome) scaffolds are left in place with their original seqid so nothing
crashes; they simply never appear in the printed seqids.txt row, and JCVI's
own Track.get_coords() already returns [None, None] (silently skipped by
draw_blocks) for any gene whose seqid isn't in the plotted seqids list.

Numbering rule: take the N largest GENE-BEARING scaffolds (by length, from
the .fai, restricted to scaffolds with >=1 annotated gene in BED -- see
"Why gene-bearing" below) as "the chromosomes". If this species' Z scaffold
(from z_chromosome.txt) is among them, it's labeled "Z" instead of consuming
a number; if it's missing from the top-N (unusual -- e.g. a very small Z),
it's still force-included and labeled "Z" (a warning is printed), UNLESS it
has zero genes (see below). W is always force-included regardless of size
rank (real W chromosomes in these snakes are often much smaller than
autosomes) and labeled "W", same zero-gene caveat. The rest of the top-N,
sorted by descending length, are numbered "1".."N" (minus 1 if Z was among
them). Multiple comma-separated Z or W scaffolds (e.g. Za/Zb, or unlocalized
W contigs) are labeled "Z","Zb","Zc"... / "W","Wb","Wc"... in the order
given.

Why gene-bearing: 04.plot_macrosynteny.sh's karyotype call uses --basepair
for true physical-length chromosome bars. In that mode, jcvi.graphics.
karyotype derives each displayed chromosome's size from the MAX GENE END
COORDINATE on that seqid (jcvi/graphics/karyotype.py's Karyotype.__init__,
the `else` branch of `if generank: ... else: sz = dict((x, max(z.end for z
in bed.sub_bed(x))) for x in seqids)`) -- it does NOT read the .fai in this
mode. A "chromosome" with zero annotated genes crashes with "max() arg is
an empty sequence". A large scaffold with no genes at all is very likely
unannotated repeat/heterochromatic sequence rather than a true chromosome
anyway, so we skip it and backfill from the next-largest gene-bearing
scaffold instead of crashing -- you still get N genuine entries whenever the
assembly has that many gene-bearing scaffolds. Known sex scaffolds
(Z/W) that turn out to have zero genes are dropped with a loud warning
rather than crashing the whole plot -- if that happens for real, the fix is
either better annotation coverage there or removing that entry until then.

Usage:
  build_plot_seqids.py SPECIES BED FAI N Z_SCAFFOLDS W_SCAFFOLDS OUT_BED

  SPECIES       species name (for warning messages only)
  BED           input .longest.bed (or .longest.flipped.bed)
  FAI           genome .fa.fai (for scaffold lengths)
  N             chromosome count (autosomes + Z) for this species
  Z_SCAFFOLDS   comma-separated Z scaffold accession(s), or "" if unknown
  W_SCAFFOLDS   comma-separated W scaffold accession(s), or "" if unknown/male
  OUT_BED       path to write the renamed plotting bed to

Prints the seqids.txt row (comma-separated labels, chromosome order) to
stdout. Prints nothing else to stdout; all diagnostics go to stderr.

---------------------------------------------------------------------------
Z-COLUMN ALIGNMENT BOOKKEEPING (added 2026-08-20)
---------------------------------------------------------------------------
In --basepair mode, jcvi.graphics.karyotype.Track lays out each species'
bars left-to-right using ratio = span/total, where `total` is the SUM of
each plotted seqid's "size" -- and in --basepair mode that size is NOT the
.fai scaffold length, it's max(gene END coordinate) on that seqid (see
Karyotype.__init__'s `else` branch: `sz = dict((x, max(z.end for z in
bed.sub_bed(x))) for x in seqids)`). So the W-gap patch already applied to
karyotype.py trivially right-aligns every species' W column (W is always
the LAST seqid in every row), but Z is NOT last (W follows it), so its
natural x-position varies per species with autosome count/total bp, and
needs a per-species DYNAMIC extra gap to line up -- unlike W's fixed 0.025.

This script now replicates JCVI's exact Track x-position formula (xstart=
0.18, xend=0.80 -- BOTH CONFIRMED CONSTANT across every species' row in
this project's own 03.run_synteny.sh, which is what makes a single shared
target alignment position meaningful at all) to compute, per species, the
natural (un-padded) x position where its first "Z"-labeled bar would start
-- using sizes computed from THIS SAME renamed OUT_BED via the identical
max-gene-end rule JCVI itself uses, not .fai lengths (those would give the
wrong ratio/total and silently produce a wrong alignment). The result is
appended (one line per species with a resolved Z) to a shared
z_alignment.tsv sidecar file living next to OUT_BED; 03.run_synteny.sh
aggregates that after all species have run into a per-species delta file
(z_gap_deltas.tsv, delta = max(z_start_x) - own z_start_x) that a
karyotype.py patch reads via the
KARYOTYPE_Z_GAP_FILE env var and adds as extra gap immediately before each
species' Z bar, mirroring the existing W-gap patch's shape.

GAP_DEFAULT below MUST match jcvi.graphics.karyotype's actual default
`gap` value (the CLI --gap default / Track gap parameter) for this to be
correct -- 0.01 is JCVI's default and is what is used here. Check your installed version with:
  grep -n "gap" "$(mamba run -p "$JCVI_ENV_PREFIX" python3 -c \
    'import jcvi.graphics.karyotype as k; print(k.__file__)')" | grep -i default
before trusting z_alignment.tsv's numbers for anything other than a sanity
check, and update GAP_DEFAULT here (and re-run) if it turns out to differ.
"""
import os
import sys

GAP_DEFAULT = 0.01     # JCVI karyotype default gap (see note above)
MAX_SEQIDS = 16
XSTART = 0.18           # confirmed constant -- 03.run_synteny.sh x_start
XEND = 0.80              # confirmed constant -- 03.run_synteny.sh x_end


def compute_sizes(out_bed_path, wanted_labels):
    """Replicate jcvi's --basepair sizing rule: size[seqid] = max BED `end`
    (3rd column) among rows on that (already-renamed) seqid. Only tracks
    labels we actually care about (wanted_labels) -- cheap single pass."""
    sizes = {lab: 0 for lab in wanted_labels}
    with open(out_bed_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            seqid = parts[0]
            if seqid not in sizes:
                continue
            try:
                end = int(parts[2])
            except ValueError:
                continue
            if end > sizes[seqid]:
                sizes[seqid] = end
    return sizes


def compute_z_start_x(generic_labels, z_labels, w_labels, sizes):
    """Returns the natural (pre-alignment-patch) x position where this
    species' FIRST Z-labeled bar starts, replicating JCVI Track's exact
    xstart/ratio/gap formula -- or None if this species has no Z label."""
    if not z_labels:
        return None
    all_labels = generic_labels + z_labels + w_labels
    nseqids = len(all_labels)
    if nseqids == 0:
        return None
    total = sum(sizes.get(lab, 0) for lab in all_labels)
    if total <= 0:
        return None
    gap = GAP_DEFAULT
    if nseqids > MAX_SEQIDS:
        gap = min(gap, gap * MAX_SEQIDS / nseqids + 0.001)
    rpad = 1 - XEND
    span = 1 - XSTART - rpad - gap * (nseqids - 1)
    ratio = span / total
    generic_total = sum(sizes.get(lab, 0) for lab in generic_labels)
    return XSTART + ratio * generic_total + gap * len(generic_labels)


def main():
    if len(sys.argv) != 8:
        sys.exit(f"Usage: {sys.argv[0]} SPECIES BED FAI N Z_SCAFFOLDS W_SCAFFOLDS OUT_BED")
    species, bed_path, fai_path, n_str, z_arg, w_arg, out_bed = sys.argv[1:8]
    n = int(n_str)

    lengths = {}
    with open(fai_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                lengths[parts[0]] = int(parts[1])

    def parse_list(arg):
        return [s.strip() for s in arg.split(",") if s.strip()]

    z_scafs = parse_list(z_arg)
    w_scafs = parse_list(w_arg)

    for s in z_scafs + w_scafs:
        if s not in lengths:
            print(f"WARNING: [{species}] sex scaffold '{s}' not found in {fai_path} -- skipped", file=sys.stderr)
    z_scafs = [s for s in z_scafs if s in lengths]
    w_scafs = [s for s in w_scafs if s in lengths]

    # Which scaffolds actually have >=1 gene in this species' bed (see the
    # "Why gene-bearing" note above -- required for --basepair mode to be
    # able to size the chromosome bar at all).
    bed_seqids_with_genes = set()
    with open(bed_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts and parts[0]:
                bed_seqids_with_genes.add(parts[0])

    for s in z_scafs + w_scafs:
        if s not in bed_seqids_with_genes:
            print(f"WARNING: [{species}] sex scaffold '{s}' has ZERO annotated genes -- "
                  f"cannot be sized/plotted in --basepair mode, dropping it entirely "
                  f"(not just unlabeled -- excluded)", file=sys.stderr)
    z_scafs = [s for s in z_scafs if s in bed_seqids_with_genes]
    w_scafs = [s for s in w_scafs if s in bed_seqids_with_genes]

    # Top-N GENE-BEARING scaffolds by length = "the chromosomes". Skip (and
    # backfill past) any large scaffold with zero genes rather than crash.
    ranked_all = sorted(lengths.items(), key=lambda kv: kv[1], reverse=True)
    n_empty_skipped = sum(1 for seqid, _ in ranked_all[:n] if seqid not in bed_seqids_with_genes)
    ranked = [(seqid, ln) for seqid, ln in ranked_all if seqid in bed_seqids_with_genes]
    top_n = [seqid for seqid, _ in ranked[:n]]
    top_n_set = set(top_n)
    if n_empty_skipped:
        print(f"NOTE: [{species}] {n_empty_skipped} scaffold(s) among the naive top-{n} "
              f"by length had zero genes -- skipped and backfilled from smaller "
              f"gene-bearing scaffolds instead", file=sys.stderr)

    for s in z_scafs:
        if s not in top_n_set:
            print(f"WARNING: [{species}] Z scaffold '{s}' is not among the top {n} "
                  f"largest gene-bearing scaffolds -- force-including it anyway, labeled Z", file=sys.stderr)

    # Build label map. Z/W scaffolds never consume a numeric slot.
    sex_set = set(z_scafs) | set(w_scafs)
    generic = [s for s in top_n if s not in sex_set]
    # generic is already in descending-length order (from `ranked`).

    label_map = {}
    for i, s in enumerate(generic, start=1):
        label_map[s] = str(i)

    def label_series(scafs, base):
        # base, then base+"_2", base+"_3", ... for extra pieces (e.g. a main
        # W scaffold plus dozens of small unlocalized W-associated contigs --
        # Malpolon_monspessulanus alone has 22). Unlimited, unlike an
        # earlier version of this that capped out at 9 lettered suffixes.
        for i, s in enumerate(scafs):
            label_map[s] = base if i == 0 else f"{base}_{i + 1}"

    label_series(z_scafs, "Z")
    label_series(w_scafs, "W")

    # Ordered seqids.txt row: numeric labels by size, then Z-series, then W-series.
    ordered_labels = [label_map[s] for s in generic]
    ordered_labels += [label_map[s] for s in z_scafs]
    ordered_labels += [label_map[s] for s in w_scafs]
    print(",".join(ordered_labels))

    # Write the renamed plotting bed -- every gene row kept, seqid renamed
    # only for kept (top-N + Z + W) scaffolds; everything else passes through
    # unchanged (harmless -- it just won't be in the seqids.txt row above, so
    # JCVI silently skips any block touching it).
    n_renamed = 0
    n_total = 0
    with open(bed_path) as fh, open(out_bed, "w") as out:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 1:
                continue
            n_total += 1
            seqid = parts[0]
            if seqid in label_map:
                parts[0] = label_map[seqid]
                n_renamed += 1
            out.write("\t".join(parts) + "\n")

    print(f"[{species}] {out_bed}: {n_renamed}/{n_total} bed rows on a kept chromosome "
          f"({len(generic)} numbered + {len(z_scafs)} Z + {len(w_scafs)} W)", file=sys.stderr)

    # ---- Z-alignment bookkeeping (see module docstring) --------------------
    generic_labels = [label_map[s] for s in generic]
    z_labels = [label_map[s] for s in z_scafs]
    w_labels = [label_map[s] for s in w_scafs]
    sizes = compute_sizes(out_bed, set(generic_labels + z_labels + w_labels))
    z_start_x = compute_z_start_x(generic_labels, z_labels, w_labels, sizes)

    align_file = os.path.join(os.path.dirname(os.path.abspath(out_bed)), "z_alignment.tsv")
    if z_start_x is not None:
        with open(align_file, "a") as af:
            af.write(f"{species}\t{z_start_x:.10f}\n")
        print(f"[{species}] natural Z-start x = {z_start_x:.6f} (appended to {align_file})",
              file=sys.stderr)
    else:
        print(f"[{species}] no Z label -- skipped Z-alignment bookkeeping", file=sys.stderr)


if __name__ == "__main__":
    main()
