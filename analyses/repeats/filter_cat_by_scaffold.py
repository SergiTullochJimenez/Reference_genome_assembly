#!/usr/bin/env python3
"""
filter_cat_by_scaffold.py -- filter a RepeatMasker .cat/.cat.gz alignment
file down to only the alignment blocks whose query scaffold is in a given
keep-list, so calcDivergenceFromAlign.pl can be run on (e.g.) just the
autosomes instead of the whole genome.

WHY THIS IS SAFE: each RepeatMasker alignment block starts with a
distinctive summary line -- four leading whitespace-separated numeric
fields (SW score, %div, %del, %ins), with the 5th field being the query
scaffold name, e.g.:

  1864 25.72 2.45 3.79 ENA_OZ184073_OZ184073.1 69 763 (310194260) rnd-1_family-362#LTR/Gypsy 1194 1879 (570) m_b1s001i0

No other line in the file (alignment body lines, the consensus "i/v/?"
diff line, or the trailing "Matrix =" / "Kimura (with divCpGMod) =" /
"Transitions / transversions =" / "Gap_init rate =" stats lines) starts
with four numeric fields like that, so this line shape unambiguously marks
the start of a new block. Everything between one such line and the next
(alignment body + its trailing stats lines, blank lines included) belongs
to that block and is kept or dropped as a whole, verbatim -- this rewrites
nothing, just removes entire blocks.

Verified against a real sample of Malpolon_monspessulanus's
.genome.fa.cat.gz (RepeatMasker 4.1.6) before being trusted here.

Usage:
  filter_cat_by_scaffold.py input.cat[.gz] scaffold_list.txt output.cat.gz

  scaffold_list.txt -- one scaffold name per line (as it appears as the
                       5th field of a block's header line, e.g. the .fai's
                       first column); blank lines and lines starting with #
                       are ignored. Typically produced by
                       build_autosome_scaffold_list.py.
  output.cat.gz     -- always written gzipped, so it's a drop-in
                       replacement for the original .cat.gz in a
                       calcDivergenceFromAlign.pl call.
"""
import sys
import re
import gzip

HEADER_RE = re.compile(r'^\s*\d+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+(\S+)\s')


def open_maybe_gzip(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def read_scaffold_list(path):
    keep = set()
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            keep.add(s)
    return keep


def main():
    if len(sys.argv) != 4:
        sys.exit(f"Usage: {sys.argv[0]} input.cat[.gz] scaffold_list.txt output.cat.gz")
    in_path, scaffold_list_path, out_path = sys.argv[1:4]

    keep = read_scaffold_list(scaffold_list_path)
    print(f"{scaffold_list_path}: {len(keep)} scaffolds to keep.", file=sys.stderr)

    n_blocks_total = 0
    n_blocks_kept = 0
    scaffolds_seen = set()
    scaffolds_seen_not_in_keep = set()

    current_block = []      # accumulated raw lines of the block being built
    current_keep = False     # whether current_block's scaffold is in `keep`

    def flush(out):
        nonlocal n_blocks_kept
        if current_block:
            if current_keep:
                out.writelines(current_block)
                n_blocks_kept += 1

    with open_maybe_gzip(in_path) as fh, gzip.open(out_path, "wt") as out:
        for line in fh:
            m = HEADER_RE.match(line)
            if m:
                # New block starting -- flush whatever we were accumulating.
                flush(out)
                n_blocks_total += 1
                scaffold = m.group(1)
                scaffolds_seen.add(scaffold)
                current_keep = scaffold in keep
                if not current_keep:
                    scaffolds_seen_not_in_keep.add(scaffold)
                current_block = [line]
            else:
                current_block.append(line)
        flush(out)  # last block in the file

    print(f"{in_path}: {n_blocks_total} alignment blocks seen, "
          f"{n_blocks_kept} kept -> {out_path}", file=sys.stderr)
    print(f"  distinct scaffolds seen: {len(scaffolds_seen)}", file=sys.stderr)
    unknown = scaffolds_seen_not_in_keep - keep
    dropped_scaffolds = scaffolds_seen - keep
    if dropped_scaffolds:
        preview = sorted(dropped_scaffolds)[:10]
        more = " ..." if len(dropped_scaffolds) > 10 else ""
        print(f"  {len(dropped_scaffolds)} scaffold(s) present in {in_path} but NOT "
              f"in {scaffold_list_path} -- all their blocks were dropped: "
              f"{preview}{more}", file=sys.stderr)
    unused_keep = keep - scaffolds_seen
    if unused_keep:
        print(f"  NOTE: {len(unused_keep)} scaffold(s) in {scaffold_list_path} had "
              f"zero alignment blocks in {in_path} (no repeats there, or a name "
              f"mismatch -- worth spot-checking a couple).", file=sys.stderr)
    if n_blocks_kept == 0:
        print(f"WARNING: 0 blocks kept -- {out_path} is effectively empty. Check "
              f"that scaffold names in {scaffold_list_path} exactly match the "
              f"5th field of {in_path}'s header lines.", file=sys.stderr)


if __name__ == "__main__":
    main()
