#!/usr/bin/env python3
"""
extract_scaffold_subset.py -- filter a BED + matching protein FASTA down to
just the entries on one or more named scaffolds.

Used for the W-chromosome linking: rather than running a full genome-wide
comparison between two non-adjacent species just to get one W-to-W edge,
this pulls out only the genes on each species' named W scaffold(s) first,
so the diamond+jcvi comparison run on the result is tiny and produces ONLY
W-vs-W synteny blocks -- nothing from autosomes or Z, which stay connected
only through the normal staircase chain.

Some assemblies split the W across a main chromosome-level scaffold PLUS
many small unlocalized/unplaced W-associated contigs (e.g. a "chromosome: W"
record plus dozens of "SUPER_W_unloc_N" contigs). Pass ALL of those as one
comma-separated list in the scaffold_names argument to pull genes from every
one of them into the same W-only subset.

Usage: extract_scaffold_subset.py input.bed input.faa scaffold_name[,scaffold_name...] \
                                   output.bed output.faa
"""
import sys


def main():
    if len(sys.argv) != 6:
        sys.exit(f"Usage: {sys.argv[0]} input.bed input.faa "
                  f"scaffold_name[,scaffold_name...] output.bed output.faa")
    in_bed, in_faa, scaffold_arg, out_bed, out_faa = sys.argv[1:6]
    scaffolds = {s.strip() for s in scaffold_arg.split(",") if s.strip()}

    keep_ids = set()
    n_bed = 0
    seen_scaffolds = set()
    with open(in_bed) as fh, open(out_bed, "w") as out:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            if parts[0] in scaffolds:
                out.write(line)
                keep_ids.add(parts[3])
                n_bed += 1
                seen_scaffolds.add(parts[0])

    n_faa = 0
    with open(in_faa) as fh, open(out_faa, "w") as out:
        write_this = False
        for line in fh:
            if line.startswith(">"):
                header_id = line[1:].split()[0]
                write_this = header_id in keep_ids
                if write_this:
                    n_faa += 1
            if write_this:
                out.write(line)

    missing = scaffolds - seen_scaffolds
    if missing:
        print(f"WARNING: {len(missing)}/{len(scaffolds)} scaffold(s) had 0 "
              f"entries in {in_bed} (no annotated genes there, or a name "
              f"typo -- case-sensitive): {sorted(missing)}", file=sys.stderr)
    if n_bed == 0:
        print(f"WARNING: none of the given scaffolds have any entries in "
              f"{in_bed}.", file=sys.stderr)
    print(f"{in_bed} scaffolds={sorted(scaffolds)}: {n_bed} BED rows, "
          f"{n_faa} proteins -> {out_bed}, {out_faa}")


if __name__ == "__main__":
    main()
