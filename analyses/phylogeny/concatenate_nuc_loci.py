#!/usr/bin/env python3
"""
concatenate_nuc_loci.py

Concatenate per-locus, aligned+trimmed nucleotide FASTA files (one file
per BUSCO locus, one sequence per species) into a single supermatrix
FASTA plus a RAxML-style partition file for IQ-TREE.

Standard-library only. Skips (with a warning, not an abort) any locus
whose taxon set doesn't match the expected species list, or whose
alignment is zero-length after trimming -- same fail-soft behaviour as
the rest of this pipeline.
"""
import argparse
import glob
import os
import sys


def read_fasta(path):
    seqs = {}
    name = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--aligned-dir", required=True,
                     help="directory of per-locus aligned+trimmed nucleotide FASTAs")
    ap.add_argument("--species-list", required=True)
    ap.add_argument("--glob-pattern", default="*.fna")
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-partition", required=True)
    args = ap.parse_args()

    with open(args.species_list) as fh:
        species = [l.strip() for l in fh if l.strip()]
    species_set = set(species)

    supermatrix = {sp: [] for sp in species}
    partitions = []
    pos = 1
    n_used, n_skipped = 0, 0

    files = sorted(glob.glob(os.path.join(args.aligned_dir, args.glob_pattern)))
    if not files:
        sys.exit(f"ERROR: no files matching {args.glob_pattern} in {args.aligned_dir}")

    for path in files:
        locus = os.path.splitext(os.path.basename(path))[0]
        seqs = read_fasta(path)

        if set(seqs) != species_set:
            missing = species_set - set(seqs)
            extra = set(seqs) - species_set
            msg = f"{locus}: taxon set mismatch"
            if missing:
                msg += f" (missing: {', '.join(sorted(missing))})"
            if extra:
                msg += f" (unexpected: {', '.join(sorted(extra))})"
            print("WARNING:", msg, "-- skipping locus", file=sys.stderr)
            n_skipped += 1
            continue

        lengths = {len(s) for s in seqs.values()}
        if len(lengths) != 1:
            print(f"WARNING: {locus}: sequences have unequal length ({lengths}) -- skipping locus",
                  file=sys.stderr)
            n_skipped += 1
            continue
        locus_len = lengths.pop()
        if locus_len == 0:
            print(f"WARNING: {locus}: zero-length alignment -- skipping locus", file=sys.stderr)
            n_skipped += 1
            continue

        for sp in species:
            supermatrix[sp].append(seqs[sp])
        partitions.append((locus, pos, pos + locus_len - 1))
        pos += locus_len
        n_used += 1

    with open(args.out_fasta, "w") as out:
        for sp in species:
            out.write(f">{sp}\n{''.join(supermatrix[sp])}\n")

    with open(args.out_partition, "w") as out:
        for locus, start, end in partitions:
            out.write(f"DNA, {locus} = {start}-{end}\n")

    total_len = pos - 1
    print(f"Concatenated {n_used} loci ({n_skipped} skipped) into supermatrix: "
          f"{len(species)} taxa x {total_len} columns", file=sys.stderr)
    print(f"Wrote {args.out_fasta} and {args.out_partition}", file=sys.stderr)


if __name__ == "__main__":
    main()
