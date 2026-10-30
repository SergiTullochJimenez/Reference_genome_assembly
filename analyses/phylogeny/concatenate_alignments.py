#!/usr/bin/env python3
"""
concatenate_alignments.py

Concatenate a directory of per-locus trimmed alignments (same species
set, one sequence per species per locus, produced by MAFFT + trimAl)
into a single supermatrix FASTA plus a RAxML/IQ-TREE-style partition
file.

Any locus that doesn't have exactly one sequence per species after
trimming (e.g. trimAl dropped a taxon entirely, or the alignment ended
up length-0) is skipped with a warning rather than silently padding
with gaps -- that kind of silent padding can quietly corrupt a
concatenated matrix, so it's safer to skip and report.
"""
import argparse
import glob
import os
import sys
from collections import defaultdict


def read_fasta(path):
    seqs = {}
    name = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
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
    ap.add_argument("--aln-dir", required=True,
                     help="Directory of trimmed per-locus alignments (*.fasta)")
    ap.add_argument("--species-list", required=True)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-partition", required=True)
    args = ap.parse_args()

    with open(args.species_list) as fh:
        species = [l.strip() for l in fh if l.strip()]
    species_set = set(species)

    aln_files = sorted(glob.glob(os.path.join(args.aln_dir, "*.fasta")))
    if not aln_files:
        sys.exit(f"No alignment files found in {args.aln_dir}")

    supermatrix = defaultdict(str)
    partitions = []
    pos = 1
    n_used = 0
    n_skipped = 0

    for path in aln_files:
        seqs = read_fasta(path)
        locus = os.path.splitext(os.path.basename(path))[0]

        if set(seqs.keys()) != species_set:
            missing = species_set - set(seqs.keys())
            extra = set(seqs.keys()) - species_set
            print(f"WARNING: skipping {locus} -- taxon set mismatch "
                  f"(missing={sorted(missing)}, unexpected={sorted(extra)})",
                  file=sys.stderr)
            n_skipped += 1
            continue

        lengths = {len(s) for s in seqs.values()}
        if len(lengths) != 1:
            print(f"WARNING: skipping {locus} -- unequal sequence lengths "
                  f"after trimming ({lengths})", file=sys.stderr)
            n_skipped += 1
            continue

        aln_len = lengths.pop()
        if aln_len == 0:
            print(f"WARNING: skipping {locus} -- zero-length alignment "
                  f"(trimAl removed everything)", file=sys.stderr)
            n_skipped += 1
            continue

        for sp in species:
            supermatrix[sp] += seqs[sp]

        partitions.append(f"LG, {locus} = {pos}-{pos + aln_len - 1}")
        pos += aln_len
        n_used += 1

    if n_used == 0:
        sys.exit("ERROR: no usable loci survived concatenation -- check trimAl output")

    with open(args.out_fasta, "w") as out:
        for sp in species:
            out.write(f">{sp}\n{supermatrix[sp]}\n")

    with open(args.out_partition, "w") as out:
        out.write("\n".join(partitions) + "\n")

    print(f"Concatenated {n_used} loci (skipped {n_skipped}), "
          f"total length {pos - 1} aa, {len(species)} taxa -> {args.out_fasta}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
