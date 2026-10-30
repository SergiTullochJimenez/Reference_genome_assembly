#!/usr/bin/env python3
"""
flip_bed.py -- mirror coordinates (and swap strand) for selected scaffolds
in a BED file, using a .fai for scaffold lengths.

Display only: so that species whose
chromosome numbering/orientation runs opposite to their neighbor's don't
show needlessly crossed synteny lines in the final plot.

For each flipped scaffold: new_start = length - old_end,
                            new_end   = length - old_start
strand is swapped (+ <-> -); scaffolds not listed pass through unchanged.

Usage: flip_bed.py input.bed input.fa.fai scaffold1,scaffold2,... output.bed
"""
import sys


def main():
    if len(sys.argv) != 5:
        sys.exit(f"Usage: {sys.argv[0]} input.bed input.fa.fai "
                  f"scaffold1,scaffold2,... output.bed")
    bed_path, fai_path, flip_list, out_path = sys.argv[1:5]
    flip_set = {s.strip() for s in flip_list.split(",") if s.strip()}

    lengths = {}
    with open(fai_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                lengths[parts[0]] = int(parts[1])

    missing = flip_set - set(lengths)
    if missing:
        print(f"WARNING: scaffold(s) not found in {fai_path}, left "
              f"unflipped: {sorted(missing)}", file=sys.stderr)

    n_flipped = 0
    with open(bed_path) as fh, open(out_path, "w") as out:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                out.write(line)
                continue
            seqid, start, end, name, score, strand = parts[:6]
            if seqid in flip_set and seqid in lengths:
                length = lengths[seqid]
                start_i, end_i = int(start), int(end)
                new_start = length - end_i
                new_end = length - start_i
                new_strand = {"+": "-", "-": "+"}.get(strand, strand)
                out.write(f"{seqid}\t{new_start}\t{new_end}\t{name}\t{score}\t{new_strand}\n")
                n_flipped += 1
            else:
                out.write(line)

    print(f"{bed_path}: flipped {n_flipped} records on scaffold(s) "
          f"{sorted(flip_set)} -> {out_path}")


if __name__ == "__main__":
    main()
