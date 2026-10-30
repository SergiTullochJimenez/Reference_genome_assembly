#!/usr/bin/env python3
"""
gtf_to_bed.py -- convert a GTF's transcript-level features into the BED6
format jcvi.compara.synteny expects (seqid, start, end, accn, score, strand).

Replaces `jcvi.formats.gff bed --type=mRNA --key=ID`, which only understands GFF3 "ID=" attributes. GTF uses
`transcript_id "..."` instead (no "ID=" key at all), so that step would
silently produce an empty/near-empty BED against GTFs. This tries GTF's transcript_id first and falls back to GFF3's ID=
so it works against either format without the caller having to know which
one they have.

Usage: gtf_to_bed.py input.gtf output.bed
"""
import sys
import re

TRANSCRIPT_FEATURES = {"transcript", "mRNA"}

ID_RE_GTF = re.compile(r'transcript_id\s+"([^"]+)"')
ID_RE_GFF3 = re.compile(r'ID=([^;]+)')


def extract_id(attrs):
    m = ID_RE_GTF.search(attrs)
    if m:
        return m.group(1)
    m = ID_RE_GFF3.search(attrs)
    if m:
        return m.group(1).strip("|")
    return None


def main():
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} input.gtf output.bed")
    in_path, out_path = sys.argv[1], sys.argv[2]

    n_total = 0
    n_written = 0
    with open(in_path) as fh, open(out_path, "w") as out:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feat_type = parts[2]
            if feat_type not in TRANSCRIPT_FEATURES:
                continue
            n_total += 1
            seqid = parts[0]
            start = int(parts[3]) - 1  # GTF is 1-based inclusive -> BED 0-based half-open
            end = int(parts[4])
            strand = parts[6] if parts[6] in ("+", "-") else "."
            tx_id = extract_id(parts[8])
            if tx_id is None:
                continue
            out.write(f"{seqid}\t{start}\t{end}\t{tx_id}\t0\t{strand}\n")
            n_written += 1

    if n_written == 0:
        print(f"WARNING: {in_path} -> 0 transcript records written "
              f"(saw {n_total} transcript/mRNA lines). Check feature-type "
              f"names and attribute format in this GTF.", file=sys.stderr)
    else:
        print(f"{in_path}: wrote {n_written}/{n_total} transcript records -> {out_path}")


if __name__ == "__main__":
    main()
