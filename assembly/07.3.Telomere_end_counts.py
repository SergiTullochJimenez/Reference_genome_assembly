"""07.3 Telomere completeness of chromosome-scale scaffolds.
Counts TTAGGG + CCCTAA copies in the first and last 10 kb of every sequence >= 5 Mb. An end is called telomeric when it
has >= 100 copies; a scaffold is T2T when both ends are telomeric.
Usage: python3 07.3.Telomere_end_counts.py <assembly.fa> <species> > <species>.telomere_end_counts.tsv
Columns: species, scaffold, length, copies in first 10 kb, copies in last 10 kb, copies in the last partial 10-kb window.
"""
import re
import sys

fa, sp = sys.argv[1], sys.argv[2]
MIN_LEN = 5_000_000
END = 10_000


def count(seq):
    return len(re.findall("TTAGGG", seq)) + len(re.findall("CCCTAA", seq))


def flush(name, buf):
    s = "".join(buf).upper()
    L = len(s)
    if L < MIN_LEN:
        return
    print(f"{sp}\t{name}\t{L}\t{count(s[:END])}\t{count(s[-END:])}\t{count(s[-(L % END or END):])}")


name, buf = None, []
with open(fa) as h:
    for line in h:
        if line.startswith(">"):
            if name:
                flush(name, buf)
            name = line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if name:
        flush(name, buf)
