#!/usr/bin/env python3
"""
remap_seqids.py -- rewrite the seqid (column 1) of a BED file using a
two-column mapping table (annotation_id<TAB>genome_id).

Why this exists: some assemblies are submitted jointly to ENA/GenBank and
NCBI RefSeq, and end up with TWO different accession systems for the same
chromosome -- e.g. the genome FASTA you have uses the ENA/GenBank accession
(">OZ187402.1 ... chromosome: 1"), while a RefSeq-sourced GTF/GFF3 for that
same assembly (look for a "#!genome-build-accession NCBI_Assembly:GCF_..."
header comment) uses RefSeq accessions instead (e.g. "NC_138548.1") in its
seqid column, for that very same chromosome. When that happens, gtf_to_bed.py
produces a BED whose seqid column matches nothing in the genome's .fai, and
every downstream step (flipping, W-scaffold extraction, DIAMOND/JCVI
synteny) silently sees ~0 matches for that species.

Where to get the mapping: NCBI publishes a sequence-report table for every
RefSeq assembly (GCF_...) that lists, per chromosome, both the RefSeq
accession and the GenBank/INSDC (ENA) accession side by side -- e.g. via
`datasets summary genome accession GCF_xxx --report sequence` (NCBI
`datasets` CLI) or the assembly's "*_assembly_report.txt" / sequence_report
on the NCBI FTP/website for that GCF accession. Build a 2-column TSV from
that: <RefSeq accession, as used in the GTF> <TAB> <GenBank/ENA accession,
as used in your genome FASTA>.

Usage: remap_seqids.py input.bed mapping.tsv output.bed
"""
import sys


def main():
    if len(sys.argv) != 4:
        sys.exit(f"Usage: {sys.argv[0]} input.bed mapping.tsv output.bed")
    in_bed, map_path, out_bed = sys.argv[1:4]

    mapping = {}
    with open(map_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2 or not parts[0].strip() or not parts[1].strip():
                continue
            mapping[parts[0].strip()] = parts[1].strip()

    n_total = 0
    n_mapped = 0
    unmapped_seen = set()
    with open(in_bed) as fh, open(out_bed, "w") as out:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            n_total += 1
            old_id = parts[0]
            new_id = mapping.get(old_id)
            if new_id is None:
                unmapped_seen.add(old_id)
                out.write(line if line.endswith("\n") else line + "\n")
                continue
            parts[0] = new_id
            out.write("\t".join(parts) + "\n")
            n_mapped += 1

    if unmapped_seen:
        preview = sorted(unmapped_seen)[:10]
        more = " ..." if len(unmapped_seen) > 10 else ""
        print(f"WARNING: {len(unmapped_seen)} distinct seqid(s) in {in_bed} had no "
              f"entry in {map_path} and were left UNCHANGED (likely still "
              f"mismatched vs the genome .fai): {preview}{more}", file=sys.stderr)
    print(f"{in_bed}: remapped {n_mapped}/{n_total} rows using {map_path} -> {out_bed}")


if __name__ == "__main__":
    main()
