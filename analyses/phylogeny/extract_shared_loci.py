#!/usr/bin/env python3
"""
extract_shared_loci.py

Find BUSCO IDs with Status == "Complete" (single-copy) in every species
listed in --species-list, then write one multi-FASTA per shared locus
(one amino-acid sequence per species, header replaced with the plain
species name) into --out-dir.

IMPORTANT: sequences are pulled from each species' own genome-mode BUSCO
run (single_copy_busco_sequences/*.faa) -- i.e. derived directly from
the genome assembly via BUSCO's internal gene prediction, NOT from the
separate GeMoMa/BRAKER3/TSEBRA annotation pipeline.

Expected layout per species (BUSCO genome-mode output, assembly/07.2.BUSCO_genome.sh):
  <busco-base>/<species>/BUSCO_<species>/run_<lineage>/full_table.tsv
  <busco-base>/<species>/BUSCO_<species>/run_<lineage>/busco_sequences/
      single_copy_busco_sequences/<busco_id>.faa

Only standard library is used (no biopython/pandas dependency).
"""
import argparse
import os
import sys


def read_full_table(path):
    """Return dict: busco_id -> status string, from a BUSCO full_table.tsv."""
    status = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            busco_id, stat = fields[0], fields[1]
            status[busco_id] = stat
    return status


def read_fasta_single(path):
    """Read a single-sequence FASTA file and return just the sequence string."""
    seq_lines = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            seq_lines.append(line.strip())
    return "".join(seq_lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--busco-base", required=True,
                     help="Base dir containing <species>/BUSCO_<species>/run_<lineage>/")
    ap.add_argument("--species-list", required=True,
                     help="Text file, one species folder name per line")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--lineage", default="squamata_odb12")
    ap.add_argument("--min-species", type=int, default=None,
                     help="Minimum species required per locus (default: all species, "
                          "i.e. strict intersection)")
    args = ap.parse_args()

    with open(args.species_list) as fh:
        species = [l.strip() for l in fh if l.strip()]

    min_species = args.min_species or len(species)
    os.makedirs(args.out_dir, exist_ok=True)

    # 1) Determine which BUSCO IDs are Complete (single-copy) per species
    complete_ids = {}
    for sp in species:
        ft = os.path.join(args.busco_base, sp, f"BUSCO_{sp}",
                           f"run_{args.lineage}", "full_table.tsv")
        if not os.path.isfile(ft):
            sys.exit(f"ERROR: missing full_table.tsv for {sp}: {ft}")
        status = read_full_table(ft)
        complete_ids[sp] = {bid for bid, st in status.items() if st == "Complete"}
        print(f"{sp}: {len(complete_ids[sp])} Complete single-copy BUSCOs", file=sys.stderr)

    # 2) Shared set across (by default) all species
    all_ids = set()
    for ids in complete_ids.values():
        all_ids |= ids

    shared = []
    for bid in sorted(all_ids):
        n = sum(1 for sp in species if bid in complete_ids[sp])
        if n >= min_species:
            shared.append(bid)

    print(f"Shared Complete single-copy BUSCOs across "
          f"{min_species}/{len(species)} species: {len(shared)}", file=sys.stderr)

    if not shared:
        sys.exit("ERROR: no shared loci found -- check --lineage / full_table.tsv parsing")

    # 3) Write one multi-FASTA per shared locus (skip any locus where the
    #    .faa is unexpectedly missing/empty for a required species)
    written = 0
    for bid in shared:
        records = []
        ok = True
        for sp in species:
            faa = os.path.join(args.busco_base, sp, f"BUSCO_{sp}",
                                f"run_{args.lineage}", "busco_sequences",
                                "single_copy_busco_sequences", f"{bid}.faa")
            if not os.path.isfile(faa):
                ok = False
                break
            seq = read_fasta_single(faa)
            if not seq:
                ok = False
                break
            records.append((sp, seq))
        if not ok or len(records) < min_species:
            continue
        out_path = os.path.join(args.out_dir, f"{bid}.faa")
        with open(out_path, "w") as out:
            for sp, seq in records:
                out.write(f">{sp}\n{seq}\n")
        written += 1

    print(f"Wrote {written} per-locus FASTA files to {args.out_dir}", file=sys.stderr)

    with open(os.path.join(args.out_dir, "_species_used.txt"), "w") as fh:
        fh.write("\n".join(species) + "\n")


if __name__ == "__main__":
    main()
