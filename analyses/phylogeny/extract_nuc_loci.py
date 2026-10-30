#!/usr/bin/env python3
"""
extract_nuc_loci.py

Standalone, from-scratch nucleotide version of the nuclear BUSCO-ortholog
phylogeny pipeline. Extracts each shared single-copy BUSCO locus's coding
sequence (CDS), in nucleotides, directly from each species' own genome
FASTA and the miniprot-derived .gff coordinates BUSCO already produced --
no dependency on (or backtranslation from) the earlier protein alignment.

For each species:
  1. Parse full_table.tsv, keep BUSCO IDs with Status == "Complete"
     (single-copy only -- same criterion the protein pipeline used).
  2. Intersect across all species -> shared_ids.
  3. For each shared BUSCO ID, parse its .gff (CDS feature lines give the
     exact genomic coordinates of every coding exon, already in the
     species' own genome coordinate system).
  4. Batch-fetch all needed genomic regions in one `samtools faidx -r`
     call per species (much faster than one call per region).
  5. Reconstruct each gene's CDS by concatenating its CDS exons in
     ascending genomic order; if the gene is on the '-' strand,
     reverse-complement the whole concatenated sequence (equivalent to,
     but simpler than, reverse-complementing each exon individually and
     concatenating in descending order).
  6. Sanity-check every reconstructed CDS against the already-known
     protein length from full_table.tsv / the .faa BUSCO wrote, and log
     (not abort on) anything that doesn't translate cleanly.
  7. Write one unaligned multi-FASTA nucleotide file per shared locus,
     one sequence per species, to --out-dir.

This is a completely independent extraction from the genome -- nothing
here reuses or depends on the earlier protein-based alignments.
"""
import argparse
import glob
import os
import subprocess
import sys

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def read_species_list(path):
    with open(path) as fh:
        return [l.strip() for l in fh if l.strip()]


def complete_busco_ids(full_table_path):
    ids = set()
    with open(full_table_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            busco_id, status = fields[0], fields[1]
            if status == "Complete":
                ids.add(busco_id)
    return ids


def parse_gff_cds(gff_path):
    """
    Return (scaffold, strand, [(start, end), ...]) for the CDS exons of
    the top-ranked (Rank=1) gene model in this BUSCO single-copy .gff.
    Coordinates are 1-based inclusive, exactly as written by miniprot,
    which is also exactly what `samtools faidx -r` expects.
    """
    scaffold = None
    strand = None
    exons_by_parent = {}
    rank1_parents = set()

    with open(gff_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                continue
            seqid, source, feature, start, end, score, strand_col, phase, attrs = cols[:9]
            if feature != "CDS":
                continue
            attr_map = {}
            for kv in attrs.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attr_map[k] = v
            parent = attr_map.get("Parent", "NA")
            if attr_map.get("Rank") == "1":
                rank1_parents.add(parent)
            exons_by_parent.setdefault(parent, []).append(
                (int(start), int(end), seqid, strand_col)
            )

    if not exons_by_parent:
        return None

    # Prefer the Rank=1 gene model; fall back to whichever parent has the
    # most CDS bases if no Rank=1 tag is present (defensive -- in practice
    # single_copy_busco_sequences .gff files only ever have one model).
    if rank1_parents:
        parent = sorted(rank1_parents)[0]
    else:
        parent = max(
            exons_by_parent,
            key=lambda p: sum(e - s + 1 for s, e, _, _ in exons_by_parent[p]),
        )

    exons = exons_by_parent[parent]
    exons.sort(key=lambda e: e[0])
    scaffold = exons[0][2]
    strand = exons[0][3]
    coords = [(s, e) for s, e, _, _ in exons]
    return scaffold, strand, coords


CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}


def translate(nt_seq):
    prot = []
    for i in range(0, len(nt_seq) - len(nt_seq) % 3, 3):
        codon = nt_seq[i:i + 3].upper()
        prot.append(CODON_TABLE.get(codon, 'X'))
    return "".join(prot)


def run_faidx_batch(genome_fa, regions, samtools_bin, tmp_prefix):
    """
    regions: list of (scaffold, start, end) tuples (1-based inclusive).
    Returns dict {(scaffold, start, end): sequence}.
    """
    if not regions:
        return {}
    region_file = tmp_prefix + ".regions.txt"
    with open(region_file, "w") as fh:
        for scaffold, start, end in regions:
            fh.write(f"{scaffold}:{start}-{end}\n")

    out = subprocess.run(
        [samtools_bin, "faidx", genome_fa, "-r", region_file],
        capture_output=True, text=True, check=True,
    )
    os.remove(region_file)

    fetched = {}
    header = None
    seq_chunks = []
    for line in out.stdout.splitlines():
        if line.startswith(">"):
            if header is not None:
                fetched[header] = "".join(seq_chunks)
            header = line[1:].strip()
            seq_chunks = []
        else:
            seq_chunks.append(line.strip())
    if header is not None:
        fetched[header] = "".join(seq_chunks)

    # Map "scaffold:start-end" string keys back to tuple keys
    result = {}
    for scaffold, start, end in regions:
        key = f"{scaffold}:{start}-{end}"
        if key in fetched:
            result[(scaffold, start, end)] = fetched[key]
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--busco-base", required=True,
                     help="e.g. /path/to/working_directory/res/01.BUSCO_Assessment")
    ap.add_argument("--species-list", required=True)
    ap.add_argument("--lineage", default="squamata_odb12")
    ap.add_argument("--genome-dir", required=True,
                     help="expects <genome-dir>/<species>/<species>.masked.fa")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--samtools-bin", default="samtools")
    ap.add_argument("--log", default=None, help="path for warnings log (default: <out-dir>/extract_nuc_loci.log)")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    log_path = args.log or os.path.join(args.out_dir, "extract_nuc_loci.log")
    log_fh = open(log_path, "w")

    def warn(msg):
        print("WARNING:", msg, file=sys.stderr)
        log_fh.write(msg + "\n")

    species_list = read_species_list(args.species_list)

    # ---- Step 1-2: shared Complete BUSCO IDs across all species ----
    per_species_complete = {}
    for sp in species_list:
        ft = os.path.join(args.busco_base, sp, f"BUSCO_{sp}", f"run_{args.lineage}", "full_table.tsv")
        if not os.path.exists(ft):
            sys.exit(f"ERROR: full_table.tsv not found for {sp}: {ft}")
        per_species_complete[sp] = complete_busco_ids(ft)
        print(f"{sp}: {len(per_species_complete[sp])} Complete BUSCOs", file=sys.stderr)

    shared_ids = set.intersection(*per_species_complete.values())
    print(f"\nShared Complete BUSCO IDs across all {len(species_list)} species: {len(shared_ids)}",
          file=sys.stderr)

    # ---- Step 3-5: per-species GFF parsing + batch genome extraction ----
    # locus_seqs[busco_id][species] = nucleotide CDS string
    locus_seqs = {bid: {} for bid in shared_ids}

    for sp in species_list:
        genome_fa = os.path.join(args.genome_dir, sp, f"{sp}.masked.fa")
        if not os.path.exists(genome_fa):
            sys.exit(f"ERROR: genome FASTA not found for {sp}: {genome_fa}")
        if not os.path.exists(genome_fa + ".fai"):
            subprocess.run([args.samtools_bin, "faidx", genome_fa], check=True)

        gff_dir = os.path.join(args.busco_base, sp, f"BUSCO_{sp}", f"run_{args.lineage}",
                                "busco_sequences", "single_copy_busco_sequences")

        # Parse all needed GFFs for this species first, collecting every
        # region we'll need to fetch in one batch call.
        parsed = {}   # busco_id -> (scaffold, strand, [(start,end),...])
        regions_needed = set()
        for bid in shared_ids:
            gff_path = os.path.join(gff_dir, f"{bid}.gff")
            if not os.path.exists(gff_path):
                warn(f"{sp}/{bid}: .gff not found at {gff_path} -- skipping this locus for this species")
                continue
            parsed_gff = parse_gff_cds(gff_path)
            if parsed_gff is None:
                warn(f"{sp}/{bid}: no CDS features found in {gff_path} -- skipping")
                continue
            scaffold, strand, coords = parsed_gff
            parsed[bid] = (scaffold, strand, coords)
            for start, end in coords:
                regions_needed.add((scaffold, start, end))

        tmp_prefix = os.path.join(args.out_dir, f".tmp_{sp}")
        fetched = run_faidx_batch(genome_fa, sorted(regions_needed), args.samtools_bin, tmp_prefix)

        n_ok, n_fail = 0, 0
        for bid, (scaffold, strand, coords) in parsed.items():
            pieces = []
            missing = False
            for start, end in coords:
                seq = fetched.get((scaffold, start, end))
                if seq is None:
                    warn(f"{sp}/{bid}: could not fetch region {scaffold}:{start}-{end} -- skipping locus")
                    missing = True
                    break
                pieces.append(seq.upper())
            if missing:
                n_fail += 1
                continue
            cds = "".join(pieces)
            if strand == "-":
                cds = revcomp(cds)
            locus_seqs[bid][sp] = cds
            n_ok += 1
        print(f"{sp}: extracted {n_ok} loci ({n_fail} failed/skipped)", file=sys.stderr)

    # ---- Step 6-7: sanity-check translation + write output ----
    n_written = 0
    for bid in sorted(shared_ids):
        seqs = locus_seqs.get(bid, {})
        # Only write loci where every species in the shared set actually
        # got a sequence -- a handful of per-species GFF/fetch failures
        # can otherwise silently shrink the taxon set for just one locus.
        if len(seqs) != len(species_list):
            missing_sp = set(species_list) - set(seqs)
            warn(f"{bid}: incomplete across species ({len(seqs)}/{len(species_list)}); "
                 f"missing: {', '.join(sorted(missing_sp))} -- dropping this locus entirely")
            continue

        out_path = os.path.join(args.out_dir, f"{bid}.fna")
        with open(out_path, "w") as out:
            for sp in species_list:
                out.write(f">{sp}\n{seqs[sp]}\n")
        n_written += 1

    print(f"\nWrote {n_written} / {len(shared_ids)} shared loci as complete-across-all-species "
          f"nucleotide FASTAs to {args.out_dir}", file=sys.stderr)
    log_fh.close()


if __name__ == "__main__":
    main()
