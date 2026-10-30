#!/usr/bin/env python3
"""
longest_isoform.py -- keep one protein (the longest) per gene, and filter
a matching BED to just those kept transcript IDs, so the protein FASTA and
BED stay 1:1 (jcvi.compara.synteny needs the BED "accn" column and the
protein FASTA headers to refer to the same transcript IDs).

Requires a gene_id <-> transcript_id mapping, taken from the same GTF used
by gtf_to_bed.py (GTF or GFF3 attributes, same ID-extraction logic).

TWO PROTEIN-HEADER CONVENTIONS HANDLED:
1. AUGUSTUS/GeMoMa-style (and NCBI-direct-accession-style): the protein
   FASTA header IS the transcript_id/ID used on the mRNA/transcript row
   (e.g. header "g1.t1" == transcript_id "g1.t1"). Matches directly.
2. Ensembl-style: the protein FASTA header is a *protein*-level stable ID
   plus version (e.g. "ENSUPKP00000001486.1"), built from the CDS row's
   `protein_id "ENSUPKP00000001486"` + `protein_version "1"` (GTF) or
   `protein_id=ENSUPKP00000001486;version=1` (GFF3) -- NOT the same string
   as transcript_id/ID at all. Confirmed against real data from
   Malpolon_monspessulanus, Natrix_maura, Hemorrhois_hippocrepis,
   Vipera_latastei, and Natrix_helvetica (all Ensembl GTF/GFF3 annotations)
   via diagnose_id_mismatch.py and check_proteins_vs_annotation.py.
   For these, this script reconstructs "{protein_id}.{version}" from each
   CDS row (linked to its transcript via the CDS row's own transcript_id
   (GTF) or Parent= (GFF3)) and uses that to resolve the protein header
   back to the transcript_id BED rows are keyed on.

Whichever convention matched, the OUTPUT protein FASTA header is always
rewritten to the transcript_id (matching the BED's accn column exactly) --
not left as the original protein-level ID -- since diamond blastp's query/
subject IDs come straight from the FASTA headers, and jcvi.compara.synteny
needs those to equal the BED accn column for anchors to resolve at all.

A protein whose header matches NEITHER convention is kept in the output
FASTA under its original header (not silently dropped), grouped under a
synthetic "__unmapped__<header>" pseudo-gene so genuinely orphaned data
isn't lost -- but note its BED row will NOT be found (the BED's accn is
never that header), so it will not participate in synteny. This is
reported explicitly, not left to be discovered as a mysterious mismatch
downstream.

Usage: longest_isoform.py input.gtf input.proteins.fa input.bed \
                           output.proteins.fa output.bed
"""
import sys
import re

TRANSCRIPT_FEATURES = {"transcript", "mRNA"}
CDS_FEATURES = {"CDS"}

TX_ID_RE_GTF = re.compile(r'transcript_id\s+"([^"]+)"')
ID_RE_GFF3 = re.compile(r'(?:^|;)ID=([^;]+)')
GENE_ID_RE_GTF = re.compile(r'gene_id\s+"([^"]+)"')
PARENT_RE_GFF3 = re.compile(r'(?:^|;)Parent=([^;]+)')

PROTEIN_ID_RE_GTF = re.compile(r'protein_id\s+"([^"]+)"')
PROTEIN_VERSION_RE_GTF = re.compile(r'protein_version\s+"([^"]+)"')
PROTEIN_ID_RE_GFF3 = re.compile(r'(?:^|;)protein_id=([^;]+)')
VERSION_RE_GFF3 = re.compile(r'(?:^|;)version=([^;]+)')


def extract(attrs, re_primary, re_fallback):
    m = re_primary.search(attrs)
    if m:
        return m.group(1)
    m = re_fallback.search(attrs)
    if m:
        return m.group(1).strip("|")
    return None


def build_maps(gtf_path):
    """Returns:
      tx_to_gene: {transcript_id_or_ID: gene_id}   -- from mRNA/transcript rows
      protein_to_tx: {"{protein_id}.{version}": transcript_id_or_ID}  -- from
        CDS rows, Ensembl-style versioned protein ID -> its transcript.
      protein_to_tx_bare: {protein_id: transcript_id_or_ID}  -- same, but
        without the version suffix, as a weaker fallback in case a proteins.fa
        was built without version suffixes.
    """
    tx_to_gene = {}
    protein_to_tx = {}
    protein_to_tx_bare = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feat, attrs = parts[2], parts[8]
            if feat in TRANSCRIPT_FEATURES:
                tx_id = extract(attrs, TX_ID_RE_GTF, ID_RE_GFF3)
                gene_id = extract(attrs, GENE_ID_RE_GTF, PARENT_RE_GFF3)
                if tx_id and gene_id:
                    tx_to_gene[tx_id] = gene_id
            elif feat in CDS_FEATURES:
                # Same linkage a CDS row has back to its own transcript:
                # GTF repeats transcript_id on every CDS line; GFF3 points
                # Parent= at the transcript's own ID= value.
                linked_tx = extract(attrs, TX_ID_RE_GTF, PARENT_RE_GFF3)
                protein_id = extract(attrs, PROTEIN_ID_RE_GTF, PROTEIN_ID_RE_GFF3)
                version = extract(attrs, PROTEIN_VERSION_RE_GTF, VERSION_RE_GFF3)
                if linked_tx and protein_id:
                    protein_to_tx_bare.setdefault(protein_id, linked_tx)
                    if version:
                        protein_to_tx.setdefault(f"{protein_id}.{version}", linked_tx)
    return tx_to_gene, protein_to_tx, protein_to_tx_bare


def resolve_header(header, tx_to_gene, protein_to_tx, protein_to_tx_bare):
    """Returns (resolved_tx_id, gene_id, via) for a protein FASTA header, or
    (None, None, None) if it matches neither convention. `via` is a short
    tag ("transcript_id", "protein_id+version", "protein_id") saying which
    rule resolved it, for reporting."""
    if header in tx_to_gene:
        return header, tx_to_gene[header], "transcript_id"
    if header in protein_to_tx:
        tx = protein_to_tx[header]
        return tx, tx_to_gene.get(tx), "protein_id+version"
    if header in protein_to_tx_bare:
        tx = protein_to_tx_bare[header]
        return tx, tx_to_gene.get(tx), "protein_id"
    return None, None, None


def read_fasta_records(fa_path):
    """Yields (header_id, [raw_body_lines_including_newline]) per record."""
    header_id = None
    body = []
    with open(fa_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header_id is not None:
                    yield header_id, body
                header_id = line[1:].split()[0]
                body = []
            else:
                body.append(line)
        if header_id is not None:
            yield header_id, body


def main():
    if len(sys.argv) != 6:
        sys.exit(f"Usage: {sys.argv[0]} input.gtf input.proteins.fa input.bed "
                  f"output.proteins.fa output.bed")
    gtf_path, in_faa, in_bed, out_faa, out_bed = sys.argv[1:6]

    tx_to_gene, protein_to_tx, protein_to_tx_bare = build_maps(gtf_path)

    # Resolve every protein header up front, and get its sequence length
    # (sum of body-line lengths) for the longest-isoform-per-gene choice.
    resolved = {}     # original_header -> (resolved_tx_id, gene_id, via, length)
    unresolved = []    # original_header -- matched neither convention
    records = {}       # original_header -> body lines (kept in memory once)
    n_via = {"transcript_id": 0, "protein_id+version": 0, "protein_id": 0}
    for header, body in read_fasta_records(in_faa):
        records[header] = body
        length = sum(len(l.strip()) for l in body)
        tx_id, gene_id, via = resolve_header(header, tx_to_gene, protein_to_tx, protein_to_tx_bare)
        if tx_id is None:
            unresolved.append(header)
            continue
        n_via[via] += 1
        if gene_id is None:
            # Resolved to a transcript_id the CDS rows know about, but that
            # transcript_id itself never appeared on an mRNA/transcript row
            # (unusual, e.g. malformed annotation) -- keep it, grouped alone.
            gene_id = f"__unmapped__{tx_id}"
        resolved[header] = (tx_id, gene_id, via, length)

    # For each gene, pick the transcript with the longest protein.
    best_per_gene = {}   # gene_id -> (header, tx_id, length)
    for header, (tx_id, gene_id, via, length) in resolved.items():
        if gene_id not in best_per_gene or length > best_per_gene[gene_id][2]:
            best_per_gene[gene_id] = (header, tx_id, length)

    # Orphan proteins (matched neither convention) are kept under their own
    # original header, each as its own pseudo-gene, so nothing is silently
    # dropped -- but they will NOT match any BED row (see docstring).
    for header in unresolved:
        best_per_gene[f"__unmapped__{header}"] = (header, None, 0)

    n_seq_written = 0
    with open(out_faa, "w") as out:
        for gene_id, (header, tx_id, length) in best_per_gene.items():
            out_header = tx_id if tx_id is not None else header
            out.write(f">{out_header}\n")
            out.writelines(records[header])
            n_seq_written += 1

    keep_tx_ids = {tx_id for (_h, tx_id, _l) in best_per_gene.values() if tx_id is not None}

    n_bed_written = 0
    with open(in_bed) as fh, open(out_bed, "w") as out:
        for line in fh:
            tx_id = line.split("\t")[3] if len(line.split("\t")) > 3 else None
            if tx_id in keep_tx_ids:
                out.write(line)
                n_bed_written += 1

    print(f"{in_faa}: {len(records)} proteins -> resolved "
          f"{n_via['transcript_id']} via transcript_id, "
          f"{n_via['protein_id+version']} via protein_id+version, "
          f"{n_via['protein_id']} via protein_id (no version) -- "
          f"{len(unresolved)} UNRESOLVED (kept under original header, will not "
          f"match any BED row). Kept {n_seq_written} proteins (longest isoform "
          f"per gene, {len(best_per_gene)} genes) -> {out_faa}")
    if n_bed_written != len(keep_tx_ids):
        print(f"WARNING: {len(keep_tx_ids)} distinct transcript IDs kept but only "
              f"matched {n_bed_written} BED rows for {in_bed} -- some kept "
              f"transcripts have no BED row (check gtf_to_bed.py output / "
              f"seqid_map coverage).", file=sys.stderr)
    if unresolved:
        print(f"WARNING: {len(unresolved)} protein header(s) matched neither "
              f"transcript_id nor protein_id -- examples: {unresolved[:10]}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
