#!/usr/bin/env python3
"""
build_partition_scaffold_list.py -- for one species, build the scaffold
keep-list for ONE of three partitions: autosomes, Z, or W. Generalizes the
earlier autosomes-only build_autosome_scaffold_list.py to also produce the
Z-only and W-only lists needed for the three-way (autosomes/Z/W) repeat
landscape + composition figures.

Reads W scaffolds from w_chromosome.txt and Z scaffolds from
z_chromosome.txt -- same two-column format (species<TAB>comma-separated
scaffold names) used throughout this project. Comment lines (#) and blank
lines are ignored in both.

PARTITION BEHAVIOUR:
  autosomes -- every scaffold in the .fai EXCEPT the W (if known) and Z (if
               known) scaffolds. If a species has NO entry in
               z_chromosome.txt yet, its Z scaffold is NOT excluded (can't
               guess it) -- reported loudly to stderr and as a comment in
               the output file, exactly as build_autosome_scaffold_list.py
               did.
  Z         -- ONLY the scaffold(s) listed for this species in
               z_chromosome.txt. If the species has no entry (Z not yet
               resolved -- currently true for the 4 new-genome males still
               waiting on macrosynteny), the output is an empty scaffold
               list with a comment explaining why, NOT an error -- the
               caller (KimuraDiv_partitions.sh) checks for this and skips
               that species/partition combination.
  W         -- ONLY the scaffold(s) listed for this species in
               w_chromosome.txt. If the species has no entry, this is
               normal/expected for a male individual (W is structurally
               absent, not "unresolved") -- output is an empty list with a
               comment saying so.

Usage:
  build_partition_scaffold_list.py genome.fa.fai w_chromosome.txt z_chromosome.txt species {autosomes,Z,W} output.txt [--min-length N]

  --min-length N  (optional, autosomes partition only) also drop any
                  scaffold shorter than N bp. Default: no length filter.
"""
import sys


def read_sex_scaffolds(path, species):
    """Returns the set of scaffold names listed for `species` in a
    w_chromosome.txt/z_chromosome.txt-style file, or None if the species
    has no entry at all (distinct from an entry with an empty list)."""
    try:
        fh = open(path)
    except FileNotFoundError:
        return None
    try:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            sp, scaf_field = parts[0].strip(), parts[1].strip()
            if sp == species:
                return {s.strip() for s in scaf_field.split(",") if s.strip()}
    finally:
        fh.close()
    return None


def read_fai_scaffolds(fai_path):
    scaffolds = []  # preserve .fai order
    lengths = {}
    with open(fai_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            scaffolds.append(parts[0])
            lengths[parts[0]] = int(parts[1])
    return scaffolds, lengths


def build_autosomes(all_scaffolds, lengths, w_set, w_known, z_set, z_known,
                     species, w_path, z_path, min_length):
    if w_known:
        print(f"  W scaffold(s) for {species}: {sorted(w_set) if w_set else '(none -- male?)'}",
              file=sys.stderr)
    else:
        print(f"  NOTE: {species} has no entry in {w_path} -- W status unresolved, "
              f"nothing excluded for W.", file=sys.stderr)

    if z_known:
        print(f"  Z scaffold(s) for {species}: {sorted(z_set) if z_set else '(none -- unexpected)'}",
              file=sys.stderr)
    else:
        print(f"  WARNING: {species} has no entry in {z_path} -- Z status UNRESOLVED. "
              f"The output list will still include the Z scaffold since it can't be "
              f"excluded without knowing its name. This is NOT a true autosomes-only "
              f"list yet.", file=sys.stderr)

    exclude = w_set | z_set
    unmatched_exclude = exclude - set(all_scaffolds)
    if unmatched_exclude:
        print(f"  WARNING: {len(unmatched_exclude)} scaffold(s) named in "
              f"{w_path}/{z_path} for {species} were NOT found in the .fai -- "
              f"typo, or wrong .fai? {sorted(unmatched_exclude)}", file=sys.stderr)

    kept = []
    n_excluded_sex = 0
    n_excluded_short = 0
    for s in all_scaffolds:
        if s in exclude:
            n_excluded_sex += 1
            continue
        if min_length and lengths[s] < min_length:
            n_excluded_short += 1
            continue
        kept.append(s)

    header_lines = [
        f"# Autosome scaffold list for {species}, built by "
        f"build_partition_scaffold_list.py.",
        f"# Excluded {n_excluded_sex} sex-scaffold(s) (W: "
        f"{'known' if w_known else 'UNKNOWN'}, Z: {'known' if z_known else 'UNKNOWN'})"
        + (f" and {n_excluded_short} scaffold(s) shorter than {min_length}bp" if min_length else "")
        + ".",
    ]
    if not z_known:
        header_lines.append(
            f"# WARNING: Z scaffold NOT excluded (unresolved for {species} as of when "
            f"this was built) -- this list is NOT true autosomes-only yet. Rerun once "
            f"{z_path} has an entry for {species}.")

    print(f"\n{species}: kept {len(kept)}/{len(all_scaffolds)} scaffolds "
          f"(excluded {n_excluded_sex} sex-scaffold(s)"
          + (f", {n_excluded_short} short scaffold(s)" if min_length else "") + ")",
          file=sys.stderr)
    if not z_known:
        print(f"*** {species}: Z UNRESOLVED -- this is autosomes+Z, not true "
              f"autosomes, until z_chromosome.txt is filled in for this species. ***",
              file=sys.stderr)

    return kept, header_lines


def build_sex_partition(all_scaffolds, sex_set, sex_known, species, label, source_path):
    """label is 'Z' or 'W'. Returns (kept_scaffolds_in_fai_order, header_lines)."""
    if not sex_known:
        reason = ("no entry in z_chromosome.txt yet -- pending (e.g. still waiting on "
                   "macrosynteny for a male individual)") if label == "Z" else \
                  ("no entry in w_chromosome.txt -- normal/expected if this is a male "
                   "individual (ZW system: only females carry a W)")
        print(f"  {species}: {label} scaffold NOT identified yet ({reason}). "
              f"Writing an EMPTY {label} scaffold list.", file=sys.stderr)
        return [], [
            f"# {label}-only scaffold list for {species}: EMPTY.",
            f"# Reason: {reason}",
        ]

    unmatched = sex_set - set(all_scaffolds)
    if unmatched:
        print(f"  WARNING: {len(unmatched)} {label} scaffold(s) named in {source_path} "
              f"for {species} were NOT found in the .fai -- typo, or wrong .fai? "
              f"{sorted(unmatched)}", file=sys.stderr)

    # Preserve .fai order among the ones that do match
    kept = [s for s in all_scaffolds if s in sex_set]
    print(f"  {species}: {label} scaffold(s) found: {kept}", file=sys.stderr)

    header_lines = [
        f"# {label}-only scaffold list for {species}, built by "
        f"build_partition_scaffold_list.py from {source_path}.",
        f"# {len(kept)}/{len(sex_set)} named {label} scaffold(s) matched the .fai.",
    ]
    return kept, header_lines


def main():
    args = sys.argv[1:]
    min_length = 0
    if "--min-length" in args:
        i = args.index("--min-length")
        min_length = int(args[i + 1])
        del args[i:i + 2]

    if len(args) != 6:
        sys.exit(f"Usage: {sys.argv[0]} genome.fa.fai w_chromosome.txt "
                  f"z_chromosome.txt species {{autosomes,Z,W}} output.txt [--min-length N]")
    fai_path, w_path, z_path, species, partition, out_path = args

    if partition not in ("autosomes", "Z", "W"):
        sys.exit(f"ERROR: partition must be one of autosomes/Z/W, got '{partition}'")

    all_scaffolds, lengths = read_fai_scaffolds(fai_path)
    print(f"{fai_path}: {len(all_scaffolds)} scaffolds.", file=sys.stderr)

    w_scaffolds = read_sex_scaffolds(w_path, species)
    z_scaffolds = read_sex_scaffolds(z_path, species)
    w_known = w_scaffolds is not None
    z_known = z_scaffolds is not None
    w_set = w_scaffolds or set()
    z_set = z_scaffolds or set()

    if partition == "autosomes":
        kept, header_lines = build_autosomes(
            all_scaffolds, lengths, w_set, w_known, z_set, z_known,
            species, w_path, z_path, min_length)
    elif partition == "Z":
        kept, header_lines = build_sex_partition(
            all_scaffolds, z_set, z_known, species, "Z", z_path)
    else:  # W
        kept, header_lines = build_sex_partition(
            all_scaffolds, w_set, w_known, species, "W", w_path)

    with open(out_path, "w") as out:
        for line in header_lines:
            out.write(line + "\n")
        for s in kept:
            out.write(s + "\n")

    print(f"{species} [{partition}]: wrote {len(kept)} scaffold(s) -> {out_path}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
