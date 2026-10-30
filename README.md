# Reference genome assembly, annotation and comparative analyses of Iberian snakes

Scripts used to assemble and annotate eight chromosome-level snake genomes from ONT, Hi-C, Illumina and RNA-seq data,
and to compare them with seven published Iberian snake genomes.

Each script starts with a short description and a `0. Configuration` block. Set the paths there (`/path/to/...`),
the species name and the number of threads; nothing else needs editing. Scripts with `#SBATCH` lines are SLURM jobs
(replace `<account>` / `<partition>` and load the tools as appropriate on your system); the others run with `bash`.

## assembly/

| Script | Step |
|---|---|
| `01.1.Concatenate.sh` | Concatenate raw read files per data type |
| `01.2.Trimming.sh` | Filtlong (ONT), fastp (Hi-C) and Arima 5-bp trimming |
| `02.Assembly.sh` | hifiasm (ONT mode, Hi-C phasing) |
| `03.1.Pre-Scaffolding.sh` | Hi-C mapping following the Arima mapping pipeline, Picard MarkDuplicates |
| `03.2.Scaffolding.sh` | YaHS |
| `04.1.Pretext.sh` | Hi-C contact map and tracks (gaps, telomeres, coverage, repeats) for PretextView |
| `04.2.Pretext_to_fasta.sh` | Apply the manual curation (agp-tpf-utils + `helpers/rapid_join2.py`) |
| `05.1.Remove_mito_scaffold.sh` | Remove the mitochondrial scaffold(s) from the nuclear assembly |
| `05.2.Decontamination.sh` | BLASTn vs nt, BlobTools, keep scaffolds by phylum |
| `06.Masking.sh` | RepeatModeler + RepeatMasker (soft-masking) |
| `07.1.Merqury_GenomeScope.sh` | QV and k-mer completeness (Merqury), GenomeScope2 |
| `07.2.BUSCO_genome.sh` | BUSCO, genome mode |
| `07.3.Telomere_end_counts.py` | Telomeric repeats at chromosome ends (T2T calls) |
| `08.Mitogenome.sh` | Mitogenome: MitoHiFi (contig mode) + Polypolish + MitoFinder |
| `helpers/` | `fasta_gaps.py`, `repetitive_redmask.py`, `rapid_join2.py`, `two_read_bam_combiner.pl`, `fetch_taxids.py` |

External: `filter_five_end.pl` from the [Arima mapping pipeline](https://github.com/ArimaGenomics/mapping_pipeline);
`pretext_to_asm.py` from [agp-tpf-utils](https://github.com/sanger-tol/agp-tpf-utils).

## annotation/

| Script | Step |
|---|---|
| `01.RNA_alignment.sh` | HISAT2 alignment of the RNA-seq reads of each tissue |
| `02.BRAKER.sh` | BRAKER3 with the merged RNA-seq alignments and protein evidence |
| `03.GeMoMa.sh` | GeMoMa from three snake reference annotations |
| `04.TSEBRA.sh` + `tsebra_config.cfg` | Combine BRAKER3 and GeMoMa with TSEBRA |
| `05.1.Protein_search.sh` | Clean predicted proteins; DIAMOND blastp vs a Colubroidea protein database |
| `05.2.InterProScan.sh` | InterProScan |
| `05.3.AGAT.sh` | Merge the functional annotation into the final GFF3/GTF |
| `06.1.BUSCO_proteins.sh`, `06.2.OMArk.sh` | Annotation completeness and consistency |

## analyses/

| Folder | Content |
|---|---|
| `phylogeny/` | BUSCO-based nuclear phylogeny: protein and nucleotide supermatrices (MAFFT, trimAl, IQ-TREE 3) and ASTRAL |
| `macrosynteny/` | MCscan/JCVI macrosynteny across 15 genomes (config tables: Z/W scaffolds, chromosome counts, flips) |
| `breakpoints/` | Chromosome rearrangements vs *Zamenis longissimus*: pipeline, breakpoint reuse, breakpoint environment, Hi-C check of fusion junctions |
| `macro_micro/` | Macro- vs microchromosomes: per-chromosome table (00.1, 00.2), gene composition, orthology and dN/dS pipeline (01), dS-controlled omega tests (02) |
| `repeats/` | Kimura divergence landscapes and repeat content per chromosome partition (autosomes, Z, W) |
| `psmc/` | Read mapping (01) and PSMC (02) |
| `heterozygosity/` | All-sites VCFs (GATK), filters, genome-wide and window heterozygosity, macro vs micro |
| `circos/` | Circos tracks and genome plots |

## Differences from the scripts as run

- Reconstructed (not saved as a script; reproduces the layout of the outputs): concatenation of the lane files of each
  tissue in `annotation/01.RNA_alignment.sh`, and the scaffold removal in `assembly/05.1.Remove_mito_scaffold.sh`
  (done by hand with the `seqkit grep` command shown).
- Added safety check: `assembly/05.2.Decontamination.sh` stops if the accession-to-taxid map is empty.
- `analyses/heterozygosity/02.Filter_and_het.sh`: the runs-of-homozygosity step of the original script was removed
  (not used).
- `analyses/macro_micro/00.1` + `00.2` rebuild the per-chromosome table from the Circos tracks of the final assemblies in
  one pass. The published table was built from an earlier assembly version and updated for two curation edits of one
  genome, averaging the 1-Mb windows of the four affected chromosomes without length weighting; the rebuilt table
  differs only in those four chromosomes, and the pooled repeat-content test gives p = 1.5e-8 instead of 1.6e-8.
- Comments were shortened, and site-specific paths, accounts and module names were replaced by placeholders.

Third-party scripts in `assembly/helpers/` keep their original author credits.
