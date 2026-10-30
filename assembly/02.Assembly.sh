#!/bin/bash
# 02 De novo assembly with hifiasm v0.25 in ONT mode, using the Hi-C reads for phasing.
# Writes the primary (consensus) contigs and both haplotypes as FASTA.
set -euo pipefail

#### 0. Configuration ####
SPECIES="Species_name"
BASE_DIR="/path/to/project"
HIFIASM="hifiasm"                    # hifiasm >= 0.25 (built with --ont support)
THREADS=24

WD="${BASE_DIR}/${SPECIES}"
indir="${WD}/res/01.1.Trimming"
outdir="${WD}/res/02.1.Assembly"
mkdir -p ${outdir}

#### 1. Assembly ####
${HIFIASM} -t${THREADS} --ont -o ${outdir}/Assembly.asm ${indir}/filtered_ONT.fq.gz \
    --h1 ${indir}/filtered_Arima_HIC_1.fq.gz --h2 ${indir}/filtered_Arima_HIC_2.fq.gz --dual-scaf -l 1

#### 2. GFA -> FASTA ####
awk '/^S/{print ">"$2"\n"$3}' ${outdir}/Assembly.asm.hic.p_ctg.gfa | fold > ${outdir}/assembly_consensus.fa
awk '/^S/{print ">"$2"\n"$3}' ${outdir}/Assembly.asm.hic.hap1.p_ctg.gfa | fold > ${outdir}/assembly_hap1.fa
awk '/^S/{print ">"$2"\n"$3}' ${outdir}/Assembly.asm.hic.hap2.p_ctg.gfa | fold > ${outdir}/assembly_hap2.fa
