#!/bin/bash
# 08 Mitochondrial genome from the ONT assembly + Illumina polishing (independent of the nuclear pipeline).
#   MitoHiFi 3.2.1 in contig mode (-c) on the hifiasm primary contigs: finds the mito contig(s), removes concatemer
#   copies, circularises, rotates to tRNA-Phe and annotates. Reference: Elaphe urartica (RefSeq NC_073569.1).
#   Polypolish 0.5.0 (2 rounds) with the first 60 M fastp-filtered Illumina pairs, mapped with bwa mem -a and filtered
#   with polypolish_insert_filter.py; depth along the sequence with samtools depth; re-annotation with MitoFinder 1.4.
# The mito scaffold is removed from the nuclear assembly in 05.1.Remove_mito_scaffold.sh.
set -euo pipefail

#### 0. Configuration ####
SPECIES="Species_name"
BASE_DIR="/path/to/project"
CONTIGS="${BASE_DIR}/${SPECIES}/res/02.1.Assembly/assembly_consensus.fa"   # hifiasm primary contigs
READS1="/path/to/fastp/${SPECIES}.1.fastq.gz"                                # fastp-filtered WGS Illumina
READS2="/path/to/fastp/${SPECIES}.2.fastq.gz"
REF_GB="/path/to/mito_reference/NC_073569.gb"                                # Elaphe urartica, RefSeq
REF_FA="/path/to/mito_reference/NC_073569.fasta"
SIF="/path/to/mitohifi.sif"                                                  # MitoHiFi 3.2.1 container
BIND_DIR="/path/to"                                                          # directory to bind into the container
MITOFINDER="/path/to/MitoFinder/mitofinder"
INSERT_FILTER="/path/to/polypolish/polypolish_insert_filter.py"
THREADS=20
NPAIRS=60000000

WD="${BASE_DIR}/${SPECIES}"
outdir="${WD}/res/08.Mitogenome"
mkdir -p ${outdir}/01.mitohifi ${outdir}/02.polish ${outdir}/03.annotation
scratch=${outdir}/02.polish/reads
mkdir -p ${scratch}
trap 'rm -rf ${scratch} ${outdir}/02.polish/*.sam' EXIT

#### 1. Reference FASTA from the GenBank file ####
if [ ! -s ${REF_FA} ]; then
    singularity exec -B ${BIND_DIR} ${SIF} python -c \
        "from Bio import SeqIO; SeqIO.convert('${REF_GB}','genbank','${REF_FA}','fasta')"
fi

#### 2. MitoHiFi on the assembly contigs ####
cd ${outdir}/01.mitohifi
if [ ! -s final_mitogenome.fasta ]; then
    singularity exec -B ${BIND_DIR} ${SIF} mitohifi.py -c ${CONTIGS} \
        -f ${REF_FA} -g ${REF_GB} -t ${THREADS} -a animal -o 2
fi
[ -s final_mitogenome.fasta ] || { echo "WARNING: MitoHiFi produced no final_mitogenome.fasta"; exit 1; }
cat contigs_stats.tsv

#### 3. Illumina polishing (Polypolish, 2 rounds) ####
cd ${outdir}/02.polish
draft=${outdir}/01.mitohifi/final_mitogenome.fasta
if [ ! -s round2.fasta ]; then
    # round 1 streams the first NPAIRS pairs into bwa and keeps mapped records; the mito-mapped reads are reused in round 2
    for round in 1 2; do
        bwa index ${draft} 2> /dev/null
        for m in 1 2; do
            R=READS${m}
            if [ ${round} -eq 1 ]; then
                { zcat ${!R} || true; } | head -n $((NPAIRS*4)) \
                    | bwa mem -t ${THREADS} -a ${draft} /dev/stdin 2> /dev/null \
                    | samtools view -h -F 4 - > aln_${m}.sam
            else
                bwa mem -t ${THREADS} -a ${draft} ${scratch}/mito_${m}.fq 2> /dev/null > aln_${m}.sam
            fi
        done
        python3 ${INSERT_FILTER} --in1 aln_1.sam --in2 aln_2.sam --out1 filt_1.sam --out2 filt_2.sam
        polypolish ${draft} filt_1.sam filt_2.sam > round${round}.fasta 2> round${round}.log
        if [ ${round} -eq 1 ]; then
            for m in 1 2; do samtools fastq -F 0x904 aln_${m}.sam > ${scratch}/mito_${m}.fq 2> /dev/null; done
            # keep only reads whose mate also mapped, so the two files stay paired
            python3 - ${scratch}/mito_1.fq ${scratch}/mito_2.fq <<'EOF'
import sys
def read(f):
    d = {}
    with open(f) as h:
        while True:
            rec = [h.readline() for _ in range(4)]
            if not rec[0]: break
            d[rec[0].split()[0].split('/')[0]] = rec
    return d
a, b = read(sys.argv[1]), read(sys.argv[2])
keep = [k for k in a if k in b]
for f, d in ((sys.argv[1], a), (sys.argv[2], b)):
    with open(f, 'w') as o:
        for k in keep: o.writelines(d[k])
print(f"mito-mapped pairs kept: {len(keep)}")
EOF
        fi
        draft=${outdir}/02.polish/round${round}.fasta
    done
    rm -f filt_*.sam
fi
sed "1s/.*/>${SPECIES}_mitochondrion/" round2.fasta > ${SPECIES}_mitogenome.fasta

#### 4. Short-read depth along the polished sequence ####
if [ ! -s coverage.txt ]; then
    bwa index ${SPECIES}_mitogenome.fasta 2> /dev/null
    bwa mem -t ${THREADS} ${SPECIES}_mitogenome.fasta ${scratch}/mito_1.fq ${scratch}/mito_2.fq 2> /dev/null \
        | samtools sort -@4 -o final.bam - && samtools index final.bam
    samtools depth -a final.bam > coverage.txt
fi
awk '{s+=$3; if(min==""||$3<min)min=$3; if($3<20)low++} END{printf "Illumina depth: mean %.1f, min %d, positions <20x: %d of %d\n", s/NR, min, low, NR}' coverage.txt

#### 5. Annotation of the polished sequence (MitoFinder) ####
cd ${outdir}/03.annotation
if [ ! -d ${SPECIES}_polished ]; then
    ${MITOFINDER} -j ${SPECIES}_polished -a ${outdir}/02.polish/${SPECIES}_mitogenome.fasta \
        -r ${REF_GB} -o 2 -p ${THREADS} --override
fi
final=$(find ${outdir}/03.annotation -type d -name "*Final_Results" | head -1)
cat ${final}/*.infos
cat ${final}/*.gb | awk '$1=="CDS"||$1=="tRNA"||$1=="rRNA"{print $1}' | sort | uniq -c
ncds=$(cat ${final}/*.gb | awk '$1=="CDS"' | wc -l)
[ "${ncds}" -eq 13 ] || echo "WARNING: ${ncds} protein-coding genes annotated (expected 13)"
