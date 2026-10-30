#!/bin/bash
# 03.1 Map the Hi-C reads to the contigs following the Arima mapping pipeline
# (https://github.com/ArimaGenomics/mapping_pipeline): reads mapped separately, 5'-end filtered, paired,
# sorted, duplicates marked with Picard MarkDuplicates v2.22, converted to BED for YaHS.
# filter_five_end.pl comes from the Arima repository; two_read_bam_combiner.pl is in helpers/.
set -euo pipefail

#### 0. Configuration ####
SPECIES="Species_name"
BASE_DIR="/path/to/project"
ARIMA_DIR="/path/to/arima_mapping_pipeline"     # folder with filter_five_end.pl
HELPERS="/path/to/this/repo/assembly/helpers"   # two_read_bam_combiner.pl
PICARD_JAR="/path/to/picard.jar"
THREADS=20

WD="${BASE_DIR}/${SPECIES}"
indir="${WD}/res/02.1.Assembly"
outbase="${WD}/res/03.1.Pre-Scaffolding"
fq1="${WD}/res/01.1.Trimming/filtered_Arima_HIC_1.fq.gz"
fq2="${WD}/res/01.1.Trimming/filtered_Arima_HIC_2.fq.gz"

for version in consensus; do         # add hap1 hap2 to scaffold the haplotypes too
    ref="${indir}/assembly_${version}.fa"
    outdir="${outbase}/${version}"
    mkdir -p "${outdir}"
    echo "Processing version: ${version}"

    #### 1. Index and map R1 / R2 separately ####
    bwa index "${ref}"
    samtools faidx "${ref}"
    bwa mem -t ${THREADS} "${ref}" "${fq1}" | samtools view -@ ${THREADS} -Sb - > "${outdir}/forward.bam"
    bwa mem -t ${THREADS} "${ref}" "${fq2}" | samtools view -@ ${THREADS} -Sb - > "${outdir}/reverse.bam"

    #### 2. 5'-end filtering and pairing ####
    samtools view -h "${outdir}/forward.bam" | perl ${ARIMA_DIR}/filter_five_end.pl | samtools view -Sb - > ${outdir}/forward_filtered.bam
    samtools view -h "${outdir}/reverse.bam" | perl ${ARIMA_DIR}/filter_five_end.pl | samtools view -Sb - > ${outdir}/reverse_filtered.bam
    perl ${HELPERS}/two_read_bam_combiner.pl ${outdir}/forward_filtered.bam ${outdir}/reverse_filtered.bam \
        | samtools view -@ ${THREADS} -Sb > ${outdir}/combined.bam
    samtools sort -@ ${THREADS} -o ${outdir}/sorted.bam ${outdir}/combined.bam

    #### 3. Mark duplicates ####
    java -jar ${PICARD_JAR} MarkDuplicates \
        I=${outdir}/sorted.bam \
        O=${outdir}/combined.dedup.bam \
        M=${outdir}/combined.metrics.txt \
        ASSUME_SORT_ORDER=coordinate \
        REMOVE_DUPLICATES=false
    samtools index ${outdir}/combined.dedup.bam

    #### 4. BED for YaHS (sorted by read name) ####
    bamToBed -i ${outdir}/combined.dedup.bam > ${outdir}/combined.dedup.bed
    sort -k 4 ${outdir}/combined.dedup.bed > ${outdir}/tmp && mv ${outdir}/tmp ${outdir}/combined.dedup.bed

    echo "Done with version: ${version}"
done
