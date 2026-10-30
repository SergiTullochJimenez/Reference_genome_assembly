#!/bin/bash
# PSMC / heterozygosity input BAMs, one per species listed in column 1 of samples_All.txt:
#   fastp v0.20 on the Illumina WGS reads (Q30, polyG/polyX trimming, adapter detection, base correction),
#   bwa mem v0.7.17 against the reference, Picard MarkDuplicates v2.22, then samtools view -q 30 -F 260
#   (mapping quality >= 30, unmapped and secondary alignments removed), index, and genome-wide mean depth
#   (samtools depth -a) written to data/coverage_v2.tsv.
# Expected layout: ${WD}/data/ILLUMINA/<Species>/ (raw reads, *_R1*/*_R2* or *_1/*_2), ${WD}/data/genomes/<Species>.fa
# Every step is skipped when its output already exists. Run: bash 01.Mapping.sh
set -euo pipefail

WD="/path/to/project/PSMCs"
samples_file="${WD}/data/samples_All.txt"   # shared with 02.PSMC.sh

THREADS_FASTP=20
THREADS_BWA=20
PICARD=/path/to/software/picard/picard.jar
SAMTOOLS=samtools

log() { echo "[$(date '+%H:%M:%S')] $*"; }

mkdir -p "${WD}/res/fastp" "${WD}/res/02.mapping" "${WD}/res/03.1.polished" \
         "${WD}/res/03.1.polished/tmp" "${WD}/res/03.2.polished" "${WD}/data"

coverage_file="${WD}/data/coverage_v2.tsv"
if [[ ! -f "${coverage_file}" ]]; then
    echo -e "Species\tMean_depth" > "${coverage_file}"
fi

find_read_files () {
    local species="$1" raw_dir="$2" mate="$3"
    {
        find "${raw_dir}" -maxdepth 1 -type f -name "*_R${mate}*.fastq.gz"
        find "${raw_dir}" -maxdepth 1 -type f -name "*_${mate}.fastq.gz"
    } | sort -u
}

while read -r SPECIES _REST; do
    [[ -z "${SPECIES:-}" || "${SPECIES}" =~ ^# ]] && continue

    log "=== ${SPECIES} ==="

    RAW_DIR="${WD}/data/ILLUMINA/${SPECIES}"
    REF="${WD}/data/genomes/${SPECIES}.fa"
    REF_GZ="${REF}.gz"

    final_bam="${WD}/res/03.2.polished/${SPECIES}_mkdup_30.bam"
    final_bai="${final_bam}.bai"

    already_done=false
    if [[ -s "${final_bam}" && -s "${final_bai}" ]]; then
        log "${SPECIES}: final filtered BAM + index already present -- skipping mapping pipeline"
        already_done=true
    fi

    if [[ "${already_done}" != "true" ]]; then

        if [[ ! -d "${RAW_DIR}" ]]; then
            log "${SPECIES}: no ${RAW_DIR} -- no Illumina data yet, skipping (not an error)"
            continue
        fi
        if ! find "${RAW_DIR}" -maxdepth 1 -type f -name "*.fastq.gz" -print -quit | grep -q .; then
            log "${SPECIES}: ${RAW_DIR} exists but has no .fastq.gz files yet -- skipping (not an error)"
            continue
        fi
        if [[ ! -f "${REF}" && ! -f "${REF_GZ}" ]]; then
            log "${SPECIES}: no reference genome found (${REF} or ${REF_GZ}) -- skipping"
            continue
        fi

        if [[ ! -s "${REF}" && -f "${REF_GZ}" ]]; then
            log "${SPECIES}: decompressing reference"
            gunzip -c "${REF_GZ}" > "${REF}"
        fi

        if [[ ! -s "${REF}.fai" ]]; then
            log "${SPECIES}: samtools faidx"
            ${SAMTOOLS} faidx "${REF}"
        else
            log "${SPECIES}: .fai already present -- skipping faidx"
        fi

        if [[ ! -s "${REF}.bwt" ]]; then
            log "${SPECIES}: bwa index"
            bwa index "${REF}"
        else
            log "${SPECIES}: bwa index already present -- skipping"
        fi

        merged_r1="${RAW_DIR}/${SPECIES}.1.fq.gz"
        merged_r2="${RAW_DIR}/${SPECIES}.2.fq.gz"

        if [[ -s "${merged_r1}" && -s "${merged_r2}" ]]; then
            log "${SPECIES}: merged FASTQs already present -- skipping merge"
        else
            r1_files=$(find_read_files "${SPECIES}" "${RAW_DIR}" 1)
            r2_files=$(find_read_files "${SPECIES}" "${RAW_DIR}" 2)

            n_r1=$(echo "${r1_files}" | grep -c . || true)
            n_r2=$(echo "${r2_files}" | grep -c . || true)

            if [[ "${n_r1}" -eq 0 || "${n_r2}" -eq 0 ]]; then
                echo "ERROR: ${SPECIES}: found ${n_r1} R1 file(s) and ${n_r2} R2 file(s) in ${RAW_DIR}" >&2
                echo "  matched by *_R1*.fastq.gz / *_R2*.fastq.gz -- check your raw filenames," >&2
                echo "  they may not follow that convention. Files present:" >&2
                find "${RAW_DIR}" -maxdepth 1 -type f -name "*.fastq.gz" >&2
                exit 1
            fi

            overlap=$(comm -12 <(echo "${r1_files}" | sort) <(echo "${r2_files}" | sort))
            if [[ -n "${overlap}" ]]; then
                echo "ERROR: ${SPECIES}: the following file(s) matched BOTH R1 and R2 patterns -- refusing to guess:" >&2
                echo "${overlap}" >&2
                exit 1
            fi

            log "${SPECIES}: merging ${n_r1} R1 file(s) + ${n_r2} R2 file(s):"
            echo "${r1_files}"  | sed 's/^/    R1: /'
            echo "${r2_files}"  | sed 's/^/    R2: /'

            echo "${r1_files}" | xargs cat > "${merged_r1}"
            echo "${r2_files}" | xargs cat > "${merged_r2}"
        fi

        fastp_r1="${WD}/res/fastp/${SPECIES}.1.fastq.gz"
        fastp_r2="${WD}/res/fastp/${SPECIES}.2.fastq.gz"

        if [[ -s "${fastp_r1}" && -s "${fastp_r2}" ]]; then
            log "${SPECIES}: fastp output already present -- skipping"
        else
            log "${SPECIES}: running fastp"
            fastp -i "${merged_r1}" \
                  -I "${merged_r2}" \
                  -o "${fastp_r1}" \
                  -O "${fastp_r2}" \
                  --detect_adapter_for_pe \
                  --adapter_sequence AGATCGGAAGAG \
                  --adapter_sequence_r2 AGATCGGAAGAG \
                  --correction --trim_poly_g --trim_poly_x -p \
                  --thread "${THREADS_FASTP}" --qualified_quality_phred 30 \
                  --json "${WD}/res/fastp/${SPECIES}.json" \
                  --html "${WD}/res/fastp/${SPECIES}.html"
        fi

        mapped_bam="${WD}/res/02.mapping/${SPECIES}.bam"

        if [[ -s "${mapped_bam}" ]]; then
            log "${SPECIES}: mapped BAM already present -- skipping bwa mem"
        else
            log "${SPECIES}: mapping reads"
            bwa mem -t "${THREADS_BWA}" -M "${REF}" \
                "${fastp_r1}" "${fastp_r2}" \
                -R "@RG\tID:${SPECIES}\tSM:${SPECIES}" | \
                ${SAMTOOLS} sort -@ 20 -o "${mapped_bam}"
        fi

        mkdup_bam="${WD}/res/03.1.polished/${SPECIES}_mkdup.bam"

        if [[ -s "${mkdup_bam}" ]]; then
            log "${SPECIES}: mkdup BAM already present -- skipping MarkDuplicates"
        else
            log "${SPECIES}: marking duplicates"
            java -jar "${PICARD}" MarkDuplicates \
                INPUT="${mapped_bam}" \
                OUTPUT="${mkdup_bam}" \
                METRICS_FILE="${WD}/res/03.1.polished/${SPECIES}.txt" \
                ASSUME_SORT_ORDER=coordinate REMOVE_DUPLICATES=TRUE \
                CREATE_INDEX=True TMP_DIR="${WD}/res/03.1.polished/tmp"
        fi

        if [[ -s "${final_bam}" ]]; then
            log "${SPECIES}: filtered BAM already present -- skipping filter step"
        else
            log "${SPECIES}: filtering BAM (q>=30)"
            ${SAMTOOLS} view -q 30 -F 260 -b \
                "${mkdup_bam}" \
                -o "${final_bam}"
        fi

        if [[ ! -s "${final_bai}" ]]; then
            ${SAMTOOLS} index "${final_bam}"
        fi

        if ! ${SAMTOOLS} quickcheck "${final_bam}"; then
            echo "ERROR: ${SPECIES}: samtools quickcheck failed on ${final_bam} -- BAM looks truncated/corrupt" >&2
            exit 1
        fi

    fi   # end of `if ! already_done` mapping-pipeline block

    if grep -q -P "^${SPECIES}\t" "${coverage_file}" 2>/dev/null; then
        log "${SPECIES}: already has a coverage_v2.tsv entry -- skipping depth calc"
    else
        log "${SPECIES}: computing whole-genome mean depth"
        mean_depth=$(${SAMTOOLS} depth -a "${final_bam}" | awk '{sum+=$3} END{print sum/NR}')
        echo -e "${SPECIES}\t${mean_depth}" >> "${coverage_file}"
        log "${SPECIES}: mean depth = ${mean_depth}"
    fi

    log "${SPECIES}: done"

done < "${samples_file}"

log "Full_script_v2.sh complete."
