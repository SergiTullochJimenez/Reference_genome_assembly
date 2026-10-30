#!/bin/bash
# PSMC (Li & Durbin 2011) per species from the BAMs of 01.Mapping.sh.
#   Regions: autosomes only -- W and Z removed by name (W_CHROM_FILE / Z_CHROM_FILE) and only the N largest
#   remaining scaffolds kept (N from CHROM_COUNT_FILE, W removed before ranking).
#   Consensus: bcftools mpileup -C50 | bcftools call -c | vcfutils.pl vcf2fq -d <d> -D <D>
#   (d = mean depth / 3, D = 2 x mean depth, per species in samples_All.txt: SAMPLE d D mu g); fq2psmcfa -q 20.
#   psmc -N25 -t25 -r4 -p "4*1+30*2+4+6+10", 20 bootstrap replicates (splitfa), psmc_plot.pl with mu and g.
# Every step writes to <file>.tmp and is renamed only on success, so reruns skip finished steps safely.
# Run: bash 02.PSMC.sh
set -euo pipefail

WD="/path/to/project"
src=/path/to/software/psmc/utils
psmc_bin="$(dirname "${src}")/psmc"

samples_file=${WD}/PSMCs/data/samples_All.txt

CHROM_COUNT_FILE=${WD}/PSMCs/data/chromosome_count.txt
Z_CHROM_FILE=${WD}/PSMCs/data/z_chromosome.txt
W_CHROM_FILE=${WD}/PSMCs/data/w_chromosome.txt

PSMC_N=25
PSMC_T=25
PSMC_R=4
PSMC_P="4*1+30*2+4+6+10"
N_BOOTSTRAP=20

MAX_THREADS=28
N_SPECIES_PARALLEL=3
N_BOOT_PARALLEL=$(( MAX_THREADS / N_SPECIES_PARALLEL ))   # 36/3 = 12

outdir1=${WD}/PSMCs/v2/01.filtered
outdir2=${WD}/PSMCs/v2/02.psmcfa
outdir3=${WD}/PSMCs/v2/03.psmc_bootstrap
outdir_plot=${WD}/PSMCs/v2/04.psmc_plot
bed_dir=${WD}/PSMCs/v2/data/autosome_beds
logdir=${WD}/PSMCs/v2/logs
mkdir -p "$outdir1" "$outdir2" "$outdir3" "$outdir_plot" "$bed_dir" "$logdir"

LOG_FILE="${logdir}/run_$(date '+%Y%m%d_%H%M%S').log"
log() { local msg="[$(date '+%H:%M:%S')] $*"; echo "$msg"; echo "$msg" >> "${LOG_FILE}"; }
trap 'log "ERROR: exited with status $? at line ${LINENO} (command: ${BASH_COMMAND})"' ERR

log "=== Definitive_script_v2.sh starting -- log file: ${LOG_FILE} ==="
log "Parallelism budget: MAX_THREADS=${MAX_THREADS}, N_SPECIES_PARALLEL=${N_SPECIES_PARALLEL}, N_BOOT_PARALLEL=${N_BOOT_PARALLEL}"
log "PSMC parameters: -N${PSMC_N} -t${PSMC_T} -r${PSMC_R} -p \"${PSMC_P}\", N_BOOTSTRAP=${N_BOOTSTRAP}"

build_autosome_bed () {
    local SAMPLE=$1
    local REF=$2
    local BED=${bed_dir}/${SAMPLE}_autosomes.bed

    if [[ -s ${BED} ]]; then
        log "${SAMPLE}: autosome BED already exists -- skipping build"
        return
    fi

    if [[ ! -s ${REF}.fai ]]; then
        log "${SAMPLE}: samtools faidx"
        samtools faidx "${REF}"
    fi

    local N
    N=$(awk -F'\t' -v sp="${SAMPLE}" '$1==sp && $1!~/^#/ {print $2; found=1} END{if(!found) exit 1}' "${CHROM_COUNT_FILE}") \
        || { echo "ERROR: ${SAMPLE}: no entry in ${CHROM_COUNT_FILE} -- add it (same values as analyses/macrosynteny/chromosome_count.txt) before running" >&2; exit 1; }

    local z_names
    z_names=$(awk -F'\t' -v sp="${SAMPLE}" '$1==sp && $1!~/^#/ {print $2}' "${Z_CHROM_FILE}")
    if [[ -z "${z_names}" ]]; then
        echo "ERROR: ${SAMPLE}: no entry in ${Z_CHROM_FILE} -- cannot safely build an autosome-only BED without knowing the Z scaffold(s)" >&2
        exit 1
    fi

    local w_names
    w_names=$(awk -F'\t' -v sp="${SAMPLE}" '$1==sp && $1!~/^#/ {print $2}' "${W_CHROM_FILE}")

    python3 - "${REF}.fai" "${N}" "${w_names}" "${z_names}" "${BED}" "${SAMPLE}" <<'PYEOF'
import sys

fai_path, n_str, w_raw, z_raw, bed_path, sample = sys.argv[1:7]
n = int(n_str)
w_set = set(w_raw.split(',')) if w_raw else set()
z_set = set(z_raw.split(',')) if z_raw else set()

rows = []
with open(fai_path) as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        name, length = parts[0], int(parts[1])
        if name in w_set:
            continue
        rows.append((name, length))

rows.sort(key=lambda r: r[1], reverse=True)
top_n = rows[:n]
kept = [r for r in top_n if r[0] not in z_set]

n_z_removed = len(top_n) - len(kept)
if n_z_removed == 0:
    sys.stderr.write(
        f"WARNING: {sample}: none of the top {n} scaffolds matched a Z name "
        f"from z_chromosome.txt ({sorted(z_set)}) -- Z may not have been "
        f"excluded (or is genuinely outside the top {n}, which would be "
        f"unusual). Check the BED before trusting this run.\n"
    )

with open(bed_path, "w") as out:
    for name, length in kept:
        out.write(f"{name}\t0\t{length}\n")

sys.stderr.write(
    f"{sample}: kept {len(kept)} autosome scaffolds "
    f"(top {n} by length after excluding {len(w_set)} W scaffold(s), "
    f"then excluding {n_z_removed} Z scaffold(s))\n"
)
PYEOF

    if [[ ! -s ${BED} ]]; then
        echo "ERROR: ${SAMPLE}: empty autosome BED -- check chromosome_count.txt/z_chromosome.txt/w_chromosome.txt entries" >&2
        exit 1
    fi
}

while read -r SAMPLE d D UVAL GVAL _REST; do
    [[ -z "${SAMPLE:-}" || "${SAMPLE}" =~ ^# ]] && continue

    if [[ -s ${outdir2}/${SAMPLE}_diploid.psmcfa && -s ${outdir3}/${SAMPLE}_bs.psmc && -s ${outdir_plot}/${SAMPLE}.plot_done ]]; then
        log "Skipping ${SAMPLE} (already fully processed under v2)"
        continue
    fi

(
    log "=== ${SAMPLE} (u=${UVAL}, g=${GVAL}) ==="

    indir=${WD}/PSMCs/res/03.2.polished
    bam=${indir}/${SAMPLE}_mkdup_30.bam
    ref=${WD}/PSMCs/data/genomes/${SAMPLE}.fa

    if [[ ! -s ${bam} ]]; then
        echo "ERROR: ${SAMPLE}: no final BAM at ${bam} -- run Full_script_v2.sh first" >&2
        exit 1
    fi

    build_autosome_bed "${SAMPLE}" "${ref}"
    bed=${bed_dir}/${SAMPLE}_autosomes.bed

    diploid_fq=${outdir1}/${SAMPLE}_diploid.fq.gz
    if [[ -s ${diploid_fq} ]]; then
        log "${SAMPLE}: diploid consensus fastq already present -- skipping"
    else
        log "${SAMPLE}: building diploid consensus fastq"
        bcftools mpileup -C50 -Ou -R "${bed}" -f "${ref}" "${bam}" | \
            bcftools call -c | \
            vcfutils.pl vcf2fq -d "${d}" -D "${D}" | \
            gzip > "${diploid_fq}.tmp"
        mv "${diploid_fq}.tmp" "${diploid_fq}"
    fi

    psmcfa=${outdir2}/${SAMPLE}_diploid.psmcfa
    if [[ -s ${psmcfa} ]]; then
        log "${SAMPLE}: psmcfa already present -- skipping fq2psmcfa"
    else
        log "${SAMPLE}: converting to psmcfa"
        "${src}/fq2psmcfa" -q 20 "${diploid_fq}" > "${psmcfa}.tmp"
        mv "${psmcfa}.tmp" "${psmcfa}"
    fi

    main_psmc=${outdir2}/${SAMPLE}.psmc
    if [[ -s ${main_psmc} ]]; then
        log "${SAMPLE}: main PSMC run already present -- skipping"
    else
        log "${SAMPLE}: running main PSMC (N${PSMC_N} t${PSMC_T} r${PSMC_R} p\"${PSMC_P}\")"
        ( cd "${outdir2}" && "${psmc_bin}" -N"${PSMC_N}" -t"${PSMC_T}" -r"${PSMC_R}" -p "${PSMC_P}" \
            -o "${main_psmc}.tmp" "$(basename "${psmcfa}")" )
        mv "${main_psmc}.tmp" "${main_psmc}"
    fi

    split_psmcfa=${outdir3}/${SAMPLE}_split.psmcfa
    if [[ ! -s ${split_psmcfa} ]]; then
        log "${SAMPLE}: splitting psmcfa for bootstrapping"
        "${src}/splitfa" "${psmcfa}" > "${split_psmcfa}.tmp"
        mv "${split_psmcfa}.tmp" "${split_psmcfa}"
    fi

    for i in $(seq 1 "${N_BOOTSTRAP}"); do
        boot_out=${outdir3}/${SAMPLE}_round_${i}.psmc
        if [[ -s ${boot_out} ]]; then
            continue
        fi
        (
            cd "${outdir3}"
            "${psmc_bin}" -N"${PSMC_N}" -b -t"${PSMC_T}" -r"${PSMC_R}" -p "${PSMC_P}" \
                -o "${boot_out}.tmp" "$(basename "${split_psmcfa}")"
            mv "${boot_out}.tmp" "${boot_out}"
        ) &
        if [[ $(jobs -r -p | wc -l) -ge ${N_BOOT_PARALLEL} ]]; then
            wait -n || true
        fi
    done
    wait

    n_boot_done=$(find "${outdir3}" -maxdepth 1 -name "${SAMPLE}_round_*.psmc" -size +0c | wc -l)
    if [[ "${n_boot_done}" -ne "${N_BOOTSTRAP}" ]]; then
        echo "ERROR: ${SAMPLE}: expected ${N_BOOTSTRAP} bootstrap replicates, found ${n_boot_done} non-empty -- check for a failed psmc -b run above" >&2
        exit 1
    fi

    combined=${outdir3}/${SAMPLE}_bs.psmc
    if [[ -s ${combined} ]]; then
        log "${SAMPLE}: combined bootstrap file already present -- skipping"
    else
        cat "${main_psmc}" "${outdir3}"/"${SAMPLE}"_round*.psmc > "${combined}.tmp"
        mv "${combined}.tmp" "${combined}"
    fi

    plot_copy=${outdir_plot}/${SAMPLE}_bs.psmc
    plot_done=${outdir_plot}/${SAMPLE}.plot_done
    if [[ -s ${plot_done} ]]; then
        log "${SAMPLE}: plot already present -- skipping"
    else
        cp "${combined}" "${plot_copy}"
        ( cd "${outdir_plot}" && "${src}/psmc_plot.pl" -u "${UVAL}" -g "${GVAL}" \
            -x 10000 -X 10000000 -Y 40 \
            "${SAMPLE}_bs.psmc" "$(basename "${plot_copy}")" )
        date '+%Y-%m-%d %H:%M:%S %Z' > "${plot_done}.tmp"
        mv "${plot_done}.tmp" "${plot_done}"
    fi

    log "${SAMPLE}: done"
) &

if [[ $(jobs -r -p | wc -l) -ge ${N_SPECIES_PARALLEL} ]]; then
    wait -n || true
fi

done < "${samples_file}"
wait

log "Definitive_script_v2.sh complete."
