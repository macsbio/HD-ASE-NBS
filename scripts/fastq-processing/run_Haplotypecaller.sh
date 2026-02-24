#!/usr/bin/env bash
set -euo pipefail

# source "$(conda info --base)/etc/profile.d/conda.sh"
# conda activate gatk_env

REF="Homo_sapiens.GRCh37.chr1-22_X.fa"
IN_DIR="../STAR-output-last-try"
OUT_DIR="genotype_per_sample_files"
THREADS=4           # threads per HC process
MAX_JOBS=4          # number of samples run in parallel

mkdir -p "${OUT_DIR}"

process_sample() {
  local idx="$1"
  local SAMPLE="SRR${idx}"
  local IN_BAM="${IN_DIR}/${SAMPLE}_RG.bam"
  local OUT_VCF="${OUT_DIR}/${SAMPLE}.raw.vcf.gz"

  if [[ ! -f "${IN_BAM}" ]]; then
    echo "[WARN] Skipping ${SAMPLE}: BAM not found: ${IN_BAM}"
    return 0
  fi

  # ensure BAM index
  [[ -f "${IN_BAM}.bai" ]] || samtools index -@ "${THREADS}" "${IN_BAM}"

  echo "[INFO] HaplotypeCaller ${SAMPLE}"
  gatk HaplotypeCaller \
    --native-pair-hmm-threads "${THREADS}" \
    -R "${REF}" \
    -I "${IN_BAM}" \
    -O "${OUT_VCF}" \
    --dont-use-soft-clipped-bases true \
    --standard-min-confidence-threshold-for-calling 30.0 \
    --pair-hmm-implementation FASTEST_AVAILABLE
}

export -f process_sample
export THREADS REF IN_DIR OUT_DIR

seq 1747192 1747211 | parallel -j "${MAX_JOBS}" --halt soon,fail=1 process_sample {}

echo "[INFO] HaplotypeCaller finished for all samples. Outputs in ${OUT_DIR}"
