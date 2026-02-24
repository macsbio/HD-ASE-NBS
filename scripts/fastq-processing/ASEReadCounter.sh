#!/usr/bin/env bash
set -euo pipefail

# Optional: activate conda env if needed
# source "$(conda info --base)/etc/profile.d/conda.sh"
# conda activate gatk_env

REF="Homo_sapiens.GRCh37.chr1-22_X.fa"
BAM_DIR="../STAR-WASP-output"
VCF_DIR="genotype_sample_files"   # per-sample VCFs from HaplotypeCaller
OUT_DIR="ASE_results"

MAX_JOBS=4   # how many samples in parallel

mkdir -p "${OUT_DIR}"

process_sample() {
  local i="$1"
  local SAMPLE="SRR${i}"

  local IN_BAM="${BAM_DIR}/${SAMPLE}_WASP_RG.bam"
  # use bgzipped per-sample VCF (random-access + indexed with IndexFeatureFile)
  local VCF="${VCF_DIR}/${SAMPLE}.biallelic.snps.vcf.gz"
  local OUT_CSV="${OUT_DIR}/${SAMPLE}_ASEresults.csv"

  if [[ ! -f "${IN_BAM}" ]]; then
    echo "[WARN] Skipping ${SAMPLE}: BAM not found: ${IN_BAM}"
    return 0
  fi

  if [[ ! -f "${VCF}" ]]; then
    echo "[WARN] Skipping ${SAMPLE}: VCF not found: ${VCF}"
    return 0
  fi

  echo "[INFO] ASEReadCounter on ${SAMPLE}"
  gatk ASEReadCounter \
    -R "${REF}" \
    -I "${IN_BAM}" \
    -V "${VCF}" \
    -O "${OUT_CSV}" \
    --min-depth-of-non-filtered-base 10 \
    --min-mapping-quality 20 \
    --min-base-quality 5 \
    --output-format CSV
}

export -f process_sample
export REF BAM_DIR VCF_DIR OUT_DIR

# Requires GNU parallel
seq 1747192 1747211 | parallel -j "${MAX_JOBS}" --halt soon,fail=1 process_sample {}

echo "[INFO] ASEReadCounter finished for all samples. Outputs in ${OUT_DIR}"
