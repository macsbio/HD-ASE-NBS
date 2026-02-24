#!/usr/bin/env bash
set -euo pipefail

IN_DIR="genotype_sample_files"       # folder with *.raw.vcf.gz
OUT_DIR="genotype_sample_files"      # can be same as IN_DIR
MAX_JOBS=4                           # how many samples in parallel

mkdir -p "${OUT_DIR}"

filter_sample() {
  local idx="$1"
  local SAMPLE="SRR${idx}"
  local IN_VCF="${IN_DIR}/${SAMPLE}.raw.vcf.gz"
  local OUT_VCF="${OUT_DIR}/${SAMPLE}.biallelic.snps.vcf.gz"

  if [[ ! -f "${IN_VCF}" ]]; then
    echo "[WARN] Skipping ${SAMPLE}: VCF not found: ${IN_VCF}"
    return 0
  fi

  echo "[INFO] Filtering to biallelic SNPs: ${SAMPLE}"
  bcftools view -m2 -M2 -v snps -Oz -o "${OUT_VCF}" "${IN_VCF}"
  tabix -f "${OUT_VCF}"
}

export -f filter_sample
export IN_DIR OUT_DIR

# run from SRR1747192 to SRR1747211
seq 1747192 1747211 | parallel -j "${MAX_JOBS}" --halt soon,fail=1 filter_sample {}

echo "[INFO] Biallelic SNP filtering finished for all samples. Outputs in ${OUT_DIR}"
