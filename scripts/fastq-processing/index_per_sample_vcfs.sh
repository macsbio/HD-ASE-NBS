#!/usr/bin/env bash
set -euo pipefail

# Directory containing per-sample VCFs
# Adjust if needed
VCF_DIR="genotype_sample_files"

cd "${VCF_DIR}"

for i in $(seq 1747192 1747211); do
  VCF="SRR${i}.biallelic.snps.vcf.gz"

  if [[ ! -f "${VCF}" ]]; then
    echo "[WARN] Missing ${VCF}, skipping"
    continue
  fi

  echo "[INFO] Indexing ${VCF}"
  gatk IndexFeatureFile -I "${VCF}"
done

echo "[INFO] Finished indexing all per-sample VCFs in ${VCF_DIR}"
