#!/usr/bin/env bash
set -euo pipefail

GENOME_DIR="GrCh37_STAR"
READ_DIR="trimmed-reads"
VCF_DIR="genotype_sample_files"
OUT_DIR="../STAR-WASP-output"

THREADS=20      # threads per STAR run (increased to 16)
MAX_JOBS=1      # run samples sequentially (no parallel execution)

mkdir -p "${OUT_DIR}"

align_sample() {
  local idx="$1"
  local SAMPLE="SRR${idx}"

  local R1="${READ_DIR}/${SAMPLE}_1_trimmed.fastq"
  local R2="${READ_DIR}/${SAMPLE}_2_trimmed.fastq"

  # VCF paths (compressed and uncompressed)
  local VCF_GZ="${VCF_DIR}/${SAMPLE}.biallelic.snps.vcf.gz"
  local VCF="${VCF_DIR}/${SAMPLE}.biallelic.snps.vcf"

  local PREFIX="${OUT_DIR}/${SAMPLE}_WASP_"

  # Check FASTQs
  if [[ ! -f "${R1}" || ! -f "${R2}" ]]; then
    echo "[WARN] Skipping ${SAMPLE}: missing FASTQ(s): ${R1} or ${R2}"
    return 0
  fi

  # Ensure we have an uncompressed VCF for STAR
  if [[ ! -f "${VCF}" ]]; then
    if [[ -f "${VCF_GZ}" ]]; then
      echo "[INFO] Decompressing ${VCF_GZ} -> ${VCF}"
      gunzip -c "${VCF_GZ}" > "${VCF}"
    fi
  fi

  # Final VCF existence check
  if [[ ! -f "${VCF}" ]]; then
    echo "[WARN] Skipping ${SAMPLE}: VCF not found (neither ${VCF} nor ${VCF_GZ})"
    return 0
  fi

  echo "[INFO] STAR WASP alignment for ${SAMPLE}"

  STAR \
    --runMode alignReads \
    --genomeDir "${GENOME_DIR}" \
    --outSAMtype BAM SortedByCoordinate \
    --readFilesIn "${R1}" "${R2}" \
    --outFileNamePrefix "${PREFIX}" \
    --twopassMode Basic \
    --runThreadN "${THREADS}" \
    --waspOutputMode SAMtag \
    --varVCFfile "${VCF}"
}

export -f align_sample
export GENOME_DIR READ_DIR VCF_DIR OUT_DIR THREADS

# Run the alignments sequentially (no parallel execution)
for idx in $(seq 1747192 1747211); do
  align_sample "${idx}"
done

echo "[INFO] STAR WASP alignment finished for all samples. Outputs in ${OUT_DIR}"
