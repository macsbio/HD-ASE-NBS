#!/usr/bin/env bash
set -euo pipefail

IN_DIR="trimmed-reads"     # where SRR*_1/2*.fastq(.gz) live
OUT_DIR="qc-fastqc"        # FastQC output directory
THREADS=16                  # per-FastQC threads

mkdir -p "${OUT_DIR}"

for i in $(seq 1747192 1747211); do
  SAMPLE="SRR${i}"

  # Match either .fastq or .fastq.gz
  R1=$(ls -1 "${IN_DIR}/${SAMPLE}_1"*.fastq* 2>/dev/null | head -n1 || true)
  R2=$(ls -1 "${IN_DIR}/${SAMPLE}_2"*.fastq* 2>/dev/null | head -n1 || true)

  if [[ -z "${R1}" || -z "${R2}" ]]; then
    echo "[WARN] Skipping ${SAMPLE}: files not found in ${IN_DIR}"
    continue
  fi

  echo "[INFO] FastQC ${SAMPLE}"
  fastqc --threads "${THREADS}" --outdir "${OUT_DIR}" "${R1}" "${R2}"
done

echo "[INFO] FastQC finished. Reports in ${OUT_DIR}"

