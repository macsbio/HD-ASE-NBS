#!/usr/bin/env bash
set -euo pipefail

# Config
IN_DIR="../fastq_files"          # where SRR*_1.fastq / SRR*_2.fastq live
OUT_DIR="trimmed-reads"          # output folder
THREADS=16

mkdir -p "${OUT_DIR}"

# Loop SRR1747192 .. SRR1747211 (inclusive)
for i in $(seq 1747192 1747211); do
  SAMPLE="SRR${i}"
  R1="${IN_DIR}/${SAMPLE}_1.fastq"
  R2="${IN_DIR}/${SAMPLE}_2.fastq"

  # Also support .fastq.gz transparently if .fastq not found
  if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    if [[ -f "${R1}.gz" && -f "${R2}.gz" ]]; then
      R1="${R1}.gz"
      R2="${R2}.gz"
    else
      echo "[WARN] Skipping ${SAMPLE}: input FASTQs not found."
      continue
    fi
  fi

  OUT_R1="${OUT_DIR}/${SAMPLE}_1_trimmed.fastq"
  OUT_R2="${OUT_DIR}/${SAMPLE}_2_trimmed.fastq"
  HTML="${OUT_DIR}/${SAMPLE}_fastp.html"
  JSON="${OUT_DIR}/${SAMPLE}_fastp.json"
  LOG="${OUT_DIR}/${SAMPLE}_fastp.log"

  # Skip if already done (both outputs exist)
  if [[ -f "$OUT_R1" && -f "$OUT_R2" ]]; then
    echo "[INFO] ${SAMPLE} already trimmed. Skipping."
    continue
  fi

  echo "[INFO] Processing ${SAMPLE} ..."
  fastp \
    -i "${R1}" -I "${R2}" \
    -o "${OUT_R1}" -O "${OUT_R2}" \
    --detect_adapter_for_pe \
    --cut_front --cut_tail --cut_right_mean_quality 20 \
    -q 20 -u 30 -n 5 -l 30 \
    --thread "${THREADS}" \
    --html "${HTML}" \
    --json "${JSON}" \
    2> "${LOG}"

  echo "[INFO] Done ${SAMPLE}"
done

echo "[INFO] All samples processed."

