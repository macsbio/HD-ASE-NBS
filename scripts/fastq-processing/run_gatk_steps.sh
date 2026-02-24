#!/usr/bin/env bash
set -euo pipefail

# init conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gatk_env

REF="Homo_sapiens.GRCh37.chr1-22_X.fa"
DICT="Homo_sapiens.GRCh37.chr1-22_X.dict"
OUT_DIR="../STAR-output-last-try"

mkdir -p "${OUT_DIR}"

##############################
# 0) PREPARE REFERENCE FILES #
##############################
if [[ ! -f "${DICT}" ]]; then
  gatk CreateSequenceDictionary \
    -R "${REF}" \
    -O "${DICT}"
fi

##############################
# Per-sample processing loop #
##############################
for i in $(seq 1747192 1747211); do
  SAMPLE="SRR${i}"

  IN_BAM="${OUT_DIR}/${SAMPLE}Aligned.sortedByCoord.out.bam"
  SPLIT_BAM="${OUT_DIR}/${SAMPLE}_split.bam"
  RG_BAM="${OUT_DIR}/${SAMPLE}_RG.bam"

  if [[ ! -f "${IN_BAM}" ]]; then
    echo "[WARN] Skipping ${SAMPLE}: input BAM not found: ${IN_BAM}"
    continue
  fi

  # temp paths on fast local storage
  tmp_in="/tmp/${SAMPLE}.bam"
  tmp_split="/tmp/${SAMPLE}_split.bam"
  tmp_rg="/tmp/${SAMPLE}_RG.bam"

  echo "[INFO] Copy to local tmp -> ${tmp_in}"
  cp "${IN_BAM}" "${tmp_in}"
  samtools index "${tmp_in}"

  echo "[INFO] SplitNCigarReads -> ${tmp_split}"
  gatk --java-options "-Xmx4G" SplitNCigarReads \
    -R "${REF}" \
    -I "${tmp_in}" \
    -O "${tmp_split}"

  echo "[INFO] AddOrReplaceReadGroups -> ${tmp_rg}"
  gatk AddOrReplaceReadGroups \
    -I "${tmp_split}" \
    -O "${tmp_rg}" \
    -RGID "${SAMPLE}" \
    -RGLB lib1 \
    -RGPL ILLUMINA \
    -RGPU unit1 \
    -RGSM "${SAMPLE}"

  echo "[INFO] samtools index -> ${tmp_rg}.bai"
  samtools index "${tmp_rg}"

  echo "[INFO] Copy back to PVC -> ${RG_BAM}"
  cp "${tmp_rg}" "${RG_BAM}"
  cp "${tmp_rg}.bai" "${RG_BAM}.bai"

  # optional cleanup
  rm -f "${tmp_in}" "${tmp_in}.bai" "${tmp_split}" "${tmp_rg}" "${tmp_rg}.bai"

  echo "[INFO] Done ${SAMPLE}"
done

echo "[INFO] All samples processed."


