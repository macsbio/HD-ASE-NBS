#!/usr/bin/env bash
set -euo pipefail

# Optional: activate conda env if needed
# source "$(conda info --base)/etc/profile.d/conda.sh"
# conda activate gatk_env

IN_DIR="../STAR-WASP-output"   # directory with *_WASP_Aligned.sortedByCoord.out.bam
OUT_DIR="../STAR-WASP-output"  # write RG BAMs in same dir
THREADS=16                      # threads for samtools index (optional)

# mkdir -p "${OUT_DIR}"

for i in $(seq 1747192 1747211); do
  SAMPLE="SRR${i}"
  IN_BAM="${IN_DIR}/${SAMPLE}_WASP_Aligned.sortedByCoord.out.bam"
  OUT_BAM="${OUT_DIR}/${SAMPLE}_WASP_RG.bam"

  if [[ ! -f "${IN_BAM}" ]]; then
    echo "[WARN] Skipping ${SAMPLE}: input BAM not found: ${IN_BAM}"
    continue
  fi

  echo "[INFO] AddOrReplaceReadGroups for ${SAMPLE}"
  gatk AddOrReplaceReadGroups \
    -I "${IN_BAM}" \
    -O "${OUT_BAM}" \
    -RGID "${SAMPLE}" \
    -RGLB "lib1" \
    -RGPL "ILLUMINA" \
    -RGPU "${SAMPLE}" \
    -RGSM "${SAMPLE}" \
    --CREATE_INDEX true

  # CREATE_INDEX=true already creates a .bai, but you can force re-index with samtools if you want:
  echo "[INFO] Indexing ${OUT_BAM}"
  samtools index -@ "${THREADS}" "${OUT_BAM}"

done

echo "[INFO] Finished AddOrReplaceReadGroups + indexing for SRR1747192–SRR1747211"
