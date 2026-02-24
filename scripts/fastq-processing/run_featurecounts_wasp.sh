#!/usr/bin/env bash
set -euo pipefail

# Optional: activate conda env if needed
# source "$(conda info --base)/etc/profile.d/conda.sh"
# conda activate gatk_env   # or env where featureCounts is installed

# Paths / parameters
BAM_DIR="../STAR-WASP-output"          # where *_WASP_RG.bam live
GTF="filtered_chr1-22_X.gtf"           # annotation used for STAR
OUT_COUNTS="WASP_featureCounts.txt"    # output matrix
THREADS=16                             # threads for featureCounts
STRAND=0                               # 0=unstranded, 1=fr-stranded, 2=rf-stranded

# Collect BAMs
BAM_LIST=()
for i in $(seq 1747192 1747211); do
  SAMPLE="SRR${i}"
  BAM="${BAM_DIR}/${SAMPLE}_WASP_RG.bam"
  if [[ -f "${BAM}" ]]; then
    echo "[INFO] Adding BAM: ${BAM}"
    BAM_LIST+=("${BAM}")
  else
    echo "[WARN] Missing BAM for ${SAMPLE}: ${BAM}"
  fi
done

if [[ ${#BAM_LIST[@]} -eq 0 ]]; then
  echo "[ERROR] No BAM files found, aborting."
  exit 1
fi

echo "[INFO] Running featureCounts on ${#BAM_LIST[@]} BAM files"

featureCounts \
  -T "${THREADS}" \
  -a "${GTF}" \
  -o "${OUT_COUNTS}" \
  -p -B -C \
  -s "${STRAND}" \
  "${BAM_LIST[@]}"

echo "[INFO] featureCounts finished. Output: ${OUT_COUNTS}"
