#!/usr/bin/env bash
set -e

GENOME_DIR="GrCh37_STAR"
OUT_DIR="../STAR-output-last-try"
THREADS=30

mkdir -p "${OUT_DIR}"

# Loop SRR1747192 .. SRR1747211
for i in $(seq 1747192 1747211); do
  SAMPLE="SRR${i}"
  R1="trimmed-reads/${SAMPLE}_1_trimmed.fastq"
  R2="trimmed-reads/${SAMPLE}_2_trimmed.fastq"

  STAR --runMode alignReads \
       --runThreadN "${THREADS}" \
       --genomeDir "${GENOME_DIR}" \
       --readFilesIn "${R1}" "${R2}" \
       --outFileNamePrefix "${OUT_DIR}/${SAMPLE}" \
       --outSAMtype BAM SortedByCoordinate \
       --twopassMode Basic
done

