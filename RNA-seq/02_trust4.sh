#!/bin/bash
set -euo pipefail

# Input:
#   paired-end FASTQ files: {sample}_trimmed_R1.fastq.gz and {sample}_trimmed_R2.fastq.gz
# Reference:
#   human_IMGT+C.fa
#
# Example:
#   N0085_trimmed_R1.fastq.gz
#   N0085_trimmed_R2.fastq.gz
#
# Usage:
#   bash 02_trust4.sh N0085

FASTQ_DIR="trimmed_fastq"
REF="human_IMGT+C.fa"
OUT_DIR="trust4_results"
THREADS=8

mkdir -p "${OUT_DIR}"

SAMPLE="$1"

R1="${FASTQ_DIR}/${SAMPLE}_trimmed_R1.fastq.gz"
R2="${FASTQ_DIR}/${SAMPLE}_trimmed_R2.fastq.gz"
OUT_PREFIX="${OUT_DIR}/${SAMPLE}_TRUST4"

if [[ ! -f "${R1}" || ! -f "${R2}" ]]; then
  echo "[ERROR] Missing FASTQ: ${R1} or ${R2}" >&2
  exit 1
fi

if [[ ! -f "${REF}" ]]; then
  echo "[ERROR] Missing reference file: ${REF}" >&2
  exit 1
fi

run-trust4 \
  -1 "${R1}" \
  -2 "${R2}" \
  -f "${REF}" \
  -o "${OUT_PREFIX}" \
  -t "${THREADS}"
