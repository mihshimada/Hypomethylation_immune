#!/bin/bash
set -euo pipefail

# Input:
#   raw FASTQ files: {sample}_1.fastq.gz and {sample}_2.fastq.gz
#
# Example:
#   N0085_1.fastq.gz
#   N0085_2.fastq.gz
#
# Usage:
#   bash 01_fastp.sh N0085

RAW_DIR="raw_fastq"
OUT_DIR="trimmed_fastq"
REPORT_DIR="fastp_reports"

mkdir -p "${OUT_DIR}" "${REPORT_DIR}"

SAMPLE="$1"

R1="${RAW_DIR}/${SAMPLE}_1.fastq.gz"
R2="${RAW_DIR}/${SAMPLE}_2.fastq.gz"

fastp \
  -i "${R1}" \
  -I "${R2}" \
  -o "${OUT_DIR}/${SAMPLE}_trimmed_R1.fastq.gz" \
  -O "${OUT_DIR}/${SAMPLE}_trimmed_R2.fastq.gz" \
  --detect_adapter_for_pe \
  --trim_poly_g \
  -q 20 \
  -u 30 \
  -l 50 \
  -w 8 \
  --html "${REPORT_DIR}/${SAMPLE}_fastp_report.html" \
  --json "${REPORT_DIR}/${SAMPLE}_fastp_report.json"
