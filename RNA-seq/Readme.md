# RNA-seq analysis scripts

This directory contains example scripts used for RNA-seq preprocessing and T cell receptor (TCR) repertoire analysis.

## Overview

The analysis consisted of the following steps:

1. Read trimming and quality filtering using fastp (v0.23.4)
2. TCR repertoire reconstruction using TRUST4
3. Clonality analysis of TRA and TRB repertoires from TRUST4 AIRR output

These scripts are simplified versions of the commands used in the study and are provided for transparency and reproducibility.

---

## Input file naming

Paired-end FASTQ files are expected to follow this naming convention:

{sample}_1.fastq.gz
{sample}_2.fastq.gz

Example:

N0085_1.fastq.gz
N0085_2.fastq.gz

---

## Scripts

### 01_fastp.sh

Performs adapter trimming and quality filtering of paired-end RNA-seq reads using fastp (v0.23.4).

Main parameters:

- --detect_adapter_for_pe
- --trim_poly_g
- quality cutoff = 20
- maximum low-quality base percentage = 30
- minimum read length = 50 bp

Example:

bash 01_fastp.sh N0085

Output:

trimmed_fastq/N0085_trimmed_R1.fastq.gz
trimmed_fastq/N0085_trimmed_R2.fastq.gz
fastp_reports/N0085_fastp_report.html
fastp_reports/N0085_fastp_report.json

---

### 02_trust4.sh

Runs TRUST4 on trimmed RNA-seq reads.

Input:

trimmed_fastq/{sample}_trimmed_R1.fastq.gz
trimmed_fastq/{sample}_trimmed_R2.fastq.gz

Example:

bash 02_trust4.sh N0085

Output:

trust4_results/{sample}_TRUST4_airr.tsv

---

### 03_trust4_clonality.py

Summarizes TRUST4 AIRR output and calculates clonality metrics for TRA and TRB repertoires.

Clonotype definitions:

AA : CDR3 amino acid sequence (junction_aa)
NT : CDR3 nucleotide sequence (junction)

Metrics calculated:

- number of clones
- top1 frequency
- top10 cumulative frequency
- clonality (1 − normalized Shannon entropy)

Example:

python 03_trust4_clonality.py trust4_results

Output:

trust4_results/trust4_clonality_TRA_TRB.csv

---

## Notes

These scripts represent simplified versions of the commands used in the study and may require adaptation depending on the computing environment.
