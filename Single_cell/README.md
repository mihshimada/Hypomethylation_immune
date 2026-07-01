# Single-cell RNA-seq and TCR analysis scripts

This directory contains example scripts used for single-cell RNA-seq (scRNA-seq) and paired single-cell TCR (scTCR) repertoire analysis of PBMCs from NT1 patients and controls.

## Overview

The analysis consisted of the following steps (run in numerical order):

1. **01_QC_filtering.R** — Per-sample quality control and filtering of the CellRanger gene-expression matrices (Seurat).
2. **02_scDblFinder.R** — Doublet detection and removal (scDblFinder), then merging into a single object.
3. **03_normalization_clustering.R** — SCTransform normalization, Harmony batch integration, clustering, and Azimuth reference-based annotation of the PBMC object.
4. **04_lineage_subset_reclustering.R** — Subsetting one lineage (T / B / Mono / NK) and re-clustering it. For the adaptive lineages, antigen-receptor V/(D)/J genes are removed from the variable features so that clusters reflect cell state rather than clonal receptor usage.
5. **05_tcell_subset_clustering.R** — Re-clustering of the T cell subset and removal of low-quality / artifact clusters.
6. **06_tcell_annotation_proportion.R** — Annotation of the T cell subsets, annotated UMAP, and comparison of subset abundance between conditions (NT1 vs Ctrl).
7. **07_tcell_tcr_repertoire.R** — Integration of paired alpha-beta TCR clonotypes (scRepertoire), per-subset TCR recovery, and clonotype overlap analysis between subsets.

## Requirements

R (>= 4.5.0), with Seurat (v5.5.0), harmony, scDblFinder, Azimuth, scRepertoire (v2.6.2), dplyr, tidyr, ggplot2, ggpubr.

These scripts are simplified versions of the commands used in the study and are provided for transparency and reproducibility.
