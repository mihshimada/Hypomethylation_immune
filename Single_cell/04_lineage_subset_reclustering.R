#!/usr/bin/env Rscript
# =============================================================================
# 04_lineage_subset_reclustering.R
#
# Lineage-level re-clustering of an annotated PBMC object.
#
# Given the Azimuth-annotated PBMC object, this script subsets one lineage
# (T / B / Mono / NK), re-runs SCTransform + Harmony + clustering on that
# lineage alone, and writes the re-clustered object plus marker tables. For
# the adaptive lineages, the antigen-receptor V/(D)/J genes are removed from
# the variable features so that clusters reflect cell state rather than clonal
# receptor usage.
#
# Usage:
#   Rscript 04_lineage_subset_reclustering.R <T|B|Mono|NK>
#
# Input : obj_pbmc_azimuth.rds   (output of 03_normalization_clustering.R)
# Output: obj_<lineage>_subset.rds
#         markers_<lineage>.csv, top_markers_<lineage>.csv
#
# Note: T-cell downstream artifact removal and annotation are handled in
#       05_tcell_subset_clustering.R and 06_tcell_annotation.R.
#
# Environment: HPC / R (Seurat v5, glmGamPoi)
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(future)
})

plan("sequential")
options(future.globals.maxSize = 50 * 1024^3)

set.seed(1234)

# ---- Command-line argument --------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript 04_lineage_subset_reclustering.R <T|B|Mono|NK>")
}
subset_name <- args[1]

# =============================================================================
# Paths
# =============================================================================
rds_dir      <- "."                                  # adjust to your layout
input_rds    <- file.path(rds_dir, "obj_pbmc_azimuth.rds")
out_base_dir <- file.path(rds_dir, "subset_azimuth")
dir.create(out_base_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Parameters
# =============================================================================
subset_npcs    <- 50
subset_pcs_use <- 1:20

# Clustering resolution per lineage
resolution_list <- list(T = 0.3, B = 0.2, Mono = 0.2, NK = 0.2)

# Azimuth predicted.celltype.l2 labels that define each lineage
subset_definitions <- list(
  T = c(
    "CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 CTL", "CD4 Proliferating",
    "CD8 Naive", "CD8 TCM", "CD8 TEM", "CD8 Proliferating",
    "Treg", "dnT", "gdT", "MAIT"
  ),
  B    = c("B naive", "B intermediate", "B memory", "Plasmablast"),
  Mono = c("CD14 Mono", "CD16 Mono", "cDC1", "cDC2", "ASDC", "pDC"),
  NK   = c("NK", "NK Proliferating", "NK_CD56bright")
)

# Antigen-receptor genes removed from variable features
#   T: alpha/beta TCR V and J genes (constant-region C genes are kept)
#   B: immunoglobulin heavy/kappa/lambda V genes
receptor_exclude_pattern <- list(
  T = "^(TRAV|TRAJ|TRBV|TRBJ)",
  B = "^(IGHV|IGKV|IGLV)"
)

if (!subset_name %in% names(subset_definitions)) {
  stop("subset_name must be one of: T, B, Mono, NK")
}

# =============================================================================
# Helper functions
# =============================================================================
exclude_receptor_genes <- function(obj, pattern, label = "receptor") {
  to_drop <- grep(pattern, rownames(obj), value = TRUE)
  cat("\n", label, " genes removed from variable features: ",
      length(to_drop), "\n", sep = "")
  VariableFeatures(obj) <- setdiff(VariableFeatures(obj), to_drop)
  cat("Variable features remaining: ", length(VariableFeatures(obj)), "\n")
  obj
}

run_marker_analysis <- function(obj, subset_name, out_dir) {
  DefaultAssay(obj) <- "RNA"
  obj <- JoinLayers(obj)
  obj <- NormalizeData(obj, verbose = FALSE)

  markers <- FindAllMarkers(obj, only.pos = TRUE,
                            min.pct = 0.25, logfc.threshold = 0.25)
  if (nrow(markers) == 0) { cat("No markers found.\n"); return(invisible(NULL)) }

  top_markers <- markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
    ungroup()

  write.csv(markers,
            file.path(out_dir, paste0("markers_", subset_name, ".csv")),
            row.names = FALSE)
  write.csv(top_markers,
            file.path(out_dir, paste0("top_markers_", subset_name, ".csv")),
            row.names = FALSE)
  cat("Saved marker tables for", subset_name, "\n")
}

# =============================================================================
# Load and subset
# =============================================================================
cat("========================================\n")
cat("Lineage subset re-clustering:", subset_name, "\n")
cat("Time:", format(Sys.time()), "\n")
cat("========================================\n\n")

if (!file.exists(input_rds)) {
  stop("Azimuth-annotated PBMC object not found: ", input_rds)
}

obj <- readRDS(input_rds)
if (!"predicted.celltype.l2" %in% colnames(obj@meta.data)) {
  stop("predicted.celltype.l2 not found. Run Azimuth (step 03) first.")
}

labels_use     <- subset_definitions[[subset_name]]
resolution_use <- resolution_list[[subset_name]]

sub_obj <- subset(obj, subset = predicted.celltype.l2 %in% labels_use)
rm(obj); gc()

cat("Cells in", subset_name, "subset:", ncol(sub_obj), "\n")
print(table(sub_obj$predicted.celltype.l2))
if (ncol(sub_obj) < 100) stop("Too few cells in subset: ", subset_name)

# Slim down to the RNA assay before recomputing SCT on the subset
DefaultAssay(sub_obj) <- "RNA"
sub_obj <- DietSeurat(sub_obj, assays = "RNA", dimreducs = NULL, graphs = NULL)
gc()

# =============================================================================
# SCTransform on the lineage subset
# =============================================================================
cat("\nRunning SCTransform...\n")
sub_obj <- SCTransform(sub_obj, vars.to.regress = "percent.mt",
                       method = "glmGamPoi", conserve.memory = TRUE,
                       verbose = TRUE)
DefaultAssay(sub_obj) <- "SCT"

# Remove antigen-receptor V/(D)/J genes from variable features (T and B only)
if (subset_name %in% names(receptor_exclude_pattern)) {
  sub_obj <- exclude_receptor_genes(
    sub_obj,
    pattern = receptor_exclude_pattern[[subset_name]],
    label   = subset_name
  )
}

# =============================================================================
# PCA -> Harmony -> clustering -> UMAP
# =============================================================================
cat("\nRunning PCA / Harmony / clustering / UMAP...\n")
sub_obj <- RunPCA(sub_obj, assay = "SCT", npcs = subset_npcs, verbose = FALSE)
sub_obj <- RunHarmony(sub_obj, group.by.vars = "batch",
                      reduction.use = "pca", dims.use = subset_pcs_use,
                      verbose = TRUE)
sub_obj <- FindNeighbors(sub_obj, reduction = "harmony", dims = subset_pcs_use)
sub_obj <- FindClusters(sub_obj, resolution = resolution_use)
sub_obj <- RunUMAP(sub_obj, reduction = "harmony", dims = subset_pcs_use,
                   seed.use = 1234)

cat("\nCluster counts:\n")
print(table(sub_obj$seurat_clusters))

# =============================================================================
# Save and run marker analysis
# =============================================================================
out_rds <- file.path(out_base_dir, paste0("obj_", subset_name, "_subset.rds"))
saveRDS(sub_obj, out_rds)
cat("\nSaved:", out_rds, "\n")

Idents(sub_obj) <- "seurat_clusters"
run_marker_analysis(sub_obj, subset_name, out_base_dir)

cat("\nDone:", subset_name, "\n")
cat("Finished at:", format(Sys.time()), "\n")
