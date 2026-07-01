#!/usr/bin/env Rscript
# =============================================================================
# 02_scDblFinder.R
# Input  : sample_list_filtered_annotated.rds (from 01_QC_filtering.R)
# Output : sample_list_singlets.rds
# Note   : scDblFinder is run per sample independently (recommended practice)
#          set.seed() is fixed for reproducibility; exact results may still
#          vary slightly across environments due to scDblFinder internals.
#          Analyses in the manuscript were run with:
#            scDblFinder v1.x.x, Seurat v5.x.x, R v4.x.x
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(dplyr)
})

set.seed(1234)

# =============================================================================
# Paths
# =============================================================================

in_rds  <- "./sample_list_filtered_annotated.rds"
out_rds <- "./sample_list_singlets.rds"

# =============================================================================
# Load
# =============================================================================

cat("=== Loading filtered sample list ===\n")
sample_list <- readRDS(in_rds)
cat("Samples:", paste(names(sample_list), collapse = ", "), "\n")

# =============================================================================
# Run scDblFinder per sample
# =============================================================================

cat("\n=== Running scDblFinder (per sample) ===\n")

singlet_list <- lapply(names(sample_list), function(sid) {

  cat("\n--- Sample:", sid, "---\n")
  obj <- sample_list[[sid]]

  # Preliminary clustering for scDblFinder (cluster-based mode)
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.3, verbose = FALSE)

  # Convert to SCE and run scDblFinder
  sce <- as.SingleCellExperiment(obj)
  sce <- scDblFinder(sce, clusters = "seurat_clusters")

  # Extract singlet barcodes
  singlet_barcodes <- colnames(sce)[sce$scDblFinder.class == "singlet"]
  cat("  Total cells:", ncol(obj), "\n")
  cat("  Singlets:", length(singlet_barcodes), "\n")
  cat("  Doublets removed:", ncol(obj) - length(singlet_barcodes), "\n")

  # Subset original Seurat object to singlets only
  obj <- subset(obj, cells = singlet_barcodes)

  # Keep RNA assay only before merge
  DefaultAssay(obj) <- "RNA"
  obj <- DietSeurat(obj, assays = "RNA", dimreducs = NULL, graphs = NULL)

  obj
})

names(singlet_list) <- names(sample_list)

# =============================================================================
# Merge all samples
# =============================================================================

cat("\n=== Merging samples ===\n")

merged_obj <- merge(
  x           = singlet_list[[1]],
  y           = singlet_list[-1],
  add.cell.ids = names(singlet_list)
)

cat("Total cells after merge:", ncol(merged_obj), "\n")
cat("Cells per sample:\n")
print(table(merged_obj$sample_id))

# =============================================================================
# Save
# =============================================================================

saveRDS(merged_obj, out_rds)
cat("\nSaved:", out_rds, "\n")
