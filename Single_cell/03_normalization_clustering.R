#!/usr/bin/env Rscript
# =============================================================================
# 03_normalization_clustering.R
# Input  : sample_list_singlets.rds (from 02_scDblFinder.R)
# Output : obj_pbmc_pre_azimuth.rds (pre-Azimuth), obj_pbmc_azimuth.rds (post-Azimuth)
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

# =============================================================================
# Paths
# =============================================================================
setwd("/path/to/SingleCell/results/output/rds_F")
in_rds  <- "./sample_list_singlets.rds"

# =============================================================================
# Parameters
# =============================================================================

NPCS        <- 50
PCS_USE     <- 1:20
RESOLUTION  <- 0.5

# =============================================================================
# Load
# =============================================================================

cat("=== Loading merged singlet object ===\n")
obj <- readRDS(in_rds)
cat("Total cells:", ncol(obj), "\n")

# =============================================================================
# SCTransform
# =============================================================================

cat("\n=== SCTransform ===\n")
obj <- SCTransform(
  obj,
  vars.to.regress  = "percent.mt",
  verbose          = TRUE,
  method           = "glmGamPoi"
)

# =============================================================================
# PCA
# =============================================================================

cat("\n=== PCA ===\n")
DefaultAssay(obj) <- "SCT"
obj <- RunPCA(obj, npcs = NPCS, verbose = FALSE)

# =============================================================================
# Harmony batch correction
# =============================================================================

cat("\n=== Harmony ===\n")
obj <- RunHarmony(
  object        = obj,
  group.by.vars = "batch",
  reduction.use = "pca",
  dims.use      = PCS_USE,
  verbose       = TRUE
)

# =============================================================================
# Clustering and UMAP
# =============================================================================

cat("\n=== Clustering and UMAP ===\n")
obj <- FindNeighbors(obj, reduction = "harmony", dims = PCS_USE, verbose = FALSE)
obj <- FindClusters(obj, resolution = RESOLUTION, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "harmony", dims = PCS_USE, seed.use = 1234)

cat("Cluster counts:\n")
print(table(obj$seurat_clusters))

# =============================================================================
# FindAllMarkers (for annotation reference)
# =============================================================================

cat("\n=== FindAllMarkers ===\n")
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
obj <- NormalizeData(obj, verbose = FALSE)

markers <- FindAllMarkers(
  obj,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  verbose         = FALSE
)

write.csv(markers, "./markers_pbmc_260528.csv", row.names = FALSE)
cat("Saved markers: ./markers_pbmc.csv\n")


# =============================================================================
# Save
# =============================================================================
obj_save <- DietSeurat(
  obj,
  graphs = NULL
)

saveRDS(obj_save, "./obj_pbmc_pre_azimuth.rds")
cat("Saved pre-Azimuth object\n")





# =============================================================================
# Cell type annotation
# Cluster-to-celltype mapping based on Azimuth
# Index corresponds to cluster number (0-based).
# =============================================================================

library(Azimuth)
obj <- RunAzimuth(obj, reference = "pbmcref")


# =============================================================================
# Save
# =============================================================================
saveRDS(obj, "./obj_pbmc_azimuth.rds")
cat("Saved Azimuth object\n")
