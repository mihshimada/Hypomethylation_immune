#!/usr/bin/env Rscript
# =============================================================================
# 05_tcell_subset_clustering.R
#
# T cell subset re-clustering and removal of low-quality / artifact clusters.
#
# Input : obj_T_subset.rds   (output of 04_lineage_subset_reclustering.R)
#         T cells subsetted from the annotated PBMC object; Harmony-integrated;
#         Azimuth label-transfer already applied. Still contains an IEG-driven
#         stress-artifact cluster and an ENSG00000280441-driven artifact cluster.
#
# Output: obj_T_clean.rds
#         Final, artifact-free T cell object used for all downstream analyses:
#         annotation, proportion testing, and TCR repertoire analysis.
#
# Strategy:
#   Two artifact clusters are removed in two sequential rounds. After each
#   removal the object is re-integrated (Harmony) and re-clustered, because
#   removing a sizeable population changes the neighbourhood graph. The cluster
#   indices removed below are specific to THIS dataset (identified from the
#   marker tables); all random steps are seeded so the indices are stable on
#   re-run. The rationale for each removal is documented inline.
#
# Environment: local (R / Seurat v5)
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

set.seed(1234)

# ---- Paths ------------------------------------------------------------------
data_dir <- "subset_azimuth"            # same directory used by step 04
in_path  <- file.path(data_dir, "obj_T_subset.rds")
out_path <- file.path(data_dir, "obj_T_clean.rds")

# ---- Parameters -------------------------------------------------------------
harmony_var <- "batch"   # variable corrected for by Harmony
n_dims      <- 20        # Harmony dimensions used downstream
cluster_res <- 0.3       # FindClusters resolution
umap_seed   <- 1234

# ---- Helper: standard re-integration + clustering ---------------------------
# Re-runs PCA -> Harmony -> neighbours -> clusters -> UMAP on the SCT assay.
# Used after each artifact-cluster removal.
recluster <- function(obj) {
  DefaultAssay(obj) <- "SCT"
  obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
  obj <- RunHarmony(obj, group.by.vars = harmony_var,
                    reduction.use = "pca", dims.use = 1:n_dims,
                    verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:n_dims)
  obj <- FindClusters(obj, resolution = cluster_res)
  obj <- RunUMAP(obj, reduction = "harmony", dims = 1:n_dims,
                 seed.use = umap_seed)
  obj
}


# =============================================================================
# Load input
# =============================================================================
obj <- readRDS(in_path)
cat("Loaded:", ncol(obj), "T cells,",
    length(levels(obj$seurat_clusters)), "clusters\n")
print(table(obj$seurat_clusters))


# =============================================================================
# Round 1 - remove the IEG (immediate-early gene) stress-artifact cluster
# =============================================================================
# In obj_T_subset.rds this is cluster 10.
# Rationale: this cluster is defined by immediate-early genes (FOS, JUN, FOSB,
# CD69), i.e. an ex-vivo stress / handling response rather than a genuine T
# cell state, and was strongly skewed toward a single sample.
# (Identified from markers_T.csv.)
ieg_cluster <- "10"

obj <- subset(obj, seurat_clusters != ieg_cluster)
obj <- recluster(obj)
cat("\n[Round 1] Removed IEG cluster", ieg_cluster, "->",
    ncol(obj), "cells,",
    length(levels(obj$seurat_clusters)), "clusters\n")


# =============================================================================
# Round 2 - remove the ENSG00000280441 processing-artifact cluster
# =============================================================================
# After round-1 re-clustering this is cluster 11.
# Rationale: this cluster's dominant marker is the uncharacterized locus
# ENSG00000280441, a reported sample-processing artifact rather than a bona
# fide cell state (Chen et al., Brief Bioinform 2025, bbaf527; bioRxiv
# 2026.01.16.699571). This cluster was considered likely to be driven by this
# artifact-associated signal and was removed before annotation.
artifact_cluster <- "11"

obj <- subset(obj, seurat_clusters != artifact_cluster)
obj <- recluster(obj)
cat("[Round 2] Removed ENSG00000280441 artifact cluster", artifact_cluster, "->",
    ncol(obj), "cells,",
    length(levels(obj$seurat_clusters)), "clusters\n")


# =============================================================================
# Sanity-check plots
# =============================================================================
p1 <- DimPlot(obj, group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
  ggtitle("T cell clusters (res = 0.3)")
p2 <- DimPlot(obj, group.by = "predicted.celltype.l2", label = TRUE, repel = TRUE) +
  ggtitle("Azimuth level 2")
print(p1 + p2)

cat("\nFinal cluster sizes:\n")
print(table(obj$seurat_clusters))


# =============================================================================
# Save
# =============================================================================
saveRDS(obj, out_path)
cat("\nSaved:", out_path, "\n")

# =============================================================================
# Marker genes for the final T cell subsets (source of Supplementary Table S17)
# =============================================================================
# FindAllMarkers is run on the cluster-level identities of the cleaned object.
obj <- NormalizeData(obj, verbose = FALSE)
Idents(obj) <- "seurat_clusters"

markers_T_clean2 <- FindAllMarkers(
  obj,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25
)

top_markers_T_clean2 <- markers_T_clean2 %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
  ungroup()

write.csv(markers_T_clean2,
          file.path(data_dir, "markers_T_clean2.csv"),
          row.names = FALSE)
write.csv(top_markers_T_clean2,
          file.path(data_dir, "top_markers_T_clean2.csv"),
          row.names = FALSE)

cat("Done:", format(Sys.time()), "\n")
