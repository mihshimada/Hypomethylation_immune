#!/usr/bin/env Rscript
# =============================================================================
# 06_tcell_annotation_proportion.R
#
# Annotate the cleaned T cell clusters, draw the annotated UMAP, and compare
# subset abundance between conditions (NT1 vs Ctrl).
#
# Abundance is expressed as a fraction of total PBMC (per sample), so that the
# values reflect each subset's representation in the whole blood compartment.
# Group comparison uses an unpaired Wilcoxon test; raw p-values are reported
# (exploratory, no multiple-testing correction).
#
# Input : obj_T_clean.rds        (output of 05_tcell_subset_clustering.R)
#         obj_pbmc_azimuth.rds    (for total PBMC counts = denominator)
# Output: obj_T_annotated.rds, annotated UMAP, proportion boxplots, stacked bars
#
# Environment: local (R / Seurat v5)
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
})

set.seed(1234)

# ---- Paths ------------------------------------------------------------------
data_dir <- "subset_azimuth"
in_path   <- file.path(data_dir, "obj_T_clean.rds")
pbmc_path <- "obj_pbmc_azimuth.rds"
out_path  <- file.path(data_dir, "obj_T_annotated.rds")


# =============================================================================
# Load and annotate
# =============================================================================
obj <- readRDS(in_path)
Idents(obj) <- "seurat_clusters"

# Cluster 2:  CD8 TEM enriched for TOX/ZEB2.
# Cluster 10: KLRF1+ CD8 TEM; T-cell identity supported by TCR repertoire (08).
# Clusters 7, 11: gamma-delta T cells (TRDC+).
annotation_map <- c(
  "0"  = "CD4 TCM_1",
  "1"  = "CD8 TEM_1",
  "2"  = "CD8 TEM_3 (TOX+/ZEB2+)",
  "3"  = "CD4 Naive",
  "4"  = "CD8 TEM_2",
  "5"  = "CD4 TCM_2",
  "6"  = "CD8 Naive",
  "7"  = "gdT (Vd2+)",
  "8"  = "MAIT",
  "9"  = "Treg",
  "10" = "CD8 TEM_4 (KLRF1+)",
  "11" = "gdT (Vd1+)"
)

lab <- unname(annotation_map[as.character(obj$seurat_clusters)])
obj$annotation <- factor(lab, levels = unname(annotation_map))
obj$sample_label <- obj$sample_id  

print(table(obj$annotation))

saveRDS(obj, out_path)
cat("Saved:", out_path, "\n")


# =============================================================================
# Annotated UMAP
# =============================================================================
p_umap <- DimPlot(obj, group.by = "annotation", label = TRUE, repel = TRUE) +
  ggtitle("T cell subsets")
print(p_umap)


# =============================================================================
# Subset abundance as a fraction of total PBMC
# =============================================================================
# Denominator: total cells per sample in the PBMC object.
obj_pbmc <- readRDS(pbmc_path)
total_pbmc <- obj_pbmc@meta.data %>%
  group_by(orig.ident) %>%
  summarise(total_pbmc = n(), .groups = "drop")
rm(obj_pbmc); gc()

# Cells per sample x subset, completing absent combinations with zero.
prop_df <- obj@meta.data %>%
  group_by(orig.ident, condition, annotation) %>%
  summarise(n = n(), .groups = "drop") %>%
  complete(nesting(orig.ident, condition), annotation,
           fill = list(n = 0)) %>%
  left_join(total_pbmc, by = "orig.ident") %>%
  mutate(prop = n / total_pbmc)

prop_df$condition <- factor(prop_df$condition, levels = c("Ctrl", "NT1"))


# =============================================================================
# Boxplots + Wilcoxon (raw p-values, per subset)
# =============================================================================
p_box <- ggplot(prop_df, aes(x = condition, y = prop, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.8) +
  facet_wrap(~ annotation, scales = "free_y", ncol = 4) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 3) +
  labs(x = NULL, y = "Proportion of total PBMC") +
  theme_bw() +
  theme(legend.position = "none")
print(p_box)
