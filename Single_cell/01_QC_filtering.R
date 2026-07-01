#!/usr/bin/env Rscript
# =============================================================================
# 01_QC_filtering.R
# Input  : Cell Ranger multi output (per_sample_outs)
# Output : sample_list_filtered_annotated.rds
# =============================================================================

library(Seurat)
library(dplyr)

set.seed(1234)

# =============================================================================
# Paths
# =============================================================================

sample_dirs <- c(
  Sample1  = "/path/to/SingleCell/results/cellranger_multi/Sample1_multi/outs/per_sample_outs/Sample1_multi/sample_filtered_feature_bc_matrix",
  Sample2  = "/path/to/SingleCell/results/cellranger_multi/Sample2_multi/outs/per_sample_outs/Sample2_multi/sample_filtered_feature_bc_matrix",
  Sample3 = "/path/to/SingleCell/results/cellranger_multi/Sample3_multi/outs/per_sample_outs/Sample3_multi/sample_filtered_feature_bc_matrix",
  Sample4 = "/path/to/SingleCell/results/cellranger_multi/Sample4_multi/outs/per_sample_outs/Sample4_multi/sample_filtered_feature_bc_matrix",
  Sample5 = "/path/to/SingleCell/results/cellranger_multi/Sample5_multi/outs/per_sample_outs/Sample5_multi/sample_filtered_feature_bc_matrix",
  Sample6 = "/path/to/SingleCell/results/cellranger_multi/Sample6_multi/outs/per_sample_outs/Sample6_multi/sample_filtered_feature_bc_matrix",
  Sample7 = "/path/to/SingleCell/results/cellranger_multi/Sample7_multi/outs/per_sample_outs/Sample7_multi/sample_filtered_feature_bc_matrix",
  Sample8 = "/path/to/SingleCell/results/cellranger_multi/Sample8_multi/outs/per_sample_outs/Sample8_multi/sample_filtered_feature_bc_matrix"
)

# Sample metadata
sample_rename_map <- c(
  Sample1  = "NT1_4",
  Sample2  = "NT1_5",
  Sample3 = "NT1_1",
  Sample4 = "NT1_2",
  Sample5 = "NT1_3",
  Sample6 = "Ctrl_1",
  Sample7 = "Ctrl_2",
  Sample8 = "Ctrl_3"
)

sample_batch_map <- c(
  NT1_4  = "NovaSeq6000",
  NT1_5  = "NovaSeq6000",
  NT1_1  = "NovaSeqX",
  NT1_2  = "NovaSeqX",
  NT1_3  = "NovaSeqX",
  Ctrl_1 = "NovaSeqX",
  Ctrl_2 = "NovaSeqX",
  Ctrl_3 = "NovaSeqX"
)

out_rds <- "./sample_list_filtered_annotated.rds"

# =============================================================================
# QC thresholds
# =============================================================================

MIN_FEATURES  <- 600
MAX_FEATURES  <- 5000
MAX_MT        <- 15    # percent

# =============================================================================
# Load and QC
# =============================================================================

cat("=== Loading samples ===\n")

sample_list <- lapply(names(sample_dirs), function(old_id) {

  cat("Reading:", old_id, "\n")

  mat <- Read10X(data.dir = sample_dirs[[old_id]])
  if (is.list(mat)) {
    mat <- mat[["Gene Expression"]]
  }

  obj <- CreateSeuratObject(
    counts      = mat,
    project     = old_id,
    min.cells   = 3,
    min.features = 200
  )

  # Mitochondrial gene percentage
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

  # Metadata
  new_id         <- sample_rename_map[[old_id]]
  obj$sample_id  <- new_id
  obj$condition  <- ifelse(grepl("^NT1", new_id), "NT1", "Ctrl")
  obj$batch      <- sample_batch_map[[new_id]]

  # QC filtering
  obj <- subset(
    obj,
    subset =
      nFeature_RNA > MIN_FEATURES &
      nFeature_RNA < MAX_FEATURES &
      percent.mt   < MAX_MT
  )

  cat("  Cells after QC:", ncol(obj), "\n")
  obj
})

names(sample_list) <- unname(sample_rename_map[names(sample_dirs)])

# =============================================================================
# Save
# =============================================================================

saveRDS(sample_list, out_rds)
cat("\nSaved:", out_rds, "\n")
cat("Cells per sample:\n")
print(sapply(sample_list, ncol))
