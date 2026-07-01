#!/usr/bin/env Rscript
# =============================================================================
# 07_tcell_tcr_repertoire.R
#
# TCR repertoire analysis of the annotated T cell subsets.
#
#   1. Integrate paired alpha-beta TCR clonotypes (10x VDJ) into the annotated
#      T cell object (scRepertoire).
#   2. Per-cluster TCR recovery rate.
#   3. Clonotype overlap between subsets (cell-level, row-normalized), split by
#        condition (NT1 vs Ctrl). cell-level = fraction of a subset's cells whose
#        clonotype is also found in another subset. Computed within sample
#        (clonotypes are not shared across donors). Splitting by condition lets a
#        disease-specific sharing pattern be told apart from the baseline.
#
# Input : obj_T_annotated.rds   (output of 06_tcell_annotation_proportion.R)
#         vdj_t/<sample>/filtered_contig_annotations.csv  (10x VDJ output)
# Output: figures (no object is re-saved)
#
# Environment: local (R / Seurat v5, scRepertoire)
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(scRepertoire)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

set.seed(1234)

# ---- Paths ------------------------------------------------------------------
data_dir <- "subset_azimuth"
in_path  <- file.path(data_dir, "obj_T_annotated.rds")
vdj_dir  <- "vdj_t"   # contains one folder per sample

obj <- readRDS(in_path)
Idents(obj) <- "annotation"


# =============================================================================
# 1. Integrate TCR (paired alpha-beta clonotypes)
# =============================================================================
# Map original sample id -> renamed label. combineTCR prefixes the renamed
# label to each barcode so it matches the Seurat barcodes (e.g. "NT1_4_<bc>").
sample_map <- c(
  Sample1  = "NT1_4",
  Sample2  = "NT1_5",
  Sample3 = "NT1_1",
  Sample4 = "NT1_2",
  Sample5 = "NT1_3",
  Sample6 = "Ctrl_1",
  Sample7 = "Ctrl_2",
  Sample8 = "Ctrl_3"
)

vdj_paths <- file.path(vdj_dir, names(sample_map),
                       "filtered_contig_annotations.csv")
names(vdj_paths) <- sample_map
stopifnot(all(file.exists(vdj_paths)))

contig_list  <- lapply(vdj_paths, read.csv)
combined_tcr <- combineTCR(contig_list, samples = names(vdj_paths))

obj <- combineExpression(
  combined_tcr, obj,
  cloneCall  = "aa",
  group.by   = "sample",   # clone frequency computed within sample
  proportion = TRUE
)

cat("TCR recovery (overall):",
    round(mean(!is.na(obj$CTaa)), 3), "\n")


# =============================================================================
# 2. Per-subset TCR recovery rate
# =============================================================================
tcr_rate <- obj@meta.data %>%
  group_by(annotation) %>%
  summarise(n_cells = n(),
            n_tcr   = sum(!is.na(CTaa)),
            rate    = round(mean(!is.na(CTaa)), 3),
            .groups = "drop")
print(tcr_rate)



# =============================================================================
# 3. Clonotype overlap between subsets (cell-level), split by condition
#    cell-level = fraction of a subset's cells whose clonotype is also found in
#    the column subset. Computed within sample (clonotypes are not shared across
#    donors) and row-normalized. NT1 and Ctrl are computed separately so that a
#    condition-specific sharing pattern (e.g. expanded clones spreading across
#    effector subsets) can be distinguished from the baseline.
# =============================================================================

# subset display order: CD4 lineage -> CD8 lineage -> others
subset_order <- c(
  "CD4 Naive", "CD4 TCM_1", "CD4 TCM_2", "Treg",                  # CD4 lineage
  "CD8 Naive", "CD8 TEM_1", "CD8 TEM_2",
  "CD8 TEM_3 (TOX+/ZEB2+)", "CD8 TEM_4 (KLRF1+)",                 # CD8 lineage
  "MAIT", "gdT (Vd2+)", "gdT (Vd1+)"                              # others
)

# cell-level overlap for one condition (within sample, row-normalized)
compute_overlap_cell <- function(meta_tcr) {
  subsets <- intersect(subset_order, unique(meta_tcr$annotation))

  clone_in   <- meta_tcr %>% distinct(sample_label, CTaa, annotation)
  clono_sets <- lapply(subsets, function(s)
    clone_in %>% dplyr::filter(annotation == s) %>% distinct(sample_label, CTaa))
  names(clono_sets) <- subsets

  mat <- matrix(NA_real_, length(subsets), length(subsets),
                dimnames = list(subsets, subsets))
  for (i in seq_along(subsets)) {
    cells_i <- meta_tcr %>% dplyr::filter(annotation == subsets[i])
    ni <- nrow(cells_i)
    if (ni == 0) next
    for (j in seq_along(subsets)) {
      shared <- cells_i %>%
        semi_join(clono_sets[[j]], by = c("sample_label", "CTaa")) %>% nrow()
      mat[i, j] <- shared / ni
    }
  }
  mat
}

mat_to_long <- function(mat, group_label) {
  long <- as.data.frame(as.table(mat))
  colnames(long) <- c("row_subset", "col_subset", "frac")
  long$group <- group_label
  long
}

meta_tcr_all <- obj@meta.data %>%
  dplyr::filter(!is.na(CTaa)) %>%
  dplyr::select(CTaa, annotation, sample_label, condition)

mat_ctrl <- compute_overlap_cell(meta_tcr_all %>% dplyr::filter(condition == "Ctrl"))
mat_nt1  <- compute_overlap_cell(meta_tcr_all %>% dplyr::filter(condition == "NT1"))

long_all <- bind_rows(
  mat_to_long(mat_ctrl, "Ctrl"),
  mat_to_long(mat_nt1,  "NT1")
)
long_all$row_subset <- factor(long_all$row_subset, levels = rev(subset_order))
long_all$col_subset <- factor(long_all$col_subset, levels = subset_order)
long_all$group      <- factor(long_all$group, levels = c("Ctrl", "NT1"))

p_overlap <- ggplot(long_all, aes(x = col_subset, y = row_subset, fill = frac)) +
  geom_tile(color = "grey92") +
  geom_text(aes(label = ifelse(!is.na(frac) & frac > 0.05, sprintf("%.2f", frac), "")),
            size = 2.2) +
  facet_wrap(~ group) +
  scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1),
                      na.value = "grey95",
                      name = "Fraction\nshared\n(row-normalized)") +
  labs(x = NULL, y = NULL,
       title = "TCR clonotype overlap between T cell subsets (cell-level)",
       subtitle = "within sample, row-normalized; split by condition") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

dir.create("figures_07", showWarnings = FALSE)
ggsave(file.path("figures_07", "tcr_overlap_celllevel_by_condition.pdf"),
       p_overlap, width = 14, height = 7)
print(p_overlap)

# Note: Ctrl has only 3 donors, so sharing estimates for small subsets are
# noisier than NT1 (5 donors). Read low Ctrl values with that caveat.


# =============================================================================
# 3b. Per-sample clonotype sharing for selected CD8-effector pairs
#     The heatmaps above pool all donors of a condition into a single value,
#     which hides between-donor variability. Here the same cell-level sharing is
#     computed per donor, so that a condition difference can be checked for
#     consistency across donors (NT1 n=5 vs Ctrl n=3) rather than being driven
#     by a few samples.
# =============================================================================

# per-sample, row->col cell-level sharing (fraction of row-subset cells in this
# sample whose clonotype is also present in the column subset of the same sample)
overlap_per_sample <- function(meta_tcr) {
  res <- list()
  for (smp in unique(meta_tcr$sample_label)) {
    m    <- meta_tcr %>% dplyr::filter(sample_label == smp)
    subs <- unique(m$annotation)
    cset <- lapply(subs, function(s) unique(m$CTaa[m$annotation == s]))
    names(cset) <- subs
    for (a in subs) for (b in subs) {
      if (a == b) next
      cells_a <- m$CTaa[m$annotation == a]
      res[[length(res) + 1]] <- data.frame(
        sample_label = smp, row_subset = a, col_subset = b,
        frac = mean(cells_a %in% cset[[b]]))
    }
  }
  bind_rows(res)
}

ov_smp <- overlap_per_sample(meta_tcr_all) %>%
  left_join(meta_tcr_all %>% distinct(sample_label, condition), by = "sample_label")

# CD8 effector pairs that stood out in the pooled heatmap (row -> col)
focus_pairs <- tribble(
  ~row_subset,              ~col_subset,
  "CD8 TEM_1",              "CD8 TEM_3 (TOX+/ZEB2+)",
  "CD8 TEM_1",              "CD8 TEM_4 (KLRF1+)",
  "CD8 TEM_3 (TOX+/ZEB2+)", "CD8 TEM_4 (KLRF1+)",
  "CD8 TEM_3 (TOX+/ZEB2+)", "CD8 TEM_1",
  "CD8 TEM_4 (KLRF1+)",     "CD8 TEM_1",
  "CD8 TEM_2",              "CD8 TEM_3 (TOX+/ZEB2+)"
)

sharing_df <- ov_smp %>%
  inner_join(focus_pairs, by = c("row_subset", "col_subset")) %>%
  mutate(pair = paste0(gsub(" \\(.*", "", row_subset), " -> ",
                       gsub(" \\(.*", "", col_subset)))
sharing_df$condition <- factor(sharing_df$condition, levels = c("Ctrl", "NT1"))

p_share <- ggplot(sharing_df, aes(x = condition, y = frac, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  facet_wrap(~ pair, scales = "free_y") +
  scale_fill_manual(values = c("Ctrl" = "#A8C8A0", "NT1" = "#2E7D32")) +
  labs(x = NULL,
       y = "Fraction of row-subset cells sharing clonotype with col-subset",
       title = "Per-sample clonotype sharing (cell-level), CD8 effector pairs",
       subtitle = "each point = one donor; within sample") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

ggsave(file.path("figures_07", "tcr_sharing_persample_CD8effector.pdf"),
       p_share, width = 9, height = 6)
print(p_share)
