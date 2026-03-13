#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(minfi)
  library(bumphunter)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript 01_bumphunter_dmr_CD4.R <Mvalue_matrix.txt> <phenotype_table.txt> <output.txt>")
}

mval_file <- args[1]
pheno_file <- args[2]
out_file <- args[3]

# Load appropriately filtered and normalized M value matrix
mval <- read.table(mval_file, check.names = FALSE)
colnames(mval) <- trimws(gsub('"', '', colnames(mval)))

m_matrix <- as.matrix(mval)
storage.mode(m_matrix) <- "double"

# Load phenotype table
pheno <- read.table(
  pheno_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  row.names = 1
)

# Example covariates used for the NT1 CD4 dataset
pheno$Group  <- as.numeric(pheno$Group1)
pheno$Group2 <- as.numeric(pheno$Group2)
pheno$Group3 <- factor(pheno$Group3)
pheno$Group4 <- factor(pheno$Group4)
pheno$Group5 <- as.numeric(pheno$Group5)

# Create GenomicRatioSet
grset <- makeGenomicRatioSetFromMatrix(
  mat = m_matrix,
  pData = pheno,
  array = "IlluminaHumanMethylationEPIC",
  annotation = "ilm10b4.hg19",
  what = "M"
)

# Design matrix
design_mat <- model.matrix(
  ~ Group + Group2 + Group3 + Group4 + Group5,
  data = as.data.frame(colData(grset))
)

# Run bumphunter
dmr <- bumphunter(
  object = grset,
  design = design_mat,
  coef = 2,
  B = 1000,
  type = "M",
  smooth = FALSE,
  maxGap = 1000,
  pickCutoff = TRUE,
  pickCutoffQ = 0.99,
  nullMethod = "bootstrap",
  verbose = TRUE
)

write.table(dmr$table, out_file, quote = FALSE, sep = "\t", row.names = FALSE)
