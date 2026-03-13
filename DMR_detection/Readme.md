# Hypomethylation-associated DMR analysis

This directory contains example scripts used to identify differentially methylated regions (DMRs) associated with continuous hypomethylation levels.

## Overview

The same analytical framework was applied to multiple datasets, including NT1 CD4, NT1 CD8, and MS.  
As a representative example, the script provided here shows the analysis used for the NT1 CD4 dataset.

The input data for this analysis is an **appropriately filtered and normalized M value matrix** derived from Illumina HumanMethylationEPIC (or 450K) array data.  
DMRs were identified using the **bumphunter** algorithm implemented in the R package *bumphunter*.

The analysis relies on the R/Bioconductor packages **minfi**, **bumphunter**, and **IlluminaHumanMethylationEPICanno.ilm10b4.hg19**.

The design matrix included the hypomethylation variable of interest and several covariates.  
The exact set of covariates may vary depending on the sample set.

---

## Script

### 01_bumphunter_dmr_CD4.R

This script identifies DMRs associated with continuous hypomethylation values using **bumphunter**.

### Input

- appropriately filtered and normalized **M value matrix**
- phenotype table containing hypomethylation values and covariates

### Main analysis steps

1. Load the M value matrix
2. Create a `GenomicRatioSet` object using the **minfi** package
3. Construct a design matrix including hypomethylation and covariates
4. Identify DMRs using **bumphunter**

### Key parameters

- `B = 1000`
- `type = "M"`
- `smooth = FALSE`
- `maxGap = 1000`
- `pickCutoff = TRUE`
- `pickCutoffQ = 0.99`
- `nullMethod = "bootstrap"`

---

## Example usage

Rscript 01_bumphunter_dmr_CD4.R Mvalue_matrix.txt phenotype_table.txt output.txt

---

## Notes

This script represents a simplified version of the commands used in the study and is provided for transparency and reproducibility.  
Minor modifications may be required depending on the computing environment or dataset.







## Phenotype table

The phenotype table contains one row per sample and includes the hypomethylation index and covariates used in the design matrix.

Example format:
ID Group1 Group2 Group3 Group4 Group5
C01_4 0.781760178 26 0 1 1253.459947
C02_4 0.779833644 25 0 2 1242.368165
C03_4 0.785280242 40 0 1 1355.575394


Column description:

| Column | Description |
|-------|-------------|
| ID | sample ID |
| Group1 | hypomethylation index |
| Group2 | age |
| Group3 | disease status |
| Group4 | sex |
| Group5 | naive T cell index |
