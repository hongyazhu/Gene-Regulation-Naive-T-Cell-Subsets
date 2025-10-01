##############################################
# RNA-seq Preprocessing and Batch Correction #
# ------------------------------------------ #
# 1. Import raw counts and metadata
# 2. Apply DESeq2 variance stabilizing transformation (VST)
# 3. Remove batch effects using ComBat (sva package)
# ------------------------------------------ #
# Author: Hongya Zhu
##############################################

# Load required libraries
library(DESeq2)
library(limma)
library(sva)

#------------------------------------------------
# Step 1: Load data
#------------------------------------------------
counts_all <- read.table("counts_raw.txt", header = TRUE)
sample_info <- read.table("sample_info.txt", header = TRUE, sep = "\t")

# Ensure factors are set correctly
factor_cols <- c("Source", "Age", "CellType", "Veteran", "Clean")
sample_info[factor_cols] <- lapply(sample_info[factor_cols], factor)

#------------------------------------------------
# Step 2: Create DESeq2 object & apply VST
#------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts_all,
  colData   = sample_info,
  design    = ~ Source + Age + CellType + Veteran + Clean
)

# Keep genes with at least 1 count across samples
dds <- dds[rowSums(counts(dds) > 0) >= 1, ]

# Variance stabilizing transformation
vsd <- vst(dds, blind = TRUE)

# Save VST-normalized counts
write.table(assay(vsd), 
            file = "counts_vst_blindT.txt", 
            quote = FALSE, sep = "\t")

#------------------------------------------------
# Step 3: Batch correction with ComBat
#------------------------------------------------
batch <- sample_info$Source
design <- model.matrix(~ Age + CellType + Veteran + Clean, data = sample_info)

counts_vst <- read.table("counts_vst_blindT.txt", sep = "\t", header = TRUE)

counts_combat <- ComBat(
  dat         = as.matrix(counts_vst),
  batch       = batch,
  mod         = design,
  par.prior   = TRUE,
  prior.plots = FALSE
)

# Save ComBat-corrected data
write.table(counts_combat, 
            file = "counts_vst_blindT_combat.txt", 
            quote = FALSE, sep = "\t")
