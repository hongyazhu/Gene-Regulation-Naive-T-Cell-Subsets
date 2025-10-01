##############################################
# RNA-seq Preprocessing and Batch Correction #
# 1. Load counts & metadata
# 2. VST normalization (DESeq2)
# 3. Batch correction (ComBat)
##############################################

library(DESeq2)
library(limma)
library(sva)

# Load data
counts_all <- read.table("counts_raw.txt", header = TRUE)
sample_info <- read.table("sample_info.txt", header = TRUE, sep = "\t")

factor_cols <- c("Source", "Age", "CellType", "Veteran", "Clean")
sample_info[factor_cols] <- lapply(sample_info[factor_cols], factor)

# DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = counts_all,
  colData   = sample_info,
  design    = ~ Source + Age + CellType + Veteran + Clean
)

# Filter genes with >=1 count
dds <- dds[rowSums(counts(dds) > 0) >= 1, ]

# VST transform
vsd <- vst(dds, blind = TRUE)
write.table(assay(vsd), "counts_vst_blindT.txt", quote = FALSE, sep = "\t")

# ComBat batch correction
batch <- sample_info$Source
design <- model.matrix(~ Age + CellType + Veteran + Clean, data = sample_info)
counts_vst <- read.table("counts_vst_blindT.txt", sep = "\t", header = TRUE)

counts_combat <- ComBat(as.matrix(counts_vst), batch, design,
                        par.prior = TRUE, prior.plots = FALSE)

write.table(counts_combat, "counts_vst_blindT_combat.txt", quote = FALSE, sep = "\t")
