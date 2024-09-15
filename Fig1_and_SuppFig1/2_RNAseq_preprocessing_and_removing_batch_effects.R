library(limma) 
library(DESeq2)
library(ggplot2)

counts_all = read.table("counts_raw.txt", header = T)
sample_info = read.table("sample_info.txt", header = T, sep = "\t")

sample_info$Source <- factor(sample_info$Source) 
sample_info$Age <- factor(sample_info$Age) 
sample_info$CellType <- factor(sample_info$CellType) 
#state <- factor(sample_info$State) 
sample_info$Veteran <- factor(sample_info$Veteran) 
sample_info$Clean <- factor(sample_info$Clean) 

# Define experimental factors and design matrix. 
batch <- factor(sample_info$Source) 
age <- factor(sample_info$Age) 
cellType <- factor(sample_info$CellType) 
#state <- factor(sample_info$State) 
veteran <- factor(sample_info$Veteran) 
clean <- factor(sample_info$Clean) 

design <- model.matrix(~age+cellType+veteran+clean)

dds <- DESeqDataSetFromMatrix(countData = counts_all,
                              colData = sample_info,
                              design= ~ Source + Age + CellType + Veteran + Clean)
keep <- rowSums(counts(dds) > 0) >= 1
dds <- dds[keep,]


# DESeq2 vst blind=T
vsd <- vst(dds, blind=T)
write.table(assay(vsd), "counts_vst_blindT.txt", quote = F, sep = "\t")


# remove batche effect with combat
library(sva)

counts = read.table("counts_vst_blindT.txt", sep = "\t", header = T)
sample_info = read.table("sample_info.txt", header = T, sep = "\t")

sample_info$Source <- factor(sample_info$Source) 
sample_info$Age <- factor(sample_info$Age) 
sample_info$CellType <- factor(sample_info$CellType) 
#state <- factor(sample_info$State) 
sample_info$Veteran <- factor(sample_info$Veteran) 
sample_info$Clean <- factor(sample_info$Clean) 

batch <- factor(sample_info$Source) 
age <- factor(sample_info$Age) 
cellType <- factor(sample_info$CellType) 
#state <- factor(sample_info$State) 
veteran <- factor(sample_info$Veteran) 
clean <- factor(sample_info$Clean) 
design <- model.matrix(~age+cellType+veteran+clean)

counts_combat = ComBat(dat=as.matrix(counts), batch=batch, mod=design, par.prior=TRUE, prior.plots=FALSE)
write.table(counts_combat, "counts_vst_blindT_combat.txt", quote = F, sep = "\t")

