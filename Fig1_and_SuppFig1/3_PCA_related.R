########################################
# RNA-seq PCA Plots & Loadings
# 1. PCA before ComBat (Supp Fig. 1A)
# 2. PCA after ComBat (Fig. 1B, Supp Figs. 1B–G)
# 3. Loadings (Fig. 1D)
# 4. Eigencorplot (Fig. 1C)
########################################

library(ggfortify)
library(matrixStats)
library(ggplot2)
library(RColorBrewer)
library(ggtext)
library(cowplot)
library(PCAtools)

# Helper: run PCA on top variable genes
run_pca <- function(counts, ntop = 500) {
  rv <- rowVars(as.matrix(counts))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  mat <- t(counts[select, ])
  prcomp(mat)
}

#----------------------------------------
# PCA before ComBat (Supp Fig. 1A)
#----------------------------------------
counts <- read.table("counts_vst_blindT.txt", sep="\t", header=TRUE)
sample_info <- read.table("sample_info.txt", sep="\t", header=TRUE)

pca <- run_pca(counts)
pca_plot <- cbind(as.data.frame(pca$x), sample_info)
pca_plot$CellType <- replace(pca_plot$CellType, pca_plot$CellType=="MP", "VM")

ggplot(pca_plot, aes(PC1, PC2, colour=Source, shape=CellType)) +
  geom_point(size=3) +
  theme_bw(base_size=14) +
  xlab("PC1 (32.56%)") + ylab("PC2 (20.5%)") +
  labs(color="Source", shape="Cell type") +
  scale_color_manual(
    limits = c('Clean Project','Mir29 Project','Smith 2018','Wissink 2015','Lin28 Project','Veteran Project'),
    labels = c('Tabilas et al. (2022)','Yee Mon et al. (2021)','Smith et al. (2018)','Wissink et al. (2015)','Lin28 Project','Veteran Project'),
    values = brewer.pal(6,"Set2")
  )
ggsave("PCA/PCA_vst_blindT_sourceCelltype.pdf", width=6.5, height=4, dpi=300)

#----------------------------------------
# PCA after ComBat (Fig. 1B)
#----------------------------------------
counts_combat <- read.table("counts_vst_blindT_combat.txt", sep="\t", header=TRUE)
pca <- run_pca(counts_combat)
pca_plot <- cbind(as.data.frame(pca$x), sample_info)
pca_plot$CellType <- replace(pca_plot$CellType, pca_plot$CellType=="MP", "VM")
pca_plot$Clean <- factor(ifelse(pca_plot$Clean=="C","Clean","Dirty"))

ggplot(pca_plot, aes(PC1, PC2, colour=CellType, shape=paste(Age,Clean))) +
  geom_point(size=3) +
  theme_bw(base_size=14) +
  xlab("PC1 (60.18%)") + ylab("PC2 (9.27%)") +
  labs(color="Cell type", shape="Age") +
  scale_color_manual(values=c('VM'='#004D7F','All'='#3D83BDFF','TN'='#56C1FF'))
ggsave("PCA/PCA_vst_blindT_combat_CelltypeCleanAge.pdf", width=6, height=4, dpi=300)

#----------------------------------------
# Loadings (Fig. 1D)
#----------------------------------------
loadings <- as.data.frame(pca$rotation[,1:2])
loadings$genename <- rownames(loadings)
loadings$PC1_info <- ifelse(loadings$PC1 > 0,"PC1+","PC1-")
write.table(loadings,"PCA/PCA_loadings12.txt",quote=FALSE,sep="\t")

ngenes <- 5
top_loadings <- rbind(head(loadings[order(loadings$PC1),],ngenes),
                      tail(loadings[order(loadings$PC1),],ngenes))
ggplot(top_loadings, aes(PC1, reorder(genename,PC1), color=PC1_info)) +
  geom_point(size=3, shape=17) +
  theme_bw() + xlab("PC1 loading") + ylab("") +
  geom_vline(xintercept=0) +
  scale_color_manual(values=c('PC1+'='#00306FFF','PC1-'='#EAD357FF'))
ggsave("PCA/PC1_loadings_top5.pdf", width=4.5, height=3, dpi=300)

#----------------------------------------
# Eigencorplot (Fig. 1C)
#----------------------------------------
metadata <- sample_info
rownames(metadata) <- gsub("-",".",metadata$Sample)
metadata$Sample <- NULL
metadata$Source <- factor(metadata$Source, levels=c("Wissink 2015","Smith 2018","Mir29 Project","Lin28 Project","Veteran Project","Clean Project"))
metadata$Age <- factor(metadata$Age, levels=c("Adult","Neo"))
metadata$CellType <- factor(metadata$CellType, levels=c("TN","All","MP"))
metadata$Veteran <- factor(metadata$Veteran, levels=c("Vet","RTE","Mat"))
metadata$Clean <- factor(metadata$Clean, levels=c("C","D"))

p <- pca(t(counts_combat), metadata=metadata)
pdf("PCA/PCA_eigencorplot_r.pdf", width=7, height=4)
eigencorplot(p,
             components=getComponents(p,1:3),
             metavars=c("CellType","Age","Clean"),
             col=brewer.pal(5,"RdYlGn"),
             corFUN="pearson",
             corUSE="pairwise.complete.obs")
dev.off()

#----------------------------------------
# Supp Figs. 1B–G: PCAs by features
#----------------------------------------
make_pc_plot <- function(df, x, y, feature, values, filename) {
  ggplot(df, aes_string(x,y,color=feature)) +
    geom_point(size=3) +
    theme_bw(base_size=14) +
    scale_color_manual(values=values) +
    ggsave(filename, width=7, height=5, dpi=300)
}

# Example: PC1 vs PC2 by CellType/Age/Clean
pc12_celltype <- ggplot(pca_plot, aes(PC1,PC2,color=CellType)) + geom_point(size=3) + theme_bw()
pc12_age      <- ggplot(pca_plot, aes(PC1,PC2,color=Age)) + geom_point(size=3) + theme_bw()
pc12_clean    <- ggplot(pca_plot, aes(PC1,PC2,color=Clean)) + geom_point(size=3) + theme_bw()

plot_grid(pc12_celltype, pc12_age, pc12_clean, ncol=2, align="hv")
ggsave("PCA/PCA_types_pc12_noshape.pdf", width=10, height=7, dpi=300)
