### Preprocessing

counts = read.table('counts_ATACseq',
                    row.names = 1, header = T)
counts = counts[, 6:ncol(counts)]
library(edgeR)
counts_filtered = counts[rowMeans(cpm(counts)) >= 1,]

librarySizes <- colSums(counts_filtered)
par(mar=c(7,4,4,4))
barplot(librarySizes, 
        names=names(librarySizes), 
        las=3, 
        main="Barplot of library sizes")

library(readxl)
sample_info = read_excel('sample_info_ATACseq.xlsx')
library(DESeq2)
library(limma) 
library(DESeq2)
library(ggplot2)

sample_info$Project <- factor(sample_info$Project) 
sample_info$Age <- factor(sample_info$Age) 
sample_info$CellType <- factor(sample_info$CellType) 
sample_info$Veteran <- factor(sample_info$Veteran) 
sample_info$Clean <- factor(sample_info$Clean) 

# Define experimental factors and design matrix. 
batch <- factor(sample_info$Project) 
age <- factor(sample_info$Age) 
cellType <- factor(sample_info$CellType) 
#state <- factor(sample_info$State) 
veteran <- factor(sample_info$Veteran) 
clean <- factor(sample_info$Clean) 

design <- model.matrix(~age+cellType+veteran+clean)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_info,
                              design= ~ Project + Age + CellType + Veteran + Clean)
keep <- rowSums(counts(dds) > 0) >= 1
dds <- dds[keep,]

# DESeq2 vst blind=T
vsd <- vst(dds, blind=T)


### PCA before batch effect removal (Supp Fig 2A)

ntop <- 500
rv <- rowVars(assay(vsd))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat <- t( assay(vsd)[select, ] )
pca <- prcomp(mat)
library(ggfortify)
autoplot(pca, data = sample_info,colour = 'Project', shape = 'CellType', size = 5) +
  theme_light(base_size = 15) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black"))
# PC1 88.57%, PC2 5.16%

autoplot(pca, x = 2, y = 3, data = sample_info,colour = 'Project', shape = 'CellType', size = 5) +
  theme_light(base_size = 15) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black"))
# PC2 5.16%, PC3 1.29%

pca_plot = as.data.frame(pca$x)
pca_plot$Sample = sample_info$Sample
pca_plot$Project = sample_info$Project
pca_plot$Age = sample_info$Age
pca_plot$CellType = sample_info$CellType
#pca_plot$CellType <- replace(pca_plot$CellType, pca_plot$CellType == 'MP', 'VM') 
pca_plot$Veteran = sample_info$Veteran
pca_plot$Clean = sample_info$Clean

library(RColorBrewer)
mycolors = c(brewer.pal(name="Set2", n = 6))[c(1:3,6)]
ggplot(pca_plot, aes(x = PC1, y= PC2)) +
  geom_point(size = 3, aes(colour = Project, shape = CellType)) +
  theme_bw(base_size = 14) +
  xlab("PC1 (88.57%)") +
  ylab("PC2 (5.16%)") +
  labs(color='Source') +
  labs(shape='Cell type') +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2)) +
  scale_color_manual(limits = c('cleanproject', 'mir29project', 'Cell2018', 'veteranproject'),
                     labels = c('Tabilas et al. (2022)', 'Yee Mon et al. (2021)', 'Smith et al. (2018)', 'Veteran Project'),
                     values = mycolors)
ggsave("PCA_vst_blindT_sourceCelltype.pdf", plot = last_plot(), device = "pdf",
       width = 6.5, height = 4, dpi = 300)


# batch effect removal using combat 
library(sva)
counts_combat = ComBat(dat=as.matrix(assay(vsd)), batch=batch, mod=design, par.prior=TRUE, prior.plots=FALSE)
library(ggfortify)
ntop <- 500
rv <- rowVars(counts_combat)
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat <- t( counts_combat[select, ] )


### PCA after batch effect removal (Fig 2A)

pca<-prcomp(mat)
rownames(pca$x) =colnames(counts_combat)
#autoplot(pca, label = T, label.size = 3.5) 

autoplot(pca, data = sample_info,colour = 'Project', shape = 'CellType', size = 5) +
  theme_light(base_size = 15) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black"))

autoplot(pca, y = 3)
autoplot(pca, y = 4)
# PC1: 81.38%
# PC2: 4.85%
# PC3: 3.12%
# PC4: 2.01%

pca_plot = as.data.frame(pca$x)
pca_plot$Sample = sample_info$Sample
pca_plot$Project = sample_info$Project
pca_plot$Age = sample_info$Age
pca_plot$CellType = sample_info$CellType
#pca_plot$CellType <- replace(pca_plot$CellType, pca_plot$CellType == 'MP', 'VM') 
pca_plot$Veteran = sample_info$Veteran
pca_plot$Clean = sample_info$Clean
pca_plot$info = paste0(pca_plot$CellType, " ", pca_plot$Age, " ", pca_plot$Clean) #, pca_plot$Veteran
pca_plot$info = factor(pca_plot$info, levels = c('TN Adult Clean', 'TN Adult Dirty', 'TN Neo Clean', 'TN Neo Dirty',
                                                 'All Adult Clean', 'All Adult Dirty', 'All Neo Clean', 'All Neo Dirty',
                                                 'VM Adult Clean', 'VM Adult Dirty', 'VM Neo Clean', 'VM Neo Dirty'))


pca_plot$age_clean = paste0(pca_plot$Age, " ", pca_plot$Clean) 
library(ggtext)
ggplot(pca_plot, aes(x = PC1, y= PC2)) +
  geom_point(size = 3, aes(colour = CellType, shape = age_clean)) +
  theme_bw(base_size = 14) +
  xlab("PC1 (81.38%)") +
  ylab("PC2 (4.85%)") +
  labs(color='Cell type') +
  labs(shape='Age') +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2)) +
  #theme(legend.box = "horizontal") +
  scale_color_manual(labels = c('VM', 'Bulk', 'TN'),
                     values = c('#004D7F', '#3D83BDFF', '#56C1FF'),
                     limits = c('VM', 'Bulk', 'TN'))+
  scale_shape_manual(limits = c('Adult Clean', 'Neo Clean', 'Adult Dirty', 'Neo Dirty'),
                     values = c(16,1,17,2))+
  theme(legend.text=element_markdown(size=12)) 
ggsave("PCA_vst_blindT_combat_CelltypeCleanAge.pdf", plot = last_plot(), device = "pdf",
       width = 6, height = 4, dpi = 300)




### PCA types (Supp Fig 2B-G)

# pc12, colored by each feature
pc12_age = ggplot(pca_plot, aes(x = PC1, y= PC2)) +
  geom_point(size = 3, aes(colour = Age)) +
  theme_bw(base_size = 14) +
  xlab("PC1 (81.38%)") +
  ylab("PC2 (4.85%)") +
  labs(color='Type') +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2)) +
  theme(legend.box = "horizontal") +
  scale_color_manual(limits = c('Neo', 'Adult'),
                     values = c('#017100', '#88FA4E')) +
  theme(legend.text=element_markdown(size=12)) +
  guides(shape = "none")

pc12_celltype = ggplot(pca_plot, aes(x = PC1, y= PC2)) +
  geom_point(size = 3, aes(colour = CellType)) +
  theme_bw(base_size = 14) +
  xlab("PC1 (81.38%)") +
  ylab("PC2 (4.85%)") +
  labs(color='Type') +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2)) +
  theme(legend.box = "horizontal") +
  scale_color_manual(limits = c('TN', 'Bulk', 'VM'),
                     values = c('#56C1FF', "#3D83BDFF", '#004D7F')) +
  theme(legend.text=element_markdown(size=12)) +
  guides(shape = "none")

pc12_clean = ggplot(pca_plot, aes(x = PC1, y= PC2)) +
  geom_point(size = 3, aes(colour = Clean)) +
  theme_bw(base_size = 14) +
  xlab("PC1 (81.38%)") +
  ylab("PC2 (4.85%)") +
  labs(color='Type') +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2)) +
  theme(legend.box = "horizontal") +
  scale_color_manual(limits = c('Dirty', 'Clean'),
                     values = c('#FF9300', '#FFFC66')) +
  theme(legend.text=element_markdown(size=12)) +
  guides(shape = "none")

library(cowplot)
pc12.combined = plot_grid(pc12_celltype, pc12_age, pc12_clean, align="hv", ncol = 2)
ggsave(
  paste0('PCA_types_pc12_noshape.pdf'),
  plot = pc12.combined,
  device = "pdf",
  width = 10,
  height = 7,
  dpi = 300
)

# pc13, colored by each feature
pc13_age = ggplot(pca_plot, aes(x = PC1, y= PC3)) +
  geom_point(size = 3, aes(colour = Age)) +
  theme_bw(base_size = 14) +
  xlab("PC1 (81.38%)") +
  ylab("PC3 (3.12%)") +
  labs(color='Type') +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2)) +
  theme(legend.box = "horizontal") +
  scale_color_manual(limits = c('Neo', 'Adult'),
                     values = c('#017100', '#88FA4E')) +
  theme(legend.text=element_markdown(size=12)) +
  guides(shape = "none")

pc13_celltype = ggplot(pca_plot, aes(x = PC1, y= PC3)) +
  geom_point(size = 3, aes(colour = CellType)) +
  theme_bw(base_size = 14) +
  xlab("PC1 (81.38%)") +
  ylab("PC3 (3.12%)") +
  labs(color='Type') +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2)) +
  theme(legend.box = "horizontal") +
  scale_color_manual(limits = c('TN', 'Bulk', 'VM'),
                     values = c('#56C1FF', "#3D83BDFF", '#004D7F')) +
  theme(legend.text=element_markdown(size=12)) +
  guides(shape = "none")

pc13_clean = ggplot(pca_plot, aes(x = PC1, y= PC3)) +
  geom_point(size = 3, aes(colour = Clean)) +
  theme_bw(base_size = 14) +
  xlab("PC1 (81.38%)") +
  ylab("PC3 (3.12%)") +
  labs(color='Type') +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2)) +
  theme(legend.box = "horizontal") +
  scale_color_manual(limits = c('Dirty', 'Clean'),
                     values = c('#FF9300', '#FFFC66')) +
  theme(legend.text=element_markdown(size=12)) +
  guides(shape = "none")

pc13.combined = plot_grid(pc13_celltype, pc13_age, pc13_clean, align="hv", ncol = 2)
ggsave(
  paste0('PCA_types_pc13_noshape.pdf'),
  plot = pc13.combined,
  device = "pdf",
  width = 10,
  height = 7,
  dpi = 300
)


### eigencorplot (Fig 2B)

library(PCAtools)
# repeat previous result
p0 <- pca(t(mat))
biplot(p0)
# add metadata
metadata = as.data.frame(sample_info)
#rownames(metadata) = gsub('-', '.', metadata$Sample)
rownames(metadata) = metadata$Sample
metadata$Sample = NULL
metadata$File_path = NULL
metadata$File_name = NULL
metadata$Project = factor(metadata$Project, levels = c("Cell2018","mir29project","veteranproject","cleanproject")) # randomly ordered
metadata$Age = factor(metadata$Age, levels = c("Adult", 'Neo')) 
metadata$CellType = factor(metadata$CellType, levels = c("TN", 'Bulk','VM')) 
metadata$Veteran = factor(metadata$Veteran, levels = c("Veteran", 'RTE','Mature')) 
metadata$Clean = factor(metadata$Clean, levels = c("Clean", 'Dirty')) 
mat_t = t(mat)
colnames(mat_t) = rownames(metadata)
p <- pca(mat_t, metadata = metadata)
metadata_tmp = metadata
colnames(metadata_tmp) = c('Project', 'VM/TN', 'Neonatal/Adult', 'Dirty/Clean', 'Veteran')
metadata_tmp = metadata_tmp[,rev(c(2,3,4,5,1))]
p <- pca(mat_t, metadata = metadata_tmp)

pdf("PCA_eigencorplot_r.pdf", width = 7, height = 4)
library(RColorBrewer)
eigencorplot(p,
             components = getComponents(p, 1:3),
             metavars = rev(c('VM/TN','Neonatal/Adult','Dirty/Clean')),
             col = brewer.pal(n = 5, name = "RdYlGn"),
             cexCorval = 1.2,
             fontCorval = 2,
             posLab = 'all',
             #rotLabX = 45,
             #scale = TRUE,
             main = bquote(PC ~ Pearson ~ r ~ metadata ~ correlates),
             #plotRsquared = TRUE,
             corFUN = 'pearson',
             corUSE = 'pairwise.complete.obs',
             corMultipleTestCorrection = 'BH',
             signifSymbols = c('****', '***', '**', '*', ''),
             signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1))
dev.off()

