# make RNAseq pca-related plots 

### before combat, show source and cell type (Supplementary Fig. 1A)

counts = read.table("counts_vst_blindT.txt", sep = "\t", header = T)
sample_info = read.table("sample_info.txt", header = T, sep = "\t")

library(ggfortify)
ntop <- 500
library(matrixStats)
rv <- rowVars(as.matrix(counts))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat <- t( counts[select, ] )

pca<-prcomp(mat)
rownames(pca$x) =colnames(counts)
autoplot(pca, label = T, label.size = 3.5) 
# PC1: 32.56%, PC2: 20.5%

pca_plot = as.data.frame(pca$x)
pca_plot$Sample = sample_info$Sample
pca_plot$Source = sample_info$Source
pca_plot$Age = sample_info$Age
pca_plot$CellType = sample_info$CellType
pca_plot$CellType <- replace(pca_plot$CellType, pca_plot$CellType == 'MP', 'VM') 
pca_plot$Veteran = sample_info$Veteran
pca_plot$Clean = sample_info$Clean

library(RColorBrewer)
mycolors = c(brewer.pal(name="Set2", n =6))
ggplot(pca_plot, aes(x = PC1, y= PC2)) +
  geom_point(size = 3, aes(colour = Source, shape = CellType)) +
  theme_bw(base_size = 14) +
  xlab("PC1 (32.56%)") +
  ylab("PC2 (20.5%)") +
  labs(color='Source') +
  labs(shape='Cell type') +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2)) +
  scale_color_manual(limits = c('Clean Project', 'Mir29 Project', 'Smith 2018', 'Wissink 2015', 'Lin28 Project', 'Veteran Project'),
                     labels = c('Tabilas et al. (2022)', 'Yee Mon et al. (2021)', 'Smith et al. (2018)', 'Wissink et al. (2015)', 'Lin28 Project', 'Veteran Project'),
                     values = mycolors)
ggsave("PCA/PCA_vst_blindT_sourceCelltype.pdf", plot = last_plot(), device = "pdf",
       width = 6.5, height = 4, dpi = 300)


### PCA after batch effect removal (Fig. 1B)

counts_combat = read.table("counts_vst_blindT_combat.txt", sep = "\t", header = T)
sample_info = read.table("sample_info.txt", header = T, sep = "\t")

library(ggfortify)
ntop <- 500
library(matrixStats)
rv <- rowVars(as.matrix(counts_combat))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat <- t( counts_combat[select, ] )

pca<-prcomp(mat)
rownames(pca$x) =colnames(counts_combat)
autoplot(pca, label = T, label.size = 3.5) 
# PC1: 60.18%, PC2: 9.27%

pca_plot = as.data.frame(pca$x)
pca_plot$Sample = sample_info$Sample
pca_plot$Source = sample_info$Source
pca_plot$Age = sample_info$Age
pca_plot$CellType = sample_info$CellType
pca_plot$CellType <- replace(pca_plot$CellType, pca_plot$CellType == 'MP', 'VM') 
pca_plot$Veteran = sample_info$Veteran
pca_plot$Clean = sample_info$Clean
pca_plot$Clean <- replace(pca_plot$Clean, pca_plot$Clean == 'C', 'Clean') 
pca_plot$Clean <- replace(pca_plot$Clean, pca_plot$Clean == 'D', 'Dirty') 
pca_plot$info = paste0(pca_plot$CellType, " ", pca_plot$Age, " ", pca_plot$Clean) #, pca_plot$Veteran
pca_plot$info = factor(pca_plot$info, levels = c('TN Adult Clean', 'TN Adult Dirty', 'TN Neo Clean', 'TN Neo Dirty',
                                                 'All Adult Clean', 'All Adult Dirty', 'All Neo Clean', 'All Neo Dirty',
                                                 'VM Adult Clean', 'VM Adult Dirty', 'VM Neo Clean', 'VM Neo Dirty'))
pca_plot$age_clean = paste0(pca_plot$Age, " ", pca_plot$Clean) 
library(ggtext)
ggplot(pca_plot, aes(x = PC1, y= PC2)) +
  geom_point(size = 3, aes(colour = CellType, shape = age_clean)) +
  theme_bw(base_size = 14) +
  xlab("PC1 (60.18%)") +
  ylab("PC2 (9.27%)") +
  labs(color='Cell type') +
  labs(shape='Age') +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2)) +
  #theme(legend.box = "horizontal") +
  scale_color_manual(labels = c('VM', 'Bulk', 'TN'),
                     values = c('#004D7F', '#3D83BDFF', '#56C1FF'),
                     limits = c('VM', 'All', 'TN'))+
  scale_shape_manual(limits = c('Adult Clean', 'Neo Clean', 'Adult Dirty', 'Neo Dirty'),
                     values = c(16,1,17,2))+
  theme(legend.text=element_markdown(size=12)) 
ggsave("PCA/PCA_vst_blindT_combat_CelltypeCleanAge.pdf", plot = last_plot(), device = "pdf",
       width = 6, height = 4, dpi = 300)


### get loadings (Fig. 1D)

loadings <- as.data.frame(pca$rotation[, seq_len(2)])
write.table(loadings, 'PCA/PCA_loadings12.txt', quote = F, sep = '\t')
loadings$genename = rownames(loadings)
for (i in 1:nrow(loadings)){
  if (loadings[i, 'PC1'] > 0){
    loadings[i, 'PC1_info'] = 'PC1+'
  } else if (loadings[i, 'PC1'] < 0){
    loadings[i, 'PC1_info'] = 'PC1-'
  }
}
loadings_sortPC1 <- loadings[order(loadings$PC1),]
ngenes = 5
loadings_sortPC1_toplot <- loadings_sortPC1[c(1:ngenes, (nrow(loadings)-ngenes+1):nrow(loadings)),]
ggplot(loadings_sortPC1_toplot, aes(y=reorder(genename, PC1), x=PC1)) +
  geom_point(stat='identity', aes(color=PC1_info), size = 3, shape = 17) +
  ylab('')+
  #ggtitle(paste0(largeclusterName, "_dRNA_top", ngenes)) + 
  #scale_color_discrete() +
  theme_bw() + 
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        axis.ticks = element_line(colour = "grey70", size = 0.2),
        panel.grid.major = element_line(colour = "grey70", size = 0.2),
        panel.grid.minor = element_blank())+
  geom_vline(aes(xintercept = 0)) +
  scale_color_manual(limits = c('PC1+', 'PC1-'), values = c('#00306FFF', '#EAD357FF'))
ggsave(filename = paste0("PCA/PC1_loadings_top", ngenes, ".pdf"), device = "pdf", height = 3, width = 4.5)


### make eigencorplot (Fig. 1C)

library(PCAtools)
# repeat previous result
p0 <- pca(t(mat))
biplot(p0)
# add metadata
metadata = sample_info
rownames(metadata) = gsub('-', '.', metadata$Sample)
metadata$Sample = NULL
metadata$Source = factor(metadata$Source, levels = c("Wissink 2015","Smith 2018","Mir29 Project","Lin28 Project","Veteran Project","Clean Project")) # randomly ordered
metadata$Age = factor(metadata$Age, levels = c("Adult", 'Neo')) 
metadata$CellType = factor(metadata$CellType, levels = c("TN", 'All','MP')) 
metadata$Veteran = factor(metadata$Veteran, levels = c("Vet", 'RTE','Mat')) 
metadata$Clean = factor(metadata$Clean, levels = c("C", 'D')) 
p <- pca(t(mat), metadata = metadata)
metadata_tmp = metadata
colnames(metadata_tmp) = c('Source', 'Neonatal/Adult', 'VM/TN', 'Veteran', 'Dirty/Clean')
metadata_tmp = metadata_tmp[,rev(c(3,2,5,4,1))]
p <- pca(t(mat), metadata = metadata_tmp)
pdf("PCA/PCA_eigencorplot_r.pdf", width = 7, height = 4)
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


### PCAs showing each type (Supp Fig. 1B-G)

# after combat


counts_combat = read.table("counts_vst_blindT_combat.txt", sep = "\t", header = T)
sample_info = read.table("sample_info.txt", header = T, sep = "\t")

library(ggfortify)
ntop <- 500
library(matrixStats)
rv <- rowVars(as.matrix(counts_combat))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat <- t( counts_combat[select, ] )

pca<-prcomp(mat)
rownames(pca$x) =colnames(counts_combat)
autoplot(pca, label = T, label.size = 3.5) 
autoplot(pca, x=1, y=3, label = T, label.size = 3.5) 
autoplot(pca, x=1, y=4, label = T, label.size = 3.5) 
# PC1: 60.18%, PC2: 9.27%
# PC3: 6.04%
# PC4: 3.99%

pca_plot = as.data.frame(pca$x)
pca_plot$Sample = sample_info$Sample
pca_plot$Source = sample_info$Source
pca_plot$Age = sample_info$Age
pca_plot$CellType = sample_info$CellType
pca_plot$CellType <- replace(pca_plot$CellType, pca_plot$CellType == 'MP', 'VM') 
pca_plot$Veteran = sample_info$Veteran
pca_plot$Clean = sample_info$Clean
pca_plot$Clean <- replace(pca_plot$Clean, pca_plot$Clean == 'C', 'Clean') 
pca_plot$Clean <- replace(pca_plot$Clean, pca_plot$Clean == 'D', 'Dirty') 
pca_plot$info = paste0(pca_plot$CellType, " ", pca_plot$Age, " ", pca_plot$Clean) #, pca_plot$Veteran
pca_plot$info = factor(pca_plot$info, levels = c('TN Adult Clean', 'TN Adult Dirty', 'TN Neo Clean', 'TN Neo Dirty',
                                                 'All Adult Clean', 'All Adult Dirty', 'All Neo Clean', 'All Neo Dirty',
                                                 'VM Adult Clean', 'VM Adult Dirty', 'VM Neo Clean', 'VM Neo Dirty'))

# pc12, colored by each feature
pc12_age = ggplot(pca_plot, aes(x = PC1, y= PC2)) +
  geom_point(size = 3, aes(colour = Age)) +
  theme_bw(base_size = 14) +
  xlab("PC1 (60.18%)") +
  ylab("PC2 (9.27%)") +
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
  xlab("PC1 (60.18%)") +
  ylab("PC2 (9.27%)") +
  labs(color='Type') +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2)) +
  theme(legend.box = "horizontal") +
  scale_color_manual(limits = c('TN', 'All', 'VM'),
                     values = c('#56C1FF', "#3D83BDFF", '#004D7F')) +
  theme(legend.text=element_markdown(size=12)) +
  guides(shape = "none")

pc12_clean = ggplot(pca_plot, aes(x = PC1, y= PC2)) +
  geom_point(size = 3, aes(colour = Clean)) +
  theme_bw(base_size = 14) +
  xlab("PC1 (60.18%)") +
  ylab("PC2 (9.27%)") +
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
  paste0('PCA/PCA_types_pc12_noshape.pdf'),
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
  xlab("PC1 (60.18%)") +
  ylab("PC3 (6.04%)") +
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
  xlab("PC1 (60.18%)") +
  ylab("PC3 (6.04%)") +
  labs(color='Type') +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2)) +
  theme(legend.box = "horizontal") +
  scale_color_manual(limits = c('TN', 'All', 'VM'),
                     values = c('#56C1FF', "#3D83BDFF", '#004D7F')) +
  theme(legend.text=element_markdown(size=12)) +
  guides(shape = "none")

pc13_clean = ggplot(pca_plot, aes(x = PC1, y= PC3)) +
  geom_point(size = 3, aes(colour = Clean)) +
  theme_bw(base_size = 14) +
  xlab("PC1 (60.18%)") +
  ylab("PC3 (6.04%)") +
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
  paste0('PCA/PCA_types_pc13_noshape.pdf'),
  plot = pc13.combined,
  device = "pdf",
  width = 10,
  height = 7,
  dpi = 300
)

