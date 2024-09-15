# Heatmap (Supp Fig. 1H-I, and Supp Tables 1-2)

counts_combat = read.table('counts_vst_blindT_combat.txt')
targ_genes = read.table('deg_padj01.txt')$V1
counts_combat_targ = counts_combat[rownames(counts_combat) %in% targ_genes,]
sample_info = read.table("sample_info.txt", header = T, sep = "\t")
sample_info$CellType <- replace(sample_info$CellType, sample_info$CellType == 'MP', 'VM') 
sample_info$Clean <- replace(sample_info$Clean, sample_info$Clean == 'C', 'Clean') 
sample_info$Clean <- replace(sample_info$Clean, sample_info$Clean == 'D', 'Dirty') 
sample_info$CellType = paste0(sample_info$CellType, " ", sample_info$Age, " ", sample_info$Clean) #, sample_info$Veteran
sample_info$CellType = factor(sample_info$CellType, levels = c('TN Adult Clean', 'TN Adult Dirty', 'TN Neo Clean', 'TN Neo Dirty',
                                                 'All Adult Clean', 'All Adult Dirty', 'All Neo Clean', 'All Neo Dirty',
                                                 'VM Adult Clean', 'VM Adult Dirty', 'VM Neo Clean', 'VM Neo Dirty'))
library(pheatmap)
save_pheatmap_pdf <- function(x, filename, width=7, height=7) { # useful function from https://stackoverflow.com/a/43051932
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
library(tidyverse)
sample_info_reorder = sample_info %>% arrange(CellType)
counts_combat_targ_reorder = counts_combat_targ %>% select(gsub('-', '.', sample_info_reorder$Sample)) 
library(RColorBrewer)
colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
breaks_5 <- c(seq(-5, 0, length.out=ceiling(length(colors)/2) + 1), 
              seq(2/length(colors), 5, length.out=floor(length(colors)/2)))

library(viridis)
cols = scales::viridis_pal(option = 'cividis')(12)
cols_anno = list(CellType = rev(c('TN Adult Clean' = cols[12],
                              'TN Adult Dirty' = cols[11],
                              'TN Neo Clean' = cols[10],
                              'TN Neo Dirty' = cols[9],
                              'All Adult Clean' = cols[8],
                              'All Adult Dirty' = cols[7],
                              'All Neo Clean' = cols[6],
                              'All Neo Dirty' = cols[5],
                              'VM Adult Clean' = cols[4],
                              'VM Adult Dirty' = cols[3],
                              'VM Neo Clean' = cols[2],
                              'VM Neo Dirty' = cols[1])))
anno = sample_info_reorder[, c('Sample', 'CellType')]
rownames(anno) = gsub('-', '.', anno$Sample)
anno$Sample = NULL
table(colnames(counts_combat_targ_reorder) == rownames(anno))

labels = c('Ccr9', 'Cxcr5', 'Cxcr3', 'Cd44', 'Gzmb', 'F2r', 'Ifit3', 'Thra') 
labs.row <- rownames(counts_combat_targ_reorder)
labs.row = ifelse(labs.row %in% labels, labs.row, '')
heatmap = pheatmap(counts_combat_targ_reorder, scale = 'row', cluster_cols=F, show_colnames = F, treeheight_row = 0, # show_rownames = F, 
                   annotation_col = anno, breaks = breaks_5, annotation_colors = cols_anno, labels_row = labs.row)
save_pheatmap_pdf(heatmap, "heatmap/heatmap_breaks5_wide_legendYlBu_someLabels.pdf", width=6, height=5)

heatmapTable = counts_combat_targ_reorder_scaled[heatmap$tree_row[["order"]],]
colnames(heatmapTable) = anno$CellType
write.table(heatmapTable, "heatmap/heatmap_breaks5_wide_legendYlBu.txt",
            quote = F, sep = '\t', row.names = T, col.names = T)


# heatmap - TFs
counts_combat = read.table('counts_vst_blindT_combat.txt')
tfs = read.table('regulators.txt')$V1
counts_combat_tf = counts_combat[rownames(counts_combat) %in% tfs,]
sample_info = read.table("sample_info.txt", header = T, sep = "\t")
sample_info$CellType <- replace(sample_info$CellType, sample_info$CellType == 'MP', 'VM') 
sample_info$Clean <- replace(sample_info$Clean, sample_info$Clean == 'C', 'Clean') 
sample_info$Clean <- replace(sample_info$Clean, sample_info$Clean == 'D', 'Dirty') 
sample_info$CellType = paste0(sample_info$CellType, " ", sample_info$Age, " ", sample_info$Clean) #, sample_info$Veteran
sample_info$CellType = factor(sample_info$CellType, levels = c('TN Adult Clean', 'TN Adult Dirty', 'TN Neo Clean', 'TN Neo Dirty',
                                                               'All Adult Clean', 'All Adult Dirty', 'All Neo Clean', 'All Neo Dirty',
                                                               'VM Adult Clean', 'VM Adult Dirty', 'VM Neo Clean', 'VM Neo Dirty'))
library(pheatmap)
library(tidyverse)
sample_info_reorder = sample_info %>% arrange(CellType)
counts_combat_tf_reorder = counts_combat_tf %>% select(gsub('-', '.', sample_info_reorder$Sample)) 
library(RColorBrewer)
colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
breaks_5 <- c(seq(-5, 0, length.out=ceiling(length(colors)/2) + 1), 
              seq(2/length(colors), 5, length.out=floor(length(colors)/2)))

library(viridis)
cols = scales::viridis_pal(option = 'cividis')(12)
cols_anno = list(CellType = rev(c('TN Adult Clean' = cols[12],
                                  'TN Adult Dirty' = cols[11],
                                  'TN Neo Clean' = cols[10],
                                  'TN Neo Dirty' = cols[9],
                                  'All Adult Clean' = cols[8],
                                  'All Adult Dirty' = cols[7],
                                  'All Neo Clean' = cols[6],
                                  'All Neo Dirty' = cols[5],
                                  'VM Adult Clean' = cols[4],
                                  'VM Adult Dirty' = cols[3],
                                  'VM Neo Clean' = cols[2],
                                  'VM Neo Dirty' = cols[1])))
anno = sample_info_reorder[, c('Sample', 'CellType')]
rownames(anno) = gsub('-', '.', anno$Sample)
anno$Sample = NULL
table(colnames(counts_combat_tf_reorder) == rownames(anno))

labels = c('Eomes', 'Foxo1', 'Lef1', 'Prdm1', 'Zbtb6') 
labs.row <- rownames(counts_combat_tf_reorder)
labs.row = ifelse(labs.row %in% labels, labs.row, '')
heatmap = pheatmap(counts_combat_tf_reorder, scale = 'row', cluster_cols=F, show_colnames = F, treeheight_row = 0, # show_rownames = F, 
                   annotation_col = anno, breaks = breaks_5, annotation_colors = cols_anno, labels_row = labs.row)
save_pheatmap_pdf(heatmap, "heatmap/heatmap_TFs_wide_legendYlBu_someLabels.pdf", width=6, height=3)

heatmapTable = counts_combat_tf_reorder_scaled[heatmap$tree_row[["order"]],]
write.table(heatmapTable, "heatmap/heatmap_TFs_wide_legendYlBu.txt",
            quote = F, sep = '\t', row.names = T, col.names = T)
