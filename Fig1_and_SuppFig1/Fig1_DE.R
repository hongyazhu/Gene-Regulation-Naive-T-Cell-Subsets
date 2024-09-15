# Differential expression analysis

library(DESeq2)
padj_cutoff = 1 # include all

# clean vs dirty

counts_dirty = read.table("Tabilas2022_RNAseq/counts_name.txt", header = T, row.names = 1)
counts_dirty = counts_dirty[,6:ncol(counts_dirty)]
coln_tmp = sub('.*CD', 'CD', colnames(counts_dirty))
coln_tmp = sub('_CKDL.*', '', coln_tmp)
targetFile = read.table('Tabilas2022_RNAseq/targetFile.txt', header = T)
targetFile$cd = sub('.ReadsPerGene.out.tab.rawCounts', '', targetFile$files)
colnames(counts_dirty) = targetFile[match(coln_tmp, targetFile$cd),]$label
counts_dirty_nobulk = counts_dirty[, !grepl("B.", colnames(counts_dirty))]

coldata_dirty_nobulk = data.frame(c("Neo", "Neo", "Neo", "Neo", "Adult", "Adult", "Adult", "Adult", "Neo", "Neo", "Neo", "Neo", "Adult", "Adult", "Adult", "Adult", "Neo", "Neo", "Neo", "Neo", "Adult", "Adult", "Adult", "Adult"),
                                  gsub("VM", "MP", substr(colnames(counts_dirty_nobulk), 1, 2)), 
                                  gsub("\\..*","", gsub("^.*?\\.","", colnames(counts_dirty_nobulk))))
colnames(coldata_dirty_nobulk) = c('age', 'type', 'clean')
rownames(coldata_dirty_nobulk) = colnames(counts_dirty_nobulk)
coldata_dirty_nobulk$age = as.factor(coldata_dirty_nobulk$age)
coldata_dirty_nobulk$type = as.factor(coldata_dirty_nobulk$type)
coldata_dirty_nobulk$clean = as.factor(coldata_dirty_nobulk$clean)

dds <- DESeqDataSetFromMatrix(countData = counts_dirty_nobulk,
                              colData = coldata_dirty_nobulk,
                              design= ~ age + type + clean)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="clean_D_vs_C")
resSig.clean_dirty <- subset(res, padj < padj_cutoff)
dirty <- as.data.frame(subset(resSig.clean_dirty, log2FoldChange > 0))
clean <- as.data.frame(subset(resSig.clean_dirty, log2FoldChange < 0))



# neo vs adult

counts_dirty = read.table("Tabilas2022_RNAseq/counts_name.txt", header = T, row.names = 1)
counts_dirty = counts_dirty[,6:ncol(counts_dirty)]
coln_tmp = sub('.*CD', 'CD', colnames(counts_dirty))
coln_tmp = sub('_CKDL.*', '', coln_tmp)
targetFile = read.table('targetFile.txt', header = T)
targetFile$cd = sub('.ReadsPerGene.out.tab.rawCounts', '', targetFile$files)
colnames(counts_dirty) = targetFile[match(coln_tmp, targetFile$cd),]$label
counts_dirty_nobulk = counts_dirty[, !grepl("B.", colnames(counts_dirty))]

coldata_dirty_nobulk = data.frame(c("Neo", "Neo", "Neo", "Neo", "Adult", "Adult", "Adult", "Adult", "Neo", "Neo", "Neo", "Neo", "Adult", "Adult", "Adult", "Adult", "Neo", "Neo", "Neo", "Neo", "Adult", "Adult", "Adult", "Adult"),
                                  gsub("VM", "MP", substr(colnames(counts_dirty_nobulk), 1, 2)), 
                                  gsub("\\..*","", gsub("^.*?\\.","", colnames(counts_dirty_nobulk))))
colnames(coldata_dirty_nobulk) = c('age', 'type', 'clean')
rownames(coldata_dirty_nobulk) = colnames(counts_dirty_nobulk)
coldata_dirty_nobulk$source = 'Dirty_Project'

counts_cell2018 = read.table("Smith2018_RNAseq/counts_name.txt", header = T, row.names = 1)
counts_cell2018 = counts_cell2018[,6:ncol(counts_cell2018)]
colnames(counts_cell2018) = c("all_neo_rep1", "all_neo_rep2", "all_adult_rep1", "all_adult_rep2",
                              "vm_neo_rep1", "vm_neo_rep2", "vm_adult_rep1", "vm_adult_rep2",
                              "tn_adult_rep1", "tn_adult_rep2", "vm_neo_5dpi_rep1", "vm_neo_5dpi_rep2",
                              "vm_adult_5dpi_rep1", "vm_adult_5dpi_rep2")

counts_cell2018 = counts_cell2018[, c(# 'all_neo_rep1', 'all_neo_rep2', "all_adult_rep1", "all_adult_rep2",
                                      "vm_neo_rep1", "vm_neo_rep2", "vm_adult_rep1", "vm_adult_rep2")]

coldata_cell2018 = data.frame(c('Neo', 'Neo', 'Adult', 'Adult'),
                              c('MP', 'MP', 'MP', 'MP'))
colnames(coldata_cell2018) = c('age', 'type')
rownames(coldata_cell2018) = colnames(counts_cell2018)
coldata_cell2018$clean = 'C'
coldata_cell2018$source = 'Cell_2018'


counts_adultNeo = cbind(counts_cell2018, counts_dirty_nobulk)
coldata_adultNeo = rbind(coldata_cell2018, coldata_dirty_nobulk)
coldata_adultNeo$age = as.factor(coldata_adultNeo$age)
coldata_adultNeo$type = as.factor(coldata_adultNeo$type)
coldata_adultNeo$clean = as.factor(coldata_adultNeo$clean)
coldata_adultNeo$source = as.factor(coldata_adultNeo$source)

dds <- DESeqDataSetFromMatrix(countData = counts_adultNeo,
                              colData = coldata_adultNeo,
                              design = ~ age + type + clean + source)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients

res.neo_adult <- results(dds, contrast=c("age", "Neo", "Adult"))
resSig.neo_adult <- subset(res.neo_adult, padj < padj_cutoff)
neo <- as.data.frame(subset(resSig.neo_adult, log2FoldChange > 0))
adult <- as.data.frame(subset(resSig.neo_adult, log2FoldChange < 0))


# tn vs vm

counts_dirty = read.table("Tabilas2022_RNAseq/counts_name.txt", header = T, row.names = 1)
counts_dirty = counts_dirty[,6:ncol(counts_dirty)]
coln_tmp = sub('.*CD', 'CD', colnames(counts_dirty))
coln_tmp = sub('_CKDL.*', '', coln_tmp)
targetFile = read.table('targetFile.txt', header = T)
targetFile$cd = sub('.ReadsPerGene.out.tab.rawCounts', '', targetFile$files)
colnames(counts_dirty) = targetFile[match(coln_tmp, targetFile$cd),]$label
counts_dirty_nobulk = counts_dirty[, !grepl("B.", colnames(counts_dirty))]

coldata_dirty_nobulk = data.frame(c("Neo", "Neo", "Neo", "Neo", "Adult", "Adult", "Adult", "Adult", "Neo", "Neo", "Neo", "Neo", "Adult", "Adult", "Adult", "Adult", "Neo", "Neo", "Neo", "Neo", "Adult", "Adult", "Adult", "Adult"),
                                  gsub("VM", "MP", substr(colnames(counts_dirty_nobulk), 1, 2)), 
                                  gsub("\\..*","", gsub("^.*?\\.","", colnames(counts_dirty_nobulk))))
colnames(coldata_dirty_nobulk) = c('age', 'type', 'clean')
rownames(coldata_dirty_nobulk) = colnames(counts_dirty_nobulk)
coldata_dirty_nobulk$source = 'Dirty_Project'
coldata_dirty_nobulk$rte = 'Mat'

counts_cell2018 = read.table("Smith2018_RNAseq/counts_name.txt", header = T, row.names = 1)
counts_cell2018 = counts_cell2018[,6:ncol(counts_cell2018)]
colnames(counts_cell2018) = c("all_neo_rep1", "all_neo_rep2", "all_adult_rep1", "all_adult_rep2",
                              "vm_neo_rep1", "vm_neo_rep2", "vm_adult_rep1", "vm_adult_rep2",
                              "tn_adult_rep1", "tn_adult_rep2", "vm_neo_5dpi_rep1", "vm_neo_5dpi_rep2",
                              "vm_adult_5dpi_rep1", "vm_adult_5dpi_rep2")

counts_cell2018 = counts_cell2018[, c(# 'all_neo_rep1', 'all_neo_rep2', "all_adult_rep1", "all_adult_rep2",
                                      "vm_neo_rep1", "vm_neo_rep2", "vm_adult_rep1", "vm_adult_rep2")]

coldata_cell2018 = data.frame(c('Neo', 'Neo', 'Adult', 'Adult'),
                              c('MP', 'MP', 'MP', 'MP'))
colnames(coldata_cell2018) = c('age', 'type')
rownames(coldata_cell2018) = colnames(counts_cell2018)
coldata_cell2018$clean = 'C'
coldata_cell2018$source = 'Cell_2018'
coldata_cell2018$rte = 'Mat'

counts_veteran = read.table("this_study_RNAseq/counts_name.txt", header = T, row.names = 1)
counts_veteran = counts_veteran[, 6:ncol(counts_veteran)]
colnames(counts_veteran)[1:8] = substr(colnames(counts_veteran)[1:8], 41, 48)
colnames(counts_veteran)[9:ncol(counts_veteran)] = substr(colnames(counts_veteran)[9:ncol(counts_veteran)], 11, 18)
colnames(counts_veteran) = c('RTE.TN1', 'RTE.MP1', 'Vet.MP2', 'Mat.MP2', 'RTE.MP2', 'Mat.TN3', 'RTE.TN4', 'RTE.MP4', 'Vet.TN1', 'RTE.TN3', 'Mat.MP3', 'RTE.MP3', 'Vet.TN4', 'Mat.TN4', 'Mat.TN1', 'Vet.MP4', 'Mat.MP4', 'Vet.MP1', 'Mat.MP1', 'Vet.TN2', 'Mat.TN2', 'RTE.TN2')

coldata_veteran = data.frame(substr(colnames(counts_veteran), 1, 3))
colnames(coldata_veteran) = c('rte')
rownames(coldata_veteran) = colnames(counts_veteran)
coldata_veteran$type = substr(colnames(counts_veteran), 5, 6)
coldata_veteran$clean = 'C'
coldata_veteran$source = 'RTE_Project'
coldata_veteran$age = 'Adult'
coldata_veteran = coldata_veteran[, c('age', 'type', 'clean', 'source', 'rte')]


counts_tnVm = do.call("cbind", list(counts_cell2018, counts_dirty_nobulk, counts_veteran))
coldata_tnVm = do.call("rbind", list(coldata_cell2018, coldata_dirty_nobulk, coldata_veteran))
coldata_tnVm$age = as.factor(coldata_tnVm$age)
coldata_tnVm$type = as.factor(coldata_tnVm$type)
coldata_tnVm$clean = as.factor(coldata_tnVm$clean)
coldata_tnVm$source = as.factor(coldata_tnVm$source)
coldata_tnVm$rte = as.factor(coldata_tnVm$rte)

dds <- DESeqDataSetFromMatrix(countData = counts_tnVm,
                              colData = coldata_tnVm,
                              design = ~ age + type + clean + source + rte)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients


res.TN_MP <- results(dds, contrast=c("type","TN","MP"))
resSig.TN_MP <- subset(res.TN_MP, padj < padj_cutoff)
TN <- as.data.frame(subset(resSig.TN_MP, log2FoldChange > 0))
MP <- as.data.frame(subset(resSig.TN_MP, log2FoldChange < 0))


# make volcano plots (Fig. 1E-G)
library(ggplot2)
library(cowplot)

resSig.TN_MP$diffexpressed <- "Not significant"
resSig.TN_MP$diffexpressed[resSig.TN_MP$padj < 0.001 & resSig.TN_MP$log2FoldChange > 0] <- "TN"
resSig.TN_MP$diffexpressed[resSig.TN_MP$padj < 0.001 & resSig.TN_MP$log2FoldChange < 0] <- "MP"
resSig.TN_MP$diffexpressed <- factor(resSig.TN_MP$diffexpressed, levels=c("TN", "MP", "Not significant"))
resSig.MP_TN = resSig.TN_MP
resSig.MP_TN$log2FoldChange_rev = -resSig.MP_TN$log2FoldChange

resSig.neo_adult$diffexpressed <- "Not significant"
resSig.neo_adult$diffexpressed[resSig.neo_adult$padj < 0.001 & resSig.neo_adult$log2FoldChange > 0] <- "Neonates"
resSig.neo_adult$diffexpressed[resSig.neo_adult$padj < 0.001 & resSig.neo_adult$log2FoldChange < 0] <- "Adult"

resSig.clean_dirty$diffexpressed <- "Not significant"
resSig.clean_dirty$diffexpressed[resSig.clean_dirty$padj < 0.001 & resSig.clean_dirty$log2FoldChange > 0] <- "Dirty"
resSig.clean_dirty$diffexpressed[resSig.clean_dirty$padj < 0.001 & resSig.clean_dirty$log2FoldChange < 0] <- "Clean"


resSig.MP_TN_df = as.data.frame(resSig.MP_TN)
min(resSig.MP_TN_df$padj[resSig.MP_TN_df$padj != 0])
min(resSig.MP_TN_df$pvalue[resSig.MP_TN_df$pvalue != 0])
resSig.MP_TN_df$padj = ifelse(resSig.MP_TN_df$padj == 0, 1e-276, resSig.MP_TN_df$padj)
resSig.MP_TN_df$name = rownames(resSig.MP_TN_df)
resSig.MP_TN_df$label = ifelse(rownames(resSig.MP_TN_df) %in% c('Il10rb', 'Cxcr5', 'Ccr2', 'Tnfsh10', 'Ifng', 'Ifngr1', 'Ifngas1', 'Tox', 'Ifngr2', 'Mtss1', 'Cxcr3', 'Ccr9', 'Ccl5', 'Ccr5', 'Tbx21', 'Eomes', 'Runx2', 'Cd44', 'Il18rap', 'Il2rb', 'Gzmm', 
                                                                'Gzma', 'Gzmb'), rownames(resSig.MP_TN_df), '')
library(ggrepel)
volp.MP_TN_label = ggplot(data=as.data.frame(resSig.MP_TN_df), aes(x=log2FoldChange_rev, y=-log10(padj), color=diffexpressed, label = label)) + 
  geom_point() +
  geom_text_repel( box.padding = 2,
                   segment.color = 'grey50',
                   max.overlaps = 10000000)+ 
  theme_bw()+
  scale_color_manual(name="", 
                     values=c('#56C1FF', '#004D7F', "#7F7F7F"), 
                     labels=c('TN', 'VM', 'Not significant')) +
  #xlim(-10, 10) +
  #ylim(NA, 250 ) +
  ylab('-log10(adj p-val)') +
  ggtitle('TN vs. VM') +
  xlab('log2FC') +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position=c(0.25, 0.85),
        legend.background = element_rect(fill="transparent", size=0.5, linetype="solid"),
        legend.key = element_rect(fill = "transparent"))
volp.MP_TN_label

resSig.neo_adult_df = as.data.frame(resSig.neo_adult)
resSig.neo_adult_df$name = rownames(resSig.neo_adult_df)
resSig.neo_adult_df$label = ifelse(rownames(resSig.neo_adult_df) %in% c('Ifngr1', 'Il2rb', 'Il10rb', 'Cxcr5', 'Ccr2', 'Runx2', 'Tnfsf10', 'Ifngas1', 'Tox', 'Ifngr2', 'Mtss1', 'Ccr9', 'Ccr5', 'Cxcr3', 'Eomes'), rownames(resSig.neo_adult_df), '')
volp.neo_adult_label = ggplot(data=as.data.frame(resSig.neo_adult_df), aes(x=log2FoldChange, y=-log10(padj), color=diffexpressed, label = label)) + 
  geom_point() + 
  geom_text_repel( box.padding = 2,
                   segment.color = 'grey50',
                   max.overlaps = 10000000)+ 
  theme_bw()+
  scale_color_manual(name="", values=c('#88FA4E', '#017100', "#7F7F7F")) +
  #xlim(-10, 10) +
  #ylim(NA, 100) +
  ylab('-log10(adj p-val)') +
  xlab('log2FC') +
  ggtitle('Adult vs. Neo')+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position=c(0.25, 0.85),
        legend.background = element_rect(fill="transparent", size=0.5, linetype="solid"),
        legend.key = element_rect(fill = "transparent"))
volp.neo_adult_label


resSig.clean_dirty_df = as.data.frame(resSig.clean_dirty)
resSig.clean_dirty_df$name = rownames(resSig.clean_dirty_df)
resSig.clean_dirty_df$label = ifelse(rownames(resSig.clean_dirty_df) %in% c('Eomes', 'Runx2', 'Ifngas1', 'Ly6e', 'Tnfsf10', 'Ccr2', 'Ifit3', 'Ifit3b', 'Cd300e', 'Cd300c2', 'Ccl5', 'Ccr5', 'Gzma','Ly6a', 'Cx3cr1', 'Klrg1', 'Gzmk'), rownames(resSig.clean_dirty_df), '')
volp.clean_dirty_label = ggplot(data=as.data.frame(resSig.clean_dirty_df), aes(x=log2FoldChange, y=-log10(padj), color=diffexpressed, label=label)) + 
  geom_point() + 
  geom_text_repel( box.padding = 2,
                   segment.color = 'grey50',
                   max.overlaps = 10000000)+ 
  theme_bw()+
  scale_color_manual(name="", values=c('#FFFC66', '#FF9300', "#7F7F7F")) +
  #xlim(-10, 10) +
  #ylim(NA, 100) +
  ylab('-log10(adj p-val)') +
  xlab('log2FC') +
  ggtitle('Clean vs. Dirty')+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position=c(0.25, 0.85),
        legend.background = element_rect(fill="transparent", size=0.5, linetype="solid"),
        legend.key = element_rect(fill = "transparent"))
volp.clean_dirty_label

volp.MP_TN_label2 = volp.MP_TN_label +  
  xlim(-10, 10)  #  ylim(NA, 20) 
volp.neo_adult_label2 = volp.neo_adult_label +  
  xlim(-10, 10) +
  ylim(NA, 150) 
volp.clean_dirty_label2 = volp.clean_dirty_label +  
  xlim(-10, 10)  +
  ylim(NA, 150) 
volp.combined_label2 = plot_grid(volp.MP_TN_label2, volp.neo_adult_label2, volp.clean_dirty_label2, align="hv", ncol = 3)
ggsave(
  paste0('/workdir/hz543/projects/Inferelator/cd8_noOurlier_allsamples/remake_figures_for_paper/DE_volp/DE_genes_nobulk_volp_labels2.pdf'),
  plot = volp.combined_label2,
  device = "pdf",
  width = 16,
  height = 5,
  dpi = 300
)

