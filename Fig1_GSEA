# GSEA plot (Fig. 1H)

library(DESeq2)
padj_cutoff = 0.001

# clean dirty

counts_dirty = read.table("/path/to/Tabilas2022_RNAseq/counts_name.txt", header = T, row.names = 1)
counts_dirty = counts_dirty[,6:ncol(counts_dirty)]
coln_tmp = sub('.*CD', 'CD', colnames(counts_dirty))
coln_tmp = sub('_CKDL.*', '', coln_tmp)
targetFile = read.table('/path/to/Tabilas2022_RNAseq/targetFile.txt', header = T)
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
res.dirty_clean <- results(dds, name="clean_D_vs_C") # changed name here!
resSig.dirty_clean <- subset(res.dirty_clean, padj < padj_cutoff)# changed name here!
dirty <- rownames(subset(resSig.dirty_clean, log2FoldChange > 0))# changed name here!
clean <- rownames(subset(resSig.dirty_clean, log2FoldChange < 0))# changed name here!



# neo vs adult

counts_dirty = read.table("/path/to/Tabilas2022_RNAseq/counts_name.txt", header = T, row.names = 1)
counts_dirty = counts_dirty[,6:ncol(counts_dirty)]
coln_tmp = sub('.*CD', 'CD', colnames(counts_dirty))
coln_tmp = sub('_CKDL.*', '', coln_tmp)
targetFile = read.table('/path/to/Tabilas2022_RNAseq/targetFile.txt', header = T)
targetFile$cd = sub('.ReadsPerGene.out.tab.rawCounts', '', targetFile$files)
colnames(counts_dirty) = targetFile[match(coln_tmp, targetFile$cd),]$label
counts_dirty_nobulk = counts_dirty[, !grepl("B.", colnames(counts_dirty))]

coldata_dirty_nobulk = data.frame(c("Neo", "Neo", "Neo", "Neo", "Adult", "Adult", "Adult", "Adult", "Neo", "Neo", "Neo", "Neo", "Adult", "Adult", "Adult", "Adult", "Neo", "Neo", "Neo", "Neo", "Adult", "Adult", "Adult", "Adult"),
                                  gsub("VM", "MP", substr(colnames(counts_dirty_nobulk), 1, 2)), 
                                  gsub("\\..*","", gsub("^.*?\\.","", colnames(counts_dirty_nobulk))))
colnames(coldata_dirty_nobulk) = c('age', 'type', 'clean')
rownames(coldata_dirty_nobulk) = colnames(counts_dirty_nobulk)
coldata_dirty_nobulk$source = 'Dirty_Project'

counts_cell2018 = read.table("/path/to/Smith2018_RNAseq/counts_name.txt", header = T, row.names = 1)
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
neo <- rownames(subset(resSig.neo_adult, log2FoldChange > 0))
adult <- rownames(subset(resSig.neo_adult, log2FoldChange < 0))





# tn vs vm


counts_dirty = read.table("/path/to/Tabilas2022_RNAseq/counts_name.txt", header = T, row.names = 1)
counts_dirty = counts_dirty[,6:ncol(counts_dirty)]
coln_tmp = sub('.*CD', 'CD', colnames(counts_dirty))
coln_tmp = sub('_CKDL.*', '', coln_tmp)
targetFile = read.table('/path/to/Tabilas2022_RNAseq/targetFile.txt', header = T)
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

counts_cell2018 = read.table("/path/to/Smith2018_RNAseq/counts_name.txt", header = T, row.names = 1)
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

counts_veteran = read.table("/workdir/hz543/data_analysis/1812_Veteran_RNAseq/counts_name.txt", header = T, row.names = 1)
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
TN <- rownames(subset(resSig.TN_MP, log2FoldChange > 0))
MP <- rownames(subset(resSig.TN_MP, log2FoldChange < 0))


# GSEA

library(ggplot2)
library(cowplot)

resSig.TN_MP$diffexpressed <- "Not DE"
resSig.TN_MP$diffexpressed[resSig.TN_MP$padj < 0.001 & resSig.TN_MP$log2FoldChange > 0] <- "TN"
resSig.TN_MP$diffexpressed[resSig.TN_MP$padj < 0.001 & resSig.TN_MP$log2FoldChange < 0] <- "MP"
resSig.TN_MP$diffexpressed <- factor(resSig.TN_MP$diffexpressed, levels=c("TN", "MP", "Not DE"))
#tn = rownames(resSig.TN_MP[resSig.TN_MP$diffexpressed == 'TN', ])
#mp = rownames(resSig.TN_MP[resSig.TN_MP$diffexpressed == 'MP', ])

resSig.neo_adult$diffexpressed <- "Not DE"
resSig.neo_adult$diffexpressed[resSig.neo_adult$padj < 0.001 & resSig.neo_adult$log2FoldChange > 0] <- "Neonates"
resSig.neo_adult$diffexpressed[resSig.neo_adult$padj < 0.001 & resSig.neo_adult$log2FoldChange < 0] <- "Adult"
#neo = rownames(resSig.neo_adult[resSig.neo_adult$diffexpressed == 'Neonates', ])
#adult = rownames(resSig.neo_adult[resSig.neo_adult$diffexpressed == 'Adult', ])

resSig.dirty_clean$diffexpressed <- "Not DE"
resSig.dirty_clean$diffexpressed[resSig.dirty_clean$padj < 0.001 & resSig.dirty_clean$log2FoldChange > 0] <- "Dirty"
resSig.dirty_clean$diffexpressed[resSig.dirty_clean$padj < 0.001 & resSig.dirty_clean$log2FoldChange < 0] <- "Clean"
#clean = rownames(resSig.dirty_clean[resSig.dirty_clean$diffexpressed == 'Clean', ])
#dirty = rownames(resSig.dirty_clean[resSig.dirty_clean$diffexpressed == 'Dirty', ])



library(fgsea)
library(qusage)
immgen_genesets = read.gmt('/workdir/hz543/data_analysis/cd8_cell_type_clusters/cluster_list_geneSetNames.gmt')
names(immgen_genesets) = c("Initial cytokine or effector response", 
                           "Preparation for cell division", 
                           "Cell cycle and division", 
                           "Naive and late memory", 
                           "Early effector late memory", 
                           "Short term effector and memory", 
                           "Memory precursor", 
                           "Naive or late effector or memory", 
                           "Short term effector or memory", 
                           "Late effector or memory" )


resSig.MP_TN = resSig.TN_MP
resSig.MP_TN$log2FoldChange_rev = -resSig.MP_TN$log2FoldChange
ranks.MP_TN = resSig.MP_TN[order(resSig.MP_TN$log2FoldChange_rev),'log2FoldChange_rev']
names(ranks.MP_TN) = rownames(resSig.MP_TN[order(resSig.MP_TN$log2FoldChange_rev),])

ranks.neo_adult = resSig.neo_adult[order(resSig.neo_adult$log2FoldChange),'log2FoldChange']
names(ranks.neo_adult) = rownames(resSig.neo_adult[order(resSig.neo_adult$log2FoldChange),])

ranks.dirty_clean = resSig.dirty_clean[order(resSig.dirty_clean$log2FoldChange),'log2FoldChange']
names(ranks.dirty_clean) = rownames(resSig.dirty_clean[order(resSig.dirty_clean$log2FoldChange),])

set.seed(42)
fgseaRes.MP_TN <- fgsea(pathways = immgen_genesets, 
                        stats    = ranks.MP_TN,
                        eps      = 0.0,
                        minSize  = 15,
                        maxSize  = 500)
set.seed(42)
fgseaRes.neo_adult <- fgsea(pathways = immgen_genesets, 
                            stats    = ranks.neo_adult,
                            eps      = 0.0,
                            minSize  = 15,
                            maxSize  = 500)
set.seed(42)
fgseaRes.dirty_clean <- fgsea(pathways = immgen_genesets, 
                              stats    = ranks.dirty_clean,
                              eps      = 0.0,
                              minSize  = 15,
                              maxSize  = 500)
fgseaRes.MP_TN[fgseaRes.MP_TN$padj < 0.05,]
fgseaRes.neo_adult[fgseaRes.neo_adult$padj < 0.05,]
fgseaRes.dirty_clean[fgseaRes.dirty_clean$padj < 0.05,]


# make dotplots
fgseaRes.MP_TN$type = 'MP_TN'
fgseaRes.neo_adult$type = 'neo_adult'
fgseaRes.dirty_clean$type = 'dirty_clean'

forplot_MP_TN = fgseaRes.MP_TN[, c('pathway', 'NES', 'padj', 'type')]
forplot_MP_TN$Direction = ifelse(forplot_MP_TN$NES < 0, 'TN', 'VM')
forplot_MP_TN$Direction = ifelse(forplot_MP_TN$padj > 0.05, 'Not significant', forplot_MP_TN$Direction)
forplot_neo_adult = fgseaRes.neo_adult[, c('pathway', 'NES', 'padj', 'type')]
forplot_neo_adult$Direction = ifelse(forplot_neo_adult$NES < 0, 'Adult', 'Neo')
forplot_neo_adult$Direction = ifelse(forplot_neo_adult$padj > 0.05, 'Not significant', forplot_neo_adult$Direction)
forplot_dirty_clean = fgseaRes.dirty_clean[, c('pathway', 'NES', 'padj', 'type')]
forplot_dirty_clean$Direction = ifelse(forplot_dirty_clean$NES < 0, 'Clean', 'Dirty')
forplot_dirty_clean$Direction = ifelse(forplot_dirty_clean$padj > 0.05, 'Not significant', forplot_dirty_clean$Direction)
forplot = rbind(forplot_MP_TN, forplot_neo_adult, forplot_dirty_clean)

forplot$pathway = factor(forplot$pathway, levels = rev(c("Cell cycle and division",
                                                         "Short term effector and memory",
                                                         'Naive or late effector or memory',
                                                         'Short term effector or memory',
                                                         'Late effector or memory',
                                                         'Initial cytokine or effector response',
                                                         'Preparation for cell division',
                                                         'Naive and late memory',
                                                         'Early effector late memory',
                                                         'Memory precursor')))
forplot = as.data.frame(forplot)
forplot$Group = NA
for (i in 1:nrow(forplot)){
  if (forplot[i, 'pathway'] %in% c("Cell cycle and division",
                                   "Short term effector and memory",
                                   'Naive or late effector or memory',
                                   'Short term effector or memory',
                                   'Late effector or memory',
                                   'Initial cytokine or effector response')){
    forplot[i, 'Group'] = 'Effector\n\ngene sets'
  } else if (forplot[i, 'pathway'] %in% c('Preparation for cell division',
                                          'Naive and late memory',
                                          'Early effector late memory',
                                          'Memory precursor')){
    forplot[i, 'Group'] = 'Memory\n\ngene sets'
  }
}
library(ggh4x)
strip <- strip_themed(background_y = elem_list_rect(fill = 'gray95'))
ggplot(forplot, aes(x = type, y = pathway)) +
  geom_point(aes(size = ifelse(-log10(padj)==0, NA, abs(-log10(padj))), color = Direction)) +
  geom_tile(fill = 'transparent', color = 'gray85')+
  ylab('') +
  xlab('') +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.ontop = TRUE,
        panel.background = element_rect(fill = "transparent", color = NA),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        #ggh4x.axis.nestline = element_line(linetype = 1),
        strip.placement = "outside", 
        strip.text.y.left = element_text(angle = 0),
        strip.background = element_rect(fill = "gray95", colour = NA),
        panel.border = element_blank()) +
  scale_colour_manual(values = c('#004D7F', '#56C1FF','#017100', '#88FA4E', '#FF9300', '#FFFC66', 'gray80'),
                      limits = c('VM', 'TN', 'Neo', 'Adult', 'Dirty', 'Clean', 'Not significant'),
                      labels = c('VM', 'TN', 'Neo', 'Adult', 'Dirty', 'Clean', 'Not significant')) + 
  scale_x_discrete(limits = c('MP_TN', 'neo_adult', 'dirty_clean'),
                   labels = c('VM, TN','Neo, Adult', 'Dirty, Clean'), position = 'top') +
  guides(size = guide_legend(title='-log10(adj p-val)')) +
  facet_grid2(Group ~ ., scales = 'free', space = 'free', switch = "y", strip = strip)
ggsave(
  paste0('/path/to/dotplot_logFC.pdf'),
  plot = last_plot(),
  device = "pdf",
  width = 7,
  height = 5,
  dpi = 300
)

# in the above plot, there are some circles that does not have lines around them. 
# to make them look better, add other not significant results that are not in the result table
add_res_neo_adult = data.frame(pathway = c('Cell cycle and division', 'Short term effector or memory', 'Preparation for cell division'),
                               NES = 0,
                               padj = 1,
                               type = 'neo_adult',
                               Direction = 'Not significant',
                               Group = c('Effector\n\ngene sets', 'Effector\n\ngene sets', 'Memory\n\ngene sets'))
add_res_dirty_clean = data.frame(pathway = c('Naive or late effector or memory', 'Short term effector or memory', 
                                             'Preparation for cell division', 'Naive and late memory', 'Memory precursor'),
                               NES = 0,
                               padj = 1,
                               type = 'dirty_clean',
                               Direction = 'Not significant',
                               Group = c('Effector\n\ngene sets', 'Effector\n\ngene sets', 
                                         'Memory\n\ngene sets', 'Memory\n\ngene sets', 'Memory\n\ngene sets'))
forplot_add_res = rbind(forplot, add_res_neo_adult, add_res_dirty_clean)
ggplot(forplot_add_res, aes(x = type, y = pathway)) +
  geom_point(aes(size = ifelse(-log10(padj)==0, NA, abs(-log10(padj))), color = Direction)) +
  geom_tile(fill = 'transparent', color = 'gray85')+
  ylab('') +
  xlab('') +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.ontop = TRUE,
        panel.background = element_rect(fill = "transparent", color = NA),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        #ggh4x.axis.nestline = element_line(linetype = 1),
        strip.placement = "outside", 
        strip.text.y.left = element_text(angle = 0),
        strip.background = element_rect(fill = "gray95", colour = NA),
        panel.border = element_blank()) +
  scale_colour_manual(values = c('#004D7F', '#56C1FF','#017100', '#88FA4E', '#FF9300', '#FFFC66', 'gray80'),
                      limits = c('VM', 'TN', 'Neo', 'Adult', 'Dirty', 'Clean', 'Not significant'),
                      labels = c('VM', 'TN', 'Neo', 'Adult', 'Dirty', 'Clean', 'Not significant')) + 
  scale_x_discrete(limits = c('MP_TN', 'neo_adult', 'dirty_clean'),
                   labels = c('VM, TN','Neo, Adult', 'Dirty, Clean'), position = 'top') +
  guides(size = guide_legend(title='-log10(adj p-val)')) +
  facet_grid2(Group ~ ., scales = 'free', space = 'free', switch = "y", strip = strip)
# + scale_size_continuous(limits = c(3, 9)) # dot size legend
ggsave(
  paste0('/path/to/dotplot_logFC_add_res.pdf'),
  plot = last_plot(),
  device = "pdf",
  width = 7,
  height = 5,
  dpi = 300
)
