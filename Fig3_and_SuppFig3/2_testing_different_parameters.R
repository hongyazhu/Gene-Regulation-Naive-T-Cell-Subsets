#### Parameters tested:
# 1. model size
# 2. bias
# 3. TFA options

### Parameters are tested with two datasets (both are not used in network prediction process):
# 1. Eomes over-expression RNA-seq (Istaces et al., 2019. GSE124913)
# 2. Foxo targets in (Webb et al., 2016)



### test model size

# Eomes over-expression RNA-seq (GSE124913). not used in network prediction process

library(caret)
library(readxl)
pool_genes = read.table('deg_padj01.txt')$V1
de = read_excel("GSE124913_RNA-Seq_CD8SP_EomesTGvsWT_DifferentialAnalysisResults.xlsx", 
                sheet = "Sheet1")
de = data.frame(de)
de_noNA = de[complete.cases(de), ] #no use
de_noNA$FDR_numeric = as.numeric(de_noNA$FDR)
de_noNA = de_noNA[complete.cases(de_noNA), ]
de_noNA_pool_genes = de_noNA[de_noNA$Gene.Name %in% pool_genes,]
de_Eomes_pos = de_noNA_pool_genes[de_noNA_pool_genes$FDR_numeric < 0.001 & de_noNA_pool_genes$Log2FC.WT.vs.EomesTG < 0, 'Gene.Name']
de_Eomes_neg = de_noNA_pool_genes[de_noNA_pool_genes$FDR_numeric < 0.001 & de_noNA_pool_genes$Log2FC.WT.vs.EomesTG > 0, 'Gene.Name']

ms10_combine = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms10_bias50_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms15_combine = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms15_bias50_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms5_combine  = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms5_bias50_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms20_combine = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms20_bias50_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)


ms15_combine_Eomes_pos = ms15_combine[ms15_combine$TF == 'Eomes' & ms15_combine$SignedQuantile == 1, 'Target']
ms15_combine_Eomes_neg = ms15_combine[ms15_combine$TF == 'Eomes' & ms15_combine$SignedQuantile == -1, 'Target']

de_Eomes_pos_list = ifelse(pool_genes %in% de_Eomes_pos, 1, 0)
ms15_combine_Eomes_pos_list = ifelse(pool_genes %in% ms15_combine_Eomes_pos, 1, 0)

xtab <- table(ms15_combine_Eomes_pos_list, de_Eomes_pos_list)
cm_ms15_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15_combine$overall['Accuracy']
cm_ms15_combine$byClass['Precision']
cm_ms15_combine$byClass['Recall']
cm_ms15_combine$byClass['F1']
cm_ms15_combine$byClass['Balanced Accuracy']


ms10_combine_Eomes_pos = ms10_combine[ms10_combine$TF == 'Eomes' & ms10_combine$SignedQuantile == 1, 'Target']
ms10_combine_Eomes_neg = ms10_combine[ms10_combine$TF == 'Eomes' & ms10_combine$SignedQuantile == -1, 'Target']

de_Eomes_pos_list = ifelse(pool_genes %in% de_Eomes_pos, 1, 0)
ms10_combine_Eomes_pos_list = ifelse(pool_genes %in% ms10_combine_Eomes_pos, 1, 0)

xtab <- table(ms10_combine_Eomes_pos_list, de_Eomes_pos_list)
cm_ms10_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10_combine$overall['Accuracy']
cm_ms10_combine$byClass['Precision']
cm_ms10_combine$byClass['Recall']
cm_ms10_combine$byClass['F1']
cm_ms10_combine$byClass['Balanced Accuracy']


ms5_combine_Eomes_pos = ms5_combine[ms5_combine$TF == 'Eomes' & ms5_combine$SignedQuantile == 1, 'Target']
ms5_combine_Eomes_neg = ms5_combine[ms5_combine$TF == 'Eomes' & ms5_combine$SignedQuantile == -1, 'Target']

de_Eomes_pos_list = ifelse(pool_genes %in% de_Eomes_pos, 1, 0)
ms5_combine_Eomes_pos_list = ifelse(pool_genes %in% ms5_combine_Eomes_pos, 1, 0)

xtab <- table(ms5_combine_Eomes_pos_list, de_Eomes_pos_list)
cm_ms5_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5_combine$overall['Accuracy']
cm_ms5_combine$byClass['Precision']
cm_ms5_combine$byClass['Recall']
cm_ms5_combine$byClass['F1']
cm_ms5_combine$byClass['Balanced Accuracy']


ms20_combine_Eomes_pos = ms20_combine[ms20_combine$TF == 'Eomes' & ms20_combine$SignedQuantile == 1, 'Target']
ms20_combine_Eomes_neg = ms20_combine[ms20_combine$TF == 'Eomes' & ms20_combine$SignedQuantile == -1, 'Target']

de_Eomes_pos_list = ifelse(pool_genes %in% de_Eomes_pos, 1, 0)
ms20_combine_Eomes_pos_list = ifelse(pool_genes %in% ms20_combine_Eomes_pos, 1, 0)

xtab <- table(ms20_combine_Eomes_pos_list, de_Eomes_pos_list)
cm_ms20_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20_combine$overall['Accuracy']
cm_ms20_combine$byClass['Precision']
cm_ms20_combine$byClass['Recall']
cm_ms20_combine$byClass['F1']
cm_ms20_combine$byClass['Balanced Accuracy']


# make plots
forplot = data.frame(matrix(ncol = 3, nrow = 4))
colnames(forplot) = c('Accuracy', 'Balanced Accuracy', 'Precision')
rownames(forplot) = c('ms5_combine', 'ms10_combine', 'ms15_combine', 'ms20_combine')
forplot['ms10_combine', 'Accuracy'] = cm_ms10_combine$overall['Accuracy']
forplot['ms15_combine', 'Accuracy'] = cm_ms15_combine$overall['Accuracy']
forplot['ms20_combine', 'Accuracy'] = cm_ms20_combine$overall['Accuracy']
forplot['ms5_combine', 'Accuracy'] = cm_ms5_combine$overall['Accuracy']

forplot['ms10_combine', 'Balanced Accuracy'] = cm_ms10_combine$byClass['Balanced Accuracy']
forplot['ms15_combine', 'Balanced Accuracy'] = cm_ms15_combine$byClass['Balanced Accuracy']
forplot['ms20_combine', 'Balanced Accuracy'] = cm_ms20_combine$byClass['Balanced Accuracy']
forplot['ms5_combine', 'Balanced Accuracy'] = cm_ms5_combine$byClass['Balanced Accuracy']

forplot['ms10_combine', 'Precision'] = cm_ms10_combine$byClass['Precision']
forplot['ms15_combine', 'Precision'] = cm_ms15_combine$byClass['Precision']
forplot['ms20_combine', 'Precision'] = cm_ms20_combine$byClass['Precision']
forplot['ms5_combine', 'Precision'] = cm_ms5_combine$byClass['Precision']


forplot$ms = c('Model size: 5', 'Model size: 10', 'Model size: 15', 'Model size: 20')
forplot$ms = factor(forplot$ms, levels = c('Model size: 5', 'Model size: 10', 'Model size: 15', 'Model size: 20'))
forplot$TFAOpt = c('comb', 'comb', 'comb', 'comb')
forplot$TFAOpt = factor(forplot$TFAOpt, levels = c('comb'))

forplot_melt = reshape2::melt(forplot)
library(ggplot2)
library(RColorBrewer)
# brewer.pal(3, 'Set3')
eomes_pos_plot = ggplot(data=forplot_melt, aes(x=variable, y=value, fill=ms)) +
  geom_bar(stat="identity", position=position_dodge(), color="black", width = 0.7)+
  #facet_wrap(~ms) + 
  #scale_fill_manual(limits = c('', 'TFmRNA', 'comb'),
  #                  labels = c('TFA', 'TFmRNA', 'combine'),
  #                  name = 'TFA options',
  #                  values = c("#8DD3C7", "#FFFFB3", "#BEBADA"))+
  theme_bw(base_size = 14) +
  xlab('') +
  ylab('Value')+ 
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle('Eomes positive targets')+                                              # Change color brewer palette
  scale_fill_brewer(palette = "RdYlBu")



# negative targets

ms15_combine_Eomes_neg = ms15_combine[ms15_combine$TF == 'Eomes' & ms15_combine$SignedQuantile == -1, 'Target']

de_Eomes_neg_list = ifelse(pool_genes %in% de_Eomes_neg, 1, 0)
ms15_combine_Eomes_neg_list = ifelse(pool_genes %in% ms15_combine_Eomes_neg, 1, 0)

xtab <- table(ms15_combine_Eomes_neg_list, de_Eomes_neg_list)
cm_ms15_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15_combine$overall['Accuracy']
cm_ms15_combine$byClass['Precision']
cm_ms15_combine$byClass['Recall']
cm_ms15_combine$byClass['F1']
cm_ms15_combine$byClass['Balanced Accuracy']


ms10_combine_Eomes_neg = ms10_combine[ms10_combine$TF == 'Eomes' & ms10_combine$SignedQuantile == -1, 'Target']

de_Eomes_neg_list = ifelse(pool_genes %in% de_Eomes_neg, 1, 0)
ms10_combine_Eomes_neg_list = ifelse(pool_genes %in% ms10_combine_Eomes_neg, 1, 0)

xtab <- table(ms10_combine_Eomes_neg_list, de_Eomes_neg_list)
cm_ms10_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10_combine$overall['Accuracy']
cm_ms10_combine$byClass['Precision']
cm_ms10_combine$byClass['Recall']
cm_ms10_combine$byClass['F1']
cm_ms10_combine$byClass['Balanced Accuracy']


ms5_combine_Eomes_neg = ms5_combine[ms5_combine$TF == 'Eomes' & ms5_combine$SignedQuantile == -1, 'Target']

de_Eomes_neg_list = ifelse(pool_genes %in% de_Eomes_neg, 1, 0)
ms5_combine_Eomes_neg_list = ifelse(pool_genes %in% ms5_combine_Eomes_neg, 1, 0)

xtab <- table(ms5_combine_Eomes_neg_list, de_Eomes_neg_list)
cm_ms5_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5_combine$overall['Accuracy']
cm_ms5_combine$byClass['Precision']
cm_ms5_combine$byClass['Recall']
cm_ms5_combine$byClass['F1']
cm_ms5_combine$byClass['Balanced Accuracy']


ms20_combine_Eomes_neg = ms20_combine[ms20_combine$TF == 'Eomes' & ms20_combine$SignedQuantile == -1, 'Target']

de_Eomes_neg_list = ifelse(pool_genes %in% de_Eomes_neg, 1, 0)
ms20_combine_Eomes_neg_list = ifelse(pool_genes %in% ms20_combine_Eomes_neg, 1, 0)

xtab <- table(ms20_combine_Eomes_neg_list, de_Eomes_neg_list)
cm_ms20_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20_combine$overall['Accuracy']
cm_ms20_combine$byClass['Precision']
cm_ms20_combine$byClass['Recall']
cm_ms20_combine$byClass['F1']
cm_ms20_combine$byClass['Balanced Accuracy']

# library(Metrics)
# Metrics::recall(ms15_TFmRNA_Eomes_neg_list, de_Eomes_neg_list)
# Metrics::precision(ms15_TFmRNA_Eomes_neg_list, de_Eomes_neg_list)
# Metrics::f1(ms15_TFmRNA_Eomes_neg_list, de_Eomes_neg_list) # wrong somehow


# make plots
forplotneg = data.frame(matrix(ncol = 3, nrow = 4))
colnames(forplotneg) = c('Accuracy', 'Balanced Accuracy', 'Precision')
rownames(forplotneg) = c('ms5_combine', 'ms10_combine', 'ms15_combine', 'ms20_combine')
forplotneg['ms10_combine', 'Accuracy'] = cm_ms10_combine$overall['Accuracy']
forplotneg['ms15_combine', 'Accuracy'] = cm_ms15_combine$overall['Accuracy']
forplotneg['ms20_combine', 'Accuracy'] = cm_ms20_combine$overall['Accuracy']
forplotneg['ms5_combine', 'Accuracy'] = cm_ms5_combine$overall['Accuracy']

forplotneg['ms10_combine', 'Balanced Accuracy'] = cm_ms10_combine$byClass['Balanced Accuracy']
forplotneg['ms15_combine', 'Balanced Accuracy'] = cm_ms15_combine$byClass['Balanced Accuracy']
forplotneg['ms20_combine', 'Balanced Accuracy'] = cm_ms20_combine$byClass['Balanced Accuracy']
forplotneg['ms5_combine', 'Balanced Accuracy'] = cm_ms5_combine$byClass['Balanced Accuracy']

forplotneg['ms10_combine', 'Precision'] = cm_ms10_combine$byClass['Precision']
forplotneg['ms15_combine', 'Precision'] = cm_ms15_combine$byClass['Precision']
forplotneg['ms20_combine', 'Precision'] = cm_ms20_combine$byClass['Precision']
forplotneg['ms5_combine', 'Precision'] = cm_ms5_combine$byClass['Precision']


forplotneg$ms = c('Model size: 5', 'Model size: 10', 'Model size: 15', 'Model size: 20')
forplotneg$ms = factor(forplotneg$ms, levels = c('Model size: 5', 'Model size: 10', 'Model size: 15', 'Model size: 20'))
forplotneg$TFAOpt = c('comb', 'comb', 'comb', 'comb')
forplotneg$TFAOpt = factor(forplotneg$TFAOpt, levels = c('comb'))

forplotneg_melt = reshape2::melt(forplotneg)
library(ggplot2)
# library(RColorBrewer)
# brewer.pal(3, 'Set3')
eomes_neg_plot = ggplot(data=forplotneg_melt, aes(x=variable, y=value, fill=ms)) +
  geom_bar(stat="identity", position=position_dodge(), color="black", width = 0.7)+
  #facet_wrap(~ms) + 
  #scale_fill_manual(limits = c('', 'TFmRNA', 'comb'),
  #                  labels = c('TFA', 'TFmRNA', 'combine'),
  #                  name = 'TFA options',
  #                  values = c("#8DD3C7", "#FFFFB3", "#BEBADA"))+
  theme_bw(base_size = 14) +
  xlab('') +
  ylab('Value')+ 
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle('Eomes negative targets')+                                              # Change color brewer palette
  scale_fill_brewer(palette = "RdYlBu")

library(cowplot)
eomes.combined = plot_grid(eomes_pos_plot, eomes_neg_plot, align="hv", ncol = 2)
ggsave(
  paste0('compare_ms_Eomes.pdf'),
  plot = eomes.combined,
  device = "pdf",
  width = 9,
  height = 4,
  dpi = 300
)

# Foxo targets in (Webb et al., 2016: Characterization of the direct targets of FOXO transcription factors throughout evolution)

pool_genes = read.table('deg_padj01.txt')$V1

library(readxl)

targets = read_excel("foxo1/ACEL-15-673-s002.xlsx", 
                     sheet = "Sheet1")
targets_cd8 = unique(targets$`CD8 T cells`)
targets_allcelltype = unique(targets$Core)

targets_cd8 = targets_cd8[!is.na(targets_cd8)]
targets_allcelltype = targets_allcelltype[!is.na(targets_allcelltype)]

targets_cd8_pool_genes = targets_cd8[targets_cd8 %in% pool_genes]
targets_allcelltype_pool_genes = targets_allcelltype[targets_allcelltype %in% pool_genes]


ms10_combine = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms10_bias50_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms15_combine = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms15_bias50_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms5_combine  = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms5_bias50_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms20_combine = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms20_bias50_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)

# both

ms15_combine_Foxo1 = ms15_combine[ms15_combine$TF == 'Foxo1', 'Target']

#table(ms15_combine_Foxo1 %in% pool_genes)
#table(targets_cd8_pool_genes %in% pool_genes)

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms15_combine_Foxo1_list = ifelse(pool_genes %in% ms15_combine_Foxo1, 1, 0)

xtab <- table(ms15_combine_Foxo1_list, targets_cd8_pool_genes_list)
#cm <- caret::confusionMatrix(xtab)
cm_ms15_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15_combine$overall['Accuracy']
cm_ms15_combine$byClass['Precision']
cm_ms15_combine$byClass['Recall']
cm_ms15_combine$byClass['F1']
cm_ms15_combine$byClass['Balanced Accuracy']


ms10_combine_Foxo1 = ms10_combine[ms10_combine$TF == 'Foxo1', 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms10_combine_Foxo1_list = ifelse(pool_genes %in% ms10_combine_Foxo1, 1, 0)

xtab <- table(ms10_combine_Foxo1_list, targets_cd8_pool_genes_list)
cm_ms10_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10_combine$overall['Accuracy']
cm_ms10_combine$byClass['Precision']
cm_ms10_combine$byClass['Recall']
cm_ms10_combine$byClass['F1']
cm_ms10_combine$byClass['Balanced Accuracy']


ms5_combine_Foxo1 = ms5_combine[ms5_combine$TF == 'Foxo1', 'Target']

#table(ms5_combine_Foxo1 %in% pool_genes)
#table(targets_cd8_pool_genes %in% pool_genes)

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms5_combine_Foxo1_list = ifelse(pool_genes %in% ms5_combine_Foxo1, 1, 0)

xtab <- table(ms5_combine_Foxo1_list, targets_cd8_pool_genes_list)
#cm <- caret::confusionMatrix(xtab)
cm_ms5_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5_combine$overall['Accuracy']
cm_ms5_combine$byClass['Precision']
cm_ms5_combine$byClass['Recall']
cm_ms5_combine$byClass['F1']
cm_ms5_combine$byClass['Balanced Accuracy']


ms20_combine_Foxo1 = ms20_combine[ms20_combine$TF == 'Foxo1', 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms20_combine_Foxo1_list = ifelse(pool_genes %in% ms20_combine_Foxo1, 1, 0)

xtab <- table(ms20_combine_Foxo1_list, targets_cd8_pool_genes_list)
cm_ms20_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20_combine$overall['Accuracy']
cm_ms20_combine$byClass['Precision']
cm_ms20_combine$byClass['Recall']
cm_ms20_combine$byClass['F1']
cm_ms20_combine$byClass['Balanced Accuracy']


# make plots
forplot = data.frame(matrix(ncol = 3, nrow = 4))
colnames(forplot) = c('Accuracy', 'Balanced Accuracy', 'Precision')
rownames(forplot) = c('ms5_combine', 'ms10_combine', 'ms15_combine', 'ms20_combine')
forplot['ms10_combine', 'Accuracy'] = cm_ms10_combine$overall['Accuracy']
forplot['ms15_combine', 'Accuracy'] = cm_ms15_combine$overall['Accuracy']
forplot['ms20_combine', 'Accuracy'] = cm_ms20_combine$overall['Accuracy']
forplot['ms5_combine', 'Accuracy'] = cm_ms5_combine$overall['Accuracy']

forplot['ms10_combine', 'Balanced Accuracy'] = cm_ms10_combine$byClass['Balanced Accuracy']
forplot['ms15_combine', 'Balanced Accuracy'] = cm_ms15_combine$byClass['Balanced Accuracy']
forplot['ms20_combine', 'Balanced Accuracy'] = cm_ms20_combine$byClass['Balanced Accuracy']
forplot['ms5_combine', 'Balanced Accuracy'] = cm_ms5_combine$byClass['Balanced Accuracy']

forplot['ms10_combine', 'Precision'] = cm_ms10_combine$byClass['Precision']
forplot['ms15_combine', 'Precision'] = cm_ms15_combine$byClass['Precision']
forplot['ms20_combine', 'Precision'] = cm_ms20_combine$byClass['Precision']
forplot['ms5_combine', 'Precision'] = cm_ms5_combine$byClass['Precision']

forplot$ms = c('Model size: 5', 'Model size: 10', 'Model size: 15', 'Model size: 20')
forplot$ms = factor(forplot$ms, levels = c('Model size: 5', 'Model size: 10', 'Model size: 15', 'Model size: 20'))
forplot$TFAOpt = c('combine', 'combine', 'combine', 'combine')
forplot$TFAOpt = factor(forplot$TFAOpt, levels = c('combine'))

forplot_melt = reshape2::melt(forplot)
library(ggplot2)
# library(RColorBrewer)
# brewer.pal(3, 'Set3')
foxo1_plot = ggplot(data=forplot_melt, aes(x=variable, y=value, fill=ms)) +
  geom_bar(stat="identity", position=position_dodge(), color="black", width = 0.7)+
  #facet_wrap(~ms) + 
  #scale_fill_manual(limits = c('TFmRNA', '', 'combine'),
  #                  labels = c('TFmRNA', 'TFA', 'combine'),
  #                  name = 'TFA options',
  #                  values = c("#8DD3C7", "#FFFFB3", "#BEBADA"))+
  theme_bw(base_size = 14) +
  xlab('') +
  ylab('Value')+ 
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle('Foxo1 targets')+                                              # Change color brewer palette
  scale_fill_brewer(palette = "RdYlBu")

ggsave(
  paste0('compare_ms_Foxo1.pdf'),
  plot = foxo1_plot,
  device = "pdf",
  width = 4,
  height = 4,
  dpi = 300
)



# positive targets
ms15_combine_Foxo1_pos = ms15_combine[ms15_combine$TF == 'Foxo1' & ms15_combine$SignedQuantile == 1, 'Target']
ms15_combine_Foxo1_neg = ms15_combine[ms15_combine$TF == 'Foxo1' & ms15_combine$SignedQuantile == -1, 'Target']

#table(ms15_combine_Foxo1_pos %in% pool_genes)
#table(ms15_combine_Foxo1_neg %in% pool_genes)
#table(targets_cd8_pool_genes %in% pool_genes)
#table(targets_cd8_pool_genes %in% pool_genes)

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms15_combine_Foxo1_pos_list = ifelse(pool_genes %in% ms15_combine_Foxo1_pos, 1, 0)

xtab <- table(ms15_combine_Foxo1_pos_list, targets_cd8_pool_genes_list)
#cm <- caret::confusionMatrix(xtab)
cm_ms15_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15_combine$overall['Accuracy']
cm_ms15_combine$byClass['Precision']
cm_ms15_combine$byClass['Recall']
cm_ms15_combine$byClass['F1']
cm_ms15_combine$byClass['Balanced Accuracy']


ms10_combine_Foxo1_pos = ms10_combine[ms10_combine$TF == 'Foxo1' & ms10_combine$SignedQuantile == 1, 'Target']
ms10_combine_Foxo1_neg = ms10_combine[ms10_combine$TF == 'Foxo1' & ms10_combine$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms10_combine_Foxo1_pos_list = ifelse(pool_genes %in% ms10_combine_Foxo1_pos, 1, 0)

xtab <- table(ms10_combine_Foxo1_pos_list, targets_cd8_pool_genes_list)
cm_ms10_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10_combine$overall['Accuracy']
cm_ms10_combine$byClass['Precision']
cm_ms10_combine$byClass['Recall']
cm_ms10_combine$byClass['F1']
cm_ms10_combine$byClass['Balanced Accuracy']


ms5_combine_Foxo1_pos = ms5_combine[ms5_combine$TF == 'Foxo1' & ms5_combine$SignedQuantile == 1, 'Target']
ms5_combine_Foxo1_neg = ms5_combine[ms5_combine$TF == 'Foxo1' & ms5_combine$SignedQuantile == -1, 'Target']

#table(ms5_combine_Foxo1_pos %in% pool_genes)
#table(ms5_combine_Foxo1_neg %in% pool_genes)
#table(targets_cd8_pool_genes %in% pool_genes)
#table(targets_cd8_pool_genes %in% pool_genes)

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms5_combine_Foxo1_pos_list = ifelse(pool_genes %in% ms5_combine_Foxo1_pos, 1, 0)

xtab <- table(ms5_combine_Foxo1_pos_list, targets_cd8_pool_genes_list)
#cm <- caret::confusionMatrix(xtab)
cm_ms5_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5_combine$overall['Accuracy']
cm_ms5_combine$byClass['Precision']
cm_ms5_combine$byClass['Recall']
cm_ms5_combine$byClass['F1']
cm_ms5_combine$byClass['Balanced Accuracy']


ms20_combine_Foxo1_pos = ms20_combine[ms20_combine$TF == 'Foxo1' & ms20_combine$SignedQuantile == 1, 'Target']
ms20_combine_Foxo1_neg = ms20_combine[ms20_combine$TF == 'Foxo1' & ms20_combine$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms20_combine_Foxo1_pos_list = ifelse(pool_genes %in% ms20_combine_Foxo1_pos, 1, 0)

xtab <- table(ms20_combine_Foxo1_pos_list, targets_cd8_pool_genes_list)
cm_ms20_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20_combine$overall['Accuracy']
cm_ms20_combine$byClass['Precision']
cm_ms20_combine$byClass['Recall']
cm_ms20_combine$byClass['F1']
cm_ms20_combine$byClass['Balanced Accuracy']


# make plots
forplot = data.frame(matrix(ncol = 3, nrow = 4))
colnames(forplot) = c('Accuracy', 'Balanced Accuracy', 'Precision')
rownames(forplot) = c('ms5_combine', 'ms10_combine', 'ms15_combine', 'ms20_combine')
forplot['ms10_combine', 'Accuracy'] = cm_ms10_combine$overall['Accuracy']
forplot['ms15_combine', 'Accuracy'] = cm_ms15_combine$overall['Accuracy']
forplot['ms20_combine', 'Accuracy'] = cm_ms20_combine$overall['Accuracy']
forplot['ms5_combine', 'Accuracy'] = cm_ms5_combine$overall['Accuracy']

forplot['ms10_combine', 'Balanced Accuracy'] = cm_ms10_combine$byClass['Balanced Accuracy']
forplot['ms15_combine', 'Balanced Accuracy'] = cm_ms15_combine$byClass['Balanced Accuracy']
forplot['ms20_combine', 'Balanced Accuracy'] = cm_ms20_combine$byClass['Balanced Accuracy']
forplot['ms5_combine', 'Balanced Accuracy'] = cm_ms5_combine$byClass['Balanced Accuracy']

forplot['ms10_combine', 'Precision'] = cm_ms10_combine$byClass['Precision']
forplot['ms15_combine', 'Precision'] = cm_ms15_combine$byClass['Precision']
forplot['ms20_combine', 'Precision'] = cm_ms20_combine$byClass['Precision']
forplot['ms5_combine', 'Precision'] = cm_ms5_combine$byClass['Precision']

forplot$ms = c('Model size: 5', 'Model size: 10', 'Model size: 15', 'Model size: 20')
forplot$ms = factor(forplot$ms, levels = c('Model size: 5', 'Model size: 10', 'Model size: 15', 'Model size: 20'))
forplot$TFAOpt = c('combine', 'combine', 'combine', 'combine')
forplot$TFAOpt = factor(forplot$TFAOpt, levels = c('combine'))

forplot_melt = reshape2::melt(forplot)
library(ggplot2)
# library(RColorBrewer)
# brewer.pal(3, 'Set3')
foxo1_pos_plot = ggplot(data=forplot_melt, aes(x=variable, y=value, fill=ms)) +
  geom_bar(stat="identity", position=position_dodge(), color="black", width = 0.7)+
  #facet_wrap(~ms) + 
  #scale_fill_manual(limits = c('TFmRNA', '', 'combine'),
  #                  labels = c('TFmRNA', 'TFA', 'combine'),
  #                  name = 'TFA options',
  #                  values = c("#8DD3C7", "#FFFFB3", "#BEBADA"))+
  theme_bw(base_size = 14) +
  xlab('') +
  ylab('Value')+ 
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle('Foxo1 positive targets')+                                              # Change color brewer palette
  scale_fill_brewer(palette = "RdYlBu")




# negative targets

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms15_combine_Foxo1_neg_list = ifelse(pool_genes %in% ms15_combine_Foxo1_neg, 1, 0)

xtab <- table(ms15_combine_Foxo1_neg_list, targets_cd8_pool_genes_list)
#cm <- caret::confusionMatrix(xtab)
cm_ms15_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15_combine$overall['Accuracy']
cm_ms15_combine$byClass['Precision']
cm_ms15_combine$byClass['Recall']
cm_ms15_combine$byClass['F1']
cm_ms15_combine$byClass['Balanced Accuracy']


ms10_combine_Foxo1_neg = ms10_combine[ms10_combine$TF == 'Foxo1' & ms10_combine$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms10_combine_Foxo1_neg_list = ifelse(pool_genes %in% ms10_combine_Foxo1_neg, 1, 0)

xtab <- table(ms10_combine_Foxo1_neg_list, targets_cd8_pool_genes_list)
cm_ms10_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10_combine$overall['Accuracy']
cm_ms10_combine$byClass['Precision']
cm_ms10_combine$byClass['Recall']
cm_ms10_combine$byClass['F1']
cm_ms10_combine$byClass['Balanced Accuracy']


targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms20_combine_Foxo1_neg_list = ifelse(pool_genes %in% ms20_combine_Foxo1_neg, 1, 0)

xtab <- table(ms20_combine_Foxo1_neg_list, targets_cd8_pool_genes_list)
#cm <- caret::confusionMatrix(xtab)
cm_ms20_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20_combine$overall['Accuracy']
cm_ms20_combine$byClass['Precision']
cm_ms20_combine$byClass['Recall']
cm_ms20_combine$byClass['F1']
cm_ms20_combine$byClass['Balanced Accuracy']


ms5_combine_Foxo1_neg = ms5_combine[ms5_combine$TF == 'Foxo1' & ms5_combine$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms5_combine_Foxo1_neg_list = ifelse(pool_genes %in% ms5_combine_Foxo1_neg, 1, 0)

xtab <- table(ms5_combine_Foxo1_neg_list, targets_cd8_pool_genes_list)
cm_ms5_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5_combine$overall['Accuracy']
cm_ms5_combine$byClass['Precision']
cm_ms5_combine$byClass['Recall']
cm_ms5_combine$byClass['F1']
cm_ms5_combine$byClass['Balanced Accuracy']



# make plots
forplotneg = data.frame(matrix(ncol = 3, nrow = 4))
colnames(forplotneg) = c('Accuracy', 'Balanced Accuracy', 'Precision')
rownames(forplotneg) = c('ms5_combine', 'ms10_combine', 'ms15_combine', 'ms20_combine')
forplotneg['ms10_combine', 'Accuracy'] = cm_ms10_combine$overall['Accuracy']
forplotneg['ms15_combine', 'Accuracy'] = cm_ms15_combine$overall['Accuracy']
forplotneg['ms5_combine', 'Accuracy'] = cm_ms5_combine$overall['Accuracy']
forplotneg['ms20_combine', 'Accuracy'] = cm_ms20_combine$overall['Accuracy']

forplotneg['ms10_combine', 'Balanced Accuracy'] = cm_ms10_combine$byClass['Balanced Accuracy']
forplotneg['ms15_combine', 'Balanced Accuracy'] = cm_ms15_combine$byClass['Balanced Accuracy']
forplotneg['ms5_combine', 'Balanced Accuracy'] = cm_ms5_combine$byClass['Balanced Accuracy']
forplotneg['ms20_combine', 'Balanced Accuracy'] = cm_ms20_combine$byClass['Balanced Accuracy']

forplotneg['ms10_combine', 'Precision'] = cm_ms10_combine$byClass['Precision']
forplotneg['ms15_combine', 'Precision'] = cm_ms15_combine$byClass['Precision']
forplotneg['ms5_combine', 'Precision'] = cm_ms5_combine$byClass['Precision']
forplotneg['ms20_combine', 'Precision'] = cm_ms20_combine$byClass['Precision']

forplotneg$ms = c('Model size: 5', 'Model size: 10', 'Model size: 15', 'Model size: 20')
forplotneg$ms = factor(forplotneg$ms, levels = c('Model size: 5', 'Model size: 10', 'Model size: 15', 'Model size: 20'))
forplotneg$TFAOpt = c('combine', 'combine', 'combine', 'combine')
forplotneg$TFAOpt = factor(forplotneg$TFAOpt, levels = c('combine'))

forplotneg_melt = reshape2::melt(forplotneg)
library(ggplot2)
# library(RColorBrewer)
# brewer.pal(3, 'Set3')
foxo1_neg_plot = ggplot(data=forplotneg_melt, aes(x=variable, y=value, fill=ms)) +
  geom_bar(stat="identity", position=position_dodge(), color="black", width = 0.7)+
  #facet_wrap(~ms) + 
  #scale_fill_manual(limits = c('TFmRNA', '', 'combine'),
  #                  labels = c('TFmRNA', 'TFA', 'combine'),
  #                  name = 'TFA options',
  #                  values = c("#8DD3C7", "#FFFFB3", "#BEBADA"))+
  theme_bw(base_size = 14) +
  xlab('') +
  ylab('Value')+ 
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle('Foxo1 negative targets')+                                              # Change color brewer palette
  scale_fill_brewer(palette = "RdYlBu")

library(cowplot)
foxo1.combined = plot_grid(foxo1_pos_plot, foxo1_neg_plot, align="hv", ncol = 2)
ggsave(
  paste0('compare_ms_Foxo1_posneg.pdf'),
  plot = foxo1.combined,
  device = "pdf",
  width = 9,
  height = 4,
  dpi = 300
)


### testing bias

# Eomes

library(caret)
library(readxl)
pool_genes = read.table('deg_padj01.txt')$V1
de = read_excel("eomes/GSE124913_RNA-Seq_CD8SP_EomesTGvsWT_DifferentialAnalysisResults.xlsx", 
                sheet = "Sheet1")
de = data.frame(de)
de_noNA = de[complete.cases(de), ] #no use
de_noNA$FDR_numeric = as.numeric(de_noNA$FDR)
de_noNA = de_noNA[complete.cases(de_noNA), ]
de_noNA_pool_genes = de_noNA[de_noNA$Gene.Name %in% pool_genes,]
de_Eomes_pos = de_noNA_pool_genes[de_noNA_pool_genes$FDR_numeric < 0.001 & de_noNA_pool_genes$Log2FC.WT.vs.EomesTG < 0, 'Gene.Name']
de_Eomes_neg = de_noNA_pool_genes[de_noNA_pool_genes$FDR_numeric < 0.001 & de_noNA_pool_genes$Log2FC.WT.vs.EomesTG > 0, 'Gene.Name']


ms10_bias50  = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms10_bias50_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms10_bias5  = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms10_bias5_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms10_bias100 = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms10_bias100_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms15_bias50  = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms15_bias50_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms15_bias5  = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms15_bias5_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms15_bias100 = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms15_bias100_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)

ms5_bias50  = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms5_bias50_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms5_bias5  = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms5_bias5_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms5_bias100 = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms5_bias100_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms20_bias50  = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms20_bias50_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms20_bias5  = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms20_bias5_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms20_bias100 = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms20_bias100_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)


ms15_bias100_Eomes_pos = ms15_bias100[ms15_bias100$TF == 'Eomes' & ms15_bias100$SignedQuantile == 1, 'Target']
ms15_bias100_Eomes_neg = ms15_bias100[ms15_bias100$TF == 'Eomes' & ms15_bias100$SignedQuantile == -1, 'Target']

#table(ms15_bias100_Eomes_pos %in% pool_genes)
#table(ms15_bias100_Eomes_neg %in% pool_genes)
#table(de_Eomes_pos %in% pool_genes)
#table(de_Eomes_neg %in% pool_genes)

de_Eomes_pos_list = ifelse(pool_genes %in% de_Eomes_pos, 1, 0)
ms15_bias100_Eomes_pos_list = ifelse(pool_genes %in% ms15_bias100_Eomes_pos, 1, 0)

xtab <- table(ms15_bias100_Eomes_pos_list, de_Eomes_pos_list)
#cm <- caret::confusionMatrix(xtab)
cm_ms15_bias100 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15_bias100$overall['Accuracy']
cm_ms15_bias100$byClass['Precision']
cm_ms15_bias100$byClass['Recall']
cm_ms15_bias100$byClass['F1']
cm_ms15_bias100$byClass['Balanced Accuracy']


ms15_bias5_Eomes_pos = ms15_bias5[ms15_bias5$TF == 'Eomes' & ms15_bias5$SignedQuantile == 1, 'Target']
ms15_bias5_Eomes_neg = ms15_bias5[ms15_bias5$TF == 'Eomes' & ms15_bias5$SignedQuantile == -1, 'Target']

de_Eomes_pos_list = ifelse(pool_genes %in% de_Eomes_pos, 1, 0)
ms15_bias5_Eomes_pos_list = ifelse(pool_genes %in% ms15_bias5_Eomes_pos, 1, 0)

xtab <- table(ms15_bias5_Eomes_pos_list, de_Eomes_pos_list)
cm_ms15_bias5 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15_bias5$overall['Accuracy']
cm_ms15_bias5$byClass['Precision']
cm_ms15_bias5$byClass['Recall']
cm_ms15_bias5$byClass['F1']
cm_ms15_bias5$byClass['Balanced Accuracy']


ms15_bias50_Eomes_pos = ms15_bias50[ms15_bias50$TF == 'Eomes' & ms15_bias50$SignedQuantile == 1, 'Target']
ms15_bias50_Eomes_neg = ms15_bias50[ms15_bias50$TF == 'Eomes' & ms15_bias50$SignedQuantile == -1, 'Target']

de_Eomes_pos_list = ifelse(pool_genes %in% de_Eomes_pos, 1, 0)
ms15_bias50_Eomes_pos_list = ifelse(pool_genes %in% ms15_bias50_Eomes_pos, 1, 0)

xtab <- table(ms15_bias50_Eomes_pos_list, de_Eomes_pos_list)
cm_ms15_bias50 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15_bias50$overall['Accuracy']
cm_ms15_bias50$byClass['Precision']
cm_ms15_bias50$byClass['Recall']
cm_ms15_bias50$byClass['F1']
cm_ms15_bias50$byClass['Balanced Accuracy']



ms10_bias100_Eomes_pos = ms10_bias100[ms10_bias100$TF == 'Eomes' & ms10_bias100$SignedQuantile == 1, 'Target']
ms10_bias100_Eomes_neg = ms10_bias100[ms10_bias100$TF == 'Eomes' & ms10_bias100$SignedQuantile == -1, 'Target']

de_Eomes_pos_list = ifelse(pool_genes %in% de_Eomes_pos, 1, 0)
ms10_bias100_Eomes_pos_list = ifelse(pool_genes %in% ms10_bias100_Eomes_pos, 1, 0)

xtab <- table(ms10_bias100_Eomes_pos_list, de_Eomes_pos_list)
cm_ms10_bias100 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10_bias100$overall['Accuracy']
cm_ms10_bias100$byClass['Precision']
cm_ms10_bias100$byClass['Recall']
cm_ms10_bias100$byClass['F1']
cm_ms10_bias100$byClass['Balanced Accuracy']


ms10_bias5_Eomes_pos = ms10_bias5[ms10_bias5$TF == 'Eomes' & ms10_bias5$SignedQuantile == 1, 'Target']
ms10_bias5_Eomes_neg = ms10_bias5[ms10_bias5$TF == 'Eomes' & ms10_bias5$SignedQuantile == -1, 'Target']

de_Eomes_pos_list = ifelse(pool_genes %in% de_Eomes_pos, 1, 0)
ms10_bias5_Eomes_pos_list = ifelse(pool_genes %in% ms10_bias5_Eomes_pos, 1, 0)

xtab <- table(ms10_bias5_Eomes_pos_list, de_Eomes_pos_list)
cm_ms10_bias5 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10_bias5$overall['Accuracy']
cm_ms10_bias5$byClass['Precision']
cm_ms10_bias5$byClass['Recall']
cm_ms10_bias5$byClass['F1']
cm_ms10_bias5$byClass['Balanced Accuracy']


ms10_bias50_Eomes_pos = ms10_bias50[ms10_bias50$TF == 'Eomes' & ms10_bias50$SignedQuantile == 1, 'Target']
ms10_bias50_Eomes_neg = ms10_bias50[ms10_bias50$TF == 'Eomes' & ms10_bias50$SignedQuantile == -1, 'Target']

de_Eomes_pos_list = ifelse(pool_genes %in% de_Eomes_pos, 1, 0)
ms10_bias50_Eomes_pos_list = ifelse(pool_genes %in% ms10_bias50_Eomes_pos, 1, 0)

xtab <- table(ms10_bias50_Eomes_pos_list, de_Eomes_pos_list)
cm_ms10_bias50 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10_bias50$overall['Accuracy']
cm_ms10_bias50$byClass['Precision']
cm_ms10_bias50$byClass['Recall']
cm_ms10_bias50$byClass['F1']
cm_ms10_bias50$byClass['Balanced Accuracy']



ms5_bias100_Eomes_pos = ms5_bias100[ms5_bias100$TF == 'Eomes' & ms5_bias100$SignedQuantile == 1, 'Target']
ms5_bias100_Eomes_neg = ms5_bias100[ms5_bias100$TF == 'Eomes' & ms5_bias100$SignedQuantile == -1, 'Target']

#table(ms5_bias100_Eomes_pos %in% pool_genes)
#table(ms5_bias100_Eomes_neg %in% pool_genes)
#table(de_Eomes_pos %in% pool_genes)
#table(de_Eomes_neg %in% pool_genes)

de_Eomes_pos_list = ifelse(pool_genes %in% de_Eomes_pos, 1, 0)
ms5_bias100_Eomes_pos_list = ifelse(pool_genes %in% ms5_bias100_Eomes_pos, 1, 0)

xtab <- table(ms5_bias100_Eomes_pos_list, de_Eomes_pos_list)
#cm <- caret::confusionMatrix(xtab)
cm_ms5_bias100 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5_bias100$overall['Accuracy']
cm_ms5_bias100$byClass['Precision']
cm_ms5_bias100$byClass['Recall']
cm_ms5_bias100$byClass['F1']
cm_ms5_bias100$byClass['Balanced Accuracy']


ms5_bias5_Eomes_pos = ms5_bias5[ms5_bias5$TF == 'Eomes' & ms5_bias5$SignedQuantile == 1, 'Target']
ms5_bias5_Eomes_neg = ms5_bias5[ms5_bias5$TF == 'Eomes' & ms5_bias5$SignedQuantile == -1, 'Target']

de_Eomes_pos_list = ifelse(pool_genes %in% de_Eomes_pos, 1, 0)
ms5_bias5_Eomes_pos_list = ifelse(pool_genes %in% ms5_bias5_Eomes_pos, 1, 0)

xtab <- table(ms5_bias5_Eomes_pos_list, de_Eomes_pos_list)
cm_ms5_bias5 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5_bias5$overall['Accuracy']
cm_ms5_bias5$byClass['Precision']
cm_ms5_bias5$byClass['Recall']
cm_ms5_bias5$byClass['F1']
cm_ms5_bias5$byClass['Balanced Accuracy']


ms5_bias50_Eomes_pos = ms5_bias50[ms5_bias50$TF == 'Eomes' & ms5_bias50$SignedQuantile == 1, 'Target']
ms5_bias50_Eomes_neg = ms5_bias50[ms5_bias50$TF == 'Eomes' & ms5_bias50$SignedQuantile == -1, 'Target']

de_Eomes_pos_list = ifelse(pool_genes %in% de_Eomes_pos, 1, 0)
ms5_bias50_Eomes_pos_list = ifelse(pool_genes %in% ms5_bias50_Eomes_pos, 1, 0)

xtab <- table(ms5_bias50_Eomes_pos_list, de_Eomes_pos_list)
cm_ms5_bias50 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5_bias50$overall['Accuracy']
cm_ms5_bias50$byClass['Precision']
cm_ms5_bias50$byClass['Recall']
cm_ms5_bias50$byClass['F1']
cm_ms5_bias50$byClass['Balanced Accuracy']



ms20_bias100_Eomes_pos = ms20_bias100[ms20_bias100$TF == 'Eomes' & ms20_bias100$SignedQuantile == 1, 'Target']
ms20_bias100_Eomes_neg = ms20_bias100[ms20_bias100$TF == 'Eomes' & ms20_bias100$SignedQuantile == -1, 'Target']

de_Eomes_pos_list = ifelse(pool_genes %in% de_Eomes_pos, 1, 0)
ms20_bias100_Eomes_pos_list = ifelse(pool_genes %in% ms20_bias100_Eomes_pos, 1, 0)

xtab <- table(ms20_bias100_Eomes_pos_list, de_Eomes_pos_list)
cm_ms20_bias100 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20_bias100$overall['Accuracy']
cm_ms20_bias100$byClass['Precision']
cm_ms20_bias100$byClass['Recall']
cm_ms20_bias100$byClass['F1']
cm_ms20_bias100$byClass['Balanced Accuracy']


ms20_bias5_Eomes_pos = ms20_bias5[ms20_bias5$TF == 'Eomes' & ms20_bias5$SignedQuantile == 1, 'Target']
ms20_bias5_Eomes_neg = ms20_bias5[ms20_bias5$TF == 'Eomes' & ms20_bias5$SignedQuantile == -1, 'Target']

de_Eomes_pos_list = ifelse(pool_genes %in% de_Eomes_pos, 1, 0)
ms20_bias5_Eomes_pos_list = ifelse(pool_genes %in% ms20_bias5_Eomes_pos, 1, 0)

xtab <- table(ms20_bias5_Eomes_pos_list, de_Eomes_pos_list)
cm_ms20_bias5 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20_bias5$overall['Accuracy']
cm_ms20_bias5$byClass['Precision']
cm_ms20_bias5$byClass['Recall']
cm_ms20_bias5$byClass['F1']
cm_ms20_bias5$byClass['Balanced Accuracy']


ms20_bias50_Eomes_pos = ms20_bias50[ms20_bias50$TF == 'Eomes' & ms20_bias50$SignedQuantile == 1, 'Target']
ms20_bias50_Eomes_neg = ms20_bias50[ms20_bias50$TF == 'Eomes' & ms20_bias50$SignedQuantile == -1, 'Target']

de_Eomes_pos_list = ifelse(pool_genes %in% de_Eomes_pos, 1, 0)
ms20_bias50_Eomes_pos_list = ifelse(pool_genes %in% ms20_bias50_Eomes_pos, 1, 0)

xtab <- table(ms20_bias50_Eomes_pos_list, de_Eomes_pos_list)
cm_ms20_bias50 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20_bias50$overall['Accuracy']
cm_ms20_bias50$byClass['Precision']
cm_ms20_bias50$byClass['Recall']
cm_ms20_bias50$byClass['F1']
cm_ms20_bias50$byClass['Balanced Accuracy']


# library(Metrics)
# Metrics::recall(ms15_bias100_Eomes_pos_list, de_Eomes_pos_list)
# Metrics::precision(ms15_bias100_Eomes_pos_list, de_Eomes_pos_list)
# Metrics::f1(ms15_bias100_Eomes_pos_list, de_Eomes_pos_list) # wrong somehow


# make plots
forplot = data.frame(matrix(ncol = 3, nrow = 12))
colnames(forplot) = c('Accuracy', 'Balanced Accuracy', 'Precision')
rownames(forplot) = c('ms5_bias5', 'ms5_bias100', 'ms5_bias50', 'ms10_bias5', 'ms10_bias100', 'ms10_bias50', 
                      'ms15_bias5', 'ms15_bias100', 'ms15_bias50', 'ms20_bias5', 'ms20_bias100', 'ms20_bias50')
forplot['ms10_bias5', 'Accuracy'] = cm_ms10_bias5$overall['Accuracy']
forplot['ms10_bias100', 'Accuracy'] = cm_ms10_bias100$overall['Accuracy']
forplot['ms10_bias50', 'Accuracy'] = cm_ms10_bias50$overall['Accuracy']
forplot['ms15_bias5', 'Accuracy'] = cm_ms15_bias5$overall['Accuracy']
forplot['ms15_bias100', 'Accuracy'] = cm_ms15_bias100$overall['Accuracy']
forplot['ms15_bias50', 'Accuracy'] = cm_ms15_bias50$overall['Accuracy']

forplot['ms10_bias5', 'Balanced Accuracy'] = cm_ms10_bias5$byClass['Balanced Accuracy']
forplot['ms10_bias100', 'Balanced Accuracy'] = cm_ms10_bias100$byClass['Balanced Accuracy']
forplot['ms10_bias50', 'Balanced Accuracy'] = cm_ms10_bias50$byClass['Balanced Accuracy']
forplot['ms15_bias5', 'Balanced Accuracy'] = cm_ms15_bias5$byClass['Balanced Accuracy']
forplot['ms15_bias100', 'Balanced Accuracy'] = cm_ms15_bias100$byClass['Balanced Accuracy']
forplot['ms15_bias50', 'Balanced Accuracy'] = cm_ms15_bias50$byClass['Balanced Accuracy']

forplot['ms10_bias5', 'Precision'] = cm_ms10_bias5$byClass['Precision']
forplot['ms10_bias100', 'Precision'] = cm_ms10_bias100$byClass['Precision']
forplot['ms10_bias50', 'Precision'] = cm_ms10_bias50$byClass['Precision']
forplot['ms15_bias5', 'Precision'] = cm_ms15_bias5$byClass['Precision']
forplot['ms15_bias100', 'Precision'] = cm_ms15_bias100$byClass['Precision']
forplot['ms15_bias50', 'Precision'] = cm_ms15_bias50$byClass['Precision']


forplot['ms20_bias5', 'Accuracy'] = cm_ms20_bias5$overall['Accuracy']
forplot['ms20_bias100', 'Accuracy'] = cm_ms20_bias100$overall['Accuracy']
forplot['ms20_bias50', 'Accuracy'] = cm_ms20_bias50$overall['Accuracy']
forplot['ms5_bias5', 'Accuracy'] = cm_ms5_bias5$overall['Accuracy']
forplot['ms5_bias100', 'Accuracy'] = cm_ms5_bias100$overall['Accuracy']
forplot['ms5_bias50', 'Accuracy'] = cm_ms5_bias50$overall['Accuracy']

forplot['ms20_bias5', 'Balanced Accuracy'] = cm_ms20_bias5$byClass['Balanced Accuracy']
forplot['ms20_bias100', 'Balanced Accuracy'] = cm_ms20_bias100$byClass['Balanced Accuracy']
forplot['ms20_bias50', 'Balanced Accuracy'] = cm_ms20_bias50$byClass['Balanced Accuracy']
forplot['ms5_bias5', 'Balanced Accuracy'] = cm_ms5_bias5$byClass['Balanced Accuracy']
forplot['ms5_bias100', 'Balanced Accuracy'] = cm_ms5_bias100$byClass['Balanced Accuracy']
forplot['ms5_bias50', 'Balanced Accuracy'] = cm_ms5_bias50$byClass['Balanced Accuracy']

forplot['ms20_bias5', 'Precision'] = cm_ms20_bias5$byClass['Precision']
forplot['ms20_bias100', 'Precision'] = cm_ms20_bias100$byClass['Precision']
forplot['ms20_bias50', 'Precision'] = cm_ms20_bias50$byClass['Precision']
forplot['ms5_bias5', 'Precision'] = cm_ms5_bias5$byClass['Precision']
forplot['ms5_bias100', 'Precision'] = cm_ms5_bias100$byClass['Precision']
forplot['ms5_bias50', 'Precision'] = cm_ms5_bias50$byClass['Precision']

forplot$ms = c('Model size: 5', 'Model size: 5', 'Model size: 5', 
               'Model size: 10', 'Model size: 10', 'Model size: 10', 
               'Model size: 15', 'Model size: 15', 'Model size: 15',
               'Model size: 20', 'Model size: 20', 'Model size: 20')
forplot$ms = factor(forplot$ms, levels = c('Model size: 5', 'Model size: 10', 'Model size: 15', 'Model size: 20'))
forplot$TFAOpt = c('bias5', 'bias100', 'bias50', 'bias5', 'bias100', 'bias50', 'bias5', 'bias100', 'bias50', 'bias5', 'bias100', 'bias50')
forplot$TFAOpt = factor(forplot$TFAOpt, levels = c('bias5', 'bias50', 'bias100'))

forplot_melt = reshape2::melt(forplot)
library(ggplot2)
# library(RColorBrewer)
# brewer.pal(3, 'Set3')
eomes_pos_plot = ggplot(data=forplot_melt, aes(x=variable, y=value, fill=TFAOpt)) +
  geom_bar(stat="identity", position=position_dodge(), color="black", width = 0.7)+
  facet_wrap(~ms) + 
  scale_fill_manual(limits = c('bias5', 'bias50', 'bias100'),
                    labels = c('bias5', 'bias50', 'bias100'),
                    name = 'Bias',
                    values = c("#8DD3C7", "#FFFFB3", "#BEBADA"))+
  theme_bw(base_size = 14) +
  xlab('') +
  ylab('Value')+ 
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle('Eomes positive targets')




# negative targets

de_Eomes_neg_list = ifelse(pool_genes %in% de_Eomes_neg, 1, 0)
ms15_bias100_Eomes_neg_list = ifelse(pool_genes %in% ms15_bias100_Eomes_neg, 1, 0)

xtab <- table(ms15_bias100_Eomes_neg_list, de_Eomes_neg_list)
#cm <- caret::confusionMatrix(xtab)
cm_ms15_bias100 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15_bias100$overall['Accuracy']
cm_ms15_bias100$byClass['Precision']
cm_ms15_bias100$byClass['Recall']
cm_ms15_bias100$byClass['F1']
cm_ms15_bias100$byClass['Balanced Accuracy']


ms15_bias5_Eomes_neg = ms15_bias5[ms15_bias5$TF == 'Eomes' & ms15_bias5$SignedQuantile == -1, 'Target']

de_Eomes_neg_list = ifelse(pool_genes %in% de_Eomes_neg, 1, 0)
ms15_bias5_Eomes_neg_list = ifelse(pool_genes %in% ms15_bias5_Eomes_neg, 1, 0)

xtab <- table(ms15_bias5_Eomes_neg_list, de_Eomes_neg_list)
cm_ms15_bias5 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15_bias5$overall['Accuracy']
cm_ms15_bias5$byClass['Precision']
cm_ms15_bias5$byClass['Recall']
cm_ms15_bias5$byClass['F1']
cm_ms15_bias5$byClass['Balanced Accuracy']


ms15_bias50_Eomes_neg = ms15_bias50[ms15_bias50$TF == 'Eomes' & ms15_bias50$SignedQuantile == -1, 'Target']

de_Eomes_neg_list = ifelse(pool_genes %in% de_Eomes_neg, 1, 0)
ms15_bias50_Eomes_neg_list = ifelse(pool_genes %in% ms15_bias50_Eomes_neg, 1, 0)

xtab <- table(ms15_bias50_Eomes_neg_list, de_Eomes_neg_list)
cm_ms15_bias50 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15_bias50$overall['Accuracy']
cm_ms15_bias50$byClass['Precision']
cm_ms15_bias50$byClass['Recall']
cm_ms15_bias50$byClass['F1']
cm_ms15_bias50$byClass['Balanced Accuracy']



ms10_bias100_Eomes_neg = ms10_bias100[ms10_bias100$TF == 'Eomes' & ms10_bias100$SignedQuantile == -1, 'Target']

de_Eomes_neg_list = ifelse(pool_genes %in% de_Eomes_neg, 1, 0)
ms10_bias100_Eomes_neg_list = ifelse(pool_genes %in% ms10_bias100_Eomes_neg, 1, 0)

xtab <- table(ms10_bias100_Eomes_neg_list, de_Eomes_neg_list)
cm_ms10_bias100 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10_bias100$overall['Accuracy']
cm_ms10_bias100$byClass['Precision']
cm_ms10_bias100$byClass['Recall']
cm_ms10_bias100$byClass['F1']
cm_ms10_bias100$byClass['Balanced Accuracy']


ms10_bias5_Eomes_neg = ms10_bias5[ms10_bias5$TF == 'Eomes' & ms10_bias5$SignedQuantile == -1, 'Target']

de_Eomes_neg_list = ifelse(pool_genes %in% de_Eomes_neg, 1, 0)
ms10_bias5_Eomes_neg_list = ifelse(pool_genes %in% ms10_bias5_Eomes_neg, 1, 0)

xtab <- table(ms10_bias5_Eomes_neg_list, de_Eomes_neg_list)
cm_ms10_bias5 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10_bias5$overall['Accuracy']
cm_ms10_bias5$byClass['Precision']
cm_ms10_bias5$byClass['Recall']
cm_ms10_bias5$byClass['F1']
cm_ms10_bias5$byClass['Balanced Accuracy']


ms10_bias50_Eomes_neg = ms10_bias50[ms10_bias50$TF == 'Eomes' & ms10_bias50$SignedQuantile == -1, 'Target']

de_Eomes_neg_list = ifelse(pool_genes %in% de_Eomes_neg, 1, 0)
ms10_bias50_Eomes_neg_list = ifelse(pool_genes %in% ms10_bias50_Eomes_neg, 1, 0)

xtab <- table(ms10_bias50_Eomes_neg_list, de_Eomes_neg_list)
cm_ms10_bias50 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10_bias50$overall['Accuracy']
cm_ms10_bias50$byClass['Precision']
cm_ms10_bias50$byClass['Recall']
cm_ms10_bias50$byClass['F1']
cm_ms10_bias50$byClass['Balanced Accuracy']



de_Eomes_neg_list = ifelse(pool_genes %in% de_Eomes_neg, 1, 0)
ms5_bias100_Eomes_neg_list = ifelse(pool_genes %in% ms5_bias100_Eomes_neg, 1, 0)

xtab <- table(ms5_bias100_Eomes_neg_list, de_Eomes_neg_list)
#cm <- caret::confusionMatrix(xtab)
cm_ms5_bias100 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5_bias100$overall['Accuracy']
cm_ms5_bias100$byClass['Precision']
cm_ms5_bias100$byClass['Recall']
cm_ms5_bias100$byClass['F1']
cm_ms5_bias100$byClass['Balanced Accuracy']


ms5_bias5_Eomes_neg = ms5_bias5[ms5_bias5$TF == 'Eomes' & ms5_bias5$SignedQuantile == -1, 'Target']

de_Eomes_neg_list = ifelse(pool_genes %in% de_Eomes_neg, 1, 0)
ms5_bias5_Eomes_neg_list = ifelse(pool_genes %in% ms5_bias5_Eomes_neg, 1, 0)

xtab <- table(ms5_bias5_Eomes_neg_list, de_Eomes_neg_list)
cm_ms5_bias5 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5_bias5$overall['Accuracy']
cm_ms5_bias5$byClass['Precision']
cm_ms5_bias5$byClass['Recall']
cm_ms5_bias5$byClass['F1']
cm_ms5_bias5$byClass['Balanced Accuracy']


ms5_bias50_Eomes_neg = ms5_bias50[ms5_bias50$TF == 'Eomes' & ms5_bias50$SignedQuantile == -1, 'Target']

de_Eomes_neg_list = ifelse(pool_genes %in% de_Eomes_neg, 1, 0)
ms5_bias50_Eomes_neg_list = ifelse(pool_genes %in% ms5_bias50_Eomes_neg, 1, 0)

xtab <- table(ms5_bias50_Eomes_neg_list, de_Eomes_neg_list)
cm_ms5_bias50 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5_bias50$overall['Accuracy']
cm_ms5_bias50$byClass['Precision']
cm_ms5_bias50$byClass['Recall']
cm_ms5_bias50$byClass['F1']
cm_ms5_bias50$byClass['Balanced Accuracy']



ms20_bias100_Eomes_neg = ms20_bias100[ms20_bias100$TF == 'Eomes' & ms20_bias100$SignedQuantile == -1, 'Target']

de_Eomes_neg_list = ifelse(pool_genes %in% de_Eomes_neg, 1, 0)
ms20_bias100_Eomes_neg_list = ifelse(pool_genes %in% ms20_bias100_Eomes_neg, 1, 0)

xtab <- table(ms20_bias100_Eomes_neg_list, de_Eomes_neg_list)
cm_ms20_bias100 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20_bias100$overall['Accuracy']
cm_ms20_bias100$byClass['Precision']
cm_ms20_bias100$byClass['Recall']
cm_ms20_bias100$byClass['F1']
cm_ms20_bias100$byClass['Balanced Accuracy']


ms20_bias5_Eomes_neg = ms20_bias5[ms20_bias5$TF == 'Eomes' & ms20_bias5$SignedQuantile == -1, 'Target']

de_Eomes_neg_list = ifelse(pool_genes %in% de_Eomes_neg, 1, 0)
ms20_bias5_Eomes_neg_list = ifelse(pool_genes %in% ms20_bias5_Eomes_neg, 1, 0)

xtab <- table(ms20_bias5_Eomes_neg_list, de_Eomes_neg_list)
cm_ms20_bias5 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20_bias5$overall['Accuracy']
cm_ms20_bias5$byClass['Precision']
cm_ms20_bias5$byClass['Recall']
cm_ms20_bias5$byClass['F1']
cm_ms20_bias5$byClass['Balanced Accuracy']


ms20_bias50_Eomes_neg = ms20_bias50[ms20_bias50$TF == 'Eomes' & ms20_bias50$SignedQuantile == -1, 'Target']

de_Eomes_neg_list = ifelse(pool_genes %in% de_Eomes_neg, 1, 0)
ms20_bias50_Eomes_neg_list = ifelse(pool_genes %in% ms20_bias50_Eomes_neg, 1, 0)

xtab <- table(ms20_bias50_Eomes_neg_list, de_Eomes_neg_list)
cm_ms20_bias50 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20_bias50$overall['Accuracy']
cm_ms20_bias50$byClass['Precision']
cm_ms20_bias50$byClass['Recall']
cm_ms20_bias50$byClass['F1']
cm_ms20_bias50$byClass['Balanced Accuracy']


# library(Metrics)
# Metrics::recall(ms15_bias100_Eomes_neg_list, de_Eomes_neg_list)
# Metrics::precision(ms15_bias100_Eomes_neg_list, de_Eomes_neg_list)
# Metrics::f1(ms15_bias100_Eomes_neg_list, de_Eomes_neg_list) # wrong somehow


# make plots
forplotneg = data.frame(matrix(ncol = 3, nrow = 12))
colnames(forplotneg) = c('Accuracy', 'Balanced Accuracy', 'Precision')
rownames(forplotneg) = c('ms5_bias5', 'ms5_bias100', 'ms5_bias50', 
                         'ms10_bias5', 'ms10_bias100', 'ms10_bias50', 
                         'ms15_bias5', 'ms15_bias100', 'ms15_bias50',
                         'ms20_bias5', 'ms20_bias100', 'ms20_bias50')
forplotneg['ms10_bias5', 'Accuracy'] = cm_ms10_bias5$overall['Accuracy']
forplotneg['ms10_bias100', 'Accuracy'] = cm_ms10_bias100$overall['Accuracy']
forplotneg['ms10_bias50', 'Accuracy'] = cm_ms10_bias50$overall['Accuracy']
forplotneg['ms15_bias5', 'Accuracy'] = cm_ms15_bias5$overall['Accuracy']
forplotneg['ms15_bias100', 'Accuracy'] = cm_ms15_bias100$overall['Accuracy']
forplotneg['ms15_bias50', 'Accuracy'] = cm_ms15_bias50$overall['Accuracy']

forplotneg['ms10_bias5', 'Balanced Accuracy'] = cm_ms10_bias5$byClass['Balanced Accuracy']
forplotneg['ms10_bias100', 'Balanced Accuracy'] = cm_ms10_bias100$byClass['Balanced Accuracy']
forplotneg['ms10_bias50', 'Balanced Accuracy'] = cm_ms10_bias50$byClass['Balanced Accuracy']
forplotneg['ms15_bias5', 'Balanced Accuracy'] = cm_ms15_bias5$byClass['Balanced Accuracy']
forplotneg['ms15_bias100', 'Balanced Accuracy'] = cm_ms15_bias100$byClass['Balanced Accuracy']
forplotneg['ms15_bias50', 'Balanced Accuracy'] = cm_ms15_bias50$byClass['Balanced Accuracy']

forplotneg['ms10_bias5', 'Precision'] = cm_ms10_bias5$byClass['Precision']
forplotneg['ms10_bias100', 'Precision'] = cm_ms10_bias100$byClass['Precision']
forplotneg['ms10_bias50', 'Precision'] = cm_ms10_bias50$byClass['Precision']
forplotneg['ms15_bias5', 'Precision'] = cm_ms15_bias5$byClass['Precision']
forplotneg['ms15_bias100', 'Precision'] = cm_ms15_bias100$byClass['Precision']
forplotneg['ms15_bias50', 'Precision'] = cm_ms15_bias50$byClass['Precision']


forplotneg['ms20_bias5', 'Accuracy'] = cm_ms20_bias5$overall['Accuracy']
forplotneg['ms20_bias100', 'Accuracy'] = cm_ms20_bias100$overall['Accuracy']
forplotneg['ms20_bias50', 'Accuracy'] = cm_ms20_bias50$overall['Accuracy']
forplotneg['ms5_bias5', 'Accuracy'] = cm_ms5_bias5$overall['Accuracy']
forplotneg['ms5_bias100', 'Accuracy'] = cm_ms5_bias100$overall['Accuracy']
forplotneg['ms5_bias50', 'Accuracy'] = cm_ms5_bias50$overall['Accuracy']

forplotneg['ms20_bias5', 'Balanced Accuracy'] = cm_ms20_bias5$byClass['Balanced Accuracy']
forplotneg['ms20_bias100', 'Balanced Accuracy'] = cm_ms20_bias100$byClass['Balanced Accuracy']
forplotneg['ms20_bias50', 'Balanced Accuracy'] = cm_ms20_bias50$byClass['Balanced Accuracy']
forplotneg['ms5_bias5', 'Balanced Accuracy'] = cm_ms5_bias5$byClass['Balanced Accuracy']
forplotneg['ms5_bias100', 'Balanced Accuracy'] = cm_ms5_bias100$byClass['Balanced Accuracy']
forplotneg['ms5_bias50', 'Balanced Accuracy'] = cm_ms5_bias50$byClass['Balanced Accuracy']

forplotneg['ms20_bias5', 'Precision'] = cm_ms20_bias5$byClass['Precision']
forplotneg['ms20_bias100', 'Precision'] = cm_ms20_bias100$byClass['Precision']
forplotneg['ms20_bias50', 'Precision'] = cm_ms20_bias50$byClass['Precision']
forplotneg['ms5_bias5', 'Precision'] = cm_ms5_bias5$byClass['Precision']
forplotneg['ms5_bias100', 'Precision'] = cm_ms5_bias100$byClass['Precision']
forplotneg['ms5_bias50', 'Precision'] = cm_ms5_bias50$byClass['Precision']

forplotneg$ms = c('Model size: 5', 'Model size: 5', 'Model size: 5', 
                  'Model size: 10', 'Model size: 10', 'Model size: 10', 
                  'Model size: 15', 'Model size: 15', 'Model size: 15',
                  'Model size: 20', 'Model size: 20', 'Model size: 20')
forplotneg$ms = factor(forplotneg$ms, levels = c('Model size: 5', 'Model size: 10', 'Model size: 15', 'Model size: 20'))
forplotneg$TFAOpt = c('bias5', 'bias100', 'bias50', 'bias5', 'bias100', 'bias50', 'bias5', 'bias100', 'bias50', 'bias5', 'bias100', 'bias50')
forplotneg$TFAOpt = factor(forplotneg$TFAOpt, levels = c('bias5', 'bias50', 'bias100'))

forplotneg_melt = reshape2::melt(forplotneg)
library(ggplot2)
# library(RColorBrewer)
# brewer.pal(3, 'Set3')
eomes_neg_plot = ggplot(data=forplotneg_melt, aes(x=variable, y=value, fill=TFAOpt)) +
  geom_bar(stat="identity", position=position_dodge(), color="black", width = 0.7)+
  facet_wrap(~ms) + 
  scale_fill_manual(limits = c('bias5', 'bias50', 'bias100'),
                    labels = c('bias5', 'bias50', 'bias100'),
                    name = 'Bias',
                    values = c("#8DD3C7", "#FFFFB3", "#BEBADA"))+
  theme_bw(base_size = 14) +
  xlab('') +
  ylab('Value')+ 
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle('Eomes negative targets')

library(cowplot)
eomes.combined = plot_grid(eomes_pos_plot, eomes_neg_plot, align="hv", ncol = 2)
ggsave(
  paste0('compare_ms_bias_Eomes_remake.pdf'),
  plot = eomes.combined,
  device = "pdf",
  width = 12,
  height = 6,
  dpi = 300
)


# Foxo1

# Foxo targets in Characterization of the direct targets of FOXO transcription factors throughout evolution

pool_genes = read.table('deg_padj01.txt')$V1

library(readxl)

targets = read_excel("foxo1/ACEL-15-673-s002.xlsx", 
                     sheet = "Sheet1")
targets_cd8 = unique(targets$`CD8 T cells`)
targets_allcelltype = unique(targets$Core)

targets_cd8 = targets_cd8[!is.na(targets_cd8)]
targets_allcelltype = targets_allcelltype[!is.na(targets_allcelltype)]

targets_cd8_pool_genes = targets_cd8[targets_cd8 %in% pool_genes]
targets_allcelltype_pool_genes = targets_allcelltype[targets_allcelltype %in% pool_genes]


ms10_bias50  = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms10_bias50_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms10_bias5  = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms10_bias5_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms10_bias100 = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms10_bias100_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms15_bias50  = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms15_bias50_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms15_bias5  = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms15_bias5_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms15_bias100 = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms15_bias100_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)

ms5_bias50  = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms5_bias50_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms5_bias5  = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms5_bias5_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms5_bias100 = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms5_bias100_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms20_bias50  = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms20_bias50_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms20_bias5  = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms20_bias5_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms20_bias100 = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms20_bias100_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)

# both

ms15_bias100_Foxo1 = ms15_bias100[ms15_bias100$TF == 'Foxo1', 'Target']

#table(ms15_bias100_Foxo1 %in% pool_genes)
#table(targets_cd8_pool_genes %in% pool_genes)

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms15_bias100_Foxo1_list = ifelse(pool_genes %in% ms15_bias100_Foxo1, 1, 0)

xtab <- table(ms15_bias100_Foxo1_list, targets_cd8_pool_genes_list)
#cm <- caret::confusionMatrix(xtab)
cm_ms15_bias100 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15_bias100$overall['Accuracy']
cm_ms15_bias100$byClass['Precision']
cm_ms15_bias100$byClass['Recall']
cm_ms15_bias100$byClass['F1']
cm_ms15_bias100$byClass['Balanced Accuracy']


ms15_bias5_Foxo1 = ms15_bias5[ms15_bias5$TF == 'Foxo1', 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms15_bias5_Foxo1_list = ifelse(pool_genes %in% ms15_bias5_Foxo1, 1, 0)

xtab <- table(ms15_bias5_Foxo1_list, targets_cd8_pool_genes_list)
cm_ms15_bias5 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15_bias5$overall['Accuracy']
cm_ms15_bias5$byClass['Precision']
cm_ms15_bias5$byClass['Recall']
cm_ms15_bias5$byClass['F1']
cm_ms15_bias5$byClass['Balanced Accuracy']


ms15_bias50_Foxo1 = ms15_bias50[ms15_bias50$TF == 'Foxo1', 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms15_bias50_Foxo1_list = ifelse(pool_genes %in% ms15_bias50_Foxo1, 1, 0)

xtab <- table(ms15_bias50_Foxo1_list, targets_cd8_pool_genes_list)
cm_ms15_bias50 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15_bias50$overall['Accuracy']
cm_ms15_bias50$byClass['Precision']
cm_ms15_bias50$byClass['Recall']
cm_ms15_bias50$byClass['F1']
cm_ms15_bias50$byClass['Balanced Accuracy']



ms10_bias100_Foxo1 = ms10_bias100[ms10_bias100$TF == 'Foxo1', 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms10_bias100_Foxo1_list = ifelse(pool_genes %in% ms10_bias100_Foxo1, 1, 0)

xtab <- table(ms10_bias100_Foxo1_list, targets_cd8_pool_genes_list)
cm_ms10_bias100 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10_bias100$overall['Accuracy']
cm_ms10_bias100$byClass['Precision']
cm_ms10_bias100$byClass['Recall']
cm_ms10_bias100$byClass['F1']
cm_ms10_bias100$byClass['Balanced Accuracy']


ms10_bias5_Foxo1 = ms10_bias5[ms10_bias5$TF == 'Foxo1', 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms10_bias5_Foxo1_list = ifelse(pool_genes %in% ms10_bias5_Foxo1, 1, 0)

xtab <- table(ms10_bias5_Foxo1_list, targets_cd8_pool_genes_list)
cm_ms10_bias5 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10_bias5$overall['Accuracy']
cm_ms10_bias5$byClass['Precision']
cm_ms10_bias5$byClass['Recall']
cm_ms10_bias5$byClass['F1']
cm_ms10_bias5$byClass['Balanced Accuracy']


ms10_bias50_Foxo1 = ms10_bias50[ms10_bias50$TF == 'Foxo1', 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms10_bias50_Foxo1_list = ifelse(pool_genes %in% ms10_bias50_Foxo1, 1, 0)

xtab <- table(ms10_bias50_Foxo1_list, targets_cd8_pool_genes_list)
cm_ms10_bias50 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10_bias50$overall['Accuracy']
cm_ms10_bias50$byClass['Precision']
cm_ms10_bias50$byClass['Recall']
cm_ms10_bias50$byClass['F1']
cm_ms10_bias50$byClass['Balanced Accuracy']


ms5_bias100_Foxo1 = ms5_bias100[ms5_bias100$TF == 'Foxo1', 'Target']

#table(ms5_bias100_Foxo1 %in% pool_genes)
#table(targets_cd8_pool_genes %in% pool_genes)

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms5_bias100_Foxo1_list = ifelse(pool_genes %in% ms5_bias100_Foxo1, 1, 0)

xtab <- table(ms5_bias100_Foxo1_list, targets_cd8_pool_genes_list)
#cm <- caret::confusionMatrix(xtab)
cm_ms5_bias100 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5_bias100$overall['Accuracy']
cm_ms5_bias100$byClass['Precision']
cm_ms5_bias100$byClass['Recall']
cm_ms5_bias100$byClass['F1']
cm_ms5_bias100$byClass['Balanced Accuracy']


ms5_bias5_Foxo1 = ms5_bias5[ms5_bias5$TF == 'Foxo1', 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms5_bias5_Foxo1_list = ifelse(pool_genes %in% ms5_bias5_Foxo1, 1, 0)

xtab <- table(ms5_bias5_Foxo1_list, targets_cd8_pool_genes_list)
cm_ms5_bias5 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5_bias5$overall['Accuracy']
cm_ms5_bias5$byClass['Precision']
cm_ms5_bias5$byClass['Recall']
cm_ms5_bias5$byClass['F1']
cm_ms5_bias5$byClass['Balanced Accuracy']


ms5_bias50_Foxo1 = ms5_bias50[ms5_bias50$TF == 'Foxo1', 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms5_bias50_Foxo1_list = ifelse(pool_genes %in% ms5_bias50_Foxo1, 1, 0)

xtab <- table(ms5_bias50_Foxo1_list, targets_cd8_pool_genes_list)
cm_ms5_bias50 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5_bias50$overall['Accuracy']
cm_ms5_bias50$byClass['Precision']
cm_ms5_bias50$byClass['Recall']
cm_ms5_bias50$byClass['F1']
cm_ms5_bias50$byClass['Balanced Accuracy']



ms20_bias100_Foxo1 = ms20_bias100[ms20_bias100$TF == 'Foxo1', 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms20_bias100_Foxo1_list = ifelse(pool_genes %in% ms20_bias100_Foxo1, 1, 0)

xtab <- table(ms20_bias100_Foxo1_list, targets_cd8_pool_genes_list)
cm_ms20_bias100 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20_bias100$overall['Accuracy']
cm_ms20_bias100$byClass['Precision']
cm_ms20_bias100$byClass['Recall']
cm_ms20_bias100$byClass['F1']
cm_ms20_bias100$byClass['Balanced Accuracy']


ms20_bias5_Foxo1 = ms20_bias5[ms20_bias5$TF == 'Foxo1', 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms20_bias5_Foxo1_list = ifelse(pool_genes %in% ms20_bias5_Foxo1, 1, 0)

xtab <- table(ms20_bias5_Foxo1_list, targets_cd8_pool_genes_list)
cm_ms20_bias5 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20_bias5$overall['Accuracy']
cm_ms20_bias5$byClass['Precision']
cm_ms20_bias5$byClass['Recall']
cm_ms20_bias5$byClass['F1']
cm_ms20_bias5$byClass['Balanced Accuracy']


ms20_bias50_Foxo1 = ms20_bias50[ms20_bias50$TF == 'Foxo1', 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms20_bias50_Foxo1_list = ifelse(pool_genes %in% ms20_bias50_Foxo1, 1, 0)

xtab <- table(ms20_bias50_Foxo1_list, targets_cd8_pool_genes_list)
cm_ms20_bias50 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20_bias50$overall['Accuracy']
cm_ms20_bias50$byClass['Precision']
cm_ms20_bias50$byClass['Recall']
cm_ms20_bias50$byClass['F1']
cm_ms20_bias50$byClass['Balanced Accuracy']


# library(Metrics)
# Metrics::recall(ms15_bias100_Foxo1_list, targets_cd8_pool_genes_list)
# Metrics::precision(ms15_bias100_Foxo1_list, targets_cd8_pool_genes_list)
# Metrics::f1(ms15_bias100_Foxo1_list, targets_cd8_pool_genes_list) # wrong somehow


# make plots
forplot = data.frame(matrix(ncol = 3, nrow = 12))
colnames(forplot) = c('Accuracy', 'Balanced Accuracy', 'Precision')
rownames(forplot) = c('ms5_bias5', 'ms5_bias100', 'ms5_bias50', 'ms10_bias5', 'ms10_bias100', 'ms10_bias50', 
                      'ms15_bias5', 'ms15_bias100', 'ms15_bias50', 'ms20_bias5', 'ms20_bias100', 'ms20_bias50')
forplot['ms10_bias5', 'Accuracy'] = cm_ms10_bias5$overall['Accuracy']
forplot['ms10_bias100', 'Accuracy'] = cm_ms10_bias100$overall['Accuracy']
forplot['ms10_bias50', 'Accuracy'] = cm_ms10_bias50$overall['Accuracy']
forplot['ms15_bias5', 'Accuracy'] = cm_ms15_bias5$overall['Accuracy']
forplot['ms15_bias100', 'Accuracy'] = cm_ms15_bias100$overall['Accuracy']
forplot['ms15_bias50', 'Accuracy'] = cm_ms15_bias50$overall['Accuracy']

forplot['ms10_bias5', 'Balanced Accuracy'] = cm_ms10_bias5$byClass['Balanced Accuracy']
forplot['ms10_bias100', 'Balanced Accuracy'] = cm_ms10_bias100$byClass['Balanced Accuracy']
forplot['ms10_bias50', 'Balanced Accuracy'] = cm_ms10_bias50$byClass['Balanced Accuracy']
forplot['ms15_bias5', 'Balanced Accuracy'] = cm_ms15_bias5$byClass['Balanced Accuracy']
forplot['ms15_bias100', 'Balanced Accuracy'] = cm_ms15_bias100$byClass['Balanced Accuracy']
forplot['ms15_bias50', 'Balanced Accuracy'] = cm_ms15_bias50$byClass['Balanced Accuracy']

forplot['ms10_bias5', 'Precision'] = cm_ms10_bias5$byClass['Precision']
forplot['ms10_bias100', 'Precision'] = cm_ms10_bias100$byClass['Precision']
forplot['ms10_bias50', 'Precision'] = cm_ms10_bias50$byClass['Precision']
forplot['ms15_bias5', 'Precision'] = cm_ms15_bias5$byClass['Precision']
forplot['ms15_bias100', 'Precision'] = cm_ms15_bias100$byClass['Precision']
forplot['ms15_bias50', 'Precision'] = cm_ms15_bias50$byClass['Precision']


forplot['ms20_bias5', 'Accuracy'] = cm_ms20_bias5$overall['Accuracy']
forplot['ms20_bias100', 'Accuracy'] = cm_ms20_bias100$overall['Accuracy']
forplot['ms20_bias50', 'Accuracy'] = cm_ms20_bias50$overall['Accuracy']
forplot['ms5_bias5', 'Accuracy'] = cm_ms5_bias5$overall['Accuracy']
forplot['ms5_bias100', 'Accuracy'] = cm_ms5_bias100$overall['Accuracy']
forplot['ms5_bias50', 'Accuracy'] = cm_ms5_bias50$overall['Accuracy']

forplot['ms20_bias5', 'Balanced Accuracy'] = cm_ms20_bias5$byClass['Balanced Accuracy']
forplot['ms20_bias100', 'Balanced Accuracy'] = cm_ms20_bias100$byClass['Balanced Accuracy']
forplot['ms20_bias50', 'Balanced Accuracy'] = cm_ms20_bias50$byClass['Balanced Accuracy']
forplot['ms5_bias5', 'Balanced Accuracy'] = cm_ms5_bias5$byClass['Balanced Accuracy']
forplot['ms5_bias100', 'Balanced Accuracy'] = cm_ms5_bias100$byClass['Balanced Accuracy']
forplot['ms5_bias50', 'Balanced Accuracy'] = cm_ms5_bias50$byClass['Balanced Accuracy']

forplot['ms20_bias5', 'Precision'] = cm_ms20_bias5$byClass['Precision']
forplot['ms20_bias100', 'Precision'] = cm_ms20_bias100$byClass['Precision']
forplot['ms20_bias50', 'Precision'] = cm_ms20_bias50$byClass['Precision']
forplot['ms5_bias5', 'Precision'] = cm_ms5_bias5$byClass['Precision']
forplot['ms5_bias100', 'Precision'] = cm_ms5_bias100$byClass['Precision']
forplot['ms5_bias50', 'Precision'] = cm_ms5_bias50$byClass['Precision']

forplot$ms = c('Model size: 5', 'Model size: 5', 'Model size: 5', 
               'Model size: 10', 'Model size: 10', 'Model size: 10', 
               'Model size: 15', 'Model size: 15', 'Model size: 15',
               'Model size: 20', 'Model size: 20', 'Model size: 20')
forplot$ms = factor(forplot$ms, levels = c('Model size: 5', 'Model size: 10', 'Model size: 15', 'Model size: 20'))
forplot$TFAOpt = c('bias5', 'bias100', 'bias50', 'bias5', 'bias100', 'bias50', 'bias5', 'bias100', 'bias50', 'bias5', 'bias100', 'bias50')
forplot$TFAOpt = factor(forplot$TFAOpt, levels = c('bias5', 'bias50', 'bias100'))

forplot_melt = reshape2::melt(forplot)
library(ggplot2)
# library(RColorBrewer)
# brewer.pal(3, 'Set3')
foxo1_plot = ggplot(data=forplot_melt, aes(x=variable, y=value, fill=TFAOpt)) +
  geom_bar(stat="identity", position=position_dodge(), color="black", width = 0.7)+
  facet_wrap(~ms) + 
  scale_fill_manual(limits = c('bias5', 'bias50', 'bias100'),
                    labels = c('bias5', 'bias50', 'bias100'),
                    name = 'Bias',
                    values = c("#8DD3C7", "#FFFFB3", "#BEBADA"))+
  theme_bw(base_size = 14) +
  xlab('') +
  ylab('Value')+ 
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle('Foxo1 targets')

ggsave(
  paste0('compare_ms_bias_Foxo1_moreMs.pdf'),
  plot = foxo1_plot,
  device = "pdf",
  width = 6,
  height = 6,
  dpi = 300
)



# positive targets
ms15_bias100_Foxo1_pos = ms15_bias100[ms15_bias100$TF == 'Foxo1' & ms15_bias100$SignedQuantile == 1, 'Target']
ms15_bias100_Foxo1_neg = ms15_bias100[ms15_bias100$TF == 'Foxo1' & ms15_bias100$SignedQuantile == -1, 'Target']

#table(ms15_bias100_Foxo1_pos %in% pool_genes)
#table(ms15_bias100_Foxo1_neg %in% pool_genes)
#table(targets_cd8_pool_genes %in% pool_genes)
#table(targets_cd8_pool_genes %in% pool_genes)

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms15_bias100_Foxo1_pos_list = ifelse(pool_genes %in% ms15_bias100_Foxo1_pos, 1, 0)

xtab <- table(ms15_bias100_Foxo1_pos_list, targets_cd8_pool_genes_list)
#cm <- caret::confusionMatrix(xtab)
cm_ms15_bias100 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15_bias100$overall['Accuracy']
cm_ms15_bias100$byClass['Precision']
cm_ms15_bias100$byClass['Recall']
cm_ms15_bias100$byClass['F1']
cm_ms15_bias100$byClass['Balanced Accuracy']


ms15_bias5_Foxo1_pos = ms15_bias5[ms15_bias5$TF == 'Foxo1' & ms15_bias5$SignedQuantile == 1, 'Target']
ms15_bias5_Foxo1_neg = ms15_bias5[ms15_bias5$TF == 'Foxo1' & ms15_bias5$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms15_bias5_Foxo1_pos_list = ifelse(pool_genes %in% ms15_bias5_Foxo1_pos, 1, 0)

xtab <- table(ms15_bias5_Foxo1_pos_list, targets_cd8_pool_genes_list)
cm_ms15_bias5 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15_bias5$overall['Accuracy']
cm_ms15_bias5$byClass['Precision']
cm_ms15_bias5$byClass['Recall']
cm_ms15_bias5$byClass['F1']
cm_ms15_bias5$byClass['Balanced Accuracy']


ms15_bias50_Foxo1_pos = ms15_bias50[ms15_bias50$TF == 'Foxo1' & ms15_bias50$SignedQuantile == 1, 'Target']
ms15_bias50_Foxo1_neg = ms15_bias50[ms15_bias50$TF == 'Foxo1' & ms15_bias50$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms15_bias50_Foxo1_pos_list = ifelse(pool_genes %in% ms15_bias50_Foxo1_pos, 1, 0)

xtab <- table(ms15_bias50_Foxo1_pos_list, targets_cd8_pool_genes_list)
cm_ms15_bias50 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15_bias50$overall['Accuracy']
cm_ms15_bias50$byClass['Precision']
cm_ms15_bias50$byClass['Recall']
cm_ms15_bias50$byClass['F1']
cm_ms15_bias50$byClass['Balanced Accuracy']



ms10_bias100_Foxo1_pos = ms10_bias100[ms10_bias100$TF == 'Foxo1' & ms10_bias100$SignedQuantile == 1, 'Target']
ms10_bias100_Foxo1_neg = ms10_bias100[ms10_bias100$TF == 'Foxo1' & ms10_bias100$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms10_bias100_Foxo1_pos_list = ifelse(pool_genes %in% ms10_bias100_Foxo1_pos, 1, 0)

xtab <- table(ms10_bias100_Foxo1_pos_list, targets_cd8_pool_genes_list)
cm_ms10_bias100 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10_bias100$overall['Accuracy']
cm_ms10_bias100$byClass['Precision']
cm_ms10_bias100$byClass['Recall']
cm_ms10_bias100$byClass['F1']
cm_ms10_bias100$byClass['Balanced Accuracy']


ms10_bias5_Foxo1_pos = ms10_bias5[ms10_bias5$TF == 'Foxo1' & ms10_bias5$SignedQuantile == 1, 'Target']
ms10_bias5_Foxo1_neg = ms10_bias5[ms10_bias5$TF == 'Foxo1' & ms10_bias5$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms10_bias5_Foxo1_pos_list = ifelse(pool_genes %in% ms10_bias5_Foxo1_pos, 1, 0)

xtab <- table(ms10_bias5_Foxo1_pos_list, targets_cd8_pool_genes_list)
cm_ms10_bias5 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10_bias5$overall['Accuracy']
cm_ms10_bias5$byClass['Precision']
cm_ms10_bias5$byClass['Recall']
cm_ms10_bias5$byClass['F1']
cm_ms10_bias5$byClass['Balanced Accuracy']


ms10_bias50_Foxo1_pos = ms10_bias50[ms10_bias50$TF == 'Foxo1' & ms10_bias50$SignedQuantile == 1, 'Target']
ms10_bias50_Foxo1_neg = ms10_bias50[ms10_bias50$TF == 'Foxo1' & ms10_bias50$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms10_bias50_Foxo1_pos_list = ifelse(pool_genes %in% ms10_bias50_Foxo1_pos, 1, 0)

xtab <- table(ms10_bias50_Foxo1_pos_list, targets_cd8_pool_genes_list)
cm_ms10_bias50 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10_bias50$overall['Accuracy']
cm_ms10_bias50$byClass['Precision']
cm_ms10_bias50$byClass['Recall']
cm_ms10_bias50$byClass['F1']
cm_ms10_bias50$byClass['Balanced Accuracy']


ms5_bias100_Foxo1_pos = ms5_bias100[ms5_bias100$TF == 'Foxo1' & ms5_bias100$SignedQuantile == 1, 'Target']
ms5_bias100_Foxo1_neg = ms5_bias100[ms5_bias100$TF == 'Foxo1' & ms5_bias100$SignedQuantile == -1, 'Target']

#table(ms5_bias100_Foxo1_pos %in% pool_genes)
#table(ms5_bias100_Foxo1_neg %in% pool_genes)
#table(targets_cd8_pool_genes %in% pool_genes)
#table(targets_cd8_pool_genes %in% pool_genes)

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms5_bias100_Foxo1_pos_list = ifelse(pool_genes %in% ms5_bias100_Foxo1_pos, 1, 0)

xtab <- table(ms5_bias100_Foxo1_pos_list, targets_cd8_pool_genes_list)
#cm <- caret::confusionMatrix(xtab)
cm_ms5_bias100 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5_bias100$overall['Accuracy']
cm_ms5_bias100$byClass['Precision']
cm_ms5_bias100$byClass['Recall']
cm_ms5_bias100$byClass['F1']
cm_ms5_bias100$byClass['Balanced Accuracy']


ms5_bias5_Foxo1_pos = ms5_bias5[ms5_bias5$TF == 'Foxo1' & ms5_bias5$SignedQuantile == 1, 'Target']
ms5_bias5_Foxo1_neg = ms5_bias5[ms5_bias5$TF == 'Foxo1' & ms5_bias5$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms5_bias5_Foxo1_pos_list = ifelse(pool_genes %in% ms5_bias5_Foxo1_pos, 1, 0)

xtab <- table(ms5_bias5_Foxo1_pos_list, targets_cd8_pool_genes_list)
cm_ms5_bias5 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5_bias5$overall['Accuracy']
cm_ms5_bias5$byClass['Precision']
cm_ms5_bias5$byClass['Recall']
cm_ms5_bias5$byClass['F1']
cm_ms5_bias5$byClass['Balanced Accuracy']


ms5_bias50_Foxo1_pos = ms5_bias50[ms5_bias50$TF == 'Foxo1' & ms5_bias50$SignedQuantile == 1, 'Target']
ms5_bias50_Foxo1_neg = ms5_bias50[ms5_bias50$TF == 'Foxo1' & ms5_bias50$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms5_bias50_Foxo1_pos_list = ifelse(pool_genes %in% ms5_bias50_Foxo1_pos, 1, 0)

xtab <- table(ms5_bias50_Foxo1_pos_list, targets_cd8_pool_genes_list)
cm_ms5_bias50 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5_bias50$overall['Accuracy']
cm_ms5_bias50$byClass['Precision']
cm_ms5_bias50$byClass['Recall']
cm_ms5_bias50$byClass['F1']
cm_ms5_bias50$byClass['Balanced Accuracy']



ms20_bias100_Foxo1_pos = ms20_bias100[ms20_bias100$TF == 'Foxo1' & ms20_bias100$SignedQuantile == 1, 'Target']
ms20_bias100_Foxo1_neg = ms20_bias100[ms20_bias100$TF == 'Foxo1' & ms20_bias100$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms20_bias100_Foxo1_pos_list = ifelse(pool_genes %in% ms20_bias100_Foxo1_pos, 1, 0)

xtab <- table(ms20_bias100_Foxo1_pos_list, targets_cd8_pool_genes_list)
cm_ms20_bias100 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20_bias100$overall['Accuracy']
cm_ms20_bias100$byClass['Precision']
cm_ms20_bias100$byClass['Recall']
cm_ms20_bias100$byClass['F1']
cm_ms20_bias100$byClass['Balanced Accuracy']


ms20_bias5_Foxo1_pos = ms20_bias5[ms20_bias5$TF == 'Foxo1' & ms20_bias5$SignedQuantile == 1, 'Target']
ms20_bias5_Foxo1_neg = ms20_bias5[ms20_bias5$TF == 'Foxo1' & ms20_bias5$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms20_bias5_Foxo1_pos_list = ifelse(pool_genes %in% ms20_bias5_Foxo1_pos, 1, 0)

xtab <- table(ms20_bias5_Foxo1_pos_list, targets_cd8_pool_genes_list)
cm_ms20_bias5 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20_bias5$overall['Accuracy']
cm_ms20_bias5$byClass['Precision']
cm_ms20_bias5$byClass['Recall']
cm_ms20_bias5$byClass['F1']
cm_ms20_bias5$byClass['Balanced Accuracy']


ms20_bias50_Foxo1_pos = ms20_bias50[ms20_bias50$TF == 'Foxo1' & ms20_bias50$SignedQuantile == 1, 'Target']
ms20_bias50_Foxo1_neg = ms20_bias50[ms20_bias50$TF == 'Foxo1' & ms20_bias50$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms20_bias50_Foxo1_pos_list = ifelse(pool_genes %in% ms20_bias50_Foxo1_pos, 1, 0)

xtab <- table(ms20_bias50_Foxo1_pos_list, targets_cd8_pool_genes_list)
cm_ms20_bias50 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20_bias50$overall['Accuracy']
cm_ms20_bias50$byClass['Precision']
cm_ms20_bias50$byClass['Recall']
cm_ms20_bias50$byClass['F1']
cm_ms20_bias50$byClass['Balanced Accuracy']



# library(Metrics)
# Metrics::recall(ms15_bias100_Foxo1_pos_list, targets_cd8_pool_genes_list)
# Metrics::precision(ms15_bias100_Foxo1_pos_list, targets_cd8_pool_genes_list)
# Metrics::f1(ms15_bias100_Foxo1_pos_list, targets_cd8_pool_genes_list) # wrong somehow


# make plots
forplot = data.frame(matrix(ncol = 3, nrow = 12))
colnames(forplot) = c('Accuracy', 'Balanced Accuracy', 'Precision')
rownames(forplot) = c('ms5_bias5', 'ms5_bias100', 'ms5_bias50', 'ms10_bias5', 'ms10_bias100', 'ms10_bias50', 
                      'ms15_bias5', 'ms15_bias100', 'ms15_bias50', 'ms20_bias5', 'ms20_bias100', 'ms20_bias50')
forplot['ms10_bias5', 'Accuracy'] = cm_ms10_bias5$overall['Accuracy']
forplot['ms10_bias100', 'Accuracy'] = cm_ms10_bias100$overall['Accuracy']
forplot['ms10_bias50', 'Accuracy'] = cm_ms10_bias50$overall['Accuracy']
forplot['ms15_bias5', 'Accuracy'] = cm_ms15_bias5$overall['Accuracy']
forplot['ms15_bias100', 'Accuracy'] = cm_ms15_bias100$overall['Accuracy']
forplot['ms15_bias50', 'Accuracy'] = cm_ms15_bias50$overall['Accuracy']

forplot['ms10_bias5', 'Balanced Accuracy'] = cm_ms10_bias5$byClass['Balanced Accuracy']
forplot['ms10_bias100', 'Balanced Accuracy'] = cm_ms10_bias100$byClass['Balanced Accuracy']
forplot['ms10_bias50', 'Balanced Accuracy'] = cm_ms10_bias50$byClass['Balanced Accuracy']
forplot['ms15_bias5', 'Balanced Accuracy'] = cm_ms15_bias5$byClass['Balanced Accuracy']
forplot['ms15_bias100', 'Balanced Accuracy'] = cm_ms15_bias100$byClass['Balanced Accuracy']
forplot['ms15_bias50', 'Balanced Accuracy'] = cm_ms15_bias50$byClass['Balanced Accuracy']

forplot['ms10_bias5', 'Precision'] = cm_ms10_bias5$byClass['Precision']
forplot['ms10_bias100', 'Precision'] = cm_ms10_bias100$byClass['Precision']
forplot['ms10_bias50', 'Precision'] = cm_ms10_bias50$byClass['Precision']
forplot['ms15_bias5', 'Precision'] = cm_ms15_bias5$byClass['Precision']
forplot['ms15_bias100', 'Precision'] = cm_ms15_bias100$byClass['Precision']
forplot['ms15_bias50', 'Precision'] = cm_ms15_bias50$byClass['Precision']


forplot['ms20_bias5', 'Accuracy'] = cm_ms20_bias5$overall['Accuracy']
forplot['ms20_bias100', 'Accuracy'] = cm_ms20_bias100$overall['Accuracy']
forplot['ms20_bias50', 'Accuracy'] = cm_ms20_bias50$overall['Accuracy']
forplot['ms5_bias5', 'Accuracy'] = cm_ms5_bias5$overall['Accuracy']
forplot['ms5_bias100', 'Accuracy'] = cm_ms5_bias100$overall['Accuracy']
forplot['ms5_bias50', 'Accuracy'] = cm_ms5_bias50$overall['Accuracy']

forplot['ms20_bias5', 'Balanced Accuracy'] = cm_ms20_bias5$byClass['Balanced Accuracy']
forplot['ms20_bias100', 'Balanced Accuracy'] = cm_ms20_bias100$byClass['Balanced Accuracy']
forplot['ms20_bias50', 'Balanced Accuracy'] = cm_ms20_bias50$byClass['Balanced Accuracy']
forplot['ms5_bias5', 'Balanced Accuracy'] = cm_ms5_bias5$byClass['Balanced Accuracy']
forplot['ms5_bias100', 'Balanced Accuracy'] = cm_ms5_bias100$byClass['Balanced Accuracy']
forplot['ms5_bias50', 'Balanced Accuracy'] = cm_ms5_bias50$byClass['Balanced Accuracy']

forplot['ms20_bias5', 'Precision'] = cm_ms20_bias5$byClass['Precision']
forplot['ms20_bias100', 'Precision'] = cm_ms20_bias100$byClass['Precision']
forplot['ms20_bias50', 'Precision'] = cm_ms20_bias50$byClass['Precision']
forplot['ms5_bias5', 'Precision'] = cm_ms5_bias5$byClass['Precision']
forplot['ms5_bias100', 'Precision'] = cm_ms5_bias100$byClass['Precision']
forplot['ms5_bias50', 'Precision'] = cm_ms5_bias50$byClass['Precision']

forplot$ms = c('Model size: 5', 'Model size: 5', 'Model size: 5', 
               'Model size: 10', 'Model size: 10', 'Model size: 10', 
               'Model size: 15', 'Model size: 15', 'Model size: 15',
               'Model size: 20', 'Model size: 20', 'Model size: 20')
forplot$ms = factor(forplot$ms, levels = c('Model size: 5', 'Model size: 10', 'Model size: 15', 'Model size: 20'))
forplot$TFAOpt = c('bias5', 'bias100', 'bias50', 'bias5', 'bias100', 'bias50', 'bias5', 'bias100', 'bias50', 'bias5', 'bias100', 'bias50')
forplot$TFAOpt = factor(forplot$TFAOpt, levels = c('bias5', 'bias50', 'bias100'))

forplot_melt = reshape2::melt(forplot)
library(ggplot2)
# library(RColorBrewer)
# brewer.pal(3, 'Set3')
foxo1_pos_plot = ggplot(data=forplot_melt, aes(x=variable, y=value, fill=TFAOpt)) +
  geom_bar(stat="identity", position=position_dodge(), color="black", width = 0.7)+
  facet_wrap(~ms) + 
  scale_fill_manual(limits = c('bias5', 'bias50', 'bias100'),
                    labels = c('bias5', 'bias50', 'bias100'),
                    name = 'Bias',
                    values = c("#8DD3C7", "#FFFFB3", "#BEBADA"))+
  theme_bw(base_size = 14) +
  xlab('') +
  ylab('Value')+ 
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle('Foxo1 positive targets')




# negative targets

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms15_bias100_Foxo1_neg_list = ifelse(pool_genes %in% ms15_bias100_Foxo1_neg, 1, 0)

xtab <- table(ms15_bias100_Foxo1_neg_list, targets_cd8_pool_genes_list)
#cm <- caret::confusionMatrix(xtab)
cm_ms15_bias100 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15_bias100$overall['Accuracy']
cm_ms15_bias100$byClass['Precision']
cm_ms15_bias100$byClass['Recall']
cm_ms15_bias100$byClass['F1']
cm_ms15_bias100$byClass['Balanced Accuracy']


ms15_bias5_Foxo1_neg = ms15_bias5[ms15_bias5$TF == 'Foxo1' & ms15_bias5$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms15_bias5_Foxo1_neg_list = ifelse(pool_genes %in% ms15_bias5_Foxo1_neg, 1, 0)

xtab <- table(ms15_bias5_Foxo1_neg_list, targets_cd8_pool_genes_list)
cm_ms15_bias5 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15_bias5$overall['Accuracy']
cm_ms15_bias5$byClass['Precision']
cm_ms15_bias5$byClass['Recall']
cm_ms15_bias5$byClass['F1']
cm_ms15_bias5$byClass['Balanced Accuracy']


ms15_bias50_Foxo1_neg = ms15_bias50[ms15_bias50$TF == 'Foxo1' & ms15_bias50$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms15_bias50_Foxo1_neg_list = ifelse(pool_genes %in% ms15_bias50_Foxo1_neg, 1, 0)

xtab <- table(ms15_bias50_Foxo1_neg_list, targets_cd8_pool_genes_list)
cm_ms15_bias50 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15_bias50$overall['Accuracy']
cm_ms15_bias50$byClass['Precision']
cm_ms15_bias50$byClass['Recall']
cm_ms15_bias50$byClass['F1']
cm_ms15_bias50$byClass['Balanced Accuracy']



ms10_bias100_Foxo1_neg = ms10_bias100[ms10_bias100$TF == 'Foxo1' & ms10_bias100$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms10_bias100_Foxo1_neg_list = ifelse(pool_genes %in% ms10_bias100_Foxo1_neg, 1, 0)

xtab <- table(ms10_bias100_Foxo1_neg_list, targets_cd8_pool_genes_list)
cm_ms10_bias100 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10_bias100$overall['Accuracy']
cm_ms10_bias100$byClass['Precision']
cm_ms10_bias100$byClass['Recall']
cm_ms10_bias100$byClass['F1']
cm_ms10_bias100$byClass['Balanced Accuracy']


ms10_bias5_Foxo1_neg = ms10_bias5[ms10_bias5$TF == 'Foxo1' & ms10_bias5$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms10_bias5_Foxo1_neg_list = ifelse(pool_genes %in% ms10_bias5_Foxo1_neg, 1, 0)

xtab <- table(ms10_bias5_Foxo1_neg_list, targets_cd8_pool_genes_list)
cm_ms10_bias5 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10_bias5$overall['Accuracy']
cm_ms10_bias5$byClass['Precision']
cm_ms10_bias5$byClass['Recall']
cm_ms10_bias5$byClass['F1']
cm_ms10_bias5$byClass['Balanced Accuracy']


ms10_bias50_Foxo1_neg = ms10_bias50[ms10_bias50$TF == 'Foxo1' & ms10_bias50$SignedQuantile == 1, 'Target']
ms10_bias50_Foxo1_neg = ms10_bias50[ms10_bias50$TF == 'Foxo1' & ms10_bias50$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms10_bias50_Foxo1_neg_list = ifelse(pool_genes %in% ms10_bias50_Foxo1_neg, 1, 0)

xtab <- table(ms10_bias50_Foxo1_neg_list, targets_cd8_pool_genes_list)
cm_ms10_bias50 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10_bias50$overall['Accuracy']
cm_ms10_bias50$byClass['Precision']
cm_ms10_bias50$byClass['Recall']
cm_ms10_bias50$byClass['F1']
cm_ms10_bias50$byClass['Balanced Accuracy']



targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms20_bias100_Foxo1_neg_list = ifelse(pool_genes %in% ms20_bias100_Foxo1_neg, 1, 0)

xtab <- table(ms20_bias100_Foxo1_neg_list, targets_cd8_pool_genes_list)
#cm <- caret::confusionMatrix(xtab)
cm_ms20_bias100 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20_bias100$overall['Accuracy']
cm_ms20_bias100$byClass['Precision']
cm_ms20_bias100$byClass['Recall']
cm_ms20_bias100$byClass['F1']
cm_ms20_bias100$byClass['Balanced Accuracy']


ms20_bias5_Foxo1_neg = ms20_bias5[ms20_bias5$TF == 'Foxo1' & ms20_bias5$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms20_bias5_Foxo1_neg_list = ifelse(pool_genes %in% ms20_bias5_Foxo1_neg, 1, 0)

xtab <- table(ms20_bias5_Foxo1_neg_list, targets_cd8_pool_genes_list)
cm_ms20_bias5 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20_bias5$overall['Accuracy']
cm_ms20_bias5$byClass['Precision']
cm_ms20_bias5$byClass['Recall']
cm_ms20_bias5$byClass['F1']
cm_ms20_bias5$byClass['Balanced Accuracy']


ms20_bias50_Foxo1_neg = ms20_bias50[ms20_bias50$TF == 'Foxo1' & ms20_bias50$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms20_bias50_Foxo1_neg_list = ifelse(pool_genes %in% ms20_bias50_Foxo1_neg, 1, 0)

xtab <- table(ms20_bias50_Foxo1_neg_list, targets_cd8_pool_genes_list)
cm_ms20_bias50 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20_bias50$overall['Accuracy']
cm_ms20_bias50$byClass['Precision']
cm_ms20_bias50$byClass['Recall']
cm_ms20_bias50$byClass['F1']
cm_ms20_bias50$byClass['Balanced Accuracy']



ms5_bias100_Foxo1_neg = ms5_bias100[ms5_bias100$TF == 'Foxo1' & ms5_bias100$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms5_bias100_Foxo1_neg_list = ifelse(pool_genes %in% ms5_bias100_Foxo1_neg, 1, 0)

xtab <- table(ms5_bias100_Foxo1_neg_list, targets_cd8_pool_genes_list)
cm_ms5_bias100 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5_bias100$overall['Accuracy']
cm_ms5_bias100$byClass['Precision']
cm_ms5_bias100$byClass['Recall']
cm_ms5_bias100$byClass['F1']
cm_ms5_bias100$byClass['Balanced Accuracy']


ms5_bias5_Foxo1_neg = ms5_bias5[ms5_bias5$TF == 'Foxo1' & ms5_bias5$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms5_bias5_Foxo1_neg_list = ifelse(pool_genes %in% ms5_bias5_Foxo1_neg, 1, 0)

xtab <- table(ms5_bias5_Foxo1_neg_list, targets_cd8_pool_genes_list)
cm_ms5_bias5 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5_bias5$overall['Accuracy']
cm_ms5_bias5$byClass['Precision']
cm_ms5_bias5$byClass['Recall']
cm_ms5_bias5$byClass['F1']
cm_ms5_bias5$byClass['Balanced Accuracy']


ms5_bias50_Foxo1_neg = ms5_bias50[ms5_bias50$TF == 'Foxo1' & ms5_bias50$SignedQuantile == 1, 'Target']
ms5_bias50_Foxo1_neg = ms5_bias50[ms5_bias50$TF == 'Foxo1' & ms5_bias50$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms5_bias50_Foxo1_neg_list = ifelse(pool_genes %in% ms5_bias50_Foxo1_neg, 1, 0)

xtab <- table(ms5_bias50_Foxo1_neg_list, targets_cd8_pool_genes_list)
cm_ms5_bias50 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5_bias50$overall['Accuracy']
cm_ms5_bias50$byClass['Precision']
cm_ms5_bias50$byClass['Recall']
cm_ms5_bias50$byClass['F1']
cm_ms5_bias50$byClass['Balanced Accuracy']

# library(Metrics)
# Metrics::recall(ms15_bias100_Foxo1_neg_list, targets_cd8_pool_genes_list)
# Metrics::precision(ms15_bias100_Foxo1_neg_list, targets_cd8_pool_genes_list)
# Metrics::f1(ms15_bias100_Foxo1_neg_list, targets_cd8_pool_genes_list) # wrong somehow


# make plots
forplotneg = data.frame(matrix(ncol = 3, nrow = 12))
colnames(forplotneg) = c('Accuracy', 'Balanced Accuracy', 'Precision')
rownames(forplotneg) = c('ms5_bias5', 'ms5_bias100', 'ms5_bias50', 'ms10_bias5', 'ms10_bias100', 'ms10_bias50', 
                         'ms15_bias5', 'ms15_bias100', 'ms15_bias50', 'ms20_bias5', 'ms20_bias100', 'ms20_bias50')
forplotneg['ms10_bias5', 'Accuracy'] = cm_ms10_bias5$overall['Accuracy']
forplotneg['ms10_bias100', 'Accuracy'] = cm_ms10_bias100$overall['Accuracy']
forplotneg['ms10_bias50', 'Accuracy'] = cm_ms10_bias50$overall['Accuracy']
forplotneg['ms15_bias5', 'Accuracy'] = cm_ms15_bias5$overall['Accuracy']
forplotneg['ms15_bias100', 'Accuracy'] = cm_ms15_bias100$overall['Accuracy']
forplotneg['ms15_bias50', 'Accuracy'] = cm_ms15_bias50$overall['Accuracy']

forplotneg['ms10_bias5', 'Balanced Accuracy'] = cm_ms10_bias5$byClass['Balanced Accuracy']
forplotneg['ms10_bias100', 'Balanced Accuracy'] = cm_ms10_bias100$byClass['Balanced Accuracy']
forplotneg['ms10_bias50', 'Balanced Accuracy'] = cm_ms10_bias50$byClass['Balanced Accuracy']
forplotneg['ms15_bias5', 'Balanced Accuracy'] = cm_ms15_bias5$byClass['Balanced Accuracy']
forplotneg['ms15_bias100', 'Balanced Accuracy'] = cm_ms15_bias100$byClass['Balanced Accuracy']
forplotneg['ms15_bias50', 'Balanced Accuracy'] = cm_ms15_bias50$byClass['Balanced Accuracy']

forplotneg['ms10_bias5', 'Precision'] = cm_ms10_bias5$byClass['Precision']
forplotneg['ms10_bias100', 'Precision'] = cm_ms10_bias100$byClass['Precision']
forplotneg['ms10_bias50', 'Precision'] = cm_ms10_bias50$byClass['Precision']
forplotneg['ms15_bias5', 'Precision'] = cm_ms15_bias5$byClass['Precision']
forplotneg['ms15_bias100', 'Precision'] = cm_ms15_bias100$byClass['Precision']
forplotneg['ms15_bias50', 'Precision'] = cm_ms15_bias50$byClass['Precision']


forplotneg['ms5_bias5', 'Accuracy'] = cm_ms5_bias5$overall['Accuracy']
forplotneg['ms5_bias100', 'Accuracy'] = cm_ms5_bias100$overall['Accuracy']
forplotneg['ms5_bias50', 'Accuracy'] = cm_ms5_bias50$overall['Accuracy']
forplotneg['ms20_bias5', 'Accuracy'] = cm_ms20_bias5$overall['Accuracy']
forplotneg['ms20_bias100', 'Accuracy'] = cm_ms20_bias100$overall['Accuracy']
forplotneg['ms20_bias50', 'Accuracy'] = cm_ms20_bias50$overall['Accuracy']

forplotneg['ms5_bias5', 'Balanced Accuracy'] = cm_ms5_bias5$byClass['Balanced Accuracy']
forplotneg['ms5_bias100', 'Balanced Accuracy'] = cm_ms5_bias100$byClass['Balanced Accuracy']
forplotneg['ms5_bias50', 'Balanced Accuracy'] = cm_ms5_bias50$byClass['Balanced Accuracy']
forplotneg['ms20_bias5', 'Balanced Accuracy'] = cm_ms20_bias5$byClass['Balanced Accuracy']
forplotneg['ms20_bias100', 'Balanced Accuracy'] = cm_ms20_bias100$byClass['Balanced Accuracy']
forplotneg['ms20_bias50', 'Balanced Accuracy'] = cm_ms20_bias50$byClass['Balanced Accuracy']

forplotneg['ms5_bias5', 'Precision'] = cm_ms5_bias5$byClass['Precision']
forplotneg['ms5_bias100', 'Precision'] = cm_ms5_bias100$byClass['Precision']
forplotneg['ms5_bias50', 'Precision'] = cm_ms5_bias50$byClass['Precision']
forplotneg['ms20_bias5', 'Precision'] = cm_ms20_bias5$byClass['Precision']
forplotneg['ms20_bias100', 'Precision'] = cm_ms20_bias100$byClass['Precision']
forplotneg['ms20_bias50', 'Precision'] = cm_ms20_bias50$byClass['Precision']

forplotneg$ms = c('Model size: 5', 'Model size: 5', 'Model size: 5', 
                  'Model size: 10', 'Model size: 10', 'Model size: 10', 
                  'Model size: 15', 'Model size: 15', 'Model size: 15',
                  'Model size: 20', 'Model size: 20', 'Model size: 20')
forplotneg$ms = factor(forplotneg$ms, levels = c('Model size: 5', 'Model size: 10', 'Model size: 15', 'Model size: 20'))
forplotneg$TFAOpt = c('bias5', 'bias100', 'bias50', 'bias5', 'bias100', 'bias50', 'bias5', 'bias100', 'bias50', 'bias5', 'bias100', 'bias50')
forplotneg$TFAOpt = factor(forplotneg$TFAOpt, levels = c('bias5', 'bias50', 'bias100'))

forplotneg_melt = reshape2::melt(forplotneg)
library(ggplot2)
# library(RColorBrewer)
# brewer.pal(3, 'Set3')
foxo1_neg_plot = ggplot(data=forplotneg_melt, aes(x=variable, y=value, fill=TFAOpt)) +
  geom_bar(stat="identity", position=position_dodge(), color="black", width = 0.7)+
  facet_wrap(~ms) + 
  scale_fill_manual(limits = c('bias5', 'bias50', 'bias100'),
                    labels = c('bias5', 'bias50', 'bias100'),
                    name = 'Bias',
                    values = c("#8DD3C7", "#FFFFB3", "#BEBADA"))+
  theme_bw(base_size = 14) +
  xlab('') +
  ylab('Value')+ 
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle('Foxo1 negative targets')

library(cowplot)
foxo1.combined = plot_grid(foxo1_pos_plot, foxo1_neg_plot, align="hv", ncol = 2)
ggsave(
  paste0('compare_ms_bias_Foxo1_posneg_remake.pdf'),
  plot = foxo1.combined,
  device = "pdf",
  width = 12,
  height = 6,
  dpi = 300
)


### testing TFA options

# Eomes

library(caret)
library(readxl)
pool_genes = read.table('deg_padj01.txt')$V1
de = read_excel("eomes/GSE124913_RNA-Seq_CD8SP_EomesTGvsWT_DifferentialAnalysisResults.xlsx", 
                sheet = "Sheet1")
de = data.frame(de)
de_noNA = de[complete.cases(de), ] #no use
de_noNA$FDR_numeric = as.numeric(de_noNA$FDR)
de_noNA = de_noNA[complete.cases(de_noNA), ]
de_noNA_pool_genes = de_noNA[de_noNA$Gene.Name %in% pool_genes,]
de_Eomes_pos = de_noNA_pool_genes[de_noNA_pool_genes$FDR_numeric < 0.001 & de_noNA_pool_genes$Log2FC.WT.vs.EomesTG < 0, 'Gene.Name']
de_Eomes_neg = de_noNA_pool_genes[de_noNA_pool_genes$FDR_numeric < 0.001 & de_noNA_pool_genes$Log2FC.WT.vs.EomesTG > 0, 'Gene.Name']


ms10_combine = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms10_bias50_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms10         = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms10_bias50_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms10_TFmRNA  = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms10_bias50_TFmRNA_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms15_combine = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms15_bias50_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms15         = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms15_bias50_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms15_TFmRNA  = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms15_bias50_TFmRNA_cut01_sharedbyMorethan2nets_sp.tsv', header = T)

ms5_combine  = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms5_bias50_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms5          = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms5_bias50_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms5_TFmRNA   = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms5_bias50_TFmRNA_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms20_combine = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms20_bias50_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms20         = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms20_bias50_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms20_TFmRNA  = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms20_bias50_TFmRNA_cut01_sharedbyMorethan2nets_sp.tsv', header = T)

ms15_TFmRNA_Eomes_pos = ms15_TFmRNA[ms15_TFmRNA$TF == 'Eomes' & ms15_TFmRNA$SignedQuantile == 1, 'Target']
ms15_TFmRNA_Eomes_neg = ms15_TFmRNA[ms15_TFmRNA$TF == 'Eomes' & ms15_TFmRNA$SignedQuantile == -1, 'Target']

#table(ms15_TFmRNA_Eomes_pos %in% pool_genes)
#table(ms15_TFmRNA_Eomes_neg %in% pool_genes)
#table(de_Eomes_pos %in% pool_genes)
#table(de_Eomes_neg %in% pool_genes)

de_Eomes_pos_list = ifelse(pool_genes %in% de_Eomes_pos, 1, 0)
ms15_TFmRNA_Eomes_pos_list = ifelse(pool_genes %in% ms15_TFmRNA_Eomes_pos, 1, 0)

xtab <- table(ms15_TFmRNA_Eomes_pos_list, de_Eomes_pos_list)
#cm <- caret::confusionMatrix(xtab)
cm_ms15_TFmRNA <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15_TFmRNA$overall['Accuracy']
cm_ms15_TFmRNA$byClass['Precision']
cm_ms15_TFmRNA$byClass['Recall']
cm_ms15_TFmRNA$byClass['F1']
cm_ms15_TFmRNA$byClass['Balanced Accuracy']


ms15_Eomes_pos = ms15[ms15$TF == 'Eomes' & ms15$SignedQuantile == 1, 'Target']
ms15_Eomes_neg = ms15[ms15$TF == 'Eomes' & ms15$SignedQuantile == -1, 'Target']

de_Eomes_pos_list = ifelse(pool_genes %in% de_Eomes_pos, 1, 0)
ms15_Eomes_pos_list = ifelse(pool_genes %in% ms15_Eomes_pos, 1, 0)

xtab <- table(ms15_Eomes_pos_list, de_Eomes_pos_list)
cm_ms15 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15$overall['Accuracy']
cm_ms15$byClass['Precision']
cm_ms15$byClass['Recall']
cm_ms15$byClass['F1']
cm_ms15$byClass['Balanced Accuracy']

ms15_combine_Eomes_pos = ms15_combine[ms15_combine$TF == 'Eomes' & ms15_combine$SignedQuantile == 1, 'Target']
ms15_combine_Eomes_neg = ms15_combine[ms15_combine$TF == 'Eomes' & ms15_combine$SignedQuantile == -1, 'Target']

de_Eomes_pos_list = ifelse(pool_genes %in% de_Eomes_pos, 1, 0)
ms15_combine_Eomes_pos_list = ifelse(pool_genes %in% ms15_combine_Eomes_pos, 1, 0)

xtab <- table(ms15_combine_Eomes_pos_list, de_Eomes_pos_list)
cm_ms15_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15_combine$overall['Accuracy']
cm_ms15_combine$byClass['Precision']
cm_ms15_combine$byClass['Recall']
cm_ms15_combine$byClass['F1']
cm_ms15_combine$byClass['Balanced Accuracy']



ms10_TFmRNA_Eomes_pos = ms10_TFmRNA[ms10_TFmRNA$TF == 'Eomes' & ms10_TFmRNA$SignedQuantile == 1, 'Target']
ms10_TFmRNA_Eomes_neg = ms10_TFmRNA[ms10_TFmRNA$TF == 'Eomes' & ms10_TFmRNA$SignedQuantile == -1, 'Target']

de_Eomes_pos_list = ifelse(pool_genes %in% de_Eomes_pos, 1, 0)
ms10_TFmRNA_Eomes_pos_list = ifelse(pool_genes %in% ms10_TFmRNA_Eomes_pos, 1, 0)

xtab <- table(ms10_TFmRNA_Eomes_pos_list, de_Eomes_pos_list)
cm_ms10_TFmRNA <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10_TFmRNA$overall['Accuracy']
cm_ms10_TFmRNA$byClass['Precision']
cm_ms10_TFmRNA$byClass['Recall']
cm_ms10_TFmRNA$byClass['F1']
cm_ms10_TFmRNA$byClass['Balanced Accuracy']


ms10_Eomes_pos = ms10[ms10$TF == 'Eomes' & ms10$SignedQuantile == 1, 'Target']
ms10_Eomes_neg = ms10[ms10$TF == 'Eomes' & ms10$SignedQuantile == -1, 'Target']

de_Eomes_pos_list = ifelse(pool_genes %in% de_Eomes_pos, 1, 0)
ms10_Eomes_pos_list = ifelse(pool_genes %in% ms10_Eomes_pos, 1, 0)

xtab <- table(ms10_Eomes_pos_list, de_Eomes_pos_list)
cm_ms10 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10$overall['Accuracy']
cm_ms10$byClass['Precision']
cm_ms10$byClass['Recall']
cm_ms10$byClass['F1']
cm_ms10$byClass['Balanced Accuracy']


ms10_combine_Eomes_pos = ms10_combine[ms10_combine$TF == 'Eomes' & ms10_combine$SignedQuantile == 1, 'Target']
ms10_combine_Eomes_neg = ms10_combine[ms10_combine$TF == 'Eomes' & ms10_combine$SignedQuantile == -1, 'Target']

de_Eomes_pos_list = ifelse(pool_genes %in% de_Eomes_pos, 1, 0)
ms10_combine_Eomes_pos_list = ifelse(pool_genes %in% ms10_combine_Eomes_pos, 1, 0)

xtab <- table(ms10_combine_Eomes_pos_list, de_Eomes_pos_list)
cm_ms10_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10_combine$overall['Accuracy']
cm_ms10_combine$byClass['Precision']
cm_ms10_combine$byClass['Recall']
cm_ms10_combine$byClass['F1']
cm_ms10_combine$byClass['Balanced Accuracy']


ms5_TFmRNA_Eomes_pos = ms5_TFmRNA[ms5_TFmRNA$TF == 'Eomes' & ms5_TFmRNA$SignedQuantile == 1, 'Target']
ms5_TFmRNA_Eomes_neg = ms5_TFmRNA[ms5_TFmRNA$TF == 'Eomes' & ms5_TFmRNA$SignedQuantile == -1, 'Target']

#table(ms5_TFmRNA_Eomes_pos %in% pool_genes)
#table(ms5_TFmRNA_Eomes_neg %in% pool_genes)
#table(de_Eomes_pos %in% pool_genes)
#table(de_Eomes_neg %in% pool_genes)

de_Eomes_pos_list = ifelse(pool_genes %in% de_Eomes_pos, 1, 0)
ms5_TFmRNA_Eomes_pos_list = ifelse(pool_genes %in% ms5_TFmRNA_Eomes_pos, 1, 0)

xtab <- table(ms5_TFmRNA_Eomes_pos_list, de_Eomes_pos_list)
#cm <- caret::confusionMatrix(xtab)
cm_ms5_TFmRNA <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5_TFmRNA$overall['Accuracy']
cm_ms5_TFmRNA$byClass['Precision']
cm_ms5_TFmRNA$byClass['Recall']
cm_ms5_TFmRNA$byClass['F1']
cm_ms5_TFmRNA$byClass['Balanced Accuracy']


ms5_Eomes_pos = ms5[ms5$TF == 'Eomes' & ms5$SignedQuantile == 1, 'Target']
ms5_Eomes_neg = ms5[ms5$TF == 'Eomes' & ms5$SignedQuantile == -1, 'Target']

de_Eomes_pos_list = ifelse(pool_genes %in% de_Eomes_pos, 1, 0)
ms5_Eomes_pos_list = ifelse(pool_genes %in% ms5_Eomes_pos, 1, 0)

xtab <- table(ms5_Eomes_pos_list, de_Eomes_pos_list)
cm_ms5 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5$overall['Accuracy']
cm_ms5$byClass['Precision']
cm_ms5$byClass['Recall']
cm_ms5$byClass['F1']
cm_ms5$byClass['Balanced Accuracy']

ms5_combine_Eomes_pos = ms5_combine[ms5_combine$TF == 'Eomes' & ms5_combine$SignedQuantile == 1, 'Target']
ms5_combine_Eomes_neg = ms5_combine[ms5_combine$TF == 'Eomes' & ms5_combine$SignedQuantile == -1, 'Target']

de_Eomes_pos_list = ifelse(pool_genes %in% de_Eomes_pos, 1, 0)
ms5_combine_Eomes_pos_list = ifelse(pool_genes %in% ms5_combine_Eomes_pos, 1, 0)

xtab <- table(ms5_combine_Eomes_pos_list, de_Eomes_pos_list)
cm_ms5_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5_combine$overall['Accuracy']
cm_ms5_combine$byClass['Precision']
cm_ms5_combine$byClass['Recall']
cm_ms5_combine$byClass['F1']
cm_ms5_combine$byClass['Balanced Accuracy']



ms20_TFmRNA_Eomes_pos = ms20_TFmRNA[ms20_TFmRNA$TF == 'Eomes' & ms20_TFmRNA$SignedQuantile == 1, 'Target']
ms20_TFmRNA_Eomes_neg = ms20_TFmRNA[ms20_TFmRNA$TF == 'Eomes' & ms20_TFmRNA$SignedQuantile == -1, 'Target']

de_Eomes_pos_list = ifelse(pool_genes %in% de_Eomes_pos, 1, 0)
ms20_TFmRNA_Eomes_pos_list = ifelse(pool_genes %in% ms20_TFmRNA_Eomes_pos, 1, 0)

xtab <- table(ms20_TFmRNA_Eomes_pos_list, de_Eomes_pos_list)
cm_ms20_TFmRNA <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20_TFmRNA$overall['Accuracy']
cm_ms20_TFmRNA$byClass['Precision']
cm_ms20_TFmRNA$byClass['Recall']
cm_ms20_TFmRNA$byClass['F1']
cm_ms20_TFmRNA$byClass['Balanced Accuracy']


ms20_Eomes_pos = ms20[ms20$TF == 'Eomes' & ms20$SignedQuantile == 1, 'Target']
ms20_Eomes_neg = ms20[ms20$TF == 'Eomes' & ms20$SignedQuantile == -1, 'Target']

de_Eomes_pos_list = ifelse(pool_genes %in% de_Eomes_pos, 1, 0)
ms20_Eomes_pos_list = ifelse(pool_genes %in% ms20_Eomes_pos, 1, 0)

xtab <- table(ms20_Eomes_pos_list, de_Eomes_pos_list)
cm_ms20 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20$overall['Accuracy']
cm_ms20$byClass['Precision']
cm_ms20$byClass['Recall']
cm_ms20$byClass['F1']
cm_ms20$byClass['Balanced Accuracy']


ms20_combine_Eomes_pos = ms20_combine[ms20_combine$TF == 'Eomes' & ms20_combine$SignedQuantile == 1, 'Target']
ms20_combine_Eomes_neg = ms20_combine[ms20_combine$TF == 'Eomes' & ms20_combine$SignedQuantile == -1, 'Target']

de_Eomes_pos_list = ifelse(pool_genes %in% de_Eomes_pos, 1, 0)
ms20_combine_Eomes_pos_list = ifelse(pool_genes %in% ms20_combine_Eomes_pos, 1, 0)

xtab <- table(ms20_combine_Eomes_pos_list, de_Eomes_pos_list)
cm_ms20_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20_combine$overall['Accuracy']
cm_ms20_combine$byClass['Precision']
cm_ms20_combine$byClass['Recall']
cm_ms20_combine$byClass['F1']
cm_ms20_combine$byClass['Balanced Accuracy']


# library(Metrics)
# Metrics::recall(ms15_TFmRNA_Eomes_pos_list, de_Eomes_pos_list)
# Metrics::precision(ms15_TFmRNA_Eomes_pos_list, de_Eomes_pos_list)
# Metrics::f1(ms15_TFmRNA_Eomes_pos_list, de_Eomes_pos_list) # wrong somehow


# make plots
forplot = data.frame(matrix(ncol = 3, nrow = 12))
colnames(forplot) = c('Accuracy', 'Balanced Accuracy', 'Precision')
rownames(forplot) = c('ms5', 'ms5_TFmRNA', 'ms5_combine', 
                      'ms10', 'ms10_TFmRNA', 'ms10_combine', 
                      'ms15', 'ms15_TFmRNA', 'ms15_combine',
                      'ms20', 'ms20_TFmRNA', 'ms20_combine')
forplot['ms10', 'Accuracy'] = cm_ms10$overall['Accuracy']
forplot['ms10_TFmRNA', 'Accuracy'] = cm_ms10_TFmRNA$overall['Accuracy']
forplot['ms10_combine', 'Accuracy'] = cm_ms10_combine$overall['Accuracy']
forplot['ms15', 'Accuracy'] = cm_ms15$overall['Accuracy']
forplot['ms15_TFmRNA', 'Accuracy'] = cm_ms15_TFmRNA$overall['Accuracy']
forplot['ms15_combine', 'Accuracy'] = cm_ms15_combine$overall['Accuracy']

forplot['ms10', 'Balanced Accuracy'] = cm_ms10$byClass['Balanced Accuracy']
forplot['ms10_TFmRNA', 'Balanced Accuracy'] = cm_ms10_TFmRNA$byClass['Balanced Accuracy']
forplot['ms10_combine', 'Balanced Accuracy'] = cm_ms10_combine$byClass['Balanced Accuracy']
forplot['ms15', 'Balanced Accuracy'] = cm_ms15$byClass['Balanced Accuracy']
forplot['ms15_TFmRNA', 'Balanced Accuracy'] = cm_ms15_TFmRNA$byClass['Balanced Accuracy']
forplot['ms15_combine', 'Balanced Accuracy'] = cm_ms15_combine$byClass['Balanced Accuracy']

forplot['ms10', 'Precision'] = cm_ms10$byClass['Precision']
forplot['ms10_TFmRNA', 'Precision'] = cm_ms10_TFmRNA$byClass['Precision']
forplot['ms10_combine', 'Precision'] = cm_ms10_combine$byClass['Precision']
forplot['ms15', 'Precision'] = cm_ms15$byClass['Precision']
forplot['ms15_TFmRNA', 'Precision'] = cm_ms15_TFmRNA$byClass['Precision']
forplot['ms15_combine', 'Precision'] = cm_ms15_combine$byClass['Precision']


forplot['ms20', 'Accuracy'] = cm_ms20$overall['Accuracy']
forplot['ms20_TFmRNA', 'Accuracy'] = cm_ms20_TFmRNA$overall['Accuracy']
forplot['ms20_combine', 'Accuracy'] = cm_ms20_combine$overall['Accuracy']
forplot['ms5', 'Accuracy'] = cm_ms5$overall['Accuracy']
forplot['ms5_TFmRNA', 'Accuracy'] = cm_ms5_TFmRNA$overall['Accuracy']
forplot['ms5_combine', 'Accuracy'] = cm_ms5_combine$overall['Accuracy']

forplot['ms20', 'Balanced Accuracy'] = cm_ms20$byClass['Balanced Accuracy']
forplot['ms20_TFmRNA', 'Balanced Accuracy'] = cm_ms20_TFmRNA$byClass['Balanced Accuracy']
forplot['ms20_combine', 'Balanced Accuracy'] = cm_ms20_combine$byClass['Balanced Accuracy']
forplot['ms5', 'Balanced Accuracy'] = cm_ms5$byClass['Balanced Accuracy']
forplot['ms5_TFmRNA', 'Balanced Accuracy'] = cm_ms5_TFmRNA$byClass['Balanced Accuracy']
forplot['ms5_combine', 'Balanced Accuracy'] = cm_ms5_combine$byClass['Balanced Accuracy']

forplot['ms20', 'Precision'] = cm_ms20$byClass['Precision']
forplot['ms20_TFmRNA', 'Precision'] = cm_ms20_TFmRNA$byClass['Precision']
forplot['ms20_combine', 'Precision'] = cm_ms20_combine$byClass['Precision']
forplot['ms5', 'Precision'] = cm_ms5$byClass['Precision']
forplot['ms5_TFmRNA', 'Precision'] = cm_ms5_TFmRNA$byClass['Precision']
forplot['ms5_combine', 'Precision'] = cm_ms5_combine$byClass['Precision']

forplot$ms = c('Model size: 5', 'Model size: 5', 'Model size: 5', 
               'Model size: 10', 'Model size: 10', 'Model size: 10', 
               'Model size: 15', 'Model size: 15', 'Model size: 15',
               'Model size: 20', 'Model size: 20', 'Model size: 20')
forplot$ms = factor(forplot$ms, levels = c('Model size: 5', 'Model size: 10', 'Model size: 15', 'Model size: 20'))
forplot$TFAOpt = c('', 'TFmRNA', 'comb', '', 'TFmRNA', 'comb', '', 'TFmRNA', 'comb', '', 'TFmRNA', 'comb')
forplot$TFAOpt = factor(forplot$TFAOpt, levels = c('', 'TFmRNA', 'comb'))

forplot_melt = reshape2::melt(forplot)
library(ggplot2)
# library(RColorBrewer)
# brewer.pal(3, 'Set3')
eomes_pos_plot = ggplot(data=forplot_melt, aes(x=variable, y=value, fill=TFAOpt)) +
  geom_bar(stat="identity", position=position_dodge(), color="black", width = 0.7)+
  facet_wrap(~ms) + 
  scale_fill_manual(limits = c('', 'TFmRNA', 'comb'),
                    labels = c('TFA', 'TFmRNA', 'combine'),
                    name = 'TFA options',
                    values = c("#8DD3C7", "#FFFFB3", "#BEBADA"))+
  theme_bw(base_size = 14) +
  xlab('') +
  ylab('Value')+ 
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle('Eomes positive targets')




# negative targets

de_Eomes_neg_list = ifelse(pool_genes %in% de_Eomes_neg, 1, 0)
ms15_TFmRNA_Eomes_neg_list = ifelse(pool_genes %in% ms15_TFmRNA_Eomes_neg, 1, 0)

xtab <- table(ms15_TFmRNA_Eomes_neg_list, de_Eomes_neg_list)
#cm <- caret::confusionMatrix(xtab)
cm_ms15_TFmRNA <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15_TFmRNA$overall['Accuracy']
cm_ms15_TFmRNA$byClass['Precision']
cm_ms15_TFmRNA$byClass['Recall']
cm_ms15_TFmRNA$byClass['F1']
cm_ms15_TFmRNA$byClass['Balanced Accuracy']


ms15_Eomes_neg = ms15[ms15$TF == 'Eomes' & ms15$SignedQuantile == -1, 'Target']

de_Eomes_neg_list = ifelse(pool_genes %in% de_Eomes_neg, 1, 0)
ms15_Eomes_neg_list = ifelse(pool_genes %in% ms15_Eomes_neg, 1, 0)

xtab <- table(ms15_Eomes_neg_list, de_Eomes_neg_list)
cm_ms15 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15$overall['Accuracy']
cm_ms15$byClass['Precision']
cm_ms15$byClass['Recall']
cm_ms15$byClass['F1']
cm_ms15$byClass['Balanced Accuracy']

ms15_combine_Eomes_neg = ms15_combine[ms15_combine$TF == 'Eomes' & ms15_combine$SignedQuantile == -1, 'Target']

de_Eomes_neg_list = ifelse(pool_genes %in% de_Eomes_neg, 1, 0)
ms15_combine_Eomes_neg_list = ifelse(pool_genes %in% ms15_combine_Eomes_neg, 1, 0)

xtab <- table(ms15_combine_Eomes_neg_list, de_Eomes_neg_list)
cm_ms15_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15_combine$overall['Accuracy']
cm_ms15_combine$byClass['Precision']
cm_ms15_combine$byClass['Recall']
cm_ms15_combine$byClass['F1']
cm_ms15_combine$byClass['Balanced Accuracy']



ms10_TFmRNA_Eomes_neg = ms10_TFmRNA[ms10_TFmRNA$TF == 'Eomes' & ms10_TFmRNA$SignedQuantile == -1, 'Target']

de_Eomes_neg_list = ifelse(pool_genes %in% de_Eomes_neg, 1, 0)
ms10_TFmRNA_Eomes_neg_list = ifelse(pool_genes %in% ms10_TFmRNA_Eomes_neg, 1, 0)

xtab <- table(ms10_TFmRNA_Eomes_neg_list, de_Eomes_neg_list)
cm_ms10_TFmRNA <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10_TFmRNA$overall['Accuracy']
cm_ms10_TFmRNA$byClass['Precision']
cm_ms10_TFmRNA$byClass['Recall']
cm_ms10_TFmRNA$byClass['F1']
cm_ms10_TFmRNA$byClass['Balanced Accuracy']


ms10_Eomes_neg = ms10[ms10$TF == 'Eomes' & ms10$SignedQuantile == 1, 'Target']
ms10_Eomes_neg = ms10[ms10$TF == 'Eomes' & ms10$SignedQuantile == -1, 'Target']

de_Eomes_neg_list = ifelse(pool_genes %in% de_Eomes_neg, 1, 0)
ms10_Eomes_neg_list = ifelse(pool_genes %in% ms10_Eomes_neg, 1, 0)

xtab <- table(ms10_Eomes_neg_list, de_Eomes_neg_list)
cm_ms10 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10$overall['Accuracy']
cm_ms10$byClass['Precision']
cm_ms10$byClass['Recall']
cm_ms10$byClass['F1']
cm_ms10$byClass['Balanced Accuracy']


ms10_combine_Eomes_neg = ms10_combine[ms10_combine$TF == 'Eomes' & ms10_combine$SignedQuantile == -1, 'Target']

de_Eomes_neg_list = ifelse(pool_genes %in% de_Eomes_neg, 1, 0)
ms10_combine_Eomes_neg_list = ifelse(pool_genes %in% ms10_combine_Eomes_neg, 1, 0)

xtab <- table(ms10_combine_Eomes_neg_list, de_Eomes_neg_list)
cm_ms10_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10_combine$overall['Accuracy']
cm_ms10_combine$byClass['Precision']
cm_ms10_combine$byClass['Recall']
cm_ms10_combine$byClass['F1']
cm_ms10_combine$byClass['Balanced Accuracy']



de_Eomes_neg_list = ifelse(pool_genes %in% de_Eomes_neg, 1, 0)
ms5_TFmRNA_Eomes_neg_list = ifelse(pool_genes %in% ms5_TFmRNA_Eomes_neg, 1, 0)

xtab <- table(ms5_TFmRNA_Eomes_neg_list, de_Eomes_neg_list)
#cm <- caret::confusionMatrix(xtab)
cm_ms5_TFmRNA <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5_TFmRNA$overall['Accuracy']
cm_ms5_TFmRNA$byClass['Precision']
cm_ms5_TFmRNA$byClass['Recall']
cm_ms5_TFmRNA$byClass['F1']
cm_ms5_TFmRNA$byClass['Balanced Accuracy']


ms5_Eomes_neg = ms5[ms5$TF == 'Eomes' & ms5$SignedQuantile == -1, 'Target']

de_Eomes_neg_list = ifelse(pool_genes %in% de_Eomes_neg, 1, 0)
ms5_Eomes_neg_list = ifelse(pool_genes %in% ms5_Eomes_neg, 1, 0)

xtab <- table(ms5_Eomes_neg_list, de_Eomes_neg_list)
cm_ms5 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5$overall['Accuracy']
cm_ms5$byClass['Precision']
cm_ms5$byClass['Recall']
cm_ms5$byClass['F1']
cm_ms5$byClass['Balanced Accuracy']

ms5_combine_Eomes_neg = ms5_combine[ms5_combine$TF == 'Eomes' & ms5_combine$SignedQuantile == -1, 'Target']

de_Eomes_neg_list = ifelse(pool_genes %in% de_Eomes_neg, 1, 0)
ms5_combine_Eomes_neg_list = ifelse(pool_genes %in% ms5_combine_Eomes_neg, 1, 0)

xtab <- table(ms5_combine_Eomes_neg_list, de_Eomes_neg_list)
cm_ms5_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5_combine$overall['Accuracy']
cm_ms5_combine$byClass['Precision']
cm_ms5_combine$byClass['Recall']
cm_ms5_combine$byClass['F1']
cm_ms5_combine$byClass['Balanced Accuracy']



ms20_TFmRNA_Eomes_neg = ms20_TFmRNA[ms20_TFmRNA$TF == 'Eomes' & ms20_TFmRNA$SignedQuantile == -1, 'Target']

de_Eomes_neg_list = ifelse(pool_genes %in% de_Eomes_neg, 1, 0)
ms20_TFmRNA_Eomes_neg_list = ifelse(pool_genes %in% ms20_TFmRNA_Eomes_neg, 1, 0)

xtab <- table(ms20_TFmRNA_Eomes_neg_list, de_Eomes_neg_list)
cm_ms20_TFmRNA <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20_TFmRNA$overall['Accuracy']
cm_ms20_TFmRNA$byClass['Precision']
cm_ms20_TFmRNA$byClass['Recall']
cm_ms20_TFmRNA$byClass['F1']
cm_ms20_TFmRNA$byClass['Balanced Accuracy']


ms20_Eomes_neg = ms20[ms20$TF == 'Eomes' & ms20$SignedQuantile == 1, 'Target']
ms20_Eomes_neg = ms20[ms20$TF == 'Eomes' & ms20$SignedQuantile == -1, 'Target']

de_Eomes_neg_list = ifelse(pool_genes %in% de_Eomes_neg, 1, 0)
ms20_Eomes_neg_list = ifelse(pool_genes %in% ms20_Eomes_neg, 1, 0)

xtab <- table(ms20_Eomes_neg_list, de_Eomes_neg_list)
cm_ms20 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20$overall['Accuracy']
cm_ms20$byClass['Precision']
cm_ms20$byClass['Recall']
cm_ms20$byClass['F1']
cm_ms20$byClass['Balanced Accuracy']


ms20_combine_Eomes_neg = ms20_combine[ms20_combine$TF == 'Eomes' & ms20_combine$SignedQuantile == -1, 'Target']

de_Eomes_neg_list = ifelse(pool_genes %in% de_Eomes_neg, 1, 0)
ms20_combine_Eomes_neg_list = ifelse(pool_genes %in% ms20_combine_Eomes_neg, 1, 0)

xtab <- table(ms20_combine_Eomes_neg_list, de_Eomes_neg_list)
cm_ms20_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20_combine$overall['Accuracy']
cm_ms20_combine$byClass['Precision']
cm_ms20_combine$byClass['Recall']
cm_ms20_combine$byClass['F1']
cm_ms20_combine$byClass['Balanced Accuracy']

# library(Metrics)
# Metrics::recall(ms15_TFmRNA_Eomes_neg_list, de_Eomes_neg_list)
# Metrics::precision(ms15_TFmRNA_Eomes_neg_list, de_Eomes_neg_list)
# Metrics::f1(ms15_TFmRNA_Eomes_neg_list, de_Eomes_neg_list) # wrong somehow


# make plots
forplotneg = data.frame(matrix(ncol = 3, nrow = 12))
colnames(forplotneg) = c('Accuracy', 'Balanced Accuracy', 'Precision')
rownames(forplotneg) = c('ms5', 'ms5_TFmRNA', 'ms5_combine', 
                         'ms10', 'ms10_TFmRNA', 'ms10_combine', 
                         'ms15', 'ms15_TFmRNA', 'ms15_combine',
                         'ms20', 'ms20_TFmRNA', 'ms20_combine')
forplotneg['ms10', 'Accuracy'] = cm_ms10$overall['Accuracy']
forplotneg['ms10_TFmRNA', 'Accuracy'] = cm_ms10_TFmRNA$overall['Accuracy']
forplotneg['ms10_combine', 'Accuracy'] = cm_ms10_combine$overall['Accuracy']
forplotneg['ms15', 'Accuracy'] = cm_ms15$overall['Accuracy']
forplotneg['ms15_TFmRNA', 'Accuracy'] = cm_ms15_TFmRNA$overall['Accuracy']
forplotneg['ms15_combine', 'Accuracy'] = cm_ms15_combine$overall['Accuracy']

forplotneg['ms10', 'Balanced Accuracy'] = cm_ms10$byClass['Balanced Accuracy']
forplotneg['ms10_TFmRNA', 'Balanced Accuracy'] = cm_ms10_TFmRNA$byClass['Balanced Accuracy']
forplotneg['ms10_combine', 'Balanced Accuracy'] = cm_ms10_combine$byClass['Balanced Accuracy']
forplotneg['ms15', 'Balanced Accuracy'] = cm_ms15$byClass['Balanced Accuracy']
forplotneg['ms15_TFmRNA', 'Balanced Accuracy'] = cm_ms15_TFmRNA$byClass['Balanced Accuracy']
forplotneg['ms15_combine', 'Balanced Accuracy'] = cm_ms15_combine$byClass['Balanced Accuracy']

forplotneg['ms10', 'Precision'] = cm_ms10$byClass['Precision']
forplotneg['ms10_TFmRNA', 'Precision'] = cm_ms10_TFmRNA$byClass['Precision']
forplotneg['ms10_combine', 'Precision'] = cm_ms10_combine$byClass['Precision']
forplotneg['ms15', 'Precision'] = cm_ms15$byClass['Precision']
forplotneg['ms15_TFmRNA', 'Precision'] = cm_ms15_TFmRNA$byClass['Precision']
forplotneg['ms15_combine', 'Precision'] = cm_ms15_combine$byClass['Precision']


forplotneg['ms20', 'Accuracy'] = cm_ms20$overall['Accuracy']
forplotneg['ms20_TFmRNA', 'Accuracy'] = cm_ms20_TFmRNA$overall['Accuracy']
forplotneg['ms20_combine', 'Accuracy'] = cm_ms20_combine$overall['Accuracy']
forplotneg['ms5', 'Accuracy'] = cm_ms5$overall['Accuracy']
forplotneg['ms5_TFmRNA', 'Accuracy'] = cm_ms5_TFmRNA$overall['Accuracy']
forplotneg['ms5_combine', 'Accuracy'] = cm_ms5_combine$overall['Accuracy']

forplotneg['ms20', 'Balanced Accuracy'] = cm_ms20$byClass['Balanced Accuracy']
forplotneg['ms20_TFmRNA', 'Balanced Accuracy'] = cm_ms20_TFmRNA$byClass['Balanced Accuracy']
forplotneg['ms20_combine', 'Balanced Accuracy'] = cm_ms20_combine$byClass['Balanced Accuracy']
forplotneg['ms5', 'Balanced Accuracy'] = cm_ms5$byClass['Balanced Accuracy']
forplotneg['ms5_TFmRNA', 'Balanced Accuracy'] = cm_ms5_TFmRNA$byClass['Balanced Accuracy']
forplotneg['ms5_combine', 'Balanced Accuracy'] = cm_ms5_combine$byClass['Balanced Accuracy']

forplotneg['ms20', 'Precision'] = cm_ms20$byClass['Precision']
forplotneg['ms20_TFmRNA', 'Precision'] = cm_ms20_TFmRNA$byClass['Precision']
forplotneg['ms20_combine', 'Precision'] = cm_ms20_combine$byClass['Precision']
forplotneg['ms5', 'Precision'] = cm_ms5$byClass['Precision']
forplotneg['ms5_TFmRNA', 'Precision'] = cm_ms5_TFmRNA$byClass['Precision']
forplotneg['ms5_combine', 'Precision'] = cm_ms5_combine$byClass['Precision']


forplotneg$ms = c('Model size: 5', 'Model size: 5', 'Model size: 5', 
               'Model size: 10', 'Model size: 10', 'Model size: 10', 
               'Model size: 15', 'Model size: 15', 'Model size: 15',
               'Model size: 20', 'Model size: 20', 'Model size: 20')
forplotneg$ms = factor(forplotneg$ms, levels = c('Model size: 5', 'Model size: 10', 'Model size: 15', 'Model size: 20'))
forplotneg$TFAOpt = c('', 'TFmRNA', 'comb', '', 'TFmRNA', 'comb', '', 'TFmRNA', 'comb', '', 'TFmRNA', 'comb')
forplotneg$TFAOpt = factor(forplotneg$TFAOpt, levels = c('', 'TFmRNA', 'comb'))

forplotneg_melt = reshape2::melt(forplotneg)
library(ggplot2)
# library(RColorBrewer)
# brewer.pal(3, 'Set3')
eomes_neg_plot = ggplot(data=forplotneg_melt, aes(x=variable, y=value, fill=TFAOpt)) +
  geom_bar(stat="identity", position=position_dodge(), color="black", width = 0.7)+
  facet_wrap(~ms) + 
  scale_fill_manual(limits = c('', 'TFmRNA', 'comb'),
                    labels = c('TFA', 'TFmRNA', 'combine'),
                    name = 'TFA options',
                    values = c("#8DD3C7", "#FFFFB3", "#BEBADA"))+
  theme_bw(base_size = 14) +
  xlab('') +
  ylab('Value')+ 
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle('Eomes negative targets')

library(cowplot)
eomes.combined = plot_grid(eomes_pos_plot, eomes_neg_plot, align="hv", ncol = 2)
ggsave(
  paste0('compare_ms_tfa_Eomes_remake.pdf'),
  plot = eomes.combined,
  device = "pdf",
  width = 11,
  height = 6,
  dpi = 300
)


# Foxo1

# Foxo targets in Characterization of the direct targets of FOXO transcription factors throughout evolution

pool_genes = read.table('deg_padj01.txt')$V1

library(readxl)

targets = read_excel("foxo1/ACEL-15-673-s002.xlsx", 
                     sheet = "Sheet1")
targets_cd8 = unique(targets$`CD8 T cells`)
targets_allcelltype = unique(targets$Core)

targets_cd8 = targets_cd8[!is.na(targets_cd8)]
targets_allcelltype = targets_allcelltype[!is.na(targets_allcelltype)]

targets_cd8_pool_genes = targets_cd8[targets_cd8 %in% pool_genes]
targets_allcelltype_pool_genes = targets_allcelltype[targets_allcelltype %in% pool_genes]


ms10_combine = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms10_bias50_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms10         = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms10_bias50_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms10_TFmRNA  = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms10_bias50_TFmRNA_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms15_combine = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms15_bias50_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms15         = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms15_bias50_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms15_TFmRNA  = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms15_bias50_TFmRNA_cut01_sharedbyMorethan2nets_sp.tsv', header = T)

ms5_combine  = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms5_bias50_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms5          = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms5_bias50_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms5_TFmRNA   = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms5_bias50_TFmRNA_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms20_combine = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms20_bias50_maxComb_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms20         = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms20_bias50_cut01_sharedbyMorethan2nets_sp.tsv', header = T)
ms20_TFmRNA  = read.table('runs/sharedbyMorethan2netsFrom5seeds/chip_ko_network_atac_ms20_bias50_TFmRNA_cut01_sharedbyMorethan2nets_sp.tsv', header = T)

# both

ms15_combine_Foxo1 = ms15_combine[ms15_combine$TF == 'Foxo1', 'Target']

#table(ms15_combine_Foxo1 %in% pool_genes)
#table(targets_cd8_pool_genes %in% pool_genes)

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms15_combine_Foxo1_list = ifelse(pool_genes %in% ms15_combine_Foxo1, 1, 0)

xtab <- table(ms15_combine_Foxo1_list, targets_cd8_pool_genes_list)
#cm <- caret::confusionMatrix(xtab)
cm_ms15_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15_combine$overall['Accuracy']
cm_ms15_combine$byClass['Precision']
cm_ms15_combine$byClass['Recall']
cm_ms15_combine$byClass['F1']
cm_ms15_combine$byClass['Balanced Accuracy']


ms15_TFmRNA_Foxo1 = ms15_TFmRNA[ms15_TFmRNA$TF == 'Foxo1', 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms15_TFmRNA_Foxo1_list = ifelse(pool_genes %in% ms15_TFmRNA_Foxo1, 1, 0)

xtab <- table(ms15_TFmRNA_Foxo1_list, targets_cd8_pool_genes_list)
cm_ms15_TFmRNA <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15_TFmRNA$overall['Accuracy']
cm_ms15_TFmRNA$byClass['Precision']
cm_ms15_TFmRNA$byClass['Recall']
cm_ms15_TFmRNA$byClass['F1']
cm_ms15_TFmRNA$byClass['Balanced Accuracy']


ms15_Foxo1 = ms15[ms15$TF == 'Foxo1', 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms15_Foxo1_list = ifelse(pool_genes %in% ms15_Foxo1, 1, 0)

xtab <- table(ms15_Foxo1_list, targets_cd8_pool_genes_list)
cm_ms15 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15$overall['Accuracy']
cm_ms15$byClass['Precision']
cm_ms15$byClass['Recall']
cm_ms15$byClass['F1']
cm_ms15$byClass['Balanced Accuracy']



ms10_combine_Foxo1 = ms10_combine[ms10_combine$TF == 'Foxo1', 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms10_combine_Foxo1_list = ifelse(pool_genes %in% ms10_combine_Foxo1, 1, 0)

xtab <- table(ms10_combine_Foxo1_list, targets_cd8_pool_genes_list)
cm_ms10_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10_combine$overall['Accuracy']
cm_ms10_combine$byClass['Precision']
cm_ms10_combine$byClass['Recall']
cm_ms10_combine$byClass['F1']
cm_ms10_combine$byClass['Balanced Accuracy']


ms10_TFmRNA_Foxo1 = ms10_TFmRNA[ms10_TFmRNA$TF == 'Foxo1', 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms10_TFmRNA_Foxo1_list = ifelse(pool_genes %in% ms10_TFmRNA_Foxo1, 1, 0)

xtab <- table(ms10_TFmRNA_Foxo1_list, targets_cd8_pool_genes_list)
cm_ms10_TFmRNA <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10_TFmRNA$overall['Accuracy']
cm_ms10_TFmRNA$byClass['Precision']
cm_ms10_TFmRNA$byClass['Recall']
cm_ms10_TFmRNA$byClass['F1']
cm_ms10_TFmRNA$byClass['Balanced Accuracy']


ms10_Foxo1 = ms10[ms10$TF == 'Foxo1', 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms10_Foxo1_list = ifelse(pool_genes %in% ms10_Foxo1, 1, 0)

xtab <- table(ms10_Foxo1_list, targets_cd8_pool_genes_list)
cm_ms10 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10$overall['Accuracy']
cm_ms10$byClass['Precision']
cm_ms10$byClass['Recall']
cm_ms10$byClass['F1']
cm_ms10$byClass['Balanced Accuracy']


ms5_combine_Foxo1 = ms5_combine[ms5_combine$TF == 'Foxo1', 'Target']

#table(ms5_combine_Foxo1 %in% pool_genes)
#table(targets_cd8_pool_genes %in% pool_genes)

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms5_combine_Foxo1_list = ifelse(pool_genes %in% ms5_combine_Foxo1, 1, 0)

xtab <- table(ms5_combine_Foxo1_list, targets_cd8_pool_genes_list)
#cm <- caret::confusionMatrix(xtab)
cm_ms5_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5_combine$overall['Accuracy']
cm_ms5_combine$byClass['Precision']
cm_ms5_combine$byClass['Recall']
cm_ms5_combine$byClass['F1']
cm_ms5_combine$byClass['Balanced Accuracy']


ms5_TFmRNA_Foxo1 = ms5_TFmRNA[ms5_TFmRNA$TF == 'Foxo1', 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms5_TFmRNA_Foxo1_list = ifelse(pool_genes %in% ms5_TFmRNA_Foxo1, 1, 0)

xtab <- table(ms5_TFmRNA_Foxo1_list, targets_cd8_pool_genes_list)
cm_ms5_TFmRNA <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5_TFmRNA$overall['Accuracy']
cm_ms5_TFmRNA$byClass['Precision']
cm_ms5_TFmRNA$byClass['Recall']
cm_ms5_TFmRNA$byClass['F1']
cm_ms5_TFmRNA$byClass['Balanced Accuracy']


ms5_Foxo1 = ms5[ms5$TF == 'Foxo1', 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms5_Foxo1_list = ifelse(pool_genes %in% ms5_Foxo1, 1, 0)

xtab <- table(ms5_Foxo1_list, targets_cd8_pool_genes_list)
cm_ms5 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5$overall['Accuracy']
cm_ms5$byClass['Precision']
cm_ms5$byClass['Recall']
cm_ms5$byClass['F1']
cm_ms5$byClass['Balanced Accuracy']



ms20_combine_Foxo1 = ms20_combine[ms20_combine$TF == 'Foxo1', 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms20_combine_Foxo1_list = ifelse(pool_genes %in% ms20_combine_Foxo1, 1, 0)

xtab <- table(ms20_combine_Foxo1_list, targets_cd8_pool_genes_list)
cm_ms20_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20_combine$overall['Accuracy']
cm_ms20_combine$byClass['Precision']
cm_ms20_combine$byClass['Recall']
cm_ms20_combine$byClass['F1']
cm_ms20_combine$byClass['Balanced Accuracy']


ms20_TFmRNA_Foxo1 = ms20_TFmRNA[ms20_TFmRNA$TF == 'Foxo1', 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms20_TFmRNA_Foxo1_list = ifelse(pool_genes %in% ms20_TFmRNA_Foxo1, 1, 0)

xtab <- table(ms20_TFmRNA_Foxo1_list, targets_cd8_pool_genes_list)
cm_ms20_TFmRNA <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20_TFmRNA$overall['Accuracy']
cm_ms20_TFmRNA$byClass['Precision']
cm_ms20_TFmRNA$byClass['Recall']
cm_ms20_TFmRNA$byClass['F1']
cm_ms20_TFmRNA$byClass['Balanced Accuracy']


ms20_Foxo1 = ms20[ms20$TF == 'Foxo1', 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms20_Foxo1_list = ifelse(pool_genes %in% ms20_Foxo1, 1, 0)

xtab <- table(ms20_Foxo1_list, targets_cd8_pool_genes_list)
cm_ms20 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20$overall['Accuracy']
cm_ms20$byClass['Precision']
cm_ms20$byClass['Recall']
cm_ms20$byClass['F1']
cm_ms20$byClass['Balanced Accuracy']


# library(Metrics)
# Metrics::recall(ms15_combine_Foxo1_list, targets_cd8_pool_genes_list)
# Metrics::precision(ms15_combine_Foxo1_list, targets_cd8_pool_genes_list)
# Metrics::f1(ms15_combine_Foxo1_list, targets_cd8_pool_genes_list) # wrong somehow


# make plots
forplot = data.frame(matrix(ncol = 3, nrow = 12))
colnames(forplot) = c('Accuracy', 'Balanced Accuracy', 'Precision')
rownames(forplot) = c('ms5_TFmRNA', 'ms5_combine', 'ms5', 'ms10_TFmRNA', 'ms10_combine', 'ms10', 
                      'ms15_TFmRNA', 'ms15_combine', 'ms15', 'ms20_TFmRNA', 'ms20_combine', 'ms20')
forplot['ms10_TFmRNA', 'Accuracy'] = cm_ms10_TFmRNA$overall['Accuracy']
forplot['ms10_combine', 'Accuracy'] = cm_ms10_combine$overall['Accuracy']
forplot['ms10', 'Accuracy'] = cm_ms10$overall['Accuracy']
forplot['ms15_TFmRNA', 'Accuracy'] = cm_ms15_TFmRNA$overall['Accuracy']
forplot['ms15_combine', 'Accuracy'] = cm_ms15_combine$overall['Accuracy']
forplot['ms15', 'Accuracy'] = cm_ms15$overall['Accuracy']

forplot['ms10_TFmRNA', 'Balanced Accuracy'] = cm_ms10_TFmRNA$byClass['Balanced Accuracy']
forplot['ms10_combine', 'Balanced Accuracy'] = cm_ms10_combine$byClass['Balanced Accuracy']
forplot['ms10', 'Balanced Accuracy'] = cm_ms10$byClass['Balanced Accuracy']
forplot['ms15_TFmRNA', 'Balanced Accuracy'] = cm_ms15_TFmRNA$byClass['Balanced Accuracy']
forplot['ms15_combine', 'Balanced Accuracy'] = cm_ms15_combine$byClass['Balanced Accuracy']
forplot['ms15', 'Balanced Accuracy'] = cm_ms15$byClass['Balanced Accuracy']

forplot['ms10_TFmRNA', 'Precision'] = cm_ms10_TFmRNA$byClass['Precision']
forplot['ms10_combine', 'Precision'] = cm_ms10_combine$byClass['Precision']
forplot['ms10', 'Precision'] = cm_ms10$byClass['Precision']
forplot['ms15_TFmRNA', 'Precision'] = cm_ms15_TFmRNA$byClass['Precision']
forplot['ms15_combine', 'Precision'] = cm_ms15_combine$byClass['Precision']
forplot['ms15', 'Precision'] = cm_ms15$byClass['Precision']

forplot['ms20_TFmRNA', 'Accuracy'] = cm_ms20_TFmRNA$overall['Accuracy']
forplot['ms20_combine', 'Accuracy'] = cm_ms20_combine$overall['Accuracy']
forplot['ms20', 'Accuracy'] = cm_ms20$overall['Accuracy']
forplot['ms5_TFmRNA', 'Accuracy'] = cm_ms5_TFmRNA$overall['Accuracy']
forplot['ms5_combine', 'Accuracy'] = cm_ms5_combine$overall['Accuracy']
forplot['ms5', 'Accuracy'] = cm_ms5$overall['Accuracy']

forplot['ms20_TFmRNA', 'Balanced Accuracy'] = cm_ms20_TFmRNA$byClass['Balanced Accuracy']
forplot['ms20_combine', 'Balanced Accuracy'] = cm_ms20_combine$byClass['Balanced Accuracy']
forplot['ms20', 'Balanced Accuracy'] = cm_ms20$byClass['Balanced Accuracy']
forplot['ms5_TFmRNA', 'Balanced Accuracy'] = cm_ms5_TFmRNA$byClass['Balanced Accuracy']
forplot['ms5_combine', 'Balanced Accuracy'] = cm_ms5_combine$byClass['Balanced Accuracy']
forplot['ms5', 'Balanced Accuracy'] = cm_ms5$byClass['Balanced Accuracy']

forplot['ms20_TFmRNA', 'Precision'] = cm_ms20_TFmRNA$byClass['Precision']
forplot['ms20_combine', 'Precision'] = cm_ms20_combine$byClass['Precision']
forplot['ms20', 'Precision'] = cm_ms20$byClass['Precision']
forplot['ms5_TFmRNA', 'Precision'] = cm_ms5_TFmRNA$byClass['Precision']
forplot['ms5_combine', 'Precision'] = cm_ms5_combine$byClass['Precision']
forplot['ms5', 'Precision'] = cm_ms5$byClass['Precision']

forplot$ms = c('Model size: 5', 'Model size: 5', 'Model size: 5', 
               'Model size: 10', 'Model size: 10', 'Model size: 10', 
               'Model size: 15', 'Model size: 15', 'Model size: 15',
               'Model size: 20', 'Model size: 20', 'Model size: 20')
forplot$ms = factor(forplot$ms, levels = c('Model size: 5', 'Model size: 10', 'Model size: 15', 'Model size: 20'))
forplot$TFAOpt = c('TFmRNA', 'combine', '', 'TFmRNA', 'combine', '', 'TFmRNA', 'combine', '', 'TFmRNA', 'combine', '')
forplot$TFAOpt = factor(forplot$TFAOpt, levels = c('TFmRNA', '', 'combine'))

forplot_melt = reshape2::melt(forplot)
library(ggplot2)
# library(RColorBrewer)
# brewer.pal(3, 'Set3')
foxo1_plot = ggplot(data=forplot_melt, aes(x=variable, y=value, fill=TFAOpt)) +
  geom_bar(stat="identity", position=position_dodge(), color="black", width = 0.7)+
  facet_wrap(~ms) + 
  scale_fill_manual(limits = c('TFmRNA', '', 'combine'),
                    labels = c('TFmRNA', 'TFA', 'combine'),
                    name = 'TFA options',
                    values = c("#8DD3C7", "#FFFFB3", "#BEBADA"))+
  theme_bw(base_size = 14) +
  xlab('') +
  ylab('Value')+ 
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle('Foxo1 targets')

ggsave(
  paste0('compare_ms_tfa_Foxo1_remake.pdf'),
  plot = foxo1_plot,
  device = "pdf",
  width = 6,
  height = 6,
  dpi = 300
)



# positive targets
ms15_combine_Foxo1_pos = ms15_combine[ms15_combine$TF == 'Foxo1' & ms15_combine$SignedQuantile == 1, 'Target']
ms15_combine_Foxo1_neg = ms15_combine[ms15_combine$TF == 'Foxo1' & ms15_combine$SignedQuantile == -1, 'Target']

#table(ms15_combine_Foxo1_pos %in% pool_genes)
#table(ms15_combine_Foxo1_neg %in% pool_genes)
#table(targets_cd8_pool_genes %in% pool_genes)
#table(targets_cd8_pool_genes %in% pool_genes)

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms15_combine_Foxo1_pos_list = ifelse(pool_genes %in% ms15_combine_Foxo1_pos, 1, 0)

xtab <- table(ms15_combine_Foxo1_pos_list, targets_cd8_pool_genes_list)
#cm <- caret::confusionMatrix(xtab)
cm_ms15_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15_combine$overall['Accuracy']
cm_ms15_combine$byClass['Precision']
cm_ms15_combine$byClass['Recall']
cm_ms15_combine$byClass['F1']
cm_ms15_combine$byClass['Balanced Accuracy']


ms15_TFmRNA_Foxo1_pos = ms15_TFmRNA[ms15_TFmRNA$TF == 'Foxo1' & ms15_TFmRNA$SignedQuantile == 1, 'Target']
ms15_TFmRNA_Foxo1_neg = ms15_TFmRNA[ms15_TFmRNA$TF == 'Foxo1' & ms15_TFmRNA$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms15_TFmRNA_Foxo1_pos_list = ifelse(pool_genes %in% ms15_TFmRNA_Foxo1_pos, 1, 0)

xtab <- table(ms15_TFmRNA_Foxo1_pos_list, targets_cd8_pool_genes_list)
cm_ms15_TFmRNA <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15_TFmRNA$overall['Accuracy']
cm_ms15_TFmRNA$byClass['Precision']
cm_ms15_TFmRNA$byClass['Recall']
cm_ms15_TFmRNA$byClass['F1']
cm_ms15_TFmRNA$byClass['Balanced Accuracy']


ms15_Foxo1_pos = ms15[ms15$TF == 'Foxo1' & ms15$SignedQuantile == 1, 'Target']
ms15_Foxo1_neg = ms15[ms15$TF == 'Foxo1' & ms15$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms15_Foxo1_pos_list = ifelse(pool_genes %in% ms15_Foxo1_pos, 1, 0)

xtab <- table(ms15_Foxo1_pos_list, targets_cd8_pool_genes_list)
cm_ms15 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15$overall['Accuracy']
cm_ms15$byClass['Precision']
cm_ms15$byClass['Recall']
cm_ms15$byClass['F1']
cm_ms15$byClass['Balanced Accuracy']



ms10_combine_Foxo1_pos = ms10_combine[ms10_combine$TF == 'Foxo1' & ms10_combine$SignedQuantile == 1, 'Target']
ms10_combine_Foxo1_neg = ms10_combine[ms10_combine$TF == 'Foxo1' & ms10_combine$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms10_combine_Foxo1_pos_list = ifelse(pool_genes %in% ms10_combine_Foxo1_pos, 1, 0)

xtab <- table(ms10_combine_Foxo1_pos_list, targets_cd8_pool_genes_list)
cm_ms10_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10_combine$overall['Accuracy']
cm_ms10_combine$byClass['Precision']
cm_ms10_combine$byClass['Recall']
cm_ms10_combine$byClass['F1']
cm_ms10_combine$byClass['Balanced Accuracy']


ms10_TFmRNA_Foxo1_pos = ms10_TFmRNA[ms10_TFmRNA$TF == 'Foxo1' & ms10_TFmRNA$SignedQuantile == 1, 'Target']
ms10_TFmRNA_Foxo1_neg = ms10_TFmRNA[ms10_TFmRNA$TF == 'Foxo1' & ms10_TFmRNA$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms10_TFmRNA_Foxo1_pos_list = ifelse(pool_genes %in% ms10_TFmRNA_Foxo1_pos, 1, 0)

xtab <- table(ms10_TFmRNA_Foxo1_pos_list, targets_cd8_pool_genes_list)
cm_ms10_TFmRNA <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10_TFmRNA$overall['Accuracy']
cm_ms10_TFmRNA$byClass['Precision']
cm_ms10_TFmRNA$byClass['Recall']
cm_ms10_TFmRNA$byClass['F1']
cm_ms10_TFmRNA$byClass['Balanced Accuracy']


ms10_Foxo1_pos = ms10[ms10$TF == 'Foxo1' & ms10$SignedQuantile == 1, 'Target']
ms10_Foxo1_neg = ms10[ms10$TF == 'Foxo1' & ms10$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms10_Foxo1_pos_list = ifelse(pool_genes %in% ms10_Foxo1_pos, 1, 0)

xtab <- table(ms10_Foxo1_pos_list, targets_cd8_pool_genes_list)
cm_ms10 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10$overall['Accuracy']
cm_ms10$byClass['Precision']
cm_ms10$byClass['Recall']
cm_ms10$byClass['F1']
cm_ms10$byClass['Balanced Accuracy']


ms5_combine_Foxo1_pos = ms5_combine[ms5_combine$TF == 'Foxo1' & ms5_combine$SignedQuantile == 1, 'Target']
ms5_combine_Foxo1_neg = ms5_combine[ms5_combine$TF == 'Foxo1' & ms5_combine$SignedQuantile == -1, 'Target']

#table(ms5_combine_Foxo1_pos %in% pool_genes)
#table(ms5_combine_Foxo1_neg %in% pool_genes)
#table(targets_cd8_pool_genes %in% pool_genes)
#table(targets_cd8_pool_genes %in% pool_genes)

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms5_combine_Foxo1_pos_list = ifelse(pool_genes %in% ms5_combine_Foxo1_pos, 1, 0)

xtab <- table(ms5_combine_Foxo1_pos_list, targets_cd8_pool_genes_list)
#cm <- caret::confusionMatrix(xtab)
cm_ms5_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5_combine$overall['Accuracy']
cm_ms5_combine$byClass['Precision']
cm_ms5_combine$byClass['Recall']
cm_ms5_combine$byClass['F1']
cm_ms5_combine$byClass['Balanced Accuracy']


ms5_TFmRNA_Foxo1_pos = ms5_TFmRNA[ms5_TFmRNA$TF == 'Foxo1' & ms5_TFmRNA$SignedQuantile == 1, 'Target']
ms5_TFmRNA_Foxo1_neg = ms5_TFmRNA[ms5_TFmRNA$TF == 'Foxo1' & ms5_TFmRNA$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms5_TFmRNA_Foxo1_pos_list = ifelse(pool_genes %in% ms5_TFmRNA_Foxo1_pos, 1, 0)

xtab <- table(ms5_TFmRNA_Foxo1_pos_list, targets_cd8_pool_genes_list)
cm_ms5_TFmRNA <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5_TFmRNA$overall['Accuracy']
cm_ms5_TFmRNA$byClass['Precision']
cm_ms5_TFmRNA$byClass['Recall']
cm_ms5_TFmRNA$byClass['F1']
cm_ms5_TFmRNA$byClass['Balanced Accuracy']


ms5_Foxo1_pos = ms5[ms5$TF == 'Foxo1' & ms5$SignedQuantile == 1, 'Target']
ms5_Foxo1_neg = ms5[ms5$TF == 'Foxo1' & ms5$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms5_Foxo1_pos_list = ifelse(pool_genes %in% ms5_Foxo1_pos, 1, 0)

xtab <- table(ms5_Foxo1_pos_list, targets_cd8_pool_genes_list)
cm_ms5 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5$overall['Accuracy']
cm_ms5$byClass['Precision']
cm_ms5$byClass['Recall']
cm_ms5$byClass['F1']
cm_ms5$byClass['Balanced Accuracy']



ms20_combine_Foxo1_pos = ms20_combine[ms20_combine$TF == 'Foxo1' & ms20_combine$SignedQuantile == 1, 'Target']
ms20_combine_Foxo1_neg = ms20_combine[ms20_combine$TF == 'Foxo1' & ms20_combine$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms20_combine_Foxo1_pos_list = ifelse(pool_genes %in% ms20_combine_Foxo1_pos, 1, 0)

xtab <- table(ms20_combine_Foxo1_pos_list, targets_cd8_pool_genes_list)
cm_ms20_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20_combine$overall['Accuracy']
cm_ms20_combine$byClass['Precision']
cm_ms20_combine$byClass['Recall']
cm_ms20_combine$byClass['F1']
cm_ms20_combine$byClass['Balanced Accuracy']


ms20_TFmRNA_Foxo1_pos = ms20_TFmRNA[ms20_TFmRNA$TF == 'Foxo1' & ms20_TFmRNA$SignedQuantile == 1, 'Target']
ms20_TFmRNA_Foxo1_neg = ms20_TFmRNA[ms20_TFmRNA$TF == 'Foxo1' & ms20_TFmRNA$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms20_TFmRNA_Foxo1_pos_list = ifelse(pool_genes %in% ms20_TFmRNA_Foxo1_pos, 1, 0)

xtab <- table(ms20_TFmRNA_Foxo1_pos_list, targets_cd8_pool_genes_list)
cm_ms20_TFmRNA <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20_TFmRNA$overall['Accuracy']
cm_ms20_TFmRNA$byClass['Precision']
cm_ms20_TFmRNA$byClass['Recall']
cm_ms20_TFmRNA$byClass['F1']
cm_ms20_TFmRNA$byClass['Balanced Accuracy']


ms20_Foxo1_pos = ms20[ms20$TF == 'Foxo1' & ms20$SignedQuantile == 1, 'Target']
ms20_Foxo1_neg = ms20[ms20$TF == 'Foxo1' & ms20$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms20_Foxo1_pos_list = ifelse(pool_genes %in% ms20_Foxo1_pos, 1, 0)

xtab <- table(ms20_Foxo1_pos_list, targets_cd8_pool_genes_list)
cm_ms20 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20$overall['Accuracy']
cm_ms20$byClass['Precision']
cm_ms20$byClass['Recall']
cm_ms20$byClass['F1']
cm_ms20$byClass['Balanced Accuracy']



# library(Metrics)
# Metrics::recall(ms15_combine_Foxo1_pos_list, targets_cd8_pool_genes_list)
# Metrics::precision(ms15_combine_Foxo1_pos_list, targets_cd8_pool_genes_list)
# Metrics::f1(ms15_combine_Foxo1_pos_list, targets_cd8_pool_genes_list) # wrong somehow


# make plots
forplot = data.frame(matrix(ncol = 3, nrow = 12))
colnames(forplot) = c('Accuracy', 'Balanced Accuracy', 'Precision')
rownames(forplot) = c('ms5_TFmRNA', 'ms5_combine', 'ms5', 'ms10_TFmRNA', 'ms10_combine', 'ms10', 
                      'ms15_TFmRNA', 'ms15_combine', 'ms15', 'ms20_TFmRNA', 'ms20_combine', 'ms20')
forplot['ms10_TFmRNA', 'Accuracy'] = cm_ms10_TFmRNA$overall['Accuracy']
forplot['ms10_combine', 'Accuracy'] = cm_ms10_combine$overall['Accuracy']
forplot['ms10', 'Accuracy'] = cm_ms10$overall['Accuracy']
forplot['ms15_TFmRNA', 'Accuracy'] = cm_ms15_TFmRNA$overall['Accuracy']
forplot['ms15_combine', 'Accuracy'] = cm_ms15_combine$overall['Accuracy']
forplot['ms15', 'Accuracy'] = cm_ms15$overall['Accuracy']

forplot['ms10_TFmRNA', 'Balanced Accuracy'] = cm_ms10_TFmRNA$byClass['Balanced Accuracy']
forplot['ms10_combine', 'Balanced Accuracy'] = cm_ms10_combine$byClass['Balanced Accuracy']
forplot['ms10', 'Balanced Accuracy'] = cm_ms10$byClass['Balanced Accuracy']
forplot['ms15_TFmRNA', 'Balanced Accuracy'] = cm_ms15_TFmRNA$byClass['Balanced Accuracy']
forplot['ms15_combine', 'Balanced Accuracy'] = cm_ms15_combine$byClass['Balanced Accuracy']
forplot['ms15', 'Balanced Accuracy'] = cm_ms15$byClass['Balanced Accuracy']

forplot['ms10_TFmRNA', 'Precision'] = cm_ms10_TFmRNA$byClass['Precision']
forplot['ms10_combine', 'Precision'] = cm_ms10_combine$byClass['Precision']
forplot['ms10', 'Precision'] = cm_ms10$byClass['Precision']
forplot['ms15_TFmRNA', 'Precision'] = cm_ms15_TFmRNA$byClass['Precision']
forplot['ms15_combine', 'Precision'] = cm_ms15_combine$byClass['Precision']
forplot['ms15', 'Precision'] = cm_ms15$byClass['Precision']


forplot['ms20_TFmRNA', 'Accuracy'] = cm_ms20_TFmRNA$overall['Accuracy']
forplot['ms20_combine', 'Accuracy'] = cm_ms20_combine$overall['Accuracy']
forplot['ms20', 'Accuracy'] = cm_ms20$overall['Accuracy']
forplot['ms5_TFmRNA', 'Accuracy'] = cm_ms5_TFmRNA$overall['Accuracy']
forplot['ms5_combine', 'Accuracy'] = cm_ms5_combine$overall['Accuracy']
forplot['ms5', 'Accuracy'] = cm_ms5$overall['Accuracy']

forplot['ms20_TFmRNA', 'Balanced Accuracy'] = cm_ms20_TFmRNA$byClass['Balanced Accuracy']
forplot['ms20_combine', 'Balanced Accuracy'] = cm_ms20_combine$byClass['Balanced Accuracy']
forplot['ms20', 'Balanced Accuracy'] = cm_ms20$byClass['Balanced Accuracy']
forplot['ms5_TFmRNA', 'Balanced Accuracy'] = cm_ms5_TFmRNA$byClass['Balanced Accuracy']
forplot['ms5_combine', 'Balanced Accuracy'] = cm_ms5_combine$byClass['Balanced Accuracy']
forplot['ms5', 'Balanced Accuracy'] = cm_ms5$byClass['Balanced Accuracy']

forplot['ms20_TFmRNA', 'Precision'] = cm_ms20_TFmRNA$byClass['Precision']
forplot['ms20_combine', 'Precision'] = cm_ms20_combine$byClass['Precision']
forplot['ms20', 'Precision'] = cm_ms20$byClass['Precision']
forplot['ms5_TFmRNA', 'Precision'] = cm_ms5_TFmRNA$byClass['Precision']
forplot['ms5_combine', 'Precision'] = cm_ms5_combine$byClass['Precision']
forplot['ms5', 'Precision'] = cm_ms5$byClass['Precision']

forplot$ms = c('Model size: 5', 'Model size: 5', 'Model size: 5', 
               'Model size: 10', 'Model size: 10', 'Model size: 10', 
               'Model size: 15', 'Model size: 15', 'Model size: 15',
               'Model size: 20', 'Model size: 20', 'Model size: 20')
forplot$ms = factor(forplot$ms, levels = c('Model size: 5', 'Model size: 10', 'Model size: 15', 'Model size: 20'))
forplot$TFAOpt = c('TFmRNA', 'combine', '', 'TFmRNA', 'combine', '', 'TFmRNA', 'combine', '', 'TFmRNA', 'combine', '')
forplot$TFAOpt = factor(forplot$TFAOpt, levels = c('TFmRNA', '', 'combine'))

forplot_melt = reshape2::melt(forplot)
library(ggplot2)
# library(RColorBrewer)
# brewer.pal(3, 'Set3')
foxo1_pos_plot = ggplot(data=forplot_melt, aes(x=variable, y=value, fill=TFAOpt)) +
  geom_bar(stat="identity", position=position_dodge(), color="black", width = 0.7)+
  facet_wrap(~ms) + 
  scale_fill_manual(limits = c('TFmRNA', '', 'combine'),
                    labels = c('TFmRNA', 'TFA', 'combine'),
                    name = 'TFA options',
                    values = c("#8DD3C7", "#FFFFB3", "#BEBADA"))+
  theme_bw(base_size = 14) +
  xlab('') +
  ylab('Value')+ 
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle('Foxo1 positive targets')




# negative targets

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms15_combine_Foxo1_neg_list = ifelse(pool_genes %in% ms15_combine_Foxo1_neg, 1, 0)

xtab <- table(ms15_combine_Foxo1_neg_list, targets_cd8_pool_genes_list)
#cm <- caret::confusionMatrix(xtab)
cm_ms15_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15_combine$overall['Accuracy']
cm_ms15_combine$byClass['Precision']
cm_ms15_combine$byClass['Recall']
cm_ms15_combine$byClass['F1']
cm_ms15_combine$byClass['Balanced Accuracy']


ms15_TFmRNA_Foxo1_neg = ms15_TFmRNA[ms15_TFmRNA$TF == 'Foxo1' & ms15_TFmRNA$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms15_TFmRNA_Foxo1_neg_list = ifelse(pool_genes %in% ms15_TFmRNA_Foxo1_neg, 1, 0)

xtab <- table(ms15_TFmRNA_Foxo1_neg_list, targets_cd8_pool_genes_list)
cm_ms15_TFmRNA <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15_TFmRNA$overall['Accuracy']
cm_ms15_TFmRNA$byClass['Precision']
cm_ms15_TFmRNA$byClass['Recall']
cm_ms15_TFmRNA$byClass['F1']
cm_ms15_TFmRNA$byClass['Balanced Accuracy']


ms15_Foxo1_neg = ms15[ms15$TF == 'Foxo1' & ms15$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms15_Foxo1_neg_list = ifelse(pool_genes %in% ms15_Foxo1_neg, 1, 0)

xtab <- table(ms15_Foxo1_neg_list, targets_cd8_pool_genes_list)
cm_ms15 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms15$overall['Accuracy']
cm_ms15$byClass['Precision']
cm_ms15$byClass['Recall']
cm_ms15$byClass['F1']
cm_ms15$byClass['Balanced Accuracy']



ms10_combine_Foxo1_neg = ms10_combine[ms10_combine$TF == 'Foxo1' & ms10_combine$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms10_combine_Foxo1_neg_list = ifelse(pool_genes %in% ms10_combine_Foxo1_neg, 1, 0)

xtab <- table(ms10_combine_Foxo1_neg_list, targets_cd8_pool_genes_list)
cm_ms10_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10_combine$overall['Accuracy']
cm_ms10_combine$byClass['Precision']
cm_ms10_combine$byClass['Recall']
cm_ms10_combine$byClass['F1']
cm_ms10_combine$byClass['Balanced Accuracy']


ms10_TFmRNA_Foxo1_neg = ms10_TFmRNA[ms10_TFmRNA$TF == 'Foxo1' & ms10_TFmRNA$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms10_TFmRNA_Foxo1_neg_list = ifelse(pool_genes %in% ms10_TFmRNA_Foxo1_neg, 1, 0)

xtab <- table(ms10_TFmRNA_Foxo1_neg_list, targets_cd8_pool_genes_list)
cm_ms10_TFmRNA <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10_TFmRNA$overall['Accuracy']
cm_ms10_TFmRNA$byClass['Precision']
cm_ms10_TFmRNA$byClass['Recall']
cm_ms10_TFmRNA$byClass['F1']
cm_ms10_TFmRNA$byClass['Balanced Accuracy']


ms10_Foxo1_neg = ms10[ms10$TF == 'Foxo1' & ms10$SignedQuantile == 1, 'Target']
ms10_Foxo1_neg = ms10[ms10$TF == 'Foxo1' & ms10$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms10_Foxo1_neg_list = ifelse(pool_genes %in% ms10_Foxo1_neg, 1, 0)

xtab <- table(ms10_Foxo1_neg_list, targets_cd8_pool_genes_list)
cm_ms10 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms10$overall['Accuracy']
cm_ms10$byClass['Precision']
cm_ms10$byClass['Recall']
cm_ms10$byClass['F1']
cm_ms10$byClass['Balanced Accuracy']



targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms20_combine_Foxo1_neg_list = ifelse(pool_genes %in% ms20_combine_Foxo1_neg, 1, 0)

xtab <- table(ms20_combine_Foxo1_neg_list, targets_cd8_pool_genes_list)
#cm <- caret::confusionMatrix(xtab)
cm_ms20_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20_combine$overall['Accuracy']
cm_ms20_combine$byClass['Precision']
cm_ms20_combine$byClass['Recall']
cm_ms20_combine$byClass['F1']
cm_ms20_combine$byClass['Balanced Accuracy']


ms20_TFmRNA_Foxo1_neg = ms20_TFmRNA[ms20_TFmRNA$TF == 'Foxo1' & ms20_TFmRNA$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms20_TFmRNA_Foxo1_neg_list = ifelse(pool_genes %in% ms20_TFmRNA_Foxo1_neg, 1, 0)

xtab <- table(ms20_TFmRNA_Foxo1_neg_list, targets_cd8_pool_genes_list)
cm_ms20_TFmRNA <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20_TFmRNA$overall['Accuracy']
cm_ms20_TFmRNA$byClass['Precision']
cm_ms20_TFmRNA$byClass['Recall']
cm_ms20_TFmRNA$byClass['F1']
cm_ms20_TFmRNA$byClass['Balanced Accuracy']


ms20_Foxo1_neg = ms20[ms20$TF == 'Foxo1' & ms20$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms20_Foxo1_neg_list = ifelse(pool_genes %in% ms20_Foxo1_neg, 1, 0)

xtab <- table(ms20_Foxo1_neg_list, targets_cd8_pool_genes_list)
cm_ms20 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms20$overall['Accuracy']
cm_ms20$byClass['Precision']
cm_ms20$byClass['Recall']
cm_ms20$byClass['F1']
cm_ms20$byClass['Balanced Accuracy']



ms5_combine_Foxo1_neg = ms5_combine[ms5_combine$TF == 'Foxo1' & ms5_combine$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms5_combine_Foxo1_neg_list = ifelse(pool_genes %in% ms5_combine_Foxo1_neg, 1, 0)

xtab <- table(ms5_combine_Foxo1_neg_list, targets_cd8_pool_genes_list)
cm_ms5_combine <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5_combine$overall['Accuracy']
cm_ms5_combine$byClass['Precision']
cm_ms5_combine$byClass['Recall']
cm_ms5_combine$byClass['F1']
cm_ms5_combine$byClass['Balanced Accuracy']


ms5_TFmRNA_Foxo1_neg = ms5_TFmRNA[ms5_TFmRNA$TF == 'Foxo1' & ms5_TFmRNA$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms5_TFmRNA_Foxo1_neg_list = ifelse(pool_genes %in% ms5_TFmRNA_Foxo1_neg, 1, 0)

xtab <- table(ms5_TFmRNA_Foxo1_neg_list, targets_cd8_pool_genes_list)
cm_ms5_TFmRNA <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5_TFmRNA$overall['Accuracy']
cm_ms5_TFmRNA$byClass['Precision']
cm_ms5_TFmRNA$byClass['Recall']
cm_ms5_TFmRNA$byClass['F1']
cm_ms5_TFmRNA$byClass['Balanced Accuracy']


ms5_Foxo1_neg = ms5[ms5$TF == 'Foxo1' & ms5$SignedQuantile == 1, 'Target']
ms5_Foxo1_neg = ms5[ms5$TF == 'Foxo1' & ms5$SignedQuantile == -1, 'Target']

targets_cd8_pool_genes_list = ifelse(pool_genes %in% targets_cd8_pool_genes, 1, 0)
ms5_Foxo1_neg_list = ifelse(pool_genes %in% ms5_Foxo1_neg, 1, 0)

xtab <- table(ms5_Foxo1_neg_list, targets_cd8_pool_genes_list)
cm_ms5 <- caret::confusionMatrix(xtab, positive = '1', mode = 'everything')
cm_ms5$overall['Accuracy']
cm_ms5$byClass['Precision']
cm_ms5$byClass['Recall']
cm_ms5$byClass['F1']
cm_ms5$byClass['Balanced Accuracy']

# library(Metrics)
# Metrics::recall(ms15_combine_Foxo1_neg_list, targets_cd8_pool_genes_list)
# Metrics::precision(ms15_combine_Foxo1_neg_list, targets_cd8_pool_genes_list)
# Metrics::f1(ms15_combine_Foxo1_neg_list, targets_cd8_pool_genes_list) # wrong somehow


# make plots
forplotneg = data.frame(matrix(ncol = 3, nrow = 12))
colnames(forplotneg) = c('Accuracy', 'Balanced Accuracy', 'Precision')
rownames(forplotneg) = c('ms5_TFmRNA', 'ms5_combine', 'ms5', 'ms10_TFmRNA', 'ms10_combine', 'ms10', 
                         'ms15_TFmRNA', 'ms15_combine', 'ms15', 'ms20_TFmRNA', 'ms20_combine', 'ms20')
forplotneg['ms10_TFmRNA', 'Accuracy'] = cm_ms10_TFmRNA$overall['Accuracy']
forplotneg['ms10_combine', 'Accuracy'] = cm_ms10_combine$overall['Accuracy']
forplotneg['ms10', 'Accuracy'] = cm_ms10$overall['Accuracy']
forplotneg['ms15_TFmRNA', 'Accuracy'] = cm_ms15_TFmRNA$overall['Accuracy']
forplotneg['ms15_combine', 'Accuracy'] = cm_ms15_combine$overall['Accuracy']
forplotneg['ms15', 'Accuracy'] = cm_ms15$overall['Accuracy']

forplotneg['ms10_TFmRNA', 'Balanced Accuracy'] = cm_ms10_TFmRNA$byClass['Balanced Accuracy']
forplotneg['ms10_combine', 'Balanced Accuracy'] = cm_ms10_combine$byClass['Balanced Accuracy']
forplotneg['ms10', 'Balanced Accuracy'] = cm_ms10$byClass['Balanced Accuracy']
forplotneg['ms15_TFmRNA', 'Balanced Accuracy'] = cm_ms15_TFmRNA$byClass['Balanced Accuracy']
forplotneg['ms15_combine', 'Balanced Accuracy'] = cm_ms15_combine$byClass['Balanced Accuracy']
forplotneg['ms15', 'Balanced Accuracy'] = cm_ms15$byClass['Balanced Accuracy']

forplotneg['ms10_TFmRNA', 'Precision'] = cm_ms10_TFmRNA$byClass['Precision']
forplotneg['ms10_combine', 'Precision'] = cm_ms10_combine$byClass['Precision']
forplotneg['ms10', 'Precision'] = cm_ms10$byClass['Precision']
forplotneg['ms15_TFmRNA', 'Precision'] = cm_ms15_TFmRNA$byClass['Precision']
forplotneg['ms15_combine', 'Precision'] = cm_ms15_combine$byClass['Precision']
forplotneg['ms15', 'Precision'] = cm_ms15$byClass['Precision']


forplotneg['ms5_TFmRNA', 'Accuracy'] = cm_ms5_TFmRNA$overall['Accuracy']
forplotneg['ms5_combine', 'Accuracy'] = cm_ms5_combine$overall['Accuracy']
forplotneg['ms5', 'Accuracy'] = cm_ms5$overall['Accuracy']
forplotneg['ms20_TFmRNA', 'Accuracy'] = cm_ms20_TFmRNA$overall['Accuracy']
forplotneg['ms20_combine', 'Accuracy'] = cm_ms20_combine$overall['Accuracy']
forplotneg['ms20', 'Accuracy'] = cm_ms20$overall['Accuracy']

forplotneg['ms5_TFmRNA', 'Balanced Accuracy'] = cm_ms5_TFmRNA$byClass['Balanced Accuracy']
forplotneg['ms5_combine', 'Balanced Accuracy'] = cm_ms5_combine$byClass['Balanced Accuracy']
forplotneg['ms5', 'Balanced Accuracy'] = cm_ms5$byClass['Balanced Accuracy']
forplotneg['ms20_TFmRNA', 'Balanced Accuracy'] = cm_ms20_TFmRNA$byClass['Balanced Accuracy']
forplotneg['ms20_combine', 'Balanced Accuracy'] = cm_ms20_combine$byClass['Balanced Accuracy']
forplotneg['ms20', 'Balanced Accuracy'] = cm_ms20$byClass['Balanced Accuracy']

forplotneg['ms5_TFmRNA', 'Precision'] = cm_ms5_TFmRNA$byClass['Precision']
forplotneg['ms5_combine', 'Precision'] = cm_ms5_combine$byClass['Precision']
forplotneg['ms5', 'Precision'] = cm_ms5$byClass['Precision']
forplotneg['ms20_TFmRNA', 'Precision'] = cm_ms20_TFmRNA$byClass['Precision']
forplotneg['ms20_combine', 'Precision'] = cm_ms20_combine$byClass['Precision']
forplotneg['ms20', 'Precision'] = cm_ms20$byClass['Precision']

forplotneg$ms = c('Model size: 5', 'Model size: 5', 'Model size: 5', 
                  'Model size: 10', 'Model size: 10', 'Model size: 10', 
                  'Model size: 15', 'Model size: 15', 'Model size: 15',
                  'Model size: 20', 'Model size: 20', 'Model size: 20')
forplotneg$ms = factor(forplotneg$ms, levels = c('Model size: 5', 'Model size: 10', 'Model size: 15', 'Model size: 20'))
forplotneg$TFAOpt = c('TFmRNA', 'combine', '', 'TFmRNA', 'combine', '', 'TFmRNA', 'combine', '', 'TFmRNA', 'combine', '')
forplotneg$TFAOpt = factor(forplotneg$TFAOpt, levels = c('TFmRNA', '', 'combine'))

forplotneg_melt = reshape2::melt(forplotneg)
library(ggplot2)
# library(RColorBrewer)
# brewer.pal(3, 'Set3')
foxo1_neg_plot = ggplot(data=forplotneg_melt, aes(x=variable, y=value, fill=TFAOpt)) +
  geom_bar(stat="identity", position=position_dodge(), color="black", width = 0.7)+
  facet_wrap(~ms) + 
  scale_fill_manual(limits = c('TFmRNA', '', 'combine'),
                    labels = c('TFmRNA', 'TFA', 'combine'),
                    name = 'TFA options',
                    values = c("#8DD3C7", "#FFFFB3", "#BEBADA"))+
  theme_bw(base_size = 14) +
  xlab('') +
  ylab('Value')+ 
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle('Foxo1 negative targets')

library(cowplot)
foxo1.combined = plot_grid(foxo1_pos_plot, foxo1_neg_plot, align="hv", ncol = 2)
ggsave(
  paste0('compare_ms_tfa_Foxo1_posneg_remake.pdf'),
  plot = foxo1.combined,
  device = "pdf",
  width = 11,
  height = 6,
  dpi = 300
)


