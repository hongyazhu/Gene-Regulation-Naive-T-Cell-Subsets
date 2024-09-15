# getting core TFs for each naive subsets using GSEA of predicted TF targets, inspired by Miraldi et al., 2018.
# a TF is defined as a core TF for a certain naive subtype (e.g. VM) if its positive targets are enriched in that subtype (VM in this case), or its negative targets are enriched in the counterpart (TN in this case). Similar for adult vs. neonatal, and clean vs. dirty.

num_nets = 2
ms = 10
bias = 50
net = read.table(paste0('chip_ko_network_atac_ms', ms, '_bias', bias, '_maxComb_cut01_sharedbyMorethan', num_nets, 'nets_sp.tsv'),
                 header = T)

### get differential expression results 
source('Fig1_and_SuppFig1/4_DE.R')

resSig.MP_TN = resSig.TN_MP
resSig.MP_TN$log2FoldChange_rev = -resSig.MP_TN$log2FoldChange
ranks.MP_TN = resSig.MP_TN[order(resSig.MP_TN$log2FoldChange_rev),'log2FoldChange_rev']
names(ranks.MP_TN) = rownames(resSig.MP_TN[order(resSig.MP_TN$log2FoldChange_rev),])

ranks.neo_adult = resSig.neo_adult[order(resSig.neo_adult$log2FoldChange),'log2FoldChange']
names(ranks.neo_adult) = rownames(resSig.neo_adult[order(resSig.neo_adult$log2FoldChange),])

ranks.dirty_clean = resSig.dirty_clean[order(resSig.dirty_clean$log2FoldChange),'log2FoldChange']
names(ranks.dirty_clean) = rownames(resSig.dirty_clean[order(resSig.dirty_clean$log2FoldChange),])

# create list of TFs' targets, pos
net_pos = net[net$SignedQuantile == 1,]
TFlist_net_pos = unique(net_pos$TF)
TFtargetlist_net_pos = list()
for (i in 1:length(TFlist_net_pos)){
  TFtargetlist_net_pos[[i]] = net_pos[net_pos$TF == TFlist_net_pos[i], 'Target']
}
names(TFtargetlist_net_pos) <- TFlist_net_pos

library(fgsea)
library(qusage)

set.seed(42)
fgseaRes.MP_TN.net_pos <- fgsea(pathways = TFtargetlist_net_pos, 
                        stats    = ranks.MP_TN,
                        eps      = 0.0,
                        minSize  = 15,
                        maxSize  = 500)
set.seed(42)
fgseaRes.neo_adult.net_pos <- fgsea(pathways = TFtargetlist_net_pos, 
                            stats    = ranks.neo_adult,
                            eps      = 0.0,
                            minSize  = 15,
                            maxSize  = 500)
set.seed(42)
fgseaRes.dirty_clean.net_pos <- fgsea(pathways = TFtargetlist_net_pos, 
                              stats    = ranks.dirty_clean,
                              eps      = 0.0,
                              minSize  = 15,
                              maxSize  = 500)

fgseaSigRes.MP_TN.net_pos = fgseaRes.MP_TN.net_pos[fgseaRes.MP_TN.net_pos$padj < 0.05,]
fgseaSigRes.neo_adult.net_pos = fgseaRes.neo_adult.net_pos[fgseaRes.neo_adult.net_pos$padj < 0.05,]
fgseaSigRes.dirty_clean.net_pos = fgseaRes.dirty_clean.net_pos[fgseaRes.dirty_clean.net_pos$padj < 0.05,]



# create list of TFs' targets, neg
net_neg = net[net$SignedQuantile == -1,]
TFlist_net_neg = unique(net_neg$TF)
TFtargetlist_net_neg = list()
for (i in 1:length(TFlist_net_neg)){
  TFtargetlist_net_neg[[i]] = net_neg[net_neg$TF == TFlist_net_neg[i], 'Target']
}
names(TFtargetlist_net_neg) <- TFlist_net_neg

set.seed(42)
fgseaRes.MP_TN.net_neg <- fgsea(pathways = TFtargetlist_net_neg, 
                                stats    = ranks.MP_TN,
                                eps      = 0.0,
                                minSize  = 15,
                                maxSize  = 500)
set.seed(42)
fgseaRes.neo_adult.net_neg <- fgsea(pathways = TFtargetlist_net_neg, 
                                    stats    = ranks.neo_adult,
                                    eps      = 0.0,
                                    minSize  = 15,
                                    maxSize  = 500)
set.seed(42)
fgseaRes.dirty_clean.net_neg <- fgsea(pathways = TFtargetlist_net_neg, 
                                      stats    = ranks.dirty_clean,
                                      eps      = 0.0,
                                      minSize  = 15,
                                      maxSize  = 500)


fgseaRes.MP_TN.net_pos$type = 'MP_TN'
fgseaRes.MP_TN.net_pos$net = 'pos'
fgseaRes.MP_TN.net_pos$Direction = ifelse(fgseaRes.MP_TN.net_pos$NES < 0, 'TN', 'VM')
fgseaRes.neo_adult.net_pos$type = 'neo_adult'
fgseaRes.neo_adult.net_pos$net = 'pos'
fgseaRes.neo_adult.net_pos$Direction = ifelse(fgseaRes.neo_adult.net_pos$NES < 0, 'Adult', 'Neo')
fgseaRes.dirty_clean.net_pos$type = 'dirty_clean'
fgseaRes.dirty_clean.net_pos$net = 'pos'
fgseaRes.dirty_clean.net_pos$Direction = ifelse(fgseaRes.dirty_clean.net_pos$NES < 0, 'Clean', 'Dirty')

fgseaRes.MP_TN.net_neg$type = 'MP_TN'
fgseaRes.MP_TN.net_neg$net = 'neg'
fgseaRes.MP_TN.net_neg$Direction = ifelse(fgseaRes.MP_TN.net_neg$NES < 0, 'TN', 'VM')
fgseaRes.neo_adult.net_neg$type = 'neo_adult'
fgseaRes.neo_adult.net_neg$net = 'neg'
fgseaRes.neo_adult.net_neg$Direction = ifelse(fgseaRes.neo_adult.net_neg$NES < 0, 'Adult', 'Neo')
fgseaRes.dirty_clean.net_neg$type = 'dirty_clean'
fgseaRes.dirty_clean.net_neg$net = 'neg'
fgseaRes.dirty_clean.net_neg$Direction = ifelse(fgseaRes.dirty_clean.net_neg$NES < 0, 'Clean', 'Dirty')

fgseaRes.net_pos = rbind(fgseaRes.MP_TN.net_pos, fgseaRes.neo_adult.net_pos, fgseaRes.dirty_clean.net_pos)
#fgseaRes.net_pos$Direction = ifelse(fgseaRes.net_pos$padj > 0.05, 'Not Significant', fgseaRes.net_pos$Direction)
fgseaRes.net_pos = na.omit(fgseaRes.net_pos)
fgseaRes.net_pos = fgseaRes.net_pos[fgseaRes.net_pos$padj != 1]

fgseaRes.net_neg = rbind(fgseaRes.MP_TN.net_neg, fgseaRes.neo_adult.net_neg, fgseaRes.dirty_clean.net_neg)
#fgseaRes.net_neg$Direction = ifelse(fgseaRes.net_neg$padj > 0.05, 'Not Significant', fgseaRes.net_neg$Direction)
fgseaRes.net_neg = na.omit(fgseaRes.net_neg)
fgseaRes.net_neg = fgseaRes.net_neg[fgseaRes.net_neg$padj != 1]


df_pos_padj = data.frame(matrix(data = NA, ncol = length(unique(fgseaRes.net_pos$pathway)), nrow = 6))
rownames(df_pos_padj) = c('TN', 'VM', 'Adult', 'Neo', 'Clean', 'Dirty')
colnames(df_pos_padj) = unique(fgseaRes.net_pos$pathway)
for (c in colnames(df_pos_padj))
  for (r in rownames(df_pos_padj)){
    tmp = fgseaRes.net_pos[fgseaRes.net_pos$pathway == c & fgseaRes.net_pos$Direction == r, ]
    if (!is.na(df_pos_padj[r, c])){
      print(r)
      print(c)
    }
    if (nrow(tmp) == 0){
      df_pos_padj[r, c] = 1
    } else {
      df_pos_padj[r, c] = tmp$padj
    }
}

df_pos_NES = data.frame(matrix(data = NA, ncol = length(unique(fgseaRes.net_pos$pathway)), nrow = 6))
rownames(df_pos_NES) = c('TN', 'VM', 'Adult', 'Neo', 'Clean', 'Dirty')
colnames(df_pos_NES) = unique(fgseaRes.net_pos$pathway)
for (c in colnames(df_pos_NES))
  for (r in rownames(df_pos_NES)){
    tmp = fgseaRes.net_pos[fgseaRes.net_pos$pathway == c & fgseaRes.net_pos$Direction == r, ]
    if (!is.na(df_pos_NES[r, c])){
      print(r)
      print(c)
    }
    if (nrow(tmp) == 0){
      df_pos_NES[r, c] = 1
    } else {
      df_pos_NES[r, c] = tmp$NES
    }
  }


df_neg_padj = data.frame(matrix(data = NA, ncol = length(unique(fgseaRes.net_neg$pathway)), nrow = 6))
rownames(df_neg_padj) = c('TN', 'VM', 'Adult', 'Neo', 'Clean', 'Dirty')
colnames(df_neg_padj) = unique(fgseaRes.net_neg$pathway)
for (c in colnames(df_neg_padj))
  for (r in rownames(df_neg_padj)){
    tmp = fgseaRes.net_neg[fgseaRes.net_neg$pathway == c & fgseaRes.net_neg$Direction == r, ]
    if (!is.na(df_neg_padj[r, c])){
      print(r)
      print(c)
    }
    if (nrow(tmp) == 0){
      df_neg_padj[r, c] = 1
    } else {
      df_neg_padj[r, c] = tmp$padj
    }
  }

df_neg_NES = data.frame(matrix(data = NA, ncol = length(unique(fgseaRes.net_neg$pathway)), nrow = 6))
rownames(df_neg_NES) = c('TN', 'VM', 'Adult', 'Neo', 'Clean', 'Dirty')
colnames(df_neg_NES) = unique(fgseaRes.net_neg$pathway)
for (c in colnames(df_neg_NES))
  for (r in rownames(df_neg_NES)){
    tmp = fgseaRes.net_neg[fgseaRes.net_neg$pathway == c & fgseaRes.net_neg$Direction == r, ]
    if (!is.na(df_neg_NES[r, c])){
      print(r)
      print(c)
    }
    if (nrow(tmp) == 0){
      df_neg_NES[r, c] = 1
    } else {
      df_neg_NES[r, c] = tmp$NES
    }
  }


write.table(df_pos_padj, 'coreTFs/GSEA_on_TF_targets/df_pos_padj.txt', quote = F, sep = '\t')
write.table(df_pos_NES, 'coreTFs/GSEA_on_TF_targets/df_pos_NES.txt', quote = F, sep = '\t')
write.table(df_neg_padj, 'coreTFs/GSEA_on_TF_targets/df_neg_padj.txt', quote = F, sep = '\t')
write.table(df_neg_NES, 'coreTFs/GSEA_on_TF_targets/df_neg_NES.txt', quote = F, sep = '\t')


# GSEA on TF targets using ImmGen gene sets (Best et al., 2013) was performed with tfTarget_GSEA.sh in https://github.com/emiraldi/infTRN_lassoStARS

# make plots for core TFs and their ImmGen enrichment results 
# Fig 4A-B show TFs with ImmGen enrichments. Supp Fig 4A-D was done similarly showing all core TFs without filtering.

# core TFs, filtered by adjusted p value immgen != 0

up = read.table('df_pos_padj.txt')
dn = read.table('df_neg_padj.txt')

up_immgen = read.table('TFtargets_GSEA/GSEA/chip_ko_network_atac_ms10_bias50_maxComb_cut01_sharedbyMorethan2nets_ImmGenGeneSets_Praw0p1_dir_wCut0p0_minSet5/ImmGenGeneSets_fdr100_up_pval.txt')
dn_immgen = read.table('TFtargets_GSEA/GSEA/chip_ko_network_atac_ms10_bias50_maxComb_cut01_sharedbyMorethan2nets_ImmGenGeneSets_Praw0p1_dir_wCut0p0_minSet5/ImmGenGeneSets_fdr100_down_pval.txt')

up_6subtypes = up[c('TN', 'VM', 'Adult', 'Neo', 'Clean', 'Dirty'), ]
dn_6subtypes = dn[c('TN', 'VM', 'Adult', 'Neo', 'Clean', 'Dirty'), ]

cutoff = 0.001
up_6subtypes_cutoff = up_6subtypes[,colSums(up_6subtypes < cutoff) > 0]
dn_6subtypes_cutoff = dn_6subtypes[,colSums(dn_6subtypes < cutoff) > 0]

up_immgen_cutoff = up_immgen[,colSums(up_immgen) > 0]
dn_immgen_cutoff = dn_immgen[,colSums(dn_immgen) > 0]
dn_immgen_cutoff['Initial_cytokine_or_effector_response',] = 1 

include = intersect(c(colnames(up_6subtypes_cutoff), colnames(dn_6subtypes_cutoff)), c(colnames(up_immgen_cutoff), colnames(dn_immgen_cutoff)))
up_6subtypes_cutoff_include = up_6subtypes_cutoff[, colnames(up_6subtypes_cutoff) %in% include]
dn_6subtypes_cutoff_include = dn_6subtypes_cutoff[, colnames(dn_6subtypes_cutoff) %in% include]
up_immgen_cutoff_include = up_immgen_cutoff[, colnames(up_immgen_cutoff) %in% include]
dn_immgen_cutoff_include = dn_immgen_cutoff[, colnames(dn_immgen_cutoff) %in% include]

# find smallers value larger than 1
min(up)
min(dn)
min(up_immgen_cutoff_include)
min(dn_immgen_cutoff_include)
min(up_immgen_cutoff_include[up_immgen_cutoff_include != 0])
min(dn_immgen_cutoff_include[dn_immgen_cutoff_include != 0])

library(dplyr)
#up_6subtypes_cutoff_include_replace0 = up_6subtypes_cutoff_include %>% mutate_all(~replace(., . == 0, 1e-10))
#dn_6subtypes_cutoff_include_replace0 = dn_6subtypes_cutoff_include %>% mutate_all(~replace(., . == 0, 1e-10))
up_immgen_cutoff_include_replace0 = up_immgen_cutoff_include %>% mutate_all(~replace(., . == 0, 1e-13))
dn_immgen_cutoff_include_replace0 = dn_immgen_cutoff_include %>% mutate_all(~replace(., . == 0, 1e-13))
### same as coreTFs plots!

# log
library(reshape2)
up_6subtypes_cutoff_include_replace0_plot = melt(as.matrix(-log10(up_6subtypes_cutoff_include)))
up_6subtypes_cutoff_include_replace0_plot$Direction = 'Up'
dn_6subtypes_cutoff_include_replace0_plot = melt(as.matrix(-log10(dn_6subtypes_cutoff_include)))
dn_6subtypes_cutoff_include_replace0_plot$Direction = 'Down'
up_immgen_cutoff_include_replace0_plot = melt(as.matrix(-log10(up_immgen_cutoff_include_replace0)))
up_immgen_cutoff_include_replace0_plot$Direction = 'Up'
dn_immgen_cutoff_include_replace0_plot = melt(as.matrix(-log10(dn_immgen_cutoff_include_replace0)))
dn_immgen_cutoff_include_replace0_plot$Direction = 'Down'

forplot = rbind(up_6subtypes_cutoff_include_replace0_plot, dn_6subtypes_cutoff_include_replace0_plot,
                up_immgen_cutoff_include_replace0_plot, dn_immgen_cutoff_include_replace0_plot)

forplot$Var1 = gsub('_', ' ', forplot$Var1)

forplot$Var1 = factor(forplot$Var1, levels = c("VM", "Neo",  "Dirty", "Adult", "TN", "Clean",
                                               "Cell cycle and division",
                                               "Short term effector and memory",
                                               'Naive or late effector or memory',
                                               'Short term effector or memory',
                                               'Late effector or memory',
                                               'Preparation for cell division',
                                               'Naive and late memory',
                                               'Early effector late memory',
                                               'Memory precursor',
                                               'Initial cytokine or effector response'))


for (i in 1:nrow(forplot)){
  if (forplot[i, 'Direction'] == 'Down')
    forplot[i, 'value'] = -forplot[i, 'value']
}
library(ggplot2)

genes_manual_include_tmp = c("Mef2c", "Tcf12", "Zfp212", "Batf3", "Nr4a2", "Jun", "Jdp2", "Rela", "Zfp768", "Nr3c1", "Fosl2", "Zfp41",   # VM
                             "Mga", "Foxj2", "Etv6", "Xbp1",    # Neo
                             "Trafd1", "Gata1", "Ikzf3", "Mxd1",   # Dirty
                             "Gata3", "Mef2a", "Atf3", "Phf21a", "Zfp637", "Gata6",   # VM Neo
                             "Irf8", "Klf1", "Rxra", "Rora", "Irf4", "Nfya", "E2f8", "Mxd3", "Pbx1", "E2f2", "Foxm1", "Myb", "E2f7", "Klf12", "E2f1", "Tfdp1", "Zbtb18", "Sp6", "Mybl2", "Zeb2", "Bhlhe40", "Tcf3", "Zfp367", "Bach2",   # VM dirty
                             "Irf7", "Stat2", "Irf9", "Stat3",  # Neo Dirty
                             "Smad3", "Runx2", "Tbx21", "Stat4", "Klf6", "Arid3a", "Eomes", "Prdm1", "Tbx6", "Nr1d2",   # Neo VM Dirty

                             "Rarg", "Ets2", "Patz1", "Tgif1", "Usf2",     # TN
                             "Dnmt1", "Ctcf",   # Adult
                             "Hsf1", "Pax5", "Klf10",  # Clean
                             "Zfp281", "Tcf7", "Irf3", "Klf11", "Thra", "Srebf2", "Zbtb6", # TN Clean
                             "Zfp422", # Adult Clean
                             "Sox4", "Foxo1", "Ikzf1", "Egr3", "Lef1"    # TN Adult Clean
)

backgroundtiles = forplot
backgroundtiles$value = 0
for (i in 1:nrow(forplot)){
  backgroundtiles[backgroundtiles$Var1 == forplot[i, 'Var1'] & backgroundtiles$Var2 == forplot[i, 'Var2'] & backgroundtiles$Direction == forplot[i, 'Direction'],'value'] = forplot[i, 'value']
}
backgroundtiles$Var2 = factor(backgroundtiles$Var2, levels = rev(genes_manual_include_tmp))
backgroundtiles_reorder = backgroundtiles[order(abs(backgroundtiles$value), decreasing = F),]

# add group label: https://stackoverflow.com/questions/23207878/ggplot2-group-x-axis-discrete-values-into-subgroups
backgroundtiles_reorder_group = backgroundtiles_reorder[backgroundtiles_reorder$Var2 %in% genes_manual_include, ]
backgroundtiles_reorder_group$Group = NA
for (i in 1:nrow(backgroundtiles_reorder_group)){
  if (backgroundtiles_reorder_group[i, 'Var2'] %in% c("Mef2c", "Tcf12", "Zfp212", "Batf3", "Nr4a2", "Jun", "Jdp2", "Rela", "Zfp768", "Nr3c1", "Fosl2", "Zfp41")){
    backgroundtiles_reorder_group[i, 'Group'] = 'VM'
  } else if (backgroundtiles_reorder_group[i, 'Var2'] %in% c("Mga", "Foxj2", "Etv6", "Xbp1")){
    backgroundtiles_reorder_group[i, 'Group'] = 'Neo'
  } else if (backgroundtiles_reorder_group[i, 'Var2'] %in% c("Trafd1", "Gata1", "Ikzf3", "Mxd1")){
    backgroundtiles_reorder_group[i, 'Group'] = 'Dirty'
  } else if (backgroundtiles_reorder_group[i, 'Var2'] %in% c("Gata3", "Mef2a", "Atf3", "Phf21a", "Zfp637", "Gata6")){
    backgroundtiles_reorder_group[i, 'Group'] = 'VM, Neo'
  } else if (backgroundtiles_reorder_group[i, 'Var2'] %in% c("Irf8", "Klf1", "Rxra", "Rora", "Irf4", "Nfya", "E2f8", "Mxd3", "Pbx1", "E2f2", "Foxm1", "Myb", "E2f7", "Klf12", "E2f1", "Tfdp1", "Zbtb18", "Sp6", "Mybl2", "Zeb2", "Bhlhe40", "Tcf3", "Zfp367", "Bach2")){
    backgroundtiles_reorder_group[i, 'Group'] = 'VM, Dirty'
  } else if (backgroundtiles_reorder_group[i, 'Var2'] %in% c("Stat3", "Irf7", "Irf9", "Stat2")){
    backgroundtiles_reorder_group[i, 'Group'] = 'Neo, Dirty'
  } else if (backgroundtiles_reorder_group[i, 'Var2'] %in% c("Smad3", "Runx2", "Tbx21", "Stat4", "Klf6", "Arid3a", "Eomes", "Prdm1", "Tbx6", "Nr1d2")){
    backgroundtiles_reorder_group[i, 'Group'] = 'VM, Neo, Dirty'
  } else if (backgroundtiles_reorder_group[i, 'Var2'] %in% c("Tgif1", "Usf2", "Patz1", "Rarg", "Ets2")){
    backgroundtiles_reorder_group[i, 'Group'] = 'TN'
  } else if (backgroundtiles_reorder_group[i, 'Var2'] %in% c("Dnmt1", "Ctcf")){
    backgroundtiles_reorder_group[i, 'Group'] = 'Adult'
  } else if (backgroundtiles_reorder_group[i, 'Var2'] %in% c("Hsf1", "Klf10", "Pax5")){
    backgroundtiles_reorder_group[i, 'Group'] = 'Clean'
  } else if (backgroundtiles_reorder_group[i, 'Var2'] %in% c("Zbtb6", "Thra", "Irf3", "Klf11", "Srebf2", "Tcf7", "Zfp281")){
    backgroundtiles_reorder_group[i, 'Group'] = 'TN, Clean'
  } else if (backgroundtiles_reorder_group[i, 'Var2'] %in% c("Zfp422")){
    backgroundtiles_reorder_group[i, 'Group'] = 'Adult, Clean'
  } else if (backgroundtiles_reorder_group[i, 'Var2'] %in% c("Sox4", "Foxo1", "Ikzf1", "Egr3", "Lef1")){
    backgroundtiles_reorder_group[i, 'Group'] = 'TN, Adult, Clean'
  } else {
    print (backgroundtiles_reorder_group[i, 'Var2'])
  }
}
backgroundtiles_reorder_group$Group = factor(backgroundtiles_reorder_group$Group, levels = rev(c('VM', 'Neo', 'Dirty', 'VM, Neo', 'VM, Dirty', 'Neo, Dirty', 'VM, Neo, Dirty', 
                                                                                                 'TN', 'Adult', 'Clean', 'TN, Clean', 'Adult, Clean', 'TN, Adult, Clean')))

backgroundtiles_reorder_group2 = backgroundtiles_reorder_group
backgroundtiles_reorder_group2$Group2 = NA
for (i in 1:nrow(backgroundtiles_reorder_group2)){
  if (backgroundtiles_reorder_group2[i, 'Var1'] %in% c("Cell cycle and division",
                                                       "Short term effector and memory",
                                                       'Naive or late effector or memory',
                                                       'Short term effector or memory',
                                                       'Late effector or memory',
                                                       'Initial cytokine or effector response')){
    backgroundtiles_reorder_group2[i, 'Group2'] = 'Effector gene sets'
  } else if (backgroundtiles_reorder_group2[i, 'Var1'] %in% c('Preparation for cell division',
                                                              'Naive and late memory',
                                                              'Early effector late memory',
                                                              'Memory precursor')){
    backgroundtiles_reorder_group2[i, 'Group2'] = 'Memory gene sets'
  }
}
backgroundtiles_reorder_group2_immgen = na.omit(backgroundtiles_reorder_group2)
backgroundtiles_reorder_group2_immgen$Var1 = factor(backgroundtiles_reorder_group2_immgen$Var1, levels = c("Cell cycle and division",
                                                                                                           "Short term effector and memory",
                                                                                                           'Naive or late effector or memory',
                                                                                                           'Short term effector or memory',
                                                                                                           'Late effector or memory',
                                                                                                           'Initial cytokine or effector response',
                                                                                                           'Preparation for cell division',
                                                                                                           'Naive and late memory',
                                                                                                           'Early effector late memory',
                                                                                                           'Memory precursor'))

# add gaps, colored
library(ggh4x)
library(colorjam)
library(jamba)
backgroundtiles_reorder_group$Group = factor(backgroundtiles_reorder_group$Group, levels = c('VM', 'Neo', 'Dirty', 'VM, Neo', 'VM, Dirty', 'Neo, Dirty', 'VM, Neo, Dirty', 
                                                                                             'TN', 'Adult', 'Clean', 'TN, Clean', 'Adult, Clean', 'TN, Adult, Clean'))
strip <- strip_themed(text_y = elem_list_text(face = rep("bold", length(unique(backgroundtiles_reorder_group$Group))),
                                              color = c('#004D7F', '#017100', '#FF9300', 
                                                        blend_colors(c('#004D7F', '#017100')), blend_colors(c('#004D7F', '#FF9300')), blend_colors(c('#017100', '#FF9300')),
                                                        blend_colors(c('#004D7F', '#017100', '#FF9300')),
                                                        '#56C1FF','#88FA4E', '#FFFC66', blend_colors(c('#56C1FF', '#FFFC66')), blend_colors(c('#88FA4E', '#FFFC66')), 
                                                        blend_colors(c('#56C1FF', '#88FA4E', '#FFFC66')))),
                      background_y = elem_list_rect(fill = 'gray95'))


backgroundtiles_reorder_group_sepsubtype = backgroundtiles_reorder_group
backgroundtiles_reorder_group_sepsubtype$subtype = NA
for (i in 1:nrow(backgroundtiles_reorder_group_sepsubtype)){
  if (backgroundtiles_reorder_group_sepsubtype[i, 'Var1'] %in% c('VM', 'TN')){
    backgroundtiles_reorder_group_sepsubtype[i, 'subtype'] = 'Phenotype'
  } else if (backgroundtiles_reorder_group_sepsubtype[i, 'Var1'] %in% c('Adult', 'Neo')){
    backgroundtiles_reorder_group_sepsubtype[i, 'subtype'] = 'Age'
  } else if (backgroundtiles_reorder_group_sepsubtype[i, 'Var1'] %in% c('Dirty', 'Clean')){
    backgroundtiles_reorder_group_sepsubtype[i, 'subtype'] = 'Microbiota'
  }
}
backgroundtiles_reorder_group_sepsubtype$subtype = factor(backgroundtiles_reorder_group_sepsubtype$subtype, levels = c('Phenotype', 'Age', 'Microbiota'))
strip <- strip_themed(text_y = elem_list_text(face = rep("bold", length(unique(backgroundtiles_reorder_group$Group))),
                                              color = c('#004D7F', '#017100', '#FF9300', 
                                                        blend_colors(c('#004D7F', '#017100')), blend_colors(c('#004D7F', '#FF9300')), blend_colors(c('#017100', '#FF9300')),
                                                        blend_colors(c('#004D7F', '#017100', '#FF9300')),
                                                        '#56C1FF','#88FA4E', '#FFFC66', blend_colors(c('#56C1FF', '#FFFC66')), blend_colors(c('#88FA4E', '#FFFC66')), 
                                                        blend_colors(c('#56C1FF', '#88FA4E', '#FFFC66')))),
                      background_y = elem_list_rect(fill = 'gray95'),
                      background_x = element_blank(),
                      text_x = element_blank())
#textlayer <- data.frame(color=c('#004D7F', '#56C1FF', '#017100','#88FA4E', '#FF9300', '#FFFC66'),
#                        type=c("VM", "TN", "Neo", "Adult",  "Dirty", "Clean"))
backgroundtiles_reorder_group_sepsubtype_subtypesonly = backgroundtiles_reorder_group_sepsubtype[backgroundtiles_reorder_group_sepsubtype$Var1 %in% c("VM", "TN", "Neo", "Adult",  "Dirty", "Clean"), ]
min(backgroundtiles_reorder_group_sepsubtype_subtypesonly$value)
#[1] -10.51353
max(backgroundtiles_reorder_group_sepsubtype_subtypesonly$value)
#[1] 8.322823
p = ggplot(data=backgroundtiles_reorder_group_sepsubtype_subtypesonly) +
  geom_tile(aes(x = Var1, y = Var2, fill = value), color = 'gray60') +
  labs(size='-log10(adj p-val)') +
  xlab('') +
  ylab('')+
  #scale_fill_gradient2(low = "#075AFF",
  #                     mid = "#FFFFCC",
  #                     high = "#FF0000") +
  scale_fill_gradientn(limits = c(-32, 32), colors = rev(hcl.colors(20, "RdBu")), name = '-log10(adj p-val)')+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.ontop = TRUE,
        panel.background = element_rect(fill = "transparent"),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        ggh4x.axis.nestline = element_line(linetype = 1),
        axis.text.x = element_text(colour = c('#004D7F', '#56C1FF','#017100', '#88FA4E', '#FF9300', '#FFFC66'),
                                   face = "italic"),
        #text = element_text(size=14)
  ) +
  #geom_text(data=textlayer, aes(colour = color, label = type, x=type), y=-1, hjust=.5, vjust = 1.1)+
  #scale_colour_manual(values=c('#004D7F', '#56C1FF', '#017100','#88FA4E', '#FF9300', '#FFFC66'),
  #                    limits=c("VM", "TN", "Neo", "Adult",  "Dirty", "Clean") )+
  scale_x_discrete(#limits = c("VM", "TN", "Neo", "Adult",  "Dirty", "Clean"),
    position = "top")  +
  #guides(y = ggh4x::guide_axis_nested(delim = "&"))+ 
  theme(strip.placement = "outside", 
        strip.text.y.left = element_text(angle = 0))+
  facet_grid2(Group~subtype, scales = 'free', space = 'free', switch = "y", strip = strip)
gt <- ggplotGrob(p)
gt$grobs[[78]]$children[[2]]$grobs[[1]]$children[[1]]$gp$col = c("#017100", "#88FA4E")
gt$grobs[[79]]$children[[2]]$grobs[[1]]$children[[1]]$gp$col = c("#FF9300", "#FFFC66")
pdf('coreTFs_sepsubtypes.pdf', 
    height = 13, width = 6)
grid::grid.draw(gt)
dev.off()


backgroundtiles_reorder_group2 = backgroundtiles_reorder_group
backgroundtiles_reorder_group2$Group2 = NA
for (i in 1:nrow(backgroundtiles_reorder_group2)){
  if (backgroundtiles_reorder_group2[i, 'Var1'] %in% c("Cell cycle and division",
                                                       "Short term effector and memory",
                                                       'Naive or late effector or memory',
                                                       'Short term effector or memory',
                                                       'Late effector or memory',
                                                       'Initial cytokine or effector response')){
    backgroundtiles_reorder_group2[i, 'Group2'] = 'Effector gene sets'
  } else if (backgroundtiles_reorder_group2[i, 'Var1'] %in% c('Preparation for cell division',
                                                              'Naive and late memory',
                                                              'Early effector late memory',
                                                              'Memory precursor')){
    backgroundtiles_reorder_group2[i, 'Group2'] = 'Memory gene sets'
  }
}
backgroundtiles_reorder_group2_immgen = na.omit(backgroundtiles_reorder_group2)
backgroundtiles_reorder_group2_immgen$Var1 = factor(backgroundtiles_reorder_group2_immgen$Var1, levels = c("Cell cycle and division",
                                                                                                           "Short term effector and memory",
                                                                                                           'Naive or late effector or memory',
                                                                                                           'Short term effector or memory',
                                                                                                           'Late effector or memory',
                                                                                                           'Initial cytokine or effector response',
                                                                                                           'Preparation for cell division',
                                                                                                           'Naive and late memory',
                                                                                                           'Early effector late memory',
                                                                                                           'Memory precursor'))
ggplot(backgroundtiles_reorder_group2_immgen, aes(x = Var1, y = Var2)) +
  geom_point(aes(size = ifelse(value==0, NA, abs(value)), color = Direction)) +
  geom_tile(fill = 'transparent', color = 'gray85') +
  labs(size='-log10(p-val)') +
  xlab('') +
  ylab('')+
  scale_colour_manual(values = c(hcl.colors(20, "PRGn")[5], hcl.colors(20, "PRGn")[20-5]),
                      limits = c('Up', 'Down'),
                      labels = c('Enriched on\npositive targets',
                                 'Enriched on\nnegative targets')) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.ontop = TRUE,
        panel.background = element_rect(fill = "transparent"),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        ggh4x.axis.nestline = element_line(linetype = 1),
        axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
        strip.background = element_rect(color = 'transparent', fill = 'gray95')#,
        #axis.text=element_blank()#, 
        #text = element_text(size=14)
  ) +
  scale_x_discrete(position = "bottom")  +
  #  guides(x = ggh4x::guide_axis_nested(delim = "&"))+ 
  scale_size(range = c(0.5,4))+
  #facet_grid( ~ , scales = 'free_x', space = "free_x")+
  facet_grid2(Group ~ Group2, scales = 'free', space = 'free', switch = "y", strip = strip)+ 
  theme(strip.placement = "outside", 
        strip.text.y.left = element_text(angle = 0), 
        strip.text.x = element_text(face = 'bold.italic'), 
        panel.border = element_blank()) 
ggsave(
  paste0('coreTFs_immgen.pdf'),
  plot = last_plot(),
  device = "pdf",
  width = 8,
  height = 14,
  dpi = 300
)



