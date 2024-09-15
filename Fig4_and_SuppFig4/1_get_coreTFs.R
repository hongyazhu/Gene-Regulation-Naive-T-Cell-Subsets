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

