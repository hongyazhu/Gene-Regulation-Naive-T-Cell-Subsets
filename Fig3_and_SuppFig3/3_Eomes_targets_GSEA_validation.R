### get Eomes targets

num_nets = 2
ms = 10
bias = 50
net = read.table(paste0('chip_ko_network_atac_ms', ms, '_bias', bias, '_maxComb_cut01_sharedbyMorethan', num_nets, 'nets_sp.tsv'),
                 header = T)

net_Eomes = net[net$TF == 'Eomes',]
write.table(net_Eomes, 'net_Eomes.tsv', sep = '\t', quote = F, row.names = F)


### make network plot for Eomes (Fig. 3B)
# similar plots generated for Foxo1, Tbx21, Bach2 (Supp Fig. 3D-F)

atac_sp = read.table("inputs/prior/output_prior/prior_atac_b_sp.tsv", header = T)

atac_sp_tmp = atac_sp
atac_sp_tmp$Regulation = paste0(atac_sp_tmp$TF, '->',  atac_sp_tmp$Target)
atac_sp_tmp_tf = atac_sp_tmp[atac_sp_tmp$TF %in% unique(net$TF), ]
atac_sp_tmp_tftarget = atac_sp_tmp_tf[atac_sp_tmp_tf$Target %in% unique(net$Target), ]
net_tmp = net
net_tmp$Regulation = paste0(net_tmp$TF, '->',  net_tmp$Target)
net_tmp$If_prior_support <- ifelse(net_tmp$Regulation %in% atac_sp_tmp_tftarget$Regulation,"1-Supported by ATAC-seq","2-Not supported by ATAC-seq")
net_tmp$Direction = ifelse(net_tmp$SignedQuantile == 1, 'Activation', 'Repression')

net_Eomes = net_tmp[net_tmp$TF == 'Eomes',]
table(net_Eomes$Direction)
#Activation Repression 
#79         23 
table(net_Eomes$If_prior_support)
#1-Supported by ATAC-seq 2-Not supported by ATAC-seq 
#92                          10 
nodes = data.frame(gene = unique(c(net_Eomes$TF, net_Eomes$Target)), type = NA)
nodes$type = ifelse(nodes$gene %in% unique(net$TF), 'TF', 'Target')
g = graph_from_data_frame(net_Eomes, directed = TRUE, vertices = nodes)
lay = create_layout(g, layout = "fr")
ggraph(g, layout = 'graphopt') + 
  geom_edge_link(aes(start_cap = label_rect(node1.name), end_cap = label_rect(node2.name), linetype = If_prior_support,color = Direction), 
                 arrow = arrow(type = "closed", length = unit(1, 'mm'))) +
  geom_node_text(aes(label = name, color = type)) +
  scale_color_manual(limits = c('TF', 'Target'),
                     values = c('#5b0e2d', '#ffa781')) +
  theme_graph() 
ggsave('Eomes_targets.svg',
       plot = last_plot(),
       device = "svg",
       width = 14,
       height = 10,
       dpi = 300
)

### GSEA of Eomes targets using ImmuneSigDB and GO Molecular Functions were performed with tfTarget_GSEA.sh from https://github.com/emiraldi/infTRN_lassoStARS

### Validating Eomes targets using an orthogonal dataset (Eomes over-expression RNAseq: GSE124913) (Fig. 3E, Supp Fig. 3G)
library(readxl)
de = read_excel("GSE124913_RNA-Seq_CD8SP_EomesTGvsWT_DifferentialAnalysisResults.xlsx", 
                sheet = "Sheet1")
de = data.frame(de)

ms = 10
bias = 50
num_nets = 2
nwFile = paste0('runs/sharedbyMorethan', num_nets, 'netsFrom5seeds/chip_ko_network_atac_ms', ms, '_bias', bias, '_maxComb_cut01_sharedbyMorethan', num_nets, 'nets_sp.tsv')
nwTab = read.table(nwFile, sep = '\t', header = T)

nw_eomes = nwTab[nwTab$TF == 'Eomes',]
nw_eomes_pos = nw_eomes[nw_eomes$SignedQuantile > 0, ]
nw_eomes_neg = nw_eomes[nw_eomes$SignedQuantile < 0, ]
nw_eomes_postargets = nw_eomes_pos$Target
nw_eomes_negtargets = nw_eomes_neg$Target

# select genes used for network inference as background
pool_genes = read.table('genes_DE/deg_padj01.txt')$V1
pool_regs = read.table('TF/regulators.txt')$V1
nwPool = unique(c(pool_genes, pool_regs))

de = de[de$Gene.Name %in% nwPool,]

de$if_target = 'Background'

postarget_rows = c()
for (i in 1:nrow(de)){
  if (de[i,'Gene.Name'] %in% nw_eomes_postargets)
    postarget_rows = c(postarget_rows, i)
}
de[postarget_rows, 'if_target'] = 'Postively regulated'

negtarget_rows = c()
for (i in 1:nrow(de)){
  if (de[i,'Gene.Name'] %in% nw_eomes_negtargets)
    negtarget_rows = c(negtarget_rows, i)
}
de[negtarget_rows, 'if_target'] = 'Negatively regulated'


de_pos = de[de$if_target == 'Postively regulated',]
de_neg = de[de$if_target == 'Negatively regulated',]
de_bg_all = de[de$if_target == 'Background',]

de_pos_genes = de_pos$Gene.Name
de_neg_genes = de_neg$Gene.Name
de_both_genes = c(de_pos_genes, de_neg_genes)
de_bg_all_genes = de_bg_all$Gene.Name


# positive targets
wilcox_test_pos <- wilcox.test(de_pos$'Log2FC.WT.vs.EomesTG', de_bg_all$'Log2FC.WT.vs.EomesTG')
p_value_pos <- wilcox_test_pos$p.value

pos_tmp = de_pos[, c('Gene.Name', 'Log2FC.WT.vs.EomesTG')]
pos_tmp$type = 'Eomes positive targets'
pos_bg_tmp = de_bg_all[, c('Gene.Name', 'Log2FC.WT.vs.EomesTG')]
pos_bg_tmp$type = 'Background genes'
pos_forplot = rbind(pos_tmp, pos_bg_tmp)
pos_forplot$Gene.Name = NULL
library(ggplot2)
ggplot(pos_forplot, aes(Log2FC.WT.vs.EomesTG, colour = type)) +
  stat_ecdf(linewidth=1)+ 
  theme_bw(base_size = 14) +
  xlim(-1, 1) +
  ylim(0, 1) + 
  xlab('Log2FC(WT/OE)') +
  ylab('Cumulative distribution') + 
  scale_color_manual(name = "",
                     values=c("#999999", "darkred"),
                     labels=c(paste0('Background genes'),
                              paste0('Eomes positive targets\n(n=', nrow(pos_tmp), ')')))+
  annotate("text", x=0.5, y=0.02, label=paste0('P-value = ', formatC(p_value_pos, format = "e", digits = 2)), size = 4)
ggsave(
  paste0('/eomes_pos_allbg.pdf'),
  plot = last_plot(),
  device = "pdf",
  width = 6,
  height = 3.5,
  dpi = 300
)

# negative targets
wilcox_test_neg <- wilcox.test(de_neg$'Log2FC.WT.vs.EomesTG', de_bg_all$'Log2FC.WT.vs.EomesTG')
p_value_neg <- wilcox_test_neg$p.value

neg_tmp = de_neg[, c('Gene.Name', 'Log2FC.WT.vs.EomesTG')]
neg_tmp$type = 'Eomes negitive targets'
neg_bg_tmp = de_bg_all[, c('Gene.Name', 'Log2FC.WT.vs.EomesTG')]
neg_bg_tmp$type = 'Background genes'
neg_forplot = rbind(neg_tmp, neg_bg_tmp)
neg_forplot$Gene.Name = NULL
library(ggplot2)
ggplot(neg_forplot, aes(Log2FC.WT.vs.EomesTG, colour = type)) +
  stat_ecdf(linewidth=1)+ 
  theme_bw(base_size = 14) +
  xlim(-1, 1) +
  ylim(0, 1) + 
  xlab('Log2FC(WT/OE)') +
  ylab('Cumulative distribution') + 
  scale_color_manual(name = "",
                     values=c("#999999", "darkred"),
                     labels=c(paste0('Background genes'),
                              paste0('Eomes negative targets\n(n=', nrow(neg_tmp), ')')))+
  annotate("text", x=0.5, y=0.02, label=paste0('P-value = ', formatC(p_value_neg, format = "e", digits = 2)), size = 4)
ggsave(
  paste0('eomes_neg_allbg.pdf'),
  plot = last_plot(),
  device = "pdf",
  width = 6,
  height = 3.5,
  dpi = 300
)


### validate Foxo1 targets using an orthogonal dataset (Foxo1 knockout RNAseq: GSE15037 profiles from GSE119085) (Supp Fig. 3H-I)

expr = read.csv('/workdir/hz543/projects/Inferelator/cd8/analyze_network/validate_identified_genes/foxo1/GSE119085_mouse-gpl1261-expr-subsetGSE15037.txt', header = T, sep = '\t')
expr = expr[, c(1,2,6,4)]
expr$'GeneName' = gsub(':.*', '', expr$Gene.Symbol..Gene.Title)
expr = expr[expr$GeneName != '---', ]
expr$fc = expr$GSM375667 - expr$GSM375666

### logfc of multiple probe sets for a single gene are taken mean
noneed_mean = expr[!expr$Gene.Symbol..Gene.Title %in% unique(expr$Gene.Symbol..Gene.Title[duplicated(expr$Gene.Symbol..Gene.Title)]),]
noneed_mean$mean_fc = noneed_mean$fc
noneed_mean_part = noneed_mean[,c(2,5,7)]

need_mean = expr[expr$Gene.Symbol..Gene.Title %in% unique(expr$Gene.Symbol..Gene.Title[duplicated(expr$Gene.Symbol..Gene.Title)]),]
need_mean_genes = unique(need_mean$Gene.Symbol..Gene.Title)
need_mean_new = data.frame(matrix(nrow = length(need_mean_genes), ncol = 3))
colnames(need_mean_new) = colnames(noneed_mean_part)

for (i in 1:length(need_mean_genes)){
  need_mean_gene = need_mean_genes[i]
  
  need_mean_gene_tab_tmp = need_mean[need_mean$Gene.Symbol..Gene.Title == need_mean_gene,]
  need_mean_new[i, 1] = unique(need_mean_gene_tab_tmp$Gene.Symbol..Gene.Title)
  need_mean_new[i, 2] = unique(need_mean_gene_tab_tmp$GeneName)
  need_mean_new[i, 3] = mean(need_mean_gene_tab_tmp$fc)
  
}

expr_m = rbind(noneed_mean_part, need_mean_new)

ms = 10
bias = 50
num_nets = 2
nwFile = paste0('/workdir/hz543/projects/Inferelator/cd8_rerun/analyze_network/runs/sharedbyMorethan', num_nets, 'netsFrom5seeds/chip_ko_network_atac_ms', ms, '_bias', bias, '_maxComb_cut01_sharedbyMorethan', num_nets, 'nets_sp.tsv')
nwTab = read.table(nwFile, sep = '\t', header = T)

nw_foxo1 = nwTab[nwTab$TF == 'Foxo1',]
nw_foxo1_pos = nw_foxo1[nw_foxo1$SignedQuantile > 0, ]
nw_foxo1_neg = nw_foxo1[nw_foxo1$SignedQuantile < 0, ]
nw_foxo1_postargets = nw_foxo1_pos$Target
nw_foxo1_negtargets = nw_foxo1_neg$Target


# select genes used for network inference as background
pool_genes = read.table('genes_DE/deg_padj01.txt')$V1
pool_regs = read.table('TF/regulators.txt')$V1
nwPool = unique(c(pool_genes, pool_regs))

included_rows = c()
for (i in 1:nrow(expr_m)){
  if (any(unlist(strsplit(expr_m[i,'GeneName'], " /// ")) %in% nwPool))
    included_rows = c(included_rows, i)
}
expr_m = expr_m[included_rows, ]

expr_m$if_target = 'Background'

postarget_rows = c()
for (i in 1:nrow(expr_m)){
  if (any(unlist(strsplit(expr_m[i,'GeneName'], " /// ")) %in% nw_foxo1_postargets))
    postarget_rows = c(postarget_rows, i)
}
expr_m[postarget_rows, 'if_target'] = 'Postively regulated'

negtarget_rows = c()
for (i in 1:nrow(expr_m)){
  if (any(unlist(strsplit(expr_m[i,'GeneName'], " /// ")) %in% nw_foxo1_negtargets))
    negtarget_rows = c(negtarget_rows, i)
}
expr_m[negtarget_rows, 'if_target'] = 'Negatively regulated'


expr_m$if_label = NA
for (i in 1:nrow(expr_m)){
  if (abs(expr_m[i, 'mean_fc']) > 2)
    expr_m[i, 'if_label'] = expr_m[i, 'GeneName']
}

expr_pos = expr_m[expr_m$if_target == 'Postively regulated',]
expr_neg = expr_m[expr_m$if_target == 'Negatively regulated',]
expr_bg_all = expr_m[expr_m$if_target == 'Background',]

exp_pos_genes = expr_pos$GeneName
exp_neg_genes = expr_neg$GeneName
exp_both_genes = c(exp_pos_genes, exp_neg_genes)


wilcox_test_pos <- wilcox.test(expr_pos$'mean_fc', expr_bg_all$'mean_fc')
p_value_pos <- wilcox_test_pos$p.value
wilcox_test_neg <- wilcox.test(expr_neg$'mean_fc', expr_bg_all$'mean_fc')
p_value_neg <- wilcox_test_neg$p.value

# positive targets

pos_tmp = expr_pos[, c('GeneName', 'mean_fc')]
pos_tmp$type = 'Foxo1 positive targets'
pos_bg_tmp = expr_bg_all[, c('GeneName', 'mean_fc')]
pos_bg_tmp$type = 'Background genes'
pos_forplot = rbind(pos_tmp, pos_bg_tmp)
pos_forplot$GeneName = NULL
library(ggplot2)
ggplot(pos_forplot, aes(mean_fc, colour = type)) +
  stat_ecdf(linewidth=1)+ 
  theme_bw(base_size = 14) +
  xlim(-1, 1) +
  ylim(0, 1) + 
  xlab('Log2FC(KO/WT)') +
  ylab('Cumulative distribution') + 
  scale_color_manual(name = "",
                     values=c("#999999", "darkblue"),
                     labels=c(paste0('Background genes'),
                              paste0('Foxo1 positive targets\n(n=', nrow(pos_tmp), ')')))+
  annotate("text", x=0.5, y=0.02, label=paste0('P-value = ', formatC(p_value_pos, format = "e", digits = 2)), size = 4)
# warning came from setting xlim and ylim
ggsave(
  paste0('foxo1_GSE119085_pos_allbg.pdf'),
  plot = last_plot(),
  device = "pdf",
  width = 6,
  height = 3.5,
  dpi = 300
)

# negative targets

neg_tmp = expr_neg[, c('GeneName', 'mean_fc')]
neg_tmp$type = 'Foxo1 negitive targets'
neg_bg_tmp = expr_bg_all[, c('GeneName', 'mean_fc')]
neg_bg_tmp$type = 'Background genes'
neg_forplot = rbind(neg_tmp, neg_bg_tmp)
neg_forplot$GeneName = NULL
library(ggplot2)
ggplot(neg_forplot, aes(mean_fc, colour = type)) +
  stat_ecdf(linewidth=1)+ 
  theme_bw(base_size = 14) +
  xlim(-1, 1) +
  ylim(0, 1) + 
  xlab('Log2FC(KO/WT)') +
  ylab('Cumulative distribution') + 
  scale_color_manual(name = "",
                     values=c("#999999", "darkblue"),
                     labels=c(paste0('Background genes'),
                              paste0('Foxo1 negative targets\n(n=', nrow(neg_tmp), ')')))+
  annotate("text", x=0.5, y=0.02, label=paste0('P-value = ', formatC(p_value_neg, format = "e", digits = 2)), size = 4)
ggsave(
  paste0('foxo1_GSE119085_neg_allbg.pdf'),
  plot = last_plot(),
  device = "pdf",
  width = 6,
  height = 3.5,
  dpi = 300
)

