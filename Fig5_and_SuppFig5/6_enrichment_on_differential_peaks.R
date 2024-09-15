# perform enrichment on genes proximal to differential peaks of Eomes in WT and Lin28b (neonatal like) cells, with ATAC-seq peaks as backrgound.
# (Fig. 5E-F)

common = unique(read.delim('MACS3_bdgdiff/EOMES_bdgdiff_WT_vs_lin28Tg_c3.0_common_peak_anno.txt', sep = '\t')$Gene.Name)
wt = unique(read.delim('MACS3_bdgdiff/EOMES_bdgdiff_WT_vs_lin28Tg_c3.0_cond1_peak_anno.txt', sep = '\t')$Gene.Name)
lin28 = unique(read.delim('MACS3_bdgdiff/EOMES_bdgdiff_WT_vs_lin28Tg_c3.0_cond2_peak_anno.txt', sep = '\t')$Gene.Name)

# background
atac = unique(read.delim('ATAC.annoPeak.txt', sep = '\t')$Gene.Name)


library(fgsea)
library(qusage)

immgen_genesets = qusage::read.gmt('cluster_list_geneSetNames.gmt')
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

h_genesets = qusage::read.gmt('MSigDB_gene_sets/mouse/mh.all.v2023.1.Mm.symbols.gmt')

# common peaks
resImmGen = fora(pathways = immgen_genesets, genes = common, universe = atac, minSize = 1, maxSize = Inf)
resImmGen = apply(resImmGen,2,as.character)
write.table(resImmGen, file = paste0('bdgdiff/ora/ora_ImmGen_common.txt'), quote = F, row.names = F, sep = '\t')
resImmGen = resImmGen[,1:5]
write.table(resImmGen, file = paste0('bdgdiff/ora/ora_ImmGen_common_noLeadingEdge.txt'), quote = F, row.names = F, sep = '\t')

resHallmark = fora(pathways = h_genesets, genes = common, universe = atac, minSize = 1, maxSize = Inf)
resHallmark = resHallmark[resHallmark$padj < 0.05]
resHallmark = apply(resHallmark,2,as.character)
write.table(resHallmark, paste0('bdgdiff/ora/ora_Hallmark_common.txt'), quote = F, row.names = F)
resHallmark = resHallmark[,1:5]
write.table(resHallmark, file = paste0('bdgdiff/ora/ora_Hallmark_common_noLeadingEdge.txt'), quote = F, row.names = F)


# peaks specific to lin28 cells
resImmGen = fora(pathways = immgen_genesets, genes = lin28, universe = atac, minSize = 1, maxSize = Inf)
resImmGen = apply(resImmGen,2,as.character)
write.table(resImmGen, file = paste0('bdgdiff/ora/ora_ImmGen_lin28.txt'), quote = F, row.names = F, sep = '\t')
resImmGen = resImmGen[,1:5]
write.table(resImmGen, file = paste0('bdgdiff/ora/ora_ImmGen_lin28_noLeadingEdge.txt'), quote = F, row.names = F, sep = '\t')

resHallmark = fora(pathways = h_genesets, genes = lin28, universe = atac, minSize = 1, maxSize = Inf)
resHallmark = resHallmark[resHallmark$padj < 0.05]
resHallmark = apply(resHallmark,2,as.character)
write.table(resHallmark, paste0('bdgdiff/ora/ora_Hallmark_lin28.txt'), quote = F, row.names = F)
resHallmark = resHallmark[,1:5]
write.table(resHallmark, file = paste0('bdgdiff/ora/ora_Hallmark_lin28_noLeadingEdge.txt'), quote = F, row.names = F)


# peaks specific to wt cells
resImmGen = fora(pathways = immgen_genesets, genes = wt, universe = atac, minSize = 1, maxSize = Inf)
resImmGen = apply(resImmGen,2,as.character)
write.table(resImmGen, file = paste0('bdgdiff/ora/ora_ImmGen_wt.txt'), quote = F, row.names = F, sep = '\t')
resImmGen = resImmGen[,1:5]
write.table(resImmGen, file = paste0('bdgdiff/ora/ora_ImmGen_wt_noLeadingEdge.txt'), quote = F, row.names = F, sep = '\t')

resHallmark = fora(pathways = h_genesets, genes = wt, universe = atac, minSize = 1, maxSize = Inf)
resHallmark = resHallmark[resHallmark$padj < 0.05]
resHallmark = apply(resHallmark,2,as.character)
write.table(resHallmark, paste0('bdgdiff/ora/ora_Hallmark_wt.txt'), quote = F, row.names = F)
resHallmark = resHallmark[,1:5]
write.table(resHallmark, file = paste0('bdgdiff/ora/ora_Hallmark_wt_noLeadingEdge.txt'), quote = F, row.names = F)


### make plots

# Hallmark plot
common_hallmark = read.table('bdgdiff/ora/ora_Hallmark_common_noLeadingEdge.txt', header = T)
wt_hallmark = read.table('bdgdiff/ora/ora_Hallmark_wt_noLeadingEdge.txt', header = T)
lin28_hallmark = read.table('bdgdiff/ora/ora_Hallmark_lin28_noLeadingEdge.txt', header = T)

common_hallmark$Type = 'Common'
wt_hallmark$Type = 'Adult'
lin28_hallmark$Type = 'Neo (Lin28)'

all_hallmark = rbind(common_hallmark[1:3,], lin28_hallmark, wt_hallmark[1:3,])

ggplot(all_hallmark, aes(x = Type, y = reorder(pathway,-log10(padj)))) +
  geom_point(aes(size = ifelse(-log10(padj)==0, NA, abs(-log10(padj))), color = Type)) +
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
        #axis.text.x=element_blank(),
        #ggh4x.axis.nestline = element_line(linetype = 1),
        panel.border = element_blank()) +
  guides(size = guide_legend(title='-log10(adj p-val)')) +
  scale_size_continuous(limits = c(2,12)) +
  scale_color_manual(values = c('#88FA4E', '#017100', '#2C403E'), 
                     limits = c('Adult', 'Neo (Lin28)', 'Common')) +
  scale_x_discrete(limits = c('Adult', 'Neo (Lin28)', 'Common'))
ggsave(paste0("bdgdiff/ora/plots/ora_Hallmark_top3.pdf"), 
       plot = last_plot(), device = "pdf", width = 7, height = 4, dpi = 200)



# ImmGen plot

common_ImmGen = read.table('bdgdiff/ora/ora_ImmGen_common_noLeadingEdge.txt', header = T, sep = '\t')
wt_ImmGen = read.table('bdgdiff/ora/ora_ImmGen_wt_noLeadingEdge.txt', header = T, sep = '\t')
lin28_ImmGen = read.table('bdgdiff/ora/ora_ImmGen_lin28_noLeadingEdge.txt', header = T, sep = '\t')

common_ImmGen$Type = 'Common'
wt_ImmGen$Type = 'Adult'
lin28_ImmGen$Type = 'Neo (Lin28)'

all_ImmGen = rbind(common_ImmGen, wt_ImmGen, lin28_ImmGen)
all_ImmGen$Significance = all_ImmGen$Type
all_ImmGen$Significance = ifelse(all_ImmGen$padj > 0.05, 'Not significant', all_ImmGen$Significance)

for (i in 1:nrow(all_ImmGen)){
  if (all_ImmGen[i, 'pathway'] %in% c("Cell cycle and division",
                                    "Short term effector and memory",
                                    'Naive or late effector or memory',
                                    'Short term effector or memory',
                                    'Late effector or memory',
                                    'Initial cytokine or effector response')){
    all_ImmGen[i, 'Group'] = 'Effector\n\ngene sets'
  } else if (all_ImmGen[i, 'pathway'] %in% c('Preparation for cell division',
                                           'Naive and late memory',
                                           'Early effector late memory',
                                           'Memory precursor')){
    all_ImmGen[i, 'Group'] = 'Memory\n\ngene sets'
  }
}
all_ImmGen$pathway = factor(all_ImmGen$pathway, levels = rev(c("Cell cycle and division",
                                                           "Short term effector and memory",
                                                           'Naive or late effector or memory',
                                                           'Short term effector or memory',
                                                           'Late effector or memory',
                                                           'Initial cytokine or effector response',
                                                           'Preparation for cell division',
                                                           'Naive and late memory',
                                                           'Early effector late memory',
                                                           'Memory precursor')))

ggplot(all_ImmGen, aes(x = Type, y = pathway)) +
  geom_point(aes(size = ifelse(-log10(padj)==0, NA, abs(-log10(padj))), color = Significance)) +
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
        strip.placement = "outside", 
        strip.text.y.left = element_text(angle = 0),
        strip.background = element_rect(fill = "gray95", colour = NA),
        panel.border = element_blank()) +
  guides(size = guide_legend(title='-log10(adj p-val)')) +
  facet_grid2(Group ~ ., scales = 'free', space = 'free', switch = "y") +
  scale_size_continuous(limits = c(0.000001,28),
                        breaks = c(5, 15, 25))+
  scale_color_manual(values = c('#88FA4E', '#017100', '#2C403E', 'gray90'), 
                     limits = c('Adult', 'Neo (Lin28)', 'Common', 'Not significant'),
                     name = 'Type') +
  scale_x_discrete(limits = c('Adult', 'Neo (Lin28)', 'Common'))
ggsave(paste0("bdgdiff/ora/plots/ora_ImmGen.pdf"), 
       plot = last_plot(), device = "pdf", width = 6.5, height = 5, dpi = 200)
