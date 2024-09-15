# Over-representation analysis (Fig 1F)

peakAnno_adult = read.csv('peak_annotation/adult_peak_anno.txt', sep = '\t')
peakAnno_neo   = read.csv('peak_annotation/neo_peak_anno.txt',   sep = '\t')
peakAnno_clean = read.csv('peak_annotation/clean_peak_anno.txt', sep = '\t')
peakAnno_dirty = read.csv('peak_annotation/dirty_peak_anno.txt', sep = '\t')
peakAnno_tn    = read.csv('peak_annotation/tn_peak_anno.txt',    sep = '\t')
peakAnno_vm    = read.csv('peak_annotation/vm_peak_anno.txt',    sep = '\t')

peakAnno_bg    = read.csv('peak_annotation/merged_peak_anno.txt',sep = '\t')

peakAnnoGenes_adult = unique(peakAnno_adult$'Gene.Name')
peakAnnoGenes_neo = unique(peakAnno_neo$'Gene.Name')
peakAnnoGenes_clean = unique(peakAnno_clean$'Gene.Name')
peakAnnoGenes_dirty = unique(peakAnno_dirty$'Gene.Name')
peakAnnoGenes_tn = unique(peakAnno_tn$'Gene.Name')
peakAnnoGenes_vm = unique(peakAnno_vm$'Gene.Name')

peakAnnoGenes_bg = unique(peakAnno_bg$'Gene.Name')

library(clusterProfiler)
immgen_genesets = read.gmt('cluster_list_geneSetNames_forClusterProfiler.gmt') # ImmGen gene sets from Best et al., 2013

ora_adult <- enricher(gene = peakAnnoGenes_adult,
                      universe = peakAnnoGenes_bg,
                      TERM2GENE=immgen_genesets,
                      qvalueCutoff = 1) # pvalueCutoff by default is 0.05, on unadjusted and adjusted pvalues
ora_neo <- enricher(gene = peakAnnoGenes_neo,
                    universe = peakAnnoGenes_bg,
                    TERM2GENE=immgen_genesets,
                    qvalueCutoff = 1)

ora_clean <- enricher(gene = peakAnnoGenes_clean,
                      universe = peakAnnoGenes_bg,
                      TERM2GENE=immgen_genesets,
                      qvalueCutoff = 1)
ora_dirty <- enricher(gene = peakAnnoGenes_dirty,
                      universe = peakAnnoGenes_bg,
                      TERM2GENE=immgen_genesets,
                      qvalueCutoff = 1)

ora_tn <- enricher(gene = peakAnnoGenes_tn,
                   universe = peakAnnoGenes_bg,
                   TERM2GENE=immgen_genesets,
                   qvalueCutoff = 1)
ora_vm <- enricher(gene = peakAnnoGenes_vm,
                   universe = peakAnnoGenes_bg,
                   TERM2GENE=immgen_genesets,
                   qvalueCutoff = 1)

oraRes_adult = ora_adult@result
oraRes_adult$Direction = ifelse(oraRes_adult$p.adjust < 0.05, 'Adult', 'Not significant')
oraRes_adult$type = 'Adult'
oraRes_neo = ora_neo@result
oraRes_neo$Direction = ifelse(oraRes_neo$p.adjust < 0.05, 'Neo', 'Not significant')
oraRes_neo$type = 'Neo'

oraRes_clean = ora_clean@result
oraRes_clean$Direction = ifelse(oraRes_clean$p.adjust < 0.05, 'Clean', 'Not significant')
oraRes_clean$type = 'Clean'
oraRes_dirty = ora_dirty@result
oraRes_dirty$Direction = ifelse(oraRes_dirty$p.adjust < 0.05, 'Dirty', 'Not significant')
oraRes_dirty$type = 'Dirty'

oraRes_tn = ora_tn@result
oraRes_tn$Direction = ifelse(oraRes_tn$p.adjust < 0.05, 'TN', 'Not significant')
oraRes_tn$type = 'TN'
oraRes_vm = ora_vm@result
oraRes_vm$Direction = ifelse(oraRes_vm$p.adjust < 0.05, 'VM', 'Not significant')
oraRes_vm$type = 'VM'

oraRes_combined = rbind(oraRes_adult, oraRes_neo, oraRes_clean, oraRes_dirty, oraRes_tn, oraRes_vm)

oraRes_combined$ID = gsub('cell cycle and division', 'Cell cycle and division', oraRes_combined$ID)
oraRes_combined$ID = gsub(',', '', oraRes_combined$ID)
oraRes_combined$ID = gsub(',', '', oraRes_combined$ID)
oraRes_combined$ID = gsub('Short-term', 'Short term', oraRes_combined$ID)
oraRes_combined$ID = gsub('short-term', 'Short term', oraRes_combined$ID)

oraRes_combined$ID = factor(oraRes_combined$ID, levels = rev(c("Cell cycle and division",
                                                         "Short term effector and memory",
                                                         'Naive or late effector or memory',
                                                         'Short term effector or memory',
                                                         'Late effector or memory',
                                                         'Initial cytokine or effector response',
                                                         'Preparation for cell division',
                                                         'Naive and late memory',
                                                         'Early effector late memory',
                                                         'Memory precursor')))
#oraRes_combined = as.data.frame(oraRes_combined)
oraRes_combined$Group = NA
for (i in 1:nrow(oraRes_combined)){
  if (oraRes_combined[i, 'ID'] %in% c("Cell cycle and division",
                                   "Short term effector and memory",
                                   'Naive or late effector or memory',
                                   'Short term effector or memory',
                                   'Late effector or memory',
                                   'Initial cytokine or effector response')){
    oraRes_combined[i, 'Group'] = 'Effector\n\ngene sets'
  } else if (oraRes_combined[i, 'ID'] %in% c('Preparation for cell division',
                                          'Naive and late memory',
                                          'Early effector late memory',
                                          'Memory precursor')){
    oraRes_combined[i, 'Group'] = 'Memory\n\ngene sets'
  }
}
library(ggh4x)
strip <- strip_themed(background_y = elem_list_rect(fill = 'gray95'))
ggplot(oraRes_combined, aes(x = type, y = ID)) +
  geom_point(aes(size = ifelse(-log10(p.adjust)==0, NA, abs(-log10(p.adjust))), color = Direction)) +
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
  scale_x_discrete(limits = c('VM', 'TN', 'Neo', 'Adult', 'Dirty', 'Clean'),
                   #labels = c('VM, TN','Neo, Adult', 'Dirty, Clean'), 
                   position = 'top') +
  guides(size = guide_legend(title='-log10(adj p-val)')) +
  facet_grid2(Group ~ ., scales = 'free', space = 'free', switch = "y", strip = strip)
ggsave(
  paste0('peak_annotation/gsa.pdf'),
  plot = last_plot(),
  device = "pdf",
  width = 7,
  height = 5,
  dpi = 300
)
