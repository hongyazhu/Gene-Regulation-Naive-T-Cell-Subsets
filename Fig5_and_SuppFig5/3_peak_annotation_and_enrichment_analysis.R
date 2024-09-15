# this code performs peak annotation and enrichment analysis of annotated genes (Fig. 5B and Supp Fig. 5A)

# peak annotation:
# cutoffs=(1e20 1e50 1e100 1e150 1e180)
# for cutoff in "${cutoffs[@]}" ; do
# annotatePeaks.pl "../peaks/peaks1e02_${cutoff}.narrowPeak" mm10 -annStats "./peaks1e02_${cutoff}.annoPeak.Stats" > "./peaks1e02_${cutoff}.annoPeak.txt" 
# done

# ATAC-seq peaks were also annotated, used as background
# annotatePeaks.pl ATACseq_peaks.bed mm10 -annStats ATAC.annoPeak.Stats > ATAC.annoPeak.txt

library(fgsea)
library(qusage)

cutoffs = c('1e150') # '1e20', '1e50', '1e100', '1e180'

immgen_genesets = qusage::read.gmt('cluster_list_geneSetNames.gmt') # from Best et al., 2013
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

reactome_genesets = qusage::read.gmt('MSigDB_gene_sets/mouse/m2.cp.reactome.v2023.1.Mm.symbols.gmt')

for (cutoff in cutoffs){
  geneannoPeak = unique(read.delim(paste0('peaks1e02_', cutoff, '.annoPeak.txt'), sep = '\t')$Gene.Name)
  geneannoPeakATAC = unique(read.delim('ATAC.annoPeak.txt', sep = '\t')$Gene.Name)
  
  resImmGen = fora(pathways = immgen_genesets, genes = geneannoPeak, universe = geneannoPeakATAC, minSize = 1, maxSize = Inf)
  resImmGen = apply(resImmGen,2,as.character)
  write.table(resImmGen, file = paste0('ora_ImmGen_', cutoff, '.txt'), quote = F, row.names = F, sep = '\t')
  resImmGen = resImmGen[,1:5]
  write.table(resImmGen, file = paste0('ora_ImmGen_', cutoff, '_noLeadingEdge.txt'), quote = F, row.names = F, sep = '\t')
  
  resReactome = fora(pathways = reactome_genesets, genes = geneannoPeak, universe = geneannoPeakATAC, minSize = 1, maxSize = Inf)
  resReactome = resReactome[resReactome$padj < 0.05]
  resReactome = apply(resReactome,2,as.character)
  write.table(resReactome, paste0('ora_Reactome_', cutoff, '.txt'), quote = F, row.names = F)
  resReactome = resReactome[,1:5]
  write.table(resReactome, file = paste0('ora_Reactome_', cutoff, '_noLeadingEdge.txt'), quote = F, row.names = F)
}


# make ora plots

library(ggplot2)

# reactome
for (cutoff in cutoffs){
  resReactome = read.table(paste0('ora_Reactome_', cutoff, '_noLeadingEdge.txt'), header = T)
  res_plot = resReactome[order(resReactome$padj),][1:10,]
  ggplot(res_plot, aes(x = 1, y = reorder(pathway,-log10(padj)))) +
    geom_point(aes(size = ifelse(-log10(padj)==0, NA, abs(-log10(padj))))) +
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
          axis.text.x=element_blank(),
          #ggh4x.axis.nestline = element_line(linetype = 1),
          panel.border = element_blank()) +
    guides(size = guide_legend(title='-log10(adj p-val)')) 
  ggsave(
    paste0('ora_Reactome_', cutoff, '.pdf'),
    plot = last_plot(),
    device = "pdf",
    width = 7,
    height = 7,
    dpi = 300
  )
}

# ImmGen
library(ggh4x)
strip <- strip_themed(background_y = elem_list_rect(fill = 'gray95'))
for (cutoff in cutoffs){
  resImmGen = read.table(paste0('ora_ImmGen_', cutoff, '_noLeadingEdge.txt'), header = T, sep = '\t')
  res_plot = resImmGen[order(resImmGen$padj),]
  res_plot$Group = NA
  for (i in 1:nrow(res_plot)){
    if (res_plot[i, 'pathway'] %in% c("Cell cycle and division",
                                     "Short term effector and memory",
                                     'Naive or late effector or memory',
                                     'Short term effector or memory',
                                     'Late effector or memory',
                                     'Initial cytokine or effector response')){
      res_plot[i, 'Group'] = 'Effector\n\ngene sets'
    } else if (res_plot[i, 'pathway'] %in% c('Preparation for cell division',
                                            'Naive and late memory',
                                            'Early effector late memory',
                                            'Memory precursor')){
      res_plot[i, 'Group'] = 'Memory\n\ngene sets'
    }
  }
  res_plot$pathway = factor(res_plot$pathway, levels = rev(c("Cell cycle and division",
                                                           "Short term effector and memory",
                                                           'Naive or late effector or memory',
                                                           'Short term effector or memory',
                                                           'Late effector or memory',
                                                           'Initial cytokine or effector response',
                                                           'Preparation for cell division',
                                                           'Naive and late memory',
                                                           'Early effector late memory',
                                                           'Memory precursor')))
  res_plot$Direction = ifelse(res_plot$padj > 0.05, 'Not significant', 'Significant')
  
  ggplot(res_plot, aes(x = 1, y = pathway)) +
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
          axis.text.x=element_blank(),
          strip.placement = "outside", 
          strip.text.y.left = element_text(angle = 0),
          strip.background = element_rect(fill = "gray95", colour = NA),
          panel.border = element_blank()) +
    guides(size = guide_legend(title='-log10(adj p-val)')) +
    facet_grid2(Group ~ ., scales = 'free', space = 'free', switch = "y") +
    scale_colour_manual(limits = c('Significant', 'Not significant'),
                        values = c('black', 'gray80'))
  ggsave(
    paste0('ora_ImmGen_', cutoff, '.pdf'),
    plot = last_plot(),
    device = "pdf",
    width = 7,
    height = 7,
    dpi = 300
  )
}


### remake plot reactome plot
cutoff = '1e150'
resReactome = read.table(paste0('ora_Reactome_', cutoff, '_noLeadingEdge.txt'), header = T)
res_plot = resReactome[order(resReactome$padj),][1:10,]
ggplot(res_plot, aes(x = 1, y = reorder(pathway,-log10(padj)))) +
  geom_point(aes(size = ifelse(-log10(padj)==0, NA, abs(-log10(padj))))) +
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
        axis.text.x=element_blank(),
        #ggh4x.axis.nestline = element_line(linetype = 1),
        panel.border = element_blank()) +
  guides(size = guide_legend(title='-log10(adj p-val)')) +
  scale_size_continuous(limits = c(3,16))
ggsave(
  paste0('ora_Reactome_', cutoff, '_remake.pdf'),
  plot = last_plot(),
  device = "pdf",
  width = 7.5,
  height = 6,
  dpi = 300
)

