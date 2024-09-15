# Validating Eomes targets in TRN with CUT&Tag data (Fig. 5C)
# Fisher's test for the overlap between the predicted network and CUT&Tag derived Eomes targets

cutoffs = c('1e20', '1e50', '1e100', '1e150', '1e180')

result = data.frame(matrix(nrow = length(cutoffs), ncol = 3))
colnames(result) = c('cutoff', 'odds_ratio', 'p_value')

for (i in 1:length(cutoffs)){
  
  cutoff = cutoffs[i]
  
  num_nets = 2
  ms = 10
  bias = 50
  net = read.table(paste0('chip_ko_network_atac_ms', ms, '_bias', bias, '_maxComb_cut01_sharedbyMorethan', num_nets, 'nets_sp.tsv'),
                   header = T)
  target_tmp = net[net$TF == 'Eomes','Target']
  
  target_bg = read.table('genes_DE/deg_padj01.txt')$V1
  nontarget_tmp = target_bg[!target_bg %in% target_tmp]
  
  anno = read.delim(paste0('peaks1e02_', cutoff, '.annoPeak.txt'))
  
  genes_all = anno$Gene.Name
  
  genes_unique = unique(genes_all)
  
  length(target_tmp)
  length(target_bg)
  
  length(intersect(target_tmp, genes_unique))
  length(intersect(target_bg, genes_unique))
  
  mat = matrix(c(length(intersect(target_tmp, genes_unique)), # targets with CUTandTag evidence
                 length(target_tmp) - length(intersect(target_tmp, genes_unique)), # targets without CUTandTag evidence
                 length(intersect(nontarget_tmp, genes_unique)), # non-targets with CUTandTag evidence
                 length(nontarget_tmp) - length(intersect(nontarget_tmp, genes_unique)) # targets without CUTandTag evidence
  ),nrow=2,ncol=2)
  
  res_fisher = fisher.test(mat)
  
  result[i, 1] = cutoff
  result[i, 2] = res_fisher$estimate
  result[i, 3] = res_fisher$p.value
  
  print(cutoff)
  print(res_fisher$p.value)
  print(res_fisher$estimate)
}


library(ggplot2)
ggplot(result, aes(x = odds_ratio, y = cutoff)) +
  geom_point(aes(size = -log10(p_value)))+
  theme_bw() +
  guides(size = guide_legend(title='-log10(p-val)')) +
  scale_size_continuous(limits = c(2,9)) +
  scale_y_discrete(limits = rev(c('1e20', '1e50', '1e100', '1e150', '1e180'))) +
  xlab('Odds ratio') +
  ylab('Cutoffs')+
  expand_limits(x = 1)+ 
  geom_vline(xintercept = 1)
ggsave(paste0("/workdir/hz543/projects/Inferelator/cd8_rerun/analyze_network/CUTandTag/Eomes_Lin28/peak_calling/different_cutoffs_after_MACS3/overlap_with_network/overlap_with_network_peaks1e02_cutoffs_withx1.pdf"), 
       plot = last_plot(), device = "pdf", width = 4, height = 3, dpi = 200)
