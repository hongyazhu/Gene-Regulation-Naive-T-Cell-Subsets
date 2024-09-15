# This code performed the following anaysis:
# 1. annotate loop ends to genes to get target genes for Micro-C derived regulations
# 2. scan motifs in each loop end to get regulators for Micro-C derived regulations 
# (linking regulators to target genes for each loop provides Micro-C derived regulations)
# 3. shuffle the TRN 10000 times
# 4. find overlaps between Micro-C derived regulations and shuffled networks
# 5. find overlap between Micro-C derived regulations and the predicted network
# 6. make plots (Fig. 3I)

### link loop ends to genes (shell)
bedmap --echo --echo-map-id-uniq --range 10000 merged_loops_open_sort.bed gencode.vM21.annotation.names.bed > merged_loops_open_geneAnnoatation_10k.bed

### motif scan adapted from priorConstruction/construct_atac_prior.R in https://github.com/emiraldi/infTRN_lassoStARS
fdr0001_geneAnnotation_10k = read.table('merged_loops_open_geneAnnoatation_10k_clean.bed')

fdr0001_loops = read.table('merged_loops_open.bedpe')
fdr0001_loops$position1 = paste0('chr', fdr0001_loops$V1, ':', 
                                 fdr0001_loops$V2, '-', 
                                 fdr0001_loops$V3)
fdr0001_loops$position2 = paste0('chr', fdr0001_loops$V4, ':', 
                                 fdr0001_loops$V5, '-', 
                                 fdr0001_loops$V6)

fdr0001_loops_geneAnnotation1 = merge(fdr0001_loops, fdr0001_geneAnnotation_10k, by.x = 'position1', by.y = 'V5')
fdr0001_loops_geneAnnotation2 = merge(fdr0001_loops, fdr0001_geneAnnotation_10k, by.x = 'position2', by.y = 'V5')

#load motifs
library(motifmatchr)
library(GenomicRanges)
library(TFBSTools)
source('infTRN_lassoStARS-master/priorConstruction/utils_prior.R') #!!! address modified
source('infTRN_lassoStARS-master/priorConstruction/utils_bedtools.R') #!!! address modified
dir_motif <- './motif_pfm'
motifs <- list()
for (ix in list.files(dir_motif,full.names=TRUE)){
  curr_motif <- as.matrix(read.delim(ix, header=FALSE, sep='\t'))
  curr_name <- tools::file_path_sans_ext(basename(ix))
  rownames(curr_motif) <- c('A','C','G','T')
  motifs[[curr_name]] <- PFMatrix(ID=curr_name, name=curr_name, bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                                  profileMatrix=100*curr_motif)
}
motifs <- do.call(PFMatrixList, motifs)
# load motif to gene mapping
file_motif_info <- './tbl_motif_2_geneName_cisbp.tsv'
motif_info <- read.delim(file_motif_info, header=FALSE, sep='\t')
colnames(motif_info) <- c('Motif','TF')

fdr0001_loops_geneAnnotation1$position1TFmotifs = NA
fdr0001_loops_geneAnnotation1$position2TFmotifs = NA
for (i in 1:nrow(fdr0001_loops_geneAnnotation1)){
  peaks_range_position1 <- GRanges(seqnames=paste0('chr', fdr0001_loops_geneAnnotation1[i,'V1.x']), 
                                   ranges=IRanges(start=fdr0001_loops_geneAnnotation1[i,'V2.x'], end=fdr0001_loops_geneAnnotation1[i,'V3.x']))
  
  motif_scan_position1 <- matchMotifs(motifs, peaks_range_position1, genome='mm10', p.cutoff=1E-5, out='positions', bg = c("A" = 0.252,"C" = 0.249, "G" = 0.249, "T" = 0.250))
  motif_scan_position1 <- as.data.frame(motif_scan_position1)[,c('seqnames', 'start', 'end', 'group_name')]
  colnames(motif_scan_position1) <- c('Chr','Start','End','Motif')
  motif_scan_position1_info = merge(motif_scan_position1, motif_info, by='Motif')[,2:5]
  
  fdr0001_loops_geneAnnotation1[i,'position1TFmotifs'] = paste(unique(motif_scan_position1_info$TF), collapse = ';')
  
  peaks_range_position2 <- GRanges(seqnames=paste0('chr', fdr0001_loops_geneAnnotation1[i,'V4.x']), 
                                   ranges=IRanges(start=fdr0001_loops_geneAnnotation1[i,'V5'], end=fdr0001_loops_geneAnnotation1[i,'V6']))
  
  motif_scan_position2 <- matchMotifs(motifs, peaks_range_position2, genome='mm10', p.cutoff=1E-5, out='positions', bg = c("A" = 0.252,"C" = 0.249, "G" = 0.249, "T" = 0.250))
  motif_scan_position2 <- as.data.frame(motif_scan_position2)[,c('seqnames', 'start', 'end', 'group_name')]
  colnames(motif_scan_position2) <- c('Chr','Start','End','Motif')
  motif_scan_position2_info = merge(motif_scan_position2, motif_info, by='Motif')[,2:5]
  
  fdr0001_loops_geneAnnotation1[i,'position2TFmotifs'] = paste(unique(motif_scan_position2_info$TF), collapse = ';')
  
}
write.table(fdr0001_loops_geneAnnotation1, 'merged_loops_open_geneAnnotation1_motifscan.bed', 
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep='\t')

fdr0001_loops_geneAnnotation2$position1TFmotifs = NA
fdr0001_loops_geneAnnotation2$position2TFmotifs = NA
for (i in 1:nrow(fdr0001_loops_geneAnnotation2)){
  peaks_range_position1 <- GRanges(seqnames=paste0('chr', fdr0001_loops_geneAnnotation2[i,'V1.x']), 
                                   ranges=IRanges(start=fdr0001_loops_geneAnnotation2[i,'V2.x'], end=fdr0001_loops_geneAnnotation2[i,'V3.x']))
  
  motif_scan_position1 <- matchMotifs(motifs, peaks_range_position1, genome='mm10', p.cutoff=1E-5, out='positions', bg = c("A" = 0.252,"C" = 0.249, "G" = 0.249, "T" = 0.250))
  motif_scan_position1 <- as.data.frame(motif_scan_position1)[,c('seqnames', 'start', 'end', 'group_name')]
  colnames(motif_scan_position1) <- c('Chr','Start','End','Motif')
  motif_scan_position1_info = merge(motif_scan_position1, motif_info, by='Motif')[,2:5]
  
  fdr0001_loops_geneAnnotation2[i,'position1TFmotifs'] = paste(unique(motif_scan_position1_info$TF), collapse = ';')
  
  peaks_range_position2 <- GRanges(seqnames=paste0('chr', fdr0001_loops_geneAnnotation2[i,'V4.x']), 
                                   ranges=IRanges(start=fdr0001_loops_geneAnnotation2[i,'V5'], end=fdr0001_loops_geneAnnotation2[i,'V6']))
  
  motif_scan_position2 <- matchMotifs(motifs, peaks_range_position2, genome='mm10', p.cutoff=1E-5, out='positions', bg = c("A" = 0.252,"C" = 0.249, "G" = 0.249, "T" = 0.250))
  motif_scan_position2 <- as.data.frame(motif_scan_position2)[,c('seqnames', 'start', 'end', 'group_name')]
  colnames(motif_scan_position2) <- c('Chr','Start','End','Motif')
  motif_scan_position2_info = merge(motif_scan_position2, motif_info, by='Motif')[,2:5]
  
  fdr0001_loops_geneAnnotation2[i,'position2TFmotifs'] = paste(unique(motif_scan_position2_info$TF), collapse = ';')
}
write.table(fdr0001_loops_geneAnnotation2, 'merged_loops_open_geneAnnotation2_motifscan.bed', 
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep='\t')

### permutation test
##### get random numbers for generating shuffled networks #####
# get network
num_nets = 2
ms = 10
bias = 50
net = read.table(paste0('chip_ko_network_atac_ms', ms, '_bias', bias, '_maxComb_cut01_sharedbyMorethan', num_nets, 'nets_sp.tsv'),
                 header = T)

require(permute)
df_shuffle = data.frame(shuffleSet(n = nrow(net), nset = 10000))
write.table(df_shuffle, 'df_shuffle_10000.txt',
            sep = '\t', quote = F, row.names = F, col.names = F)

##### get no. regulations overlapped with shuffled networks #####
num_nets = 2
ms = 10
bias = 50
net = read.table(paste0('chip_ko_network_atac_ms', ms, '_bias', bias, '_maxComb_cut01_sharedbyMorethan', num_nets, 'nets_sp.tsv'),
                 header = T)

df_shuffle = read.table('df_shuffle_10000.txt')

fdr0001_loops_geneAnnotation1 = read.table('merged_loops_open_geneAnnotation1_motifscan.bed',
                                           header = T, sep = '\t')
fdr0001_loops_geneAnnotation2 = read.table('merged_loops_open_geneAnnotation2_motifscan.bed',
                                           header = T, sep = '\t')

net_tmp_fit_HiC_nrow_both = c()
for (j in 1:10000){
  if (j %% 1000 == 0)
    print(j)
  net_tmp = net
  net_target_shuffle_tmp = net$Target[unlist(df_shuffle[j, ], use.names = FALSE)]
  net_tmp$Target = net_target_shuffle_tmp
  # check overlap with atac (expected)
  net_tmp$Regulation = paste0(net_tmp$TF, '->',  net_tmp$Target)
  #net_tmp$If_prior_support <- ifelse(net_tmp$Regulation %in% atac_sp_tmp_tftarget$Regulation,"1-Supported by ATAC-seq","2-Not supported by ATAC-seq")
  #table(net_tmp$If_prior_support)
  # check overlap with HiC
  regulations_fit_HiC_both = c()
  for (i in 1:nrow(fdr0001_loops_geneAnnotation1)){
    tflist_tmp_position1 = strsplit(fdr0001_loops_geneAnnotation1$position1TFmotifs[i],";")[[1]]
    tflist_tmp_position2 = strsplit(fdr0001_loops_geneAnnotation1$position2TFmotifs[i],";")[[1]]
    tflist_tmp = unique(c(tflist_tmp_position1, tflist_tmp_position2))
    tflist_tmp_netRegulators = tflist_tmp[tflist_tmp %in% unique(net_tmp$TF)]
    
    regulations_tmp = as.vector(outer(tflist_tmp_netRegulators, unlist(strsplit(fdr0001_loops_geneAnnotation1[i,'V4.y'],";",fixed=T)), paste, sep="->")) 
    regulations_fit_tmp = regulations_tmp[regulations_tmp %in% net_tmp$Regulation]
    regulations_fit_HiC_both = c(regulations_fit_HiC_both,regulations_fit_tmp)
  }
  for (i in 1:nrow(fdr0001_loops_geneAnnotation2)){
    tflist_tmp_position1 = strsplit(fdr0001_loops_geneAnnotation2$position1TFmotifs[i],";")[[1]]
    tflist_tmp_position2 = strsplit(fdr0001_loops_geneAnnotation2$position2TFmotifs[i],";")[[1]]
    tflist_tmp = unique(c(tflist_tmp_position1, tflist_tmp_position2))
    tflist_tmp_netRegulators = tflist_tmp[tflist_tmp %in% unique(net_tmp$TF)]
    
    regulations_tmp = as.vector(outer(tflist_tmp_netRegulators, unlist(strsplit(fdr0001_loops_geneAnnotation1[i,'V4.y'],";",fixed=T)), paste, sep="->")) 
    #regulations_tmp = paste0(tflist_tmp_netRegulators, '->', fdr0001_loops_geneAnnotation1[i,'V4.y'])
    #regulations_tmp = paste0(tflist_tmp_netRegulators, '->', unlist(strsplit(fdr0001_loops_geneAnnotation1[i,'V4.y'],";",fixed=T)))
    regulations_fit_tmp = regulations_tmp[regulations_tmp %in% net_tmp$Regulation]
    regulations_fit_HiC_both = c(regulations_fit_HiC_both,regulations_fit_tmp)
  }
  net_tmp_fit_HiC_nrow_both = c(net_tmp_fit_HiC_nrow_both, length(unique(regulations_fit_HiC_both)))
}
write.table(net_tmp_fit_HiC_nrow_both, 'perm_n10000_both.txt',
            sep = '\t', quote = F, row.names = F, col.names = F)

##### get no. regulations overlapped with real network #####
num_nets = 2
ms = 10
bias = 50
net = read.table(paste0('chip_ko_network_atac_ms', ms, '_bias', bias, '_maxComb_cut01_sharedbyMorethan', num_nets, 'nets_sp.tsv'),
                 header = T)
net_tmp = net

fdr0001_loops_geneAnnotation1 = read.table('merged_loops_open_geneAnnotation1_motifscan.bed',
                                           header = T, sep = '\t')
fdr0001_loops_geneAnnotation2 = read.table('merged_loops_open_geneAnnotation2_motifscan.bed',
                                           header = T, sep = '\t')


net_tmp$Regulation = paste0(net_tmp$TF, '->',  net_tmp$Target)
regulations_fit_HiC_both = c()
for (i in 1:nrow(fdr0001_loops_geneAnnotation1)){
  tflist_tmp_position1 = strsplit(fdr0001_loops_geneAnnotation1$position1TFmotifs[i],";")[[1]]
  tflist_tmp_position2 = strsplit(fdr0001_loops_geneAnnotation1$position2TFmotifs[i],";")[[1]]
  tflist_tmp = unique(c(tflist_tmp_position1, tflist_tmp_position2))
  tflist_tmp_netRegulators = tflist_tmp[tflist_tmp %in% unique(net_tmp$TF)]
  
  regulations_tmp = as.vector(outer(tflist_tmp_netRegulators, unlist(strsplit(fdr0001_loops_geneAnnotation1[i,'V4.y'],";",fixed=T)), paste, sep="->")) 
  #regulations_tmp = paste0(tflist_tmp_netRegulators, '->', fdr0001_loops_geneAnnotation1[i,'V4.y'])
  #regulations_tmp = paste0(tflist_tmp_netRegulators, '->', unlist(strsplit(fdr0001_loops_geneAnnotation1[i,'V4.y'],";",fixed=T)))
  regulations_fit_tmp = regulations_tmp[regulations_tmp %in% net_tmp$Regulation]
  regulations_fit_HiC_both = c(regulations_fit_HiC_both,regulations_fit_tmp)
}
for (i in 1:nrow(fdr0001_loops_geneAnnotation2)){
  tflist_tmp_position1 = strsplit(fdr0001_loops_geneAnnotation2$position1TFmotifs[i],";")[[1]]
  tflist_tmp_position2 = strsplit(fdr0001_loops_geneAnnotation2$position2TFmotifs[i],";")[[1]]
  tflist_tmp = unique(c(tflist_tmp_position1, tflist_tmp_position2))
  tflist_tmp_netRegulators = tflist_tmp[tflist_tmp %in% unique(net_tmp$TF)]
  
  regulations_tmp = as.vector(outer(tflist_tmp_netRegulators, unlist(strsplit(fdr0001_loops_geneAnnotation1[i,'V4.y'],";",fixed=T)), paste, sep="->")) 
  #regulations_tmp = paste0(tflist_tmp_netRegulators, '->', fdr0001_loops_geneAnnotation1[i,'V4.y'])
  #regulations_tmp = paste0(tflist_tmp_netRegulators, '->', unlist(strsplit(fdr0001_loops_geneAnnotation1[i,'V4.y'],";",fixed=T)))
  regulations_fit_tmp = regulations_tmp[regulations_tmp %in% net_tmp$Regulation]
  regulations_fit_HiC_both = c(regulations_fit_HiC_both,regulations_fit_tmp)
}
length(unique(regulations_fit_HiC_both)) 
print('length(unique(regulations_fit_HiC_both))')
print(length(unique(regulations_fit_HiC_both)))

##### make plots #####
perm_n10000 = read.table('perm_n10000_both.txt')$V1
perm_n10000 = as.data.frame(perm_n10000)
library(ggplot2)
net_reg_both = 4087 #length(unique(regulations_fit_HiC_both))
ggplot(perm_n10000, aes(x=perm_n10000)) + 
  geom_histogram(binwidth=5, aes(y = after_stat(density)))+ 
  geom_vline(aes(xintercept=mean(perm_n10000)), color="red", linewidth=1)+
  geom_vline(aes(xintercept=mean(perm_n10000)+sd(perm_n10000)), color="red", linewidth=1)+
  geom_vline(aes(xintercept=mean(perm_n10000)-sd(perm_n10000)), color="red", linewidth=1)+
  geom_vline(aes(xintercept=net_reg_both), color="blue", linewidth=1)+
  theme_bw()+
  ggtitle('Permutation: n = 10000 (union of proximal and distal)')+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab('No. supported regulations') +
  ylab('Density') +
  stat_function( fun = dnorm, args = list(mean = mean(perm_n10000$perm_n10000), sd = sd(perm_n10000$perm_n10000)), lwd = 0.5, col = 'red') +
  geom_text(aes(label=paste0('Mean: ', mean(perm_n10000)),y=0.005,x=3500),vjust=-1,col='red',size=5) +
  geom_text(aes(label=paste0('Standard deviation: ', format(round(sd(perm_n10000), 2), nsmall = 2)),y=0.0025,x=3500),vjust=-1,col='red',size=5) +
  geom_text(aes(label=paste0('Network: ', net_reg_both),y=0.005,x=3800),vjust=-1,col='blue',size=5)+
  geom_text(aes(label=paste0(format(round((net_reg_both-mean(perm_n10000))/sd(perm_n10000), 2), nsmall = 2), ' SD from the mean'),y=0,x=3800),vjust=-1,col='black',size=5)
ggsave(
  paste0('perm_n10000_both_plot.pdf'),
  plot = last_plot(),
  device = "pdf",
  width = 8,
  height = 3,
  dpi = 300
)
