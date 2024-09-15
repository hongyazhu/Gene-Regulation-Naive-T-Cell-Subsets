# get count table and PCA (Fig. 5D)

# get count table (shell script)
# cutoffs=(1e150) # 1e20 1e50 1e100 1e180
# for cutoff in "${cutoffs[@]}" ; do
# awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print $4, $1, $2+1, $3, "."}' "../peaks/peaks1e02_${cutoff}.narrowPeak" > "./peaks1e02_${cutoff}.saf"
# featureCounts -T 8 -p -a "./peaks1e02_${cutoff}.saf" -F SAF -o "./peaks1e02_${cutoff}_readCountInPeaks.txt" dedup-BAMS/*"EOMES"*.DEDUP.bam dedup-BAMS/*IgG*.DEDUP.bam &> "./peaks1e02_${cutoff}_readCountInPeaks.log"
# done

# PCA 
library(edgeR)
library(ggfortify)
percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}

cutoffs = c('1e150') # '1e20', '1e50', '1e100', '1e180'

for (cutoff in cutoffs){
  peakCounts = read.table(paste0('peaks1e02_', cutoff, '_readCountInPeaks.txt'), header = T)
  rownames(peakCounts) = peakCounts$Geneid
  peakCounts = peakCounts[,7:ncol(peakCounts)]
  colnames(peakCounts) = gsub('...from_Ya.dedup.BAMS.', '', colnames(peakCounts))
  colnames(peakCounts) = gsub('.DEDUP.bam', '', colnames(peakCounts))
  colnames(peakCounts) = gsub('_mCD8', '', colnames(peakCounts))
  
  peakCpm = cpm(peakCounts)
  pca <- prcomp(t(peakCpm))
  #rownames(pca$x) =colnames(peakCpm)
  autoplot(pca, label = T, label.size = 3.5) 
  
  # remove IgGs
  peakCpm_Eomes = peakCpm[,1:3]
  pca <- prcomp(t(peakCpm_Eomes))
  autoplot(pca, label = T, label.size = 3.5) 
  
  sample_info = data.frame(sample_names = colnames(peakCpm_Eomes), sample_types = c('Lin28Tg', 'WT', 'WT'))
  
  pca_plot = as.data.frame(pca$x)
  pca_plot_merged = merge(pca_plot, sample_info, by.x = 0, by.y = 'sample_names')
  
  ggplot(pca_plot_merged, aes(x = PC1, y= PC2)) +
    geom_point(size = 5, aes(colour = sample_types)) +
    theme_classic(base_size=12, base_family='ArialMT') +
    xlab(paste0("PC1(", percent(summary(pca)$importance[2,][1]), ")")) +
    ylab(paste0("PC2(", percent(summary(pca)$importance[2,][2]), "%)")) +
    labs(color='Group') +
    scale_color_manual(values = c('#88FA4E', '#017100'),
                       limits = c('WT', 'Lin28Tg')) +
    #labs(shape = '') +
    ggtitle('Principle Component Analysis (PC1&2)') +
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))
  #scale_color_manual(values = c("#0000CD", "#228B22", "#FFA500"))
  ggsave(paste0("plots/pca_peaks1e02_", cutoff, ".pdf"), plot = last_plot(), device = "pdf", width = 5, height = 3, dpi = 200)
  
}
