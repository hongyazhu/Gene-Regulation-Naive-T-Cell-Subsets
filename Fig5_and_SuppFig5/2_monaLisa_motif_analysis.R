# monaLisa analysis to identify enriched motifs, using ATAC-seq peaks as background (Fig. 5A)

library(tidyverse)
library(ggplot2)
library(monaLisa)
library(GenomicRanges)
library(TFBSTools)
library(JASPAR2020)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BiocParallel)
library(SummarizedExperiment)

pwms <- getMatrixSet(JASPAR2020,opts = list(matrixtype = "PWM",collection="CORE",tax_group = "vertebrates"))
pwms_TBox <- getMatrixSet(JASPAR2020,opts = list(matrixtype = "PWM",collection="CORE",tax_group = "vertebrates",class="T-Box factors"))

# Single set motif enrichment analysis: comparing a set of sequences to a suitable background (ATAC-seq peaks)
ATACpeaks = read.table('ATACseq_peaks.bed', 
                       col.names = c("chr","start","end"))
ATACpeaks.gr<-makeGRangesFromDataFrame(ATACpeaks,keep.extra.columns=T, seqinfo=Seqinfo(genome="mm10"))


cutoff = '1e150'
file_peaks <- paste0('peaks/peaks1e02_', cutoff, '.narrowPeak')
peaks = read.table(file_peaks, header = FALSE, 
                   col.names = c("chr","start","end","Peak_ID","score","strand","fold_change","log_pvalue","log_qvalue","Summit_distance"))
peaks_summit<-mutate(peaks, start=start+Summit_distance,end=start+1) %>%
  dplyr::select(!Summit_distance)
peaks_summit.gr<-makeGRangesFromDataFrame(peaks_summit,keep.extra.columns=T, seqinfo=Seqinfo(genome="mm10"))
peaks_500bpPeak.gr<-trim(resize(peaks_summit.gr,width = 500, fix = "center"))
peaks_500bpPeak.seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, peaks_500bpPeak.gr)

peaks_singleset_CORE_motifEnr <- calcBinnedMotifEnrR(seqs = peaks_500bpPeak.seq,
                                                     pwmL = pwms,
                                                     background = "genome",
                                                     genome = BSgenome.Mmusculus.UCSC.mm10,
                                                     genome.regions = ATACpeaks.gr, # sample from ATACseq peaks
                                                     #genome.oversample = 2, 
                                                     BPPARAM = BiocParallel::MulticoreParam(32),
                                                     verbose = TRUE)
sel <- apply(assay(peaks_singleset_CORE_motifEnr, "negLog10Padj"), 1, function(x) max(abs(x), 0, na.rm = TRUE)) > 15
sel2 <- apply(assay(peaks_singleset_CORE_motifEnr, "log2enr"), 1, function(x) max(x, -Inf, na.rm = TRUE)) > 0.7
#sel2 <- apply(assay(peaks_singleset_CORE_motifEnr, "log2enr"), 1, function(x) max(abs(x), 0, na.rm = TRUE)) > 0.7
sum(sel & sel2)
peaks_singleset_CORE_motifEnr_sel<-peaks_singleset_CORE_motifEnr[sel & sel2,]

SimMatSel <- motifSimilarity(rowData(peaks_singleset_CORE_motifEnr_sel)$motif.pfm)
hcl <- hclust(as.dist(1 - SimMatSel), method = "average")
pdf(paste0('monalisa_q001_', cutoff, '_singleset_allCoreMotifs_hcl_q15enr7.pdf'), width = 8, height = 9)
plotMotifHeatmaps(x = peaks_singleset_CORE_motifEnr_sel, which.plots = c("log2enr", "negLog10Padj"),  width = 2.0, cluster = hcl, maxEnr = 1, maxSig = 50, show_seqlogo = TRUE, width.seqlogo = 1.2)
dev.off()
