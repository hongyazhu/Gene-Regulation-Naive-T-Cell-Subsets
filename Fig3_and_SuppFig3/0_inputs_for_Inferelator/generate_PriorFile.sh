### Prior is conprised of several components:
# 1. ChIP-seq
# 2. KO RNA-seq
# 3. A published network based on EP-pairs
# 4. ATAC-seq
# Final step is to combine all these priors

### 1. ChIP-seq - Runx3 and Tcf1
# 1.1 ChIP-seq of Runx3 in naive CD8+ T cells (Lotem et al., 2013)
prefetch $(<SraAccList.txt)
fastq-dump --outdir fastq_data --split-3 $(<SraAccList.txt)
fastqc -t 3 -f fastq -o qc fastq_data/*
# bowtie2 alignment
cd fastq_data
for f in `ls -1 *.fastq.gz | sed 's/.fastq.gz//' `
do
bowtie2 -x mouse/mm10_bowtie2/mm10 -U ${f}.fastq.gz -p 4 -S ../bowtie2/${f}.sam
done
cd bowtie2
for i in $(ls *.sam | sed 's/.sam//'); do
samtools view -h -@ 5 -q 1 ${i}.sam | samtools sort -@ 5 -O bam -o ../bam_sort/${i}_filtered_sort.bam
gzip ${i}.sam
picard MarkDuplicates I=../bam_sort/${i}_filtered_sort.bam O=../picard/${i}_filtered_sort_rmDup.bam M=../picard/${i}_filtered_sort.dups.txt REMOVE_DUPLICATES=true 2> ../picard/${i}_filtered_sort.dups.log
done

cd picard
samtools merge runx3_merged.bam SRR955620_SRR955621_filtered_sort_rmDup.bam SRR955622_filtered_sort_rmDup.bam
samtools merge ctrl_merged.bam SRR955623_filtered_sort_rmDup.bam SRR955624_filtered_sort_rmDup.bam

macs2 callpeak -t runx3_merged.bam -c ctrl_merged.bam -f BAM -g mm -n runx3 -p 0.01 --outdir ../macs2_p01 2> ../macs2_p01/runx3.log
cd macs2_p01
cut -f 1-6 runx3_peaks.narrowPeak > runx3_peaks.bed


# 1.2 ChIP-seq of Tcf1 in naive CD8+ T cells (Steinke et al., 2014)
prefetch $(<SraAccList.txt)
fastq-dump --outdir fastq_data --split-3 $(<SraAccList.txt)
fastqc -t 3 -f fastq -o qc fastq_data/* # quality is high, no adaptors, no need for trimming
# bowtie2 alignment
cd fastq_data
for f in `ls -1 *.fastq.gz | sed 's/.fastq.gz//' `
do
bowtie2 -x /home/hz543/data/mouse/mm10_bowtie2/mm10 -U ${f}.fastq.gz -p 4 -S ../bowtie2/${f}.sam
done
cd bowtie2
for i in $(ls *.sam | sed 's/.sam//'); do
samtools view -h -@ 5 -q 1 ${i}.sam | samtools sort -@ 5 -O bam -o ../bam_sort/${i}_filtered_sort.bam
gzip ${i}.sam
picard MarkDuplicates I=../bam_sort/${i}_filtered_sort.bam O=../picard/${i}_filtered_sort_rmDup.bam M=../picard/${i}_filtered_sort.dups.txt REMOVE_DUPLICATES=true 2> ../picard/${i}_filtered_sort.dups.log
done

cd picard
macs2 callpeak -t SRR1024054_filtered_sort_rmDup.bam -c SRR1024055_filtered_sort_rmDup.bam -f BAM -g mm -n tcf1 -p 0.01 --outdir ../macs2_p01 2> ../macs2_p01/tcf1.log
cd macs2_p01
cut -f 1-6 tcf1_peaks.narrowPeak > tcf1_peaks.bed

# Priors for each of the ChIP-seq data generated using construct_atac_prior.R from https://github.com/emiraldi/infTRN_lassoStARS

# combine the priors form two ChIP-seq profiles
cat runx3_chip_q_sp.tsvtcf1_chip_q_sp.tsv | awk '!x[$0]++' > prior/chip_nofli1_sp.tsv


### 2. KO - Tcf1 and Bach2

# 2.1 Tcf1 KO (He et al., 2016)
prefetch $(<SraAccList.txt)
fastq-dump --outdir fastq_data --split-3 $(<SraAccList.txt)
fastqc -t 4 -f fastq -o qc fastq_data/* # quality is high (above 20), no adaptors, no need for trimming
# alignment
cd fastq_data
for f in `ls -1 *.fastq.gz | sed 's/.fastq.gz//' `
do
hisat2 -x mouse/mm10/genome -U ${f}.fastq.gz -S ../hisat2/${f}.bam -p 8 2> ../hisat2/${f}_summary.txt
done
featureCounts -a mouse/gencode_annotation/gencode.vM21.annotation.names.saf -F 'SAF' -s 0 -Q 50 -o counts_name.txt hisat2/*.bam 2> featureCounts_name.log

# Differential expression (in R)
library(DESeq2)
library(ggplot2)
counts = read.table("counts_name.txt", header = T, row.names = 1)
counts = counts[,6:ncol(counts)]
colnames(counts) = c("ctrl_rep1", "ctrl_rep2", "ko_rep1", "ko_rep2")

coldata = data.frame(c('ctrl', 'ctrl', 'ko', 'ko'))
colnames(coldata) = c('type')
rownames(coldata) = colnames(colnames)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design= ~ type)
keep <- rowSums(counts(dds) > 0) >= 1
dds <- dds[keep,]

dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, name="type_ko_vs_ctrl")
resSig <- subset(res, padj < 0.05)

sp = data.frame(Target = rownames(resSig))
sp$TF = 'Tcf7'
sp$Weight = resSig$log2FoldChange
sp <- sp[c("TF", "Target", "Weight")]
write.table(sp, "tcf1_sp.tsv", quote = F, sep = "\t", row.names = F)

# 2.2 Bach2 KO (Roychoudhuri et al., 2016)
for f in `ls *_[12].fastq.gz | sed 's/_[12].fastq.gz//' | uniq `
do
hisat2 -x mouse/mm10/genome -1 ${f}_1.fastq.gz -2 ${f}_2.fastq.gz -S ../hisat2/${f}.bam -p 8 2> ../hisat2/${f}_summary.txt
done
featureCounts -p -a mouse/gencode_annotation/gencode.vM21.annotation.names.saf -F 'SAF' -s 0 -Q 50 -o counts_name.txt hisat2/*.bam 2> featureCounts_name.log

# Differential expression (in R)
library(DESeq2)
library(ggplot2)
counts = read.table("counts_name.txt", header = T, row.names = 1)
counts = counts[,6:ncol(counts)]
colnames(counts) = c("ctrl_rep1", "ctrl_rep2", "ko_rep1", "ko_rep2")

coldata = data.frame(c('ctrl', 'ctrl', 'ko', 'ko'))
colnames(coldata) = c('type')
rownames(coldata) = colnames(colnames)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design= ~ type)
keep <- rowSums(counts(dds) > 0) >= 1
dds <- dds[keep,]

dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="type_ko_vs_ctrl")
resSig <- subset(res, padj < 0.05)

sp = data.frame(Target = rownames(resSig))
sp$TF = 'Bach2'
sp$Weight = resSig$log2FoldChange
sp <- sp[c("TF", "Target", "Weight")]
write.table(sp, "bach2_sp.tsv", quote = F, sep = "\t", row.names = F)

# combine two priors from TF KO RNA-seq
cat tcf1_sp.tsv bach2_sp.tsv | awk '!x[$0]++' > ko_sp.tsv


### 3. Tn network (He et al., 2016, constructed networks in naive, effector and memory CD8+ T cells by EP pairs) (in R)
library(readxl)
Tn_network = read_excel("TN_network/1-s2.0-S1074761316304812-mmc8.xlsx", sheet = "Tn_TF targets")
GS = Tn_network[, c('TF', 'Target gene')]
GS$Weight = 1
colnames(GS) = c("TF", "Target", "Weight")
# TF is all uppercase. make it to only the first letter is uppercase
firstup <- function(x) { # https://stackoverflow.com/a/18509816
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
GS$TF = firstup(tolower(GS$TF))
write.table(GS, "network_sp.tsv", quote = F, sep = "\t", row.names = F)


### 4. ATAC-seq prior
# ATAC-seq peaks generated as in Fig2_and_SuppFig2/1_ATACseq_pipeline_description and combined into one bed file.
# ATAC-seq prior is constructed with construct_atac_prior.R from https://github.com/emiraldi/infTRN_lassoStARS


### Final step: combine ChIP-seq, KO RNA-seq, Tn network, ATAC-seq priors (in R)
chip_sp = read.table("chip_nofli1_sp.tsv", header = T)
atac_sp = read.table("prior_atac_b_sp.tsv", header = T)
ko_sp = read.table("ko_sp.tsv", header = T)
network_sp = read.table("network_sp.tsv", header = T)

chip_atac_sp = rbind(chip_sp, atac_sp)
chip_ko_atac_sp = rbind(ko_sp, chip_atac_sp)
chip_ko_network_atac_sp = rbind(network_sp, chip_ko_atac_sp)
chip_ko_network_atac = data.frame(matrix(0, ncol = length(unique(chip_ko_network_atac_sp$TF)), nrow = length(unique(chip_ko_network_atac_sp$Target))))
colnames(chip_ko_network_atac) = unique(chip_ko_network_atac_sp$TF)
rownames(chip_ko_network_atac) = unique(chip_ko_network_atac_sp$Target)

for (i in 1:nrow(chip_ko_atac_sp)){
    chip_ko_network_atac[as.character(chip_ko_network_atac_sp[i, 'Target']), as.character(chip_ko_network_atac_sp[i, 'TF'])] = 1
}
write.table(chip_ko_network_atac, "prior/priors_more_than_atac_nofli1/chip_ko_network_atac.tsv", quote = F, sep = "\t", col.names=NA)
