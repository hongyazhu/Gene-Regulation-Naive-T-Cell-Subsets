# generate list of target genes for Inferelator, which are genes differentially expressed with adjusted p-value < 0.01 in at least one comparison (TN vs. VM, adult vs. neo,  clean vs. dirty ...).

library(DESeq2)
padj_cutoff = 0.01


# adult vs. neo

counts_cell2018 = read.table("Smith2018_RNAseq/counts_name.txt", header = T, row.names = 1)
counts_cell2018 = counts_cell2018[,6:ncol(counts_cell2018)]
colnames(counts_cell2018) = c("all_neo_rep1", "all_neo_rep2", "all_adult_rep1", "all_adult_rep2",
                              "vm_neo_rep1", "vm_neo_rep2", "vm_adult_rep1", "vm_adult_rep2",
                              "tn_adult_rep1", "tn_adult_rep2", "vm_neo_5dpi_rep1", "vm_neo_5dpi_rep2",
                              "vm_adult_5dpi_rep1", "vm_adult_5dpi_rep2")
counts_cell2018_aduNeo = counts_cell2018[, c(# 'all_neo_rep1', 'all_neo_rep2', "all_adult_rep1", "all_adult_rep2",
                              "vm_neo_rep1", "vm_neo_rep2", "vm_adult_rep1", "vm_adult_rep2")]

coldata_cell2018_aduNeo = data.frame(c(#'neo', 'neo', 'adult', 'adult', 
                                       'neo', 'neo', 'adult', 'adult')#,
                                     #c('all', 'all', 'all', 'all', 
                                    #   'vm', 'vm', 'vm', 'vm')
                                     )
colnames(coldata_cell2018_aduNeo) = c('age') # , 'type'
rownames(coldata_cell2018_aduNeo) = colnames(counts_cell2018_aduNeo)

dds <- DESeqDataSetFromMatrix(countData = counts_cell2018_aduNeo,
                              colData = coldata_cell2018_aduNeo,
                              design= ~ age) #  + type
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="age_neo_vs_adult")

resSig <- subset(res, padj < padj_cutoff)
deg_cell2018_aduNeo = rownames(resSig)



counts_mir29 = read.table("YeeMon2021_RNAseq/counts_name.txt", header = T, row.names = 1)
counts_mir29 = counts_mir29[,6:ncol(counts_mir29)]
colnames(counts_mir29) = c("NeoM29-3-VM", 
                     "Adult-1-TN", "Adult-1-VM", "NeoNC-1-TN", "NeoNC-1-VM", "NeoM29-1-TN", "NeoM29-1-VM", 
                     "Adult-2-TN", "Adult-2-VM", "NeoNC-2-TN", "NeoNC-2-VM", "NeoM29-2-TN", "NeoM29-2-VM",
                     "Adult-3-TN", "Adult-3-VM", "NeoNC-3-TN", "NeoNC-3-VM", "NeoM29-3-TN")
counts_mir29 = counts_mir29[,c(c(2:18), 1)] # reorder columns
counts_mir29_noM29 = counts_mir29[, !grepl("M29", colnames(counts_mir29))]


coldata_mir29_noM29 = data.frame(c('adult', 'adult', 'neo', 'neo', 'adult', 'adult', 'neo', 'neo', 'adult', 'adult', 'neo', 'neo'),
                                 rep(c('tn', 'vm'), 6))
colnames(coldata_mir29_noM29) = c('age', 'type')
rownames(coldata_mir29_noM29) = colnames(counts_mir29_noM29)

dds <- DESeqDataSetFromMatrix(countData = counts_mir29_noM29,
                              colData = coldata_mir29_noM29,
                              design= ~ age + type)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="age_neo_vs_adult")

resSig <- subset(res, padj < padj_cutoff)
deg_mir29_aduNeo = rownames(resSig)




counts_dirty = read.table("Tabilas2022_RNAseq/counts_name.txt", header = T, row.names = 1)
counts_dirty = counts_dirty[,6:ncol(counts_dirty)]
coln_tmp = sub('.*CD', 'CD', colnames(counts_dirty))
coln_tmp = sub('_CKDL.*', '', coln_tmp)
targetFile = read.table('Tabilas2022_RNAseq/targetFile.txt', header = T)
targetFile$cd = sub('.ReadsPerGene.out.tab.rawCounts', '', targetFile$files)
colnames(counts_dirty) = targetFile[match(coln_tmp, targetFile$cd),]$label
counts_dirty_nobulk = counts_dirty[, !grepl("B.", colnames(counts_dirty))]

coldata_dirty_nobulk = data.frame(c("Neo", "Neo", "Neo", "Neo", "Adult", "Adult", "Adult", "Adult", "Neo", "Neo", "Neo", "Neo", "Adult", "Adult", "Adult", "Adult", "Neo", "Neo", "Neo", "Neo", "Adult", "Adult", "Adult", "Adult"),
                                  tolower(gsub("VM", "MP", substr(colnames(counts_dirty_nobulk), 1, 2))), 
                                  gsub("\\..*","", gsub("^.*?\\.","", colnames(counts_dirty_nobulk))))
colnames(coldata_dirty_nobulk) = c('age', 'type', 'clean')
rownames(coldata_dirty_nobulk) = colnames(counts_dirty_nobulk)


dds <- DESeqDataSetFromMatrix(countData = counts_dirty_nobulk,
                              colData = coldata_dirty_nobulk,
                              design= ~ clean + type + age)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="age_Neo_vs_Adult")

resSig <- subset(res, padj < padj_cutoff)
deg_dirty_aduNeo = rownames(resSig)





# TN vs VM

counts_cell2018 = read.table("Smith2018_RNAseq/counts_name.txt", header = T, row.names = 1)
counts_cell2018 = counts_cell2018[,6:ncol(counts_cell2018)]
colnames(counts_cell2018) = c("all_neo_rep1", "all_neo_rep2", "all_adult_rep1", "all_adult_rep2",
                              "vm_neo_rep1", "vm_neo_rep2", "vm_adult_rep1", "vm_adult_rep2",
                              "tn_adult_rep1", "tn_adult_rep2", "vm_neo_5dpi_rep1", "vm_neo_5dpi_rep2",
                              "vm_adult_5dpi_rep1", "vm_adult_5dpi_rep2")
counts_cell2018_tnVm = counts_cell2018[, c("vm_adult_rep1", "vm_adult_rep2", "tn_adult_rep1", "tn_adult_rep2")]
coldata_cell2018_tnVm = data.frame(c('vm', 'vm', 'tn', 'tn'))
colnames(coldata_cell2018_tnVm) = c('type')
rownames(coldata_cell2018_tnVm) = colnames(counts_cell2018_tnVm)

dds <- DESeqDataSetFromMatrix(countData = counts_cell2018_tnVm,
                              colData = coldata_cell2018_tnVm,
                              design= ~ type)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="type_vm_vs_tn")

resSig <- subset(res, padj < padj_cutoff)
deg_cell2018_tnVm = rownames(resSig)




counts_mir29 = read.table("YeeMon2021_RNAseq/counts_name.txt", header = T, row.names = 1)
counts_mir29 = counts_mir29[,6:ncol(counts_mir29)]
colnames(counts_mir29) = c("NeoM29-3-VM", 
                     "Adult-1-TN", "Adult-1-VM", "NeoNC-1-TN", "NeoNC-1-VM", "NeoM29-1-TN", "NeoM29-1-VM", 
                     "Adult-2-TN", "Adult-2-VM", "NeoNC-2-TN", "NeoNC-2-VM", "NeoM29-2-TN", "NeoM29-2-VM",
                     "Adult-3-TN", "Adult-3-VM", "NeoNC-3-TN", "NeoNC-3-VM", "NeoM29-3-TN")
counts_mir29 = counts_mir29[,c(c(2:18), 1)] # reorder columns
counts_mir29_noM29 = counts_mir29[, !grepl("M29", colnames(counts_mir29))]

coldata_mir29_noM29 = data.frame(c('adult', 'adult', 'neo', 'neo', 'adult', 'adult', 'neo', 'neo', 'adult', 'adult', 'neo', 'neo'),
                                 rep(c('tn', 'vm'), 6))
colnames(coldata_mir29_noM29) = c('age', 'type')
rownames(coldata_mir29_noM29) = colnames(counts_mir29_noM29)

dds <- DESeqDataSetFromMatrix(countData = counts_mir29_noM29,
                              colData = coldata_mir29_noM29,
                              design= ~ age + type)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="type_vm_vs_tn")

resSig <- subset(res, padj < padj_cutoff)
deg_mir29_tnVM = rownames(resSig)




counts_lin28 = read.table("Lin28_RNAseq/counts_name.txt", header = T, row.names = 1)
counts_lin28 = counts_lin28[,6:ncol(counts_lin28)]
colnames(counts_lin28) = c("WT-TNa", "WT-MPa", "WT-IMPa", "Lin28b-TNa", "Lin28b-MPa", "Lin28b-IMPa",
                           "WT-TNb", "WT-MPb", "WT-IMPb", "Lin28b-TNb", "Lin28b-MPb", "Lin28b-IMPb",
                           "WT-TNc", "WT-MPc", "WT-IMPc", "Lin28b-TNc", "Lin28b-MPc", "Lin28b-IMPc",
                           "WT-TNd", "WT-MPd", "WT-IMPd", "Lin28b-TNd", "Lin28b-MPd", "Lin28b-IMPd") 
counts_lin28_noIMP = counts_lin28[, !grepl("IMP", colnames(counts_lin28))]
counts_lin28_noLin28 = counts_lin28_noIMP[, !grepl("Lin28", colnames(counts_lin28_noIMP))]
counts_lin28_noOutlier = counts_lin28_noLin28[, !grepl("WT-TNa", colnames(counts_lin28_noLin28))]
# counts_lin28_noOutlier = counts_lin28_noOutlier[, !grepl("Lin28b-MPc", colnames(counts_lin28_noOutlier))]

coldata_lin28_noOutlier = data.frame(c(rep(c('vm', 'tn'), 3), "vm"))
colnames(coldata_lin28_noOutlier) = 'type'
rownames(coldata_lin28_noOutlier) = colnames(counts_lin28_noOutlier)

dds <- DESeqDataSetFromMatrix(countData = counts_lin28_noOutlier,
                              colData = coldata_lin28_noOutlier,
                              design= ~ type)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="type_vm_vs_tn")

resSig <- subset(res, padj < padj_cutoff)
deg_lin28_tnVm = rownames(resSig)




counts_veteran = read.table("RNAseq/counts_name.txt", header = T, row.names = 1)
counts_veteran = counts_veteran[, 6:ncol(counts_veteran)]
colnames(counts_veteran)[1:8] = substr(colnames(counts_veteran)[1:8], 41, 48)
colnames(counts_veteran)[9:ncol(counts_veteran)] = substr(colnames(counts_veteran)[9:ncol(counts_veteran)], 11, 18)
colnames(counts_veteran) = c('RTE.TN1', 'RTE.MP1', 'Vet.MP2', 'Mat.MP2', 'RTE.MP2', 'Mat.TN3', 'RTE.TN4', 'RTE.MP4', 'Vet.TN1', 'RTE.TN3', 'Mat.MP3', 'RTE.MP3', 'Vet.TN4', 'Mat.TN4', 'Mat.TN1', 'Vet.MP4', 'Mat.MP4', 'Vet.MP1', 'Mat.MP1', 'Vet.TN2', 'Mat.TN2', 'RTE.TN2')
coldata_veteran = data.frame(gsub("mp", "vm", tolower(substr(gsub("^.*\\.", "", colnames(counts_veteran)), 1, 2))),
                             substr(colnames(counts_veteran), 1, 3))
colnames(coldata_veteran) = c('type', 'rte')
rownames(coldata_veteran) = colnames(counts_veteran)

dds <- DESeqDataSetFromMatrix(countData = counts_veteran,
                              colData = coldata_veteran,
                              design= ~ rte + type)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="type_vm_vs_tn")

resSig <- subset(res, padj < padj_cutoff)
deg_veteran_tnVM = rownames(resSig)


counts_dirty = read.table("Tabilas2022_RNAseq/counts_name.txt", header = T, row.names = 1)
counts_dirty = counts_dirty[,6:ncol(counts_dirty)]
coln_tmp = sub('.*CD', 'CD', colnames(counts_dirty))
coln_tmp = sub('_CKDL.*', '', coln_tmp)
targetFile = read.table('Tabilas2022_RNAseq/targetFile.txt', header = T)
targetFile$cd = sub('.ReadsPerGene.out.tab.rawCounts', '', targetFile$files)
colnames(counts_dirty) = targetFile[match(coln_tmp, targetFile$cd),]$label
counts_dirty_nobulk = counts_dirty[, !grepl("B.", colnames(counts_dirty))]

coldata_dirty_nobulk = data.frame(c("Neo", "Neo", "Neo", "Neo", "Adult", "Adult", "Adult", "Adult", "Neo", "Neo", "Neo", "Neo", "Adult", "Adult", "Adult", "Adult", "Neo", "Neo", "Neo", "Neo", "Adult", "Adult", "Adult", "Adult"),
                                  tolower(gsub("VM", "MP", substr(colnames(counts_dirty_nobulk), 1, 2))), 
                                  gsub("\\..*","", gsub("^.*?\\.","", colnames(counts_dirty_nobulk))))
colnames(coldata_dirty_nobulk) = c('age', 'type', 'clean')
rownames(coldata_dirty_nobulk) = colnames(counts_dirty_nobulk)


dds <- DESeqDataSetFromMatrix(countData = counts_dirty_nobulk,
                              colData = coldata_dirty_nobulk,
                              design= ~ clean + type + age)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="type_tn_vs_mp")

resSig <- subset(res, padj < padj_cutoff)
deg_dirty_tnVm = rownames(resSig)



# clean vs dirty

counts_dirty = read.table("Tabilas2022_RNAseq/counts_name.txt", header = T, row.names = 1)
counts_dirty = counts_dirty[,6:ncol(counts_dirty)]
coln_tmp = sub('.*CD', 'CD', colnames(counts_dirty))
coln_tmp = sub('_CKDL.*', '', coln_tmp)
targetFile = read.table('Tabilas2022_RNAseq/targetFile.txt', header = T)
targetFile$cd = sub('.ReadsPerGene.out.tab.rawCounts', '', targetFile$files)
colnames(counts_dirty) = targetFile[match(coln_tmp, targetFile$cd),]$label
counts_dirty_nobulk = counts_dirty[, !grepl("B.", colnames(counts_dirty))]

coldata_dirty_nobulk = data.frame(c("Neo", "Neo", "Neo", "Neo", "Adult", "Adult", "Adult", "Adult", "Neo", "Neo", "Neo", "Neo", "Adult", "Adult", "Adult", "Adult", "Neo", "Neo", "Neo", "Neo", "Adult", "Adult", "Adult", "Adult"),
                                  tolower(gsub("VM", "MP", substr(colnames(counts_dirty_nobulk), 1, 2))), 
                                  gsub("\\..*","", gsub("^.*?\\.","", colnames(counts_dirty_nobulk))))
colnames(coldata_dirty_nobulk) = c('age', 'type', 'clean')
rownames(coldata_dirty_nobulk) = colnames(counts_dirty_nobulk)


dds <- DESeqDataSetFromMatrix(countData = counts_dirty_nobulk,
                              colData = coldata_dirty_nobulk,
                              design= ~ clean + type + age)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="clean_D_vs_C")

resSig <- subset(res, padj < padj_cutoff)
deg_dirty_Cd = rownames(resSig)



# RTE vs Mat vs Vet

counts_veteran = read.table("RNAseq/counts_name.txt", header = T, row.names = 1)
counts_veteran = counts_veteran[, 6:ncol(counts_veteran)]
colnames(counts_veteran)[1:8] = substr(colnames(counts_veteran)[1:8], 41, 48)
colnames(counts_veteran)[9:ncol(counts_veteran)] = substr(colnames(counts_veteran)[9:ncol(counts_veteran)], 11, 18)
colnames(counts_veteran) = c('RTE.TN1', 'RTE.MP1', 'Vet.MP2', 'Mat.MP2', 'RTE.MP2', 'Mat.TN3', 'RTE.TN4', 'RTE.MP4', 'Vet.TN1', 'RTE.TN3', 'Mat.MP3', 'RTE.MP3', 'Vet.TN4', 'Mat.TN4', 'Mat.TN1', 'Vet.MP4', 'Mat.MP4', 'Vet.MP1', 'Mat.MP1', 'Vet.TN2', 'Mat.TN2', 'RTE.TN2')
coldata_veteran = data.frame(tolower(substr(gsub("^.*\\.", "", colnames(counts_veteran)), 1, 2)),
                             substr(colnames(counts_veteran), 1, 3))
colnames(coldata_veteran) = c('type', 'rte')
rownames(coldata_veteran) = colnames(counts_veteran)

dds <- DESeqDataSetFromMatrix(countData = counts_veteran,
                              colData = coldata_veteran,
                              design= ~ rte + type)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients


res_RTE_vs_Mat <- results(dds, name="rte_RTE_vs_Mat")
resSig_RTE_vs_Mat <- subset(res_RTE_vs_Mat, padj < padj_cutoff)
deg_RTE_vs_Mat = rownames(resSig_RTE_vs_Mat)

res_Vet_vs_Mat <- results(dds, name="rte_Vet_vs_Mat")
resSig_Vet_vs_Mat <- subset(res_Vet_vs_Mat, padj < padj_cutoff)
deg_Vet_vs_Mat = rownames(resSig_Vet_vs_Mat)




deg = unique(c(deg_cell2018_aduNeo, deg_mir29_aduNeo, deg_dirty_aduNeo, 
               deg_cell2018_tnVm, deg_mir29_tnVM, deg_lin28_tnVm, deg_veteran_tnVM, deg_dirty_tnVm, 
               deg_RTE_vs_Mat, deg_Vet_vs_Mat,
               deg_dirty_Cd))
               
write.table(deg, "input/genes_DE/deg_padj01.txt", row.names = FALSE, col.names=FALSE, sep = "\t", quote = F)
