### get differentially accessible peaks

counts = read.table('counts_ATACseq',
                    row.names = 1, header = T)
counts = counts[, 6:ncol(counts)]

library(readxl)
sample_info = read_excel('sample_info_ATACseq.xlsx')
gsub('/|-', '.', paste0('X', sample_info$File_path, sample_info$File_name)) == colnames(counts) # all true
colnames(counts) = sample_info$Sample

counts_cell2018 = counts[, sample_info$Project == 'Cell2018']
counts_mir29project = counts[, sample_info$Project == 'mir29project']
counts_cleanproject = counts[, sample_info$Project == 'cleanproject']
counts_veteranproject = counts[, sample_info$Project == 'veteranproject']

library(DESeq2)

padj_cutoff = 0.05


# clean dirty

counts_cleanproject_nobulk = counts_cleanproject[, !grepl("bulk_", colnames(counts_cleanproject))]

coldata_cleanproject_nobulk = data.frame(rep(c('Neo', 'Neo', 'Neo', 'Neo', 'Adult', 'Adult', 'Adult', 'Adult'), 3),
                                  toupper(substr(colnames(counts_cleanproject_nobulk), 1, 2)), 
                                  substr(colnames(counts_cleanproject_nobulk), 4, 8))
colnames(coldata_cleanproject_nobulk) = c('age', 'type', 'clean')
rownames(coldata_cleanproject_nobulk) = colnames(counts_cleanproject_nobulk)
coldata_cleanproject_nobulk$age = as.factor(coldata_cleanproject_nobulk$age)
coldata_cleanproject_nobulk$type = as.factor(coldata_cleanproject_nobulk$type)
coldata_cleanproject_nobulk$clean = as.factor(coldata_cleanproject_nobulk$clean)

dds <- DESeqDataSetFromMatrix(countData = counts_cleanproject_nobulk,
                              colData = coldata_cleanproject_nobulk,
                              design= ~ age + type + clean)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res.dirty_clean <- results(dds, name="clean_dirty_vs_clean")
dirty <- as.data.frame(subset(res.dirty_clean, log2FoldChange > 0 & padj < 0.05))
clean <- as.data.frame(subset(res.dirty_clean, log2FoldChange < 0 & padj < 0.05))
dirty_clean_insig = as.data.frame(subset(res.dirty_clean, padj > 0.05))


# adult neo

counts_cleanproject_nobulk = counts_cleanproject[, !grepl("bulk_", colnames(counts_cleanproject))]

coldata_cleanproject_nobulk = data.frame(rep(c('Neo', 'Neo', 'Neo', 'Neo', 'Adult', 'Adult', 'Adult', 'Adult'), 3),
                                         toupper(substr(colnames(counts_cleanproject_nobulk), 1, 2)), 
                                         substr(colnames(counts_cleanproject_nobulk), 4, 8))
colnames(coldata_cleanproject_nobulk) = c('age', 'type', 'clean')
rownames(coldata_cleanproject_nobulk) = colnames(counts_cleanproject_nobulk)
coldata_cleanproject_nobulk$age = as.factor(coldata_cleanproject_nobulk$age)
coldata_cleanproject_nobulk$type = as.factor(coldata_cleanproject_nobulk$type)
coldata_cleanproject_nobulk$clean = as.factor(coldata_cleanproject_nobulk$clean)
coldata_cleanproject_nobulk$source = 'Dirty_Project'


counts_cell2018 = counts_cell2018[, c('VM_1dayA', 'VM_1dayB', 'VM_28dayA', 'VM_28dayB')]

coldata_cell2018 = data.frame(c('Neo', 'Neo', 'Adult', 'Adult'),
                              c('VM', 'VM', 'VM', 'VM'))
colnames(coldata_cell2018) = c('age', 'type')
rownames(coldata_cell2018) = colnames(counts_cell2018)
coldata_cell2018$clean = 'clean'
coldata_cell2018$source = 'Cell_2018'

counts_adultNeo = cbind(counts_cell2018, counts_cleanproject_nobulk)
coldata_adultNeo = rbind(coldata_cell2018, coldata_cleanproject_nobulk)
coldata_adultNeo$age = as.factor(coldata_adultNeo$age)
coldata_adultNeo$type = as.factor(coldata_adultNeo$type)
coldata_adultNeo$clean = as.factor(coldata_adultNeo$clean)
coldata_adultNeo$source = as.factor(coldata_adultNeo$source)

dds <- DESeqDataSetFromMatrix(countData = counts_adultNeo,
                              colData = coldata_adultNeo,
                              design = ~ age + type + clean + source)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients

res.neo_adult <- results(dds, contrast=c("age", "Neo", "Adult"))
neo <- as.data.frame(subset(res.neo_adult, log2FoldChange > 0 & padj < 0.05))
adult <- as.data.frame(subset(res.neo_adult, log2FoldChange < 0 & padj < 0.05))
neo_adult_insig = as.data.frame(subset(res.neo_adult, padj > 0.05))


# tn vs vm

counts_cleanproject_nobulk = counts_cleanproject[, !grepl("bulk_", colnames(counts_cleanproject))]

coldata_cleanproject_nobulk = data.frame(rep(c('Neo', 'Neo', 'Neo', 'Neo', 'Adult', 'Adult', 'Adult', 'Adult'), 3),
                                         toupper(substr(colnames(counts_cleanproject_nobulk), 1, 2)), 
                                         substr(colnames(counts_cleanproject_nobulk), 4, 8))
colnames(coldata_cleanproject_nobulk) = c('age', 'type', 'clean')
rownames(coldata_cleanproject_nobulk) = colnames(counts_cleanproject_nobulk)
coldata_cleanproject_nobulk$age = as.factor(coldata_cleanproject_nobulk$age)
coldata_cleanproject_nobulk$type = as.factor(coldata_cleanproject_nobulk$type)
coldata_cleanproject_nobulk$clean = as.factor(coldata_cleanproject_nobulk$clean)
coldata_cleanproject_nobulk$source = 'Dirty_Project'
coldata_cleanproject_nobulk$rte = 'mature'

counts_cell2018 = counts_cell2018[, c('VM_1dayA', 'VM_1dayB', 'VM_28dayA', 'VM_28dayB')]
coldata_cell2018 = data.frame(c('Neo', 'Neo', 'Adult', 'Adult'),
                              c('VM', 'VM', 'VM', 'VM'))
colnames(coldata_cell2018) = c('age', 'type')
rownames(coldata_cell2018) = colnames(counts_cell2018)
coldata_cell2018$clean = 'clean'
coldata_cell2018$source = 'Cell_2018'
coldata_cell2018$rte = 'mature'

coldata_veteran = data.frame(substr(colnames(counts_veteranproject), 1, nchar(colnames(counts_veteranproject))-2))
colnames(coldata_veteran) = c('tmp')
rownames(coldata_veteran) = colnames(counts_veteranproject)
coldata_veteran$type = substr(colnames(counts_veteranproject), 1, 2)
coldata_veteran$rte = substr(coldata_veteran$tmp, 4, nchar(coldata_veteran$tmp))
coldata_veteran$type = as.factor(coldata_veteran$type)
coldata_veteran$rte = as.factor(coldata_veteran$rte)
coldata_veteran$tmp = NULL
coldata_veteran$clean = 'clean'
coldata_veteran$source = 'RTE_Project'
coldata_veteran$age = 'Adult'
coldata_veteran = coldata_veteran[, c('age', 'type', 'clean', 'source', 'rte')]


counts_tnVm = do.call("cbind", list(counts_cell2018, counts_cleanproject_nobulk, counts_veteranproject))
coldata_tnVm = do.call("rbind", list(coldata_cell2018, coldata_cleanproject_nobulk, coldata_veteran))
coldata_tnVm$age = as.factor(coldata_tnVm$age)
coldata_tnVm$type = as.factor(coldata_tnVm$type)
coldata_tnVm$clean = as.factor(coldata_tnVm$clean)
coldata_tnVm$source = as.factor(coldata_tnVm$source)
coldata_tnVm$rte = as.factor(coldata_tnVm$rte)

dds <- DESeqDataSetFromMatrix(countData = counts_tnVm,
                              colData = coldata_tnVm,
                              design = ~ age + type + clean + source + rte)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients


res.VM_TN <- results(dds, contrast=c("type","VM","TN"))
vm <- as.data.frame(subset(res.VM_TN, log2FoldChange > 0 & padj < 0.05))
tn <- as.data.frame(subset(res.VM_TN, log2FoldChange < 0 & padj < 0.05))
vm_tn_insig = as.data.frame(subset(res.VM_TN, padj > 0.05))



write.table(dirty, file = "DE/res/dirty.txt", row.names = T, col.names = F, quote = F, sep = "\t")
write.table(clean, file = "DE/res/clean.txt", row.names = T, col.names = F, quote = F, sep = "\t")
write.table(dirty_clean_insig, file = "DE/res/insig_dirty_clean.txt", row.names = T, col.names = F, quote = F, sep = "\t")

write.table(neo, file = "DE/res/neo.txt", row.names = T, col.names = F, quote = F, sep = "\t")
write.table(adult, file = "DE/res/adult.txt", row.names = T, col.names = F, quote = F, sep = "\t")
write.table(neo_adult_insig, file = "DE/res/insig_neo_adult.txt", row.names = T, col.names = F, quote = F, sep = "\t")

write.table(vm, file = "DE/res/vm.txt", row.names = T, col.names = F, quote = F, sep = "\t")
write.table(tn, file = "DE/res/tn.txt", row.names = T, col.names = F, quote = F, sep = "\t")
write.table(vm_tn_insig, file = "DE/res/insig_vm_tn.txt", row.names = T, col.names = F, quote = F, sep = "\t")
