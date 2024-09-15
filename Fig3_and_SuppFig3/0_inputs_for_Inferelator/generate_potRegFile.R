# get potRegFile by intersect a list of regulators with genes of interest in naive CD8+ T cell subsets

# Curated list of regulators from Supplemental_Table_S3.xlsx is from Miraldi et al. (2018) "Leveraging chromatin accessibility for transcriptional regulatory network inference in T Helper 17 Cells"
library(readxl)
TF_tab = read_excel("TF/Supplemental_Table_S3.xlsx", sheet = "Updated TF List")
TF_list = TF_tab$`Updated TF List (1577 putative TFs)`

# interset with genes of interest in naive CD8+ T cell subsets
DEgenes_wFC = read.csv("inputs/genes_DE/deg_padj01.txt", col.names = F)
colnames(DEgenes_wFC) = 'DEgenes_wFC'
DEgenes_wFC = DEgenes_wFC$DEgenes_wFC

regulators = intersect(DEgenes_wFC, TF_list)
write.table(regulators, "TF/regulators.txt", quote = F, row.names = F, col.names = F)

