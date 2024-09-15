### get Eomes targets

num_nets = 2
ms = 10
bias = 50
net = read.table(paste0('chip_ko_network_atac_ms', ms, '_bias', bias, '_maxComb_cut01_sharedbyMorethan', num_nets, 'nets_sp.tsv'),
                 header = T)

net_Eomes = net[net$TF == 'Eomes',]
write.table(net_Eomes, 'net_Eomes.tsv', sep = '\t', quote = F, row.names = F)


### make network plot for Eomes (Fig. 3B)

atac_sp = read.table("inputs/prior/output_prior/prior_atac_b_sp.tsv", header = T)

atac_sp_tmp = atac_sp
atac_sp_tmp$Regulation = paste0(atac_sp_tmp$TF, '->',  atac_sp_tmp$Target)
atac_sp_tmp_tf = atac_sp_tmp[atac_sp_tmp$TF %in% unique(net$TF), ]
atac_sp_tmp_tftarget = atac_sp_tmp_tf[atac_sp_tmp_tf$Target %in% unique(net$Target), ]
net_tmp = net
net_tmp$Regulation = paste0(net_tmp$TF, '->',  net_tmp$Target)
net_tmp$If_prior_support <- ifelse(net_tmp$Regulation %in% atac_sp_tmp_tftarget$Regulation,"1-Supported by ATAC-seq","2-Not supported by ATAC-seq")
net_tmp$Direction = ifelse(net_tmp$SignedQuantile == 1, 'Activation', 'Repression')

net_Eomes = net_tmp[net_tmp$TF == 'Eomes',]
table(net_Eomes$Direction)
#Activation Repression 
#79         23 
table(net_Eomes$If_prior_support)
#1-Supported by ATAC-seq 2-Not supported by ATAC-seq 
#92                          10 
nodes = data.frame(gene = unique(c(net_Eomes$TF, net_Eomes$Target)), type = NA)
nodes$type = ifelse(nodes$gene %in% unique(net$TF), 'TF', 'Target')
g = graph_from_data_frame(net_Eomes, directed = TRUE, vertices = nodes)
lay = create_layout(g, layout = "fr")
ggraph(g, layout = 'graphopt') + 
  geom_edge_link(aes(start_cap = label_rect(node1.name), end_cap = label_rect(node2.name), linetype = If_prior_support,color = Direction), 
                 arrow = arrow(type = "closed", length = unit(1, 'mm'))) +
  geom_node_text(aes(label = name, color = type)) +
  scale_color_manual(limits = c('TF', 'Target'),
                     values = c('#5b0e2d', '#ffa781')) +
  theme_graph() 
ggsave('Eomes_targets.svg',
       plot = last_plot(),
       device = "svg",
       width = 14,
       height = 10,
       dpi = 300
)

# GSEA of Eomes targets using ImmuneSigDB and GO Molecular Functions were performed with tfTarget_GSEA.sh from https://github.com/emiraldi/infTRN_lassoStARS
