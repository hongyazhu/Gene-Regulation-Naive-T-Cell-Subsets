# combining TFA and TF mRNA TRNs was done with combine_Th17_TRNs.m from https://github.com/emiraldi/infTRN_lassoStARS
# filtering was done with filter_Th17_TRNs_by_pcorr.sh from https://github.com/emiraldi/infTRN_lassoStARS

# I ran inferelator 5 times with 5 different seeds. Here to find TF-gene regulations shared by more than 2 runs (3 or more)
sp = read.table('chip_ko_network_atac_ms10_bias50_maxComb/chip_ko_network_atac_ms10_bias50_maxComb_cut01_sp.tsv', sep = '\t', header = T)
sp_seed7 = read.table('chip_ko_network_atac_ms10_bias50_seed26_maxComb/chip_ko_network_atac_ms10_bias50_seed26_maxComb_cut01_sp.tsv', sep = '\t', header = T)
sp_seed99 = read.table('chip_ko_network_atac_ms10_bias50_seed57_maxComb/chip_ko_network_atac_ms10_bias50_seed57_maxComb_cut01_sp.tsv', sep = '\t', header = T)
sp_seed26 = read.table('chip_ko_network_atac_ms10_bias50_seed7_maxComb/chip_ko_network_atac_ms10_bias50_seed7_maxComb_cut01_sp.tsv', sep = '\t', header = T)
sp_seed57 = read.table('chip_ko_network_atac_ms10_bias50_seed99_maxComb/chip_ko_network_atac_ms10_bias50_seed99_maxComb_cut01_sp.tsv', sep = '\t', header = T)

sp$direction = NA
for (i in 1:nrow(sp)){
  if (sp[i,'SignedQuantile'] > 0){
    sp[i,'direction'] = 'pos'
  } else if (sp[i,'SignedQuantile'] < 0){
    sp[i,'direction'] = 'neg'
  }
}
sp$'TF->Target' = paste(sp$TF, sp$Target, sep = '->')
sp$'TF->Target,direction' = paste(sp$'TF->Target', sp$direction, sep = ',')

sp_seed7$direction = NA
for (i in 1:nrow(sp_seed7)){
  if (sp_seed7[i,'SignedQuantile'] > 0){
    sp_seed7[i,'direction'] = 'pos'
  } else if (sp_seed7[i,'SignedQuantile'] < 0){
    sp_seed7[i,'direction'] = 'neg'
  }
}
sp_seed7$'TF->Target' = paste(sp_seed7$TF, sp_seed7$Target, sep = '->')
sp_seed7$'TF->Target,direction' = paste(sp_seed7$'TF->Target', sp_seed7$direction, sep = ',')

sp_seed99$direction = NA
for (i in 1:nrow(sp_seed99)){
  if (sp_seed99[i,'SignedQuantile'] > 0){
    sp_seed99[i,'direction'] = 'pos'
  } else if (sp_seed99[i,'SignedQuantile'] < 0){
    sp_seed99[i,'direction'] = 'neg'
  }
}
sp_seed99$'TF->Target' = paste(sp_seed99$TF, sp_seed99$Target, sep = '->')
sp_seed99$'TF->Target,direction' = paste(sp_seed99$'TF->Target', sp_seed99$direction, sep = ',')

sp_seed26$direction = NA
for (i in 1:nrow(sp_seed26)){
  if (sp_seed26[i,'SignedQuantile'] > 0){
    sp_seed26[i,'direction'] = 'pos'
  } else if (sp_seed26[i,'SignedQuantile'] < 0){
    sp_seed26[i,'direction'] = 'neg'
  }
}
sp_seed26$'TF->Target' = paste(sp_seed26$TF, sp_seed26$Target, sep = '->')
sp_seed26$'TF->Target,direction' = paste(sp_seed26$'TF->Target', sp_seed26$direction, sep = ',')

sp_seed57$direction = NA
for (i in 1:nrow(sp_seed57)){
  if (sp_seed57[i,'SignedQuantile'] > 0){
    sp_seed57[i,'direction'] = 'pos'
  } else if (sp_seed57[i,'SignedQuantile'] < 0){
    sp_seed57[i,'direction'] = 'neg'
  }
}
sp_seed57$'TF->Target' = paste(sp_seed57$TF, sp_seed57$Target, sep = '->')
sp_seed57$'TF->Target,direction' = paste(sp_seed57$'TF->Target', sp_seed57$direction, sep = ',')

#table(sort(table(c(sp$`TF->Target,direction`, sp_seed7$`TF->Target,direction`, sp_seed99$`TF->Target,direction`, sp_seed26$`TF->Target,direction`, sp_seed57$`TF->Target,direction`)), decreasing = T) == 5)
table(sort(table(c(sp$`TF->Target,direction`, sp_seed7$`TF->Target,direction`, sp_seed99$`TF->Target,direction`, sp_seed26$`TF->Target,direction`, sp_seed57$`TF->Target,direction`)), decreasing = T) > 4)
# 36008
table(sort(table(c(sp$`TF->Target,direction`, sp_seed7$`TF->Target,direction`, sp_seed99$`TF->Target,direction`, sp_seed26$`TF->Target,direction`, sp_seed57$`TF->Target,direction`)), decreasing = T) > 3)
# 50306
table(sort(table(c(sp$`TF->Target,direction`, sp_seed7$`TF->Target,direction`, sp_seed99$`TF->Target,direction`, sp_seed26$`TF->Target,direction`, sp_seed57$`TF->Target,direction`)), decreasing = T) > 2)
# 67620
table(sort(table(c(sp$`TF->Target,direction`, sp_seed7$`TF->Target,direction`, sp_seed99$`TF->Target,direction`, sp_seed26$`TF->Target,direction`, sp_seed57$`TF->Target,direction`)), decreasing = T) > 1)
# 95734
table(sort(table(c(sp$`TF->Target,direction`, sp_seed7$`TF->Target,direction`, sp_seed99$`TF->Target,direction`, sp_seed26$`TF->Target,direction`, sp_seed57$`TF->Target,direction`)), decreasing = T) > 0)
# 167859 in total



res_table = table(c(sp$`TF->Target,direction`, sp_seed7$`TF->Target,direction`, sp_seed99$`TF->Target,direction`, sp_seed26$`TF->Target,direction`, sp_seed57$`TF->Target,direction`))
#sp_all = merge(sp, sp_seed7, sp_seed99, sp_seed26, sp_seed57, by = 'TF->Target,direction', all = T)

MyMerge <- function(x, y){
  df <- merge(x, y, by= "TF->Target,direction", all = TRUE)
  return(df)
}
sp_all <- Reduce(MyMerge, list(sp, sp_seed7, sp_seed99, sp_seed26, sp_seed57))

# num = 4

for (num in 2){ # 1:4
  
  sp_interest_names = names(res_table[res_table > num])
  
  sp_interest = sp_all[sp_all$`TF->Target,direction` %in% sp_interest_names, ]
  
  sp_interest_less = sp_interest[, c('TF->Target,direction', 'TF.x')]
  sp_interest_less$TF = gsub('->.*', '', sp_interest_less$`TF->Target,direction`)
  sp_interest_less$Target = gsub(',.*', '', gsub('.*->', '', sp_interest_less$`TF->Target,direction`))
  sp_interest_less$direction = gsub('.*,', '', sp_interest_less$`TF->Target,direction`)
  sp_interest_less$SignedQuantile = NA
  for (i in 1:nrow(sp_interest_less)){
    if (sp_interest_less[i,'direction'] == 'pos'){
      sp_interest_less[i,'SignedQuantile'] = 1
    } else if (sp_interest_less[i,'direction'] == 'neg'){
      sp_interest_less[i,'SignedQuantile'] = -1
    }
  }
  sp_interest_format = sp_interest_less[, c('TF', 'Target', 'SignedQuantile')]
  print(num)
  print(nrow(sp_interest_format))
  write.table(sp_interest_format, paste0('runs/sharedbyMorethan', num, 'netsFrom5seeds/chip_ko_network_atac_ms10_bias50_maxComb_cut01_sharedbyMorethan', num, 'nets_sp.tsv'), sep = '\t', quote = F, row.names = F)  
  
  
  notsp = data.frame(matrix(0, ncol = length(unique(sp_interest_format$TF)), nrow = length(unique(sp_interest_format$Target))))
  colnames(notsp) = unique(sp_interest_format$TF)
  rownames(notsp) = unique(sp_interest_format$Target)
  
  for (i in 1:nrow(sp_interest_format))
    if (sp_interest_format[i, 'SignedQuantile'] > 0){
      notsp[sp_interest_format[i,'Target'], sp_interest_format[i,'TF']] = 1
    } else if (sp_interest_format[i, 'SignedQuantile'] < 0){
      notsp[sp_interest_format[i,'Target'], sp_interest_format[i,'TF']] = -1
    }
  
  notsp_new <- notsp[ order(row.names(notsp)), ]
  notsp_new <- notsp[ , order(colnames(notsp))]
  
  
  write.table(notsp_new, paste0('runs/sharedbyMorethan', num, 'netsFrom5seeds/chip_ko_network_atac_ms10_bias50_maxComb_cut01_sharedbyMorethan', num, 'nets.tsv'), sep = '\t', quote = F, col.names=NA)
}

