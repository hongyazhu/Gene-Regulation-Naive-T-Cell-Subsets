# run bagfoot analysis (Fig 2H-J)

# prep

library(bagfoot)

bamfile1 = "adult_sort.bam" 
cc1 = countReadsBAM(bamfile1);  # counts the number of cuts in the sequence file.
cutcountfile1 = makeCutCountBAM(bamfile1, refgenome = "mm10");   # generate a BedGraph file with DNase Cleavages counts
MAPPABILITY_FILES_DIRECTORY_MM10<<-'/home/hz543/ggbackup/workdir_210227/share/mappibility_35';   
tabMappability10 = MakeBiasCorrectionTableBAM(bamfile= bamfile1,
   outfile="Hexamer_adult_mm10_withMap.txt", 
   refgenome="mm10", 
   np=6,
   mapdir = MAPPABILITY_FILES_DIRECTORY_MM10,
   atac = F);

bamfile1 = "neo_sort.bam" 
cc1 = countReadsBAM(bamfile1);  # counts the number of cuts in the sequence file.
cutcountfile1 = makeCutCountBAM(bamfile1, refgenome = "mm10");   # generate a BedGraph file with DNase Cleavages counts
tabMappability10 = MakeBiasCorrectionTableBAM(bamfile= bamfile1,
   outfile="Hexamer_neo_mm10_withMap.txt", 
   refgenome="mm10", 
   np=6,
   mapdir = MAPPABILITY_FILES_DIRECTORY_MM10,
   atac = F);

bamfile1 = "clean_sort.bam" 
cc1 = countReadsBAM(bamfile1);  # counts the number of cuts in the sequence file.
cutcountfile1 = makeCutCountBAM(bamfile1, refgenome = "mm10");   # generate a BedGraph file with DNase Cleavages counts
tabMappability10 = MakeBiasCorrectionTableBAM(bamfile= bamfile1,
   outfile="Hexamer_clean_mm10_withMap.txt", 
   refgenome="mm10", 
   np=6,
   mapdir = MAPPABILITY_FILES_DIRECTORY_MM10,
   atac = F);

bamfile1 = "dirty_sort.bam" 
cc1 = countReadsBAM(bamfile1);  # counts the number of cuts in the sequence file.
cutcountfile1 = makeCutCountBAM(bamfile1, refgenome = "mm10");   # generate a BedGraph file with DNase Cleavages counts
tabMappability10 = MakeBiasCorrectionTableBAM(bamfile= bamfile1,
   outfile="Hexamer_dirty_mm10_withMap.txt", 
   refgenome="mm10", 
   np=6,
   mapdir = MAPPABILITY_FILES_DIRECTORY_MM10,
   atac = F);

bamfile1 = "tn_sort.bam"
cc1 = countReadsBAM(bamfile1);  # counts the number of cuts in the sequence file.
cutcountfile1 = makeCutCountBAM(bamfile1, refgenome = "mm10");   # generate a BedGraph file with DNase Cleavages counts
tabMappability10 = MakeBiasCorrectionTableBAM(bamfile= bamfile1,
   outfile="Hexamer_tn_mm10_withMap.txt", 
   refgenome="mm10", 
   np=6,
   mapdir = MAPPABILITY_FILES_DIRECTORY_MM10,
   atac = F);

bamfile1 ="vm_sort.bam"
cc1 = countReadsBAM(bamfile1);  # counts the number of cuts in the sequence file.
cutcountfile1 = makeCutCountBAM(bamfile1, refgenome = "mm10");   # generate a BedGraph file with DNase Cleavages counts
tabMappability10 = MakeBiasCorrectionTableBAM(bamfile= bamfile1,
   outfile="Hexamer_vm_mm10_withMap.txt", 
   refgenome="mm10", 
   np=6,
   mapdir = MAPPABILITY_FILES_DIRECTORY_MM10,
   atac = F);


# run

library('bagfoot')

calcCutCountInHotsplot <- function(combinedhotspotfile, bedgraphFile) {
	tmpBed = paste0(sample(1:10000000, 1), ".bed")
	if (file.exists(tmpBed)) {
		tmpBed = paste0(sample(1:10000000, 1), ".bed")
		if (file.exists(tmpBed)) {
			tmpBed = paste0(sample(1:10000000, 1), ".bed")
		}
	}
	system(paste0(" perl -ne 's/,/\t/g; s/\\\"//g; print;' ", combinedhotspotfile, " | grep -v \"ID\" | cut -f2-4 | sort -k1,1 -k2,2n > ", tmpBed ))
	count <- system(paste0("less ", bedgraphFile, " | grep -v \"desc\" | perl -ne 's/ /\t/g; print;' | sort -k1,1 -k2,2n | bedtools intersect -a ", tmpBed, " -b stdin -wa -wb | cut -f7 | perl -ne '$s+=$_; END{print \"$s\n\";}'"), intern=T)
	file.remove(tmpBed)
	return(as.numeric(as.character(count)))
}


prepNRunBaGFoot <- function (hotspotfile1, hotspotfile2, alias1, alias2) {  
	print("# Combining hotspot files...")
	combinedhotspotfile = combineTwoHotspots( hotspotfile1, hotspotfile2, alias1, alias2)

	bedgraph1 = paste0("../cutcount/", alias1, "_sort_AC_cutcount.bgr.gz") 
	bedgraph2 = paste0("../cutcount/", alias2, "_sort_AC_cutcount.bgr.gz") 

	print("# Counting number of cuts in combined hotspots...")
	rCount1 = calcCutCountInHotsplot(combinedhotspotfile, bedgraph1)
	rCount2 = calcCutCountInHotsplot(combinedhotspotfile, bedgraph2)
	print(paste0("Count1: ", rCount1, "  Count2: ", rCount2))


	print("# Mergning Hex files...")
	# Merging two Hex files (withMap).
	# The merging code from https://www.biostars.org/p/279020/
	hex.cond1 = read.table(paste0("../hexamer/Hexamer_", alias1, "_mm10_withMap.txt"), header=T, row.names=1)
	hex.cond2 = read.table(paste0("../hexamer/Hexamer_", alias2, "_mm10_withMap.txt"), header=T, row.names=1)
	hex.both = data.frame(ObCuts = hex.cond1$ObCuts + hex.cond2$ObCuts, GenomicPositionCount = hex.cond1$GenomicPositionCount, ObCutRatio=1, GPRatio = hex.cond1$GPRatio, CorrectionFactor = 1, row.names=row.names(hex.cond1))
	hex.both$ObCutRatio = hex.both$ObCuts / sum(hex.both[ ! grepl("other", row.names(hex.both)), ]$ObCuts)
	hex.both$CorrectionFactor = hex.both$GPRatio / hex.both$ObCutRatio
	hex.both[ grepl("other", row.names(hex.both)), ]$ObCutRatio = NA
	hex.both[ grepl("other", row.names(hex.both)), ]$CorrectionFactor = 1
	write.table(hex.both, paste0("Hexamer_", alias1, "_", alias2, "_mm10_withMap.txt"))


	print("# Setting-up options for the analysis...")
	cond1.count <- CutCount(file = bedgraph1,   # BedGraph file of cut counts (see "bigfoot_prep_example.R")
		                                                                count = rCount1,     # Found using featureCount    # Number of cuts within the hotspots
		                                                                name = alias1);      # Data Name 

	cond2.count <- CutCount(file = bedgraph2, # BedGraph file of cut counts (see "bigfoot_prep_example.R")
		                                                                count = rCount2,     # Number of cuts within the hotspots
		                                                                name = alias2);      # Data Name    

	motifdb_mm10 <- MotifDB(motiflistfile = 'bagfoot/motiflist_mm10.txt',   # This file contains the list of FIMO outputs file 
	     directory='bagfoot/mm10_fimo_m/');   # Directory of FIMO outputs 


	################ cond1 VS cond2 ##############
	gfootoption_cond1_cond2 <- GFootOption(biasfile = paste0("Hexamer_", alias1, "_", alias2, "_mm10_withMap.txt"),   # Hexamer bias frequency table (see "bigfoot_prep_example.R")
		                                        hexamer_pattern= "../nuccode/nuccode_mm10_6mer_{chr}.dat");                 # nuc. code files generated by MakeBiasCorrectionTableBAM 


	# Definition of BiGFoot analysis 
	GFoot_cond1_cond2 = GFoot(control = cond1.count,               ##  Control Data defined above
		                                        treatment = cond2.count,              ##  Treatment Data defined above
		                                        sitefile = paste0("pooled_", alias1, "_", alias2, "_hotspot.csv"),  ##  Sites file: Peak File.  This file must have "chr", "st", "ed" in the header.
		                                        motifDB = motifdb_mm10,                 ##  Information about FIMO motifs
		                                        gfootoption= gfootoption_cond1_cond2,   ##  Options as defined above
		                                        outputdir = paste0("OUTPUT_", alias1, "_vs_", alias2),  ## OUTPUT directory
		                                        cachedir = paste0("CACHE_", alias1, "_vs_", alias2) );   ## CACHE directory 

	print("# Running the analysis...")	
	GFoot_cond1_cond2 <- run(GFoot_cond1_cond2, graphout=T, yrange=c(-2.2,1.0), mc.cores=5.0);   # if graphout is T, aggregation & log-ratio plots are to be generated.  
	# As an output, the output file 'fed_fasted_On_pooled_fed_fasted_merged_hotspot_chr7_footprint_depth_table.csv' will be generated.

	print("# Generating final plots...")
	outfilename = paste0(alias1, "_", alias2, "_On_", "pooled_", alias1, "_", alias2, "_hotspot", "_footprint_depth_table.csv")
	dat= read.table(outfilename);
#	gen_bagplot(dat, dataname1=alias1, dataname2=alias2, factor=1.5);   # generates a bagplot from the calculated footprinting depths
    gen_bagplot_chisq(dat, dataname1=alias1, dataname2=alias2, factor=1.5);   #
}


prepNRunBaGFoot('adult.csv', 'neo.csv', 'adult', 'neo')
prepNRunBaGFoot('clean.csv', 'dirty.csv', 'clean', 'dirty')
prepNRunBaGFoot('tn.csv', 'vm.csv', tn, vm)


# remake figures

bagfootRes_dirty_clean = read.table("bagfoot/res_dirty_clean/bagplot_cutcount_diff_total_footprinting_depth_dirty-clean_qvalue_bagplot_output.csv", 
                                 header = T,
                                 sep = ",")
colnames(bagfootRes_dirty_clean)[3:4] = c("x", "y")
# prove plot should look same as automated generated one
library(aplpack)
bagplot(bagfootRes_dirty_clean$x, bagfootRes_dirty_clean$y, factor = 1.5,
        show.looppoints = FALSE, show.bagpoints = FALSE,
        show.baghull = TRUE, show.loophull = TRUE,
        show.whiskers = F, col.looppoint = "#000000",
        cex = 2, pch = 0, cex.lab = 1.5, cex.axis = 1.5,
        lwd = 3)

bagfootRes_dirty_clean_xy = bagfootRes_dirty_clean[,c("x", "y")]
rownames(bagfootRes_dirty_clean_xy) = bagfootRes_dirty_clean_xy$name


library(aplpack)
bag <- compute.bagplot(bagfootRes_dirty_clean_xy, approx.limit = nrow(bagfootRes_dirty_clean_xy),
                       factor = 1.5)

hull.loop <- data.frame(x = bag$hull.loop[,1], y = bag$hull.loop[,2])
hull.bag <- data.frame(x = bag$hull.bag[,1], y = bag$hull.bag[,2])
pxy.outlier <- data.frame(x = bag$pxy.outlier[,1], y = bag$pxy.outlier[,2])

findOutliersIndex <- function(data, outlier) {
  cols <- colnames(outlier)
  dat2 <- data[, cols]
  nr = nrow(outlier)
  rg = 1:nrow(dat2)
  matchingIdx = unique(sort(unlist(sapply(1:nr, function(ii) rg[(outlier[ii,
                                                                         2] == dat2[, 2]) & (outlier[ii, 1] == dat2[, 1])]))))
  if (length(matchingIdx) != nr) {
    stop("cannot locate outliers in the data.")
  }
  matchingIdx
}
nametoshow = rep("", nrow(bagfootRes_dirty_clean))
if (!is.null(pxy.outlier)) {
  outlieridx = findOutliersIndex(bagfootRes_dirty_clean, pxy.outlier)
  nametoshow[outlieridx] = as.character(bagfootRes_dirty_clean$name[outlieridx])
}
bagfootRes_dirty_clean$nametoshow = nametoshow                                         
# Finish the ggplot command

for (i in 1:nrow(bagfootRes_dirty_clean)){
  if (bagfootRes_dirty_clean[i,]$x * bagfootRes_dirty_clean[i,]$y < 0){
    bagfootRes_dirty_clean[i,]$nametoshow = ""
  }
}                     

head(bagfootRes_dirty_clean)

# print selected labels only
names_in_pre = c("EOMES", "TBX21", "Smad2::Smad3", "IRF2", "IRF3", "IRF4", "IRF7", "IRF9", "FOS", "FOS::JUN", "JUNB", "JUND","BATF", "BATF3", 
                 "EGR1", "EGR3", "SP4")
bagfootRes_dirty_clean$nametoshow_pre = bagfootRes_dirty_clean$nametoshow
`%notin%` <- Negate(`%in%`)

for (i in 1:nrow(bagfootRes_dirty_clean)){
  if (bagfootRes_dirty_clean[i, 'nametoshow_pre'] %notin% names_in_pre)
    bagfootRes_dirty_clean[i, 'nametoshow_pre'] = ""
}

library(ggrepel)
library(ggtext)
ggplot(bagfootRes_dirty_clean, aes(x = x,  y = y)) +
  geom_polygon(data = hull.loop, fill = "gray", alpha = 0.35) +
  geom_polygon(data = hull.bag, fill = "gray", alpha = 0.35) +
  geom_point(data = pxy.outlier, col = "black", pch = 16, cex = 1.5) +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        #text = element_text(family="ArialMT", size=25),
        axis.title.x = element_markdown(family="ArialMT")
        ) + 
  geom_text_repel(label=bagfootRes_dirty_clean$nametoshow_pre, max.overlaps = 100, 
                  min.segment.length = unit(0, 'lines'))+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  labs(
    x = "\u0394FA<br><span style='color:#FF9300'>Dirty</span> - <span style='color:#FFFC66'>Clean</span>",
    y = "\u0394FPD"
  )
ggsave(
  "bagfoot/remake_plots/bagfootRes_dirty_clean.pdf",
  plot = last_plot(),
  device = cairo_pdf,
  width = 10,
  height = 10,
  dpi = 300
)





# neo adult

bagfootRes_neo_adult = read.table("bagfoot/res_neo_adult/bagplot_cutcount_diff_total_footprinting_depth_neo-adult_qvalue_bagplot_output.csv", 
                                  header = T,
                                  sep = ",")
colnames(bagfootRes_neo_adult)[3:4] = c("x", "y")
# prove plot should look same as automated generated one
library(aplpack)
bagplot(bagfootRes_neo_adult$x, bagfootRes_neo_adult$y, factor = 1.5,
        show.looppoints = FALSE, show.bagpoints = FALSE,
        show.baghull = TRUE, show.loophull = TRUE,
        show.whiskers = F, col.looppoint = "#000000",
        cex = 2, pch = 0, cex.lab = 1.5, cex.axis = 1.5,
        lwd = 3)

bagfootRes_neo_adult_xy = bagfootRes_neo_adult[,c("x", "y")]
rownames(bagfootRes_neo_adult_xy) = bagfootRes_neo_adult_xy$name


library(aplpack)
bag <- compute.bagplot(bagfootRes_neo_adult_xy, approx.limit = nrow(bagfootRes_neo_adult_xy),
                       factor = 1.5)

hull.loop <- data.frame(x = bag$hull.loop[,1], y = bag$hull.loop[,2])
hull.bag <- data.frame(x = bag$hull.bag[,1], y = bag$hull.bag[,2])
pxy.outlier <- data.frame(x = bag$pxy.outlier[,1], y = bag$pxy.outlier[,2])

findOutliersIndex <- function(data, outlier) {
  cols <- colnames(outlier)
  dat2 <- data[, cols]
  nr = nrow(outlier)
  rg = 1:nrow(dat2)
  matchingIdx = unique(sort(unlist(sapply(1:nr, function(ii) rg[(outlier[ii,
                                                                         2] == dat2[, 2]) & (outlier[ii, 1] == dat2[, 1])]))))
  if (length(matchingIdx) != nr) {
    stop("cannot locate outliers in the data.")
  }
  matchingIdx
}
nametoshow = rep("", nrow(bagfootRes_neo_adult))
if (!is.null(pxy.outlier)) {
  outlieridx = findOutliersIndex(bagfootRes_neo_adult, pxy.outlier)
  nametoshow[outlieridx] = as.character(bagfootRes_neo_adult$name[outlieridx])
}
bagfootRes_neo_adult$nametoshow = nametoshow                                         
# Finish the ggplot command

for (i in 1:nrow(bagfootRes_neo_adult)){
  if (bagfootRes_neo_adult[i,]$x * bagfootRes_neo_adult[i,]$y < 0){
    bagfootRes_neo_adult[i,]$nametoshow = ""
  }
}                     

head(bagfootRes_neo_adult)

# print selected labels only
names_in_pre = c("EOMES", "MGA", "TBX21", "TBX6", "RUNX2", "RUNX3", "BATF", "BATF3", "FOS", "FOS::JUN", "JUNB", "JUND", "BACH1", "BACH2",
                 "EGR1", "EGR3", "SP4")
bagfootRes_neo_adult$nametoshow_pre = bagfootRes_neo_adult$nametoshow
`%notin%` <- Negate(`%in%`)

for (i in 1:nrow(bagfootRes_neo_adult)){
  if (bagfootRes_neo_adult[i, 'nametoshow_pre'] %notin% names_in_pre)
    bagfootRes_neo_adult[i, 'nametoshow_pre'] = ""
}

library(ggrepel)
library(ggtext)
ggplot(bagfootRes_neo_adult, aes(x = x,  y = y)) +
  geom_polygon(data = hull.loop, fill = "gray", alpha = 0.35) +
  geom_polygon(data = hull.bag, fill = "gray", alpha = 0.35) +
  geom_point(data = pxy.outlier, col = "black", pch = 16, cex = 1.5) +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        #text = element_text(family="ArialMT", size=25),
        axis.title.x = element_markdown(family="ArialMT")
  ) + 
  geom_text_repel(label=bagfootRes_neo_adult$nametoshow_pre, max.overlaps = 100, 
                  min.segment.length = unit(0, 'lines'))+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  labs(
    x = "\u0394FA<br><span style='color:#017100'>Neonatal</span> - <span style='color:#88FA4E'>Adult</span>",
    y = "\u0394FPD"
  )
ggsave(
  "bagfoot/remake_plots/bagfootRes_neo_adult.pdf",
  plot = last_plot(),
  device = cairo_pdf,
  width = 10,
  height = 10,
  dpi = 300
)




# VM TN

bagfootRes_vm_tn = read.table("bagfoot/res_vm_tn/bagplot_cutcount_diff_total_footprinting_depth_vm-tn_qvalue_bagplot_output.csv", 
                              header = T,
                              sep = ",")
colnames(bagfootRes_vm_tn)[3:4] = c("x", "y")
# prove plot should look same as automated generated one
library(aplpack)
bagplot(bagfootRes_vm_tn$x, bagfootRes_vm_tn$y, factor = 1.5,
        show.looppoints = FALSE, show.bagpoints = FALSE,
        show.baghull = TRUE, show.loophull = TRUE,
        show.whiskers = F, col.looppoint = "#000000",
        cex = 2, pch = 0, cex.lab = 1.5, cex.axis = 1.5,
        lwd = 3)

bagfootRes_vm_tn_xy = bagfootRes_vm_tn[,c("x", "y")]
rownames(bagfootRes_vm_tn_xy) = bagfootRes_vm_tn_xy$name


library(aplpack)
bag <- compute.bagplot(bagfootRes_vm_tn_xy, approx.limit = nrow(bagfootRes_vm_tn_xy),
                       factor = 1.5)

hull.loop <- data.frame(x = bag$hull.loop[,1], y = bag$hull.loop[,2])
hull.bag <- data.frame(x = bag$hull.bag[,1], y = bag$hull.bag[,2])
pxy.outlier <- data.frame(x = bag$pxy.outlier[,1], y = bag$pxy.outlier[,2])

findOutliersIndex <- function(data, outlier) {
  cols <- colnames(outlier)
  dat2 <- data[, cols]
  nr = nrow(outlier)
  rg = 1:nrow(dat2)
  matchingIdx = unique(sort(unlist(sapply(1:nr, function(ii) rg[(outlier[ii,
                                                                         2] == dat2[, 2]) & (outlier[ii, 1] == dat2[, 1])]))))
  if (length(matchingIdx) != nr) {
    stop("cannot locate outliers in the data.")
  }
  matchingIdx
}
nametoshow = rep("", nrow(bagfootRes_vm_tn))
if (!is.null(pxy.outlier)) {
  outlieridx = findOutliersIndex(bagfootRes_vm_tn, pxy.outlier)
  nametoshow[outlieridx] = as.character(bagfootRes_vm_tn$name[outlieridx])
}
bagfootRes_vm_tn$nametoshow = nametoshow                                         
# Finish the ggplot command

for (i in 1:nrow(bagfootRes_vm_tn)){
  if (bagfootRes_vm_tn[i,]$x * bagfootRes_vm_tn[i,]$y < 0){
    bagfootRes_vm_tn[i,]$nametoshow = ""
  }
}                     

head(bagfootRes_vm_tn)

# print selected labels only
names_in_pre = c("EOMES", "MGA", "TBX21", "TBX20", "TBX6", "RUNX2", "RUNX3", "BATF", "BATF3", "FOS", "FOS::JUN", "JUNB", "JUND", "BACH1", "BACH2",
                 "EGR1", "EGR3", "SP4", "SP3", "SP1", 'KLF2', 'Klf12')
bagfootRes_vm_tn$nametoshow_pre = bagfootRes_vm_tn$nametoshow
`%notin%` <- Negate(`%in%`)

for (i in 1:nrow(bagfootRes_vm_tn)){
  if (bagfootRes_vm_tn[i, 'nametoshow_pre'] %notin% names_in_pre)
    bagfootRes_vm_tn[i, 'nametoshow_pre'] = ""
}

library(ggrepel)
library(ggtext)
ggplot(bagfootRes_vm_tn, aes(x = x,  y = y)) +
  geom_polygon(data = hull.loop, fill = "gray", alpha = 0.35) +
  geom_polygon(data = hull.bag, fill = "gray", alpha = 0.35) +
  geom_point(data = pxy.outlier, col = "black", pch = 16, cex = 1.5) +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        #text = element_text(family="ArialMT", size=25),
        axis.title.x = element_markdown(family="ArialMT")
  ) + 
  geom_text_repel(label=bagfootRes_vm_tn$nametoshow_pre, max.overlaps = 100, 
                  min.segment.length = unit(0, 'lines'))+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  labs(
    x = "\u0394FA<br><span style='color:#004D7F'>VM</span> - <span style='color:#56C1FF'>TN</span>",
    y = "\u0394FPD"
  )
ggsave(
  "bagfoot/remake_plots/bagfootRes_vm_tn.pdf",
  plot = last_plot(),
  device = cairo_pdf,
  width = 10,
  height = 10,
  dpi = 300
)
