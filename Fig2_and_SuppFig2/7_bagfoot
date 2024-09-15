# run bagfoot analysis

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

