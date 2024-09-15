# to identify differential binding of Eomes between WT and Lin28b cells (which is neonatal like)

macs3 callpeak -B -t ./dedup-BAMS/WT*EOMES*.DEDUP.bam -c ./dedup-BAMS/WT*IgG*.DEDUP.bam -n WT_EOMES -f BAMPE -g mm -q 0.01 --fe-cutoff 10 --nolambda --outdir ../differential/MACS3_bdgdiff &>> ../differential/EOMES_MACS_bdgdiff.log
macs3 callpeak -B -t ./dedup-BAMS/lin28Tg*EOMES*.DEDUP.bam -c ./dedup-BAMS/lin28Tg*IgG*.DEDUP.bam -n lin28Tg_EOMES -f BAMPE -g mm -q 0.01 --fe-cutoff 10 --nolambda --outdir ../differential/MACS3_bdgdiff &>> ../differential/EOMES_MACS_bdgdiff.log

cd ../differential/MACS3_bdgdiff

# get d1 and d2 for bdgdiff:
egrep "fragments after filtering in treatment|fragments after filtering in control" WT_EOMES_peaks.xls
# fragments after filtering in treatment: 4578659
# fragments after filtering in control: 3141374
egrep "fragments after filtering in treatment|fragments after filtering in control" lin28Tg_EOMES_peaks.xls
# fragments after filtering in treatment: 2150122
# fragments after filtering in control: 1929811

macs3 bdgdiff --t1 WT_EOMES_treat_pileup.bdg --c1 WT_EOMES_control_lambda.bdg --t2 lin28Tg_EOMES_treat_pileup.bdg --c2 lin28Tg_EOMES_control_lambda.bdg --d1 3141374 --d2 1929811 --o-prefix EOMES_bdgdiff_WT_vs_lin28Tg --outdir .  &>> ../EOMES_MACS_bdgdiff.log

# add chr before chromosome
find -type f \( -name "*.bed" -o -name "*.narrowPeak" -o -name "*.bdg" \) -exec sed -i 's/^/chr/' {} \;
find -type f -name "*.bed" -exec sed -i 's/chrtrack/#track/' {} \;

# annotate to closest genes
for i in $(ls EOMES_bdgdiff*.bed | sed 's/.bed//'); do
annotatePeaks.pl ${i}.bed mm10 > ${i}_peak_anno.txt 2> ${i}_peak_anno.log
done
