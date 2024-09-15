cd DE/res/bed/
for i in $(ls | sed 's/.bed//'); do
annotatePeaks.pl ${i}.bed mm10 > peak_annotation/${i}_peak_anno.txt
done

# annoate all peaks to provide a background list for enrichment analysis
awk 'BEGIN {FS="\t"; OFS="\t"} {print $2, $3, $4, $1, 0, $5}' counts_ATACseq | tail -n +2 > counts_ATACseq.bed
annotatePeaks.pl counts_ATACseq.bed mm10 > peak_annotation/merged_peak_anno.txt
