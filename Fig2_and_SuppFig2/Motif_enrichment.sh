# In preparation for motif emrichment, convert differential accessible peaks results to bed
cd DE/res/
for i in $(ls *.txt| sed 's/.txt//'); do
awk -F'\t' '{print $1}' ${i}.txt | awk -F'\\:' '$1=$1' OFS="\t" | awk -F'\\-' '$1=$1' OFS="\t" | paste --delimiters='\t' - ${i}.txt | awk '{$6="." ; print ;}'| sed 's/ /\t/g' | cut -f 1-6 > bed/${i}.bed
done


### motif enrichment with AME

cd bed
for i in $(ls | sed 's/.bed//'); do
bedtools getfasta -fi mm10_no_alt_analysis_set_ENCODE.fasta -bed ${i}.bed -name -fo fastas/${i}.fasta
done

cd DE/motif_enrichment/fastas/
fasta-subsample insig_neo_adult.fasta 4982 > subsample/sample_insig_neo_adult.fasta
fasta-subsample insig_dirty_clean.fasta 3973 > subsample/sample_insig_dirty_clean.fasta
fasta-subsample insig_vm_tn.fasta 12063 > subsample/sample_insig_vm_tn.fasta
# these numbers are the mean number of peaks between the two conditions

cd DE/motif_enrichment/fastas/
ame --verbose 1 --o ../neo --scoring avg --method fisher --hit-lo-fraction 0.25 --evalue-report-threshold 10.0 --control subsample/sample_insig_neo_adult.fasta neo.fasta JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt

ame --verbose 1 --o ../adult --scoring avg --method fisher --hit-lo-fraction 0.25 --evalue-report-threshold 10.0 --control subsample/sample_insig_neo_adult.fasta adult.fasta JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt

ame --verbose 1 --o ../clean --scoring avg --method fisher --hit-lo-fraction 0.25 --evalue-report-threshold 10.0 --control subsample/sample_insig_dirty_clean.fasta clean.fasta JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt

ame --verbose 1 --o ../dirty --scoring avg --method fisher --hit-lo-fraction 0.25 --evalue-report-threshold 10.0 --control subsample/sample_insig_dirty_clean.fasta dirty.fasta JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt

ame --verbose 1 --o ../tn --scoring avg --method fisher --hit-lo-fraction 0.25 --evalue-report-threshold 10.0 --control subsample/sample_insig_vm_tn.fasta tn.fasta JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt

ame --verbose 1 --o ../vm --scoring avg --method fisher --hit-lo-fraction 0.25 --evalue-report-threshold 10.0 --control subsample/sample_insig_vm_tn.fasta vm.fasta JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt
