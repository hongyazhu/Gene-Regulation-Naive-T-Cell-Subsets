# Example_pipeline_processing_RNAseq

# qc
fastqc -t 10 -f fastq -o qc fastqs/*.fastq.gz

# cut adapters
trim_galore --fastqc_args "--outdir qc_cutadapt" --output_dir trim_galore --cores 2 fastqs/*.fastq.gz

# mapping
cd trim_galore
for f in `ls -1 *_trimmed.fq.gz | sed 's/_trimmed.fq.gz//' `
do
hisat2 -x /path/to/mouse/mm10/genome -U ${f}_trimmed.fq.gz -S ../hisat2/${f}.bam -p 8 2> ../hisat2/${f}_summary.txt
done

# counting reads
cd ..
featureCounts -a /path/to/mouse/gencode_annotation/gencode.vM21.annotation.ids.saf -F 'SAF' -s 0 -Q 50 -o counts.txt hisat2/*.bam > featureCounts.log

### counting reads with gene names 
featureCounts -a /path/to/mouse/gencode_annotation/gencode.vM21.annotation.names.saf -F 'SAF' -s 0 -Q 50 -o counts_name.txt hisat2/*.bam 2> featureCounts_name.log
