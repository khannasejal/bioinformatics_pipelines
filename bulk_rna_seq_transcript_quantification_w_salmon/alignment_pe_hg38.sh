#/bin/bash

raw_dir="/mnt/e/tanmay_lele_19thmay/RNAseq/mdamb231/fastq" #has the unzipped raw fastq files
trim_dir="/mnt/e/tanmay_lele_19thmay/RNAseq/mdamb231/trimmed"
salmon_index="/mnt/e/tanmay_lele_19thmay/RNAseq/mdamb231/hg38_refseq_fasta/rna/hg38_full_salmon_rna"
counts_dir="/mnt/e/tanmay_lele_19thmay/RNAseq/mdamb231/quants"
transcript_gene_map="/mnt/e/tanmay_lele_19thmay/RNAseq/mdamb231/hg38_refseq_fasta/transcript_gene_map_salmon.txt"
num_t=10

for r1 in "$raw_dir"/*_R1_001.fastq.gz; do
    baseprefix=$(basename "$r1" | sed 's/_R1_001.fastq.gz//')
    r2="$raw_dir/${baseprefix}_R2_001.fastq.gz" #change this to trim directory if needs to be run again
    echo "Processing $r1 and $r2"
    fastp -i "$r1" -I "$r2" -w $num_t -o "$trim_dir/${baseprefix}_R1_trim.fastq.gz" -O "$trim_dir/${baseprefix}_R2_trim.fastq.gz" 
    echo "Trimming done for $r1 and $r2"
    #rm "$raw_fastq_dir/${baseprefix}_R1_001.fastq.gz" "$raw_fastq_dir/${baseprefix}_R2_001.fastq.gz"  #remove the raw fastq files  
    salmon quant -i "$salmon_index" -l A -1 "$trim_dir/${baseprefix}_R1_trim.fastq.gz" -2 "$trim_dir/${baseprefix}_R2_trim.fastq.gz" -p $num_t -g "$transcript_gene_map" -o "$counts_dir/${baseprefix}_counts"
    echo "Salmon quantification done $r1 and $r2"
done
