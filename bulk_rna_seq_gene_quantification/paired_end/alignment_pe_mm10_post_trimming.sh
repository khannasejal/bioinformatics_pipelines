#/bin/bash
srr_file="SRR_Acc_List_original.txt"
raw_dir="/mnt/e/tanmay_lele_19thmay/fastq" #has the unzipped raw fastq files
hisat2_index="/mnt/e/datasets/mm10_hisat2_index/mm10"
trim_dir="/mnt/e/tanmay_lele_19thmay/trimmed"
num_t=10
bam_dir="/mnt/e/tanmay_lele_19thmay/bam_bai_folder"
gtf_dir="/mnt/e/datasets/gencode.vM10.annotation.gtf"
fc_dir="/mnt/e/tanmay_lele_19thmay/featurecounts"
srr=SRR11569251
#while read -r srr; do
echo "Processing $srr"
fasterq-dump "$srr" -e "$num_t" -O "$raw_dir"
fastp -i "$raw_dir/${srr}_1.fastq" -I "$raw_dir/${srr}_2.fastq" -o "$trim_dir/${srr}_1_trim.fastq" -O "$trim_dir/${srr}_2_trim.fastq" -w $num_t
gzip "$trim_dir/${srr}_1_trim.fastq"
gzip "$trim_dir/${srr}_2_trim.fastq"
rm "$raw_dir/${srr}_1.fastq" "$raw_dir/${srr}_2.fastq"  #remove the raw fastq files
echo "Trimming done and files zipped! $srr" 
echo "$trim_dir/${srr}_1_trim.fastq.gz"
hisat2 -x "$hisat2_index" -1 "$trim_dir/${srr}_1_trim.fastq.gz" -2 "$trim_dir/${srr}_2_trim.fastq.gz" -p "$num_t" -S "$bam_dir/${srr}.sam" 
echo "Alignment done! $srr"
samtools view -@ "$num_t" -bS "$bam_dir/${srr}.sam" | samtools sort -n -@ "$num_t" -o "$bam_dir/${srr}_nsorted.bam" #indexing not required for name sorted bam files
echo "sam to bam conversion and name sorting done! $srr"
rm "$bam_dir/${srr}.sam"
featureCounts -a "$gtf_dir" -o "$fc_dir/${srr}_counts.txt" -p --countReadPairs -s 0 -C -T "$num_t" "$bam_dir/${srr}_nsorted.bam" #the library is unstranded; for paired-end include -p, -CountReadPairs and -C to exclude the chimeric reads; both multimapping and multi-overlapping reads are excluded
echo "reads mapped to gene features using featureCounts! $srr"
#done < SRR_Acc_List_original.txt
