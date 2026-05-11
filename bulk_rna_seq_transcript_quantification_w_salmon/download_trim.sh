#/bin/bash
srr_file="SRR_Acc_List_original.txt"
raw_dir="/mnt/e/tanmay_lele_19thmay/RNAseq/mdamb231/fastq" #has the unzipped raw fastq files
trim_dir="/mnt/e/tanmay_lele_19thmay/RNAseq/mdamb231/trimmed"
qc_dir="/mnt/e/tanmay_lele_19thmay/RNAseq/mdamb231/fastqc_reports"
prefetch_dir="/mnt/e/tanmay_lele_19thmay/RNAseq/mdamb231/sra"
num_t=10
while read -r srr; do
    echo "Processing $srr"
    sra_path="$prefetch_dir/$srr/$srr.sra"
    fasterq-dump "$sra_path" -e "$num_t" -O "$raw_dir"
    fastp -i "$raw_dir/${srr}_1.fastq" -I "$raw_dir/${srr}_2.fastq" -o "$trim_dir/${srr}_1_trim.fastq" -O "$trim_dir/${srr}_2_trim.fastq" -w "$num_t"
    gzip "$trim_dir/${srr}_1_trim.fastq"
    gzip "$trim_dir/${srr}_2_trim.fastq"
    rm "$raw_dir/${srr}_1.fastq" "$raw_dir/${srr}_2.fastq"
    echo "Trimming done and files zipped! $srr"
    fastqc "$trim_dir/${srr}_1_trim.fastq.gz" -o "$qc_dir" -t "$num_t" --noextract
    fastqc "$trim_dir/${srr}_2_trim.fastq.gz" -o "$qc_dir" -t "$num_t" --noextract
    echo "Trimmed FastQC reports generated for $srr."
done < "$srr_file"