#/bin/bash
#the script is for trimming the paired-end reads; include fastqc as well for next run
search_dir=18-159
output_dir=18-159/trimmed
num_t=10

for file in "$search_dir"/*_R1.fastq.gz; do
    fullfilename=$(basename "$file")
    baseprefix=$(echo "$fullfilename" | sed 's/_R1\.fastq\.gz//')
    r1="$search_dir/${baseprefix}_R1.fastq.gz"
    r2="$search_dir/${baseprefix}_R2.fastq.gz"
    echo "Processing $r1 and $r2"
    fastp -i "$r1" -I "$r2" -o "$output_dir/${baseprefix}_R1_trim.fastq.gz" -O "$output_dir/${baseprefix}_R2_trim.fastq.gz" -w $num_t
    echo "Finished: $baseprefix"
done