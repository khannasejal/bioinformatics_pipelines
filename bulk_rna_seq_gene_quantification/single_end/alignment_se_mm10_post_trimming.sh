#/bin/bash

raw_dir=15-010
hisat2_index=mm10_hisat2_index/mm10
trim_dir=15-010/trimmed
num_t=10
bam_dir=15-010/bam_bai_folder
gtf_dir=gencode.vM10.annotation.gtf
fc_dir=15-010/featurecounts

for r1 in "$raw_dir"/*_R1_001.fastq.gz; do
    baseprefix=$(basename "$r1" | sed 's/_R1_001.fastq.gz//')
    echo "Processing $r1"

    fastp -i "$r1" -o "$trim_dir/${baseprefix}_trim.fastq.gz" -w $num_t
    echo "Trimming done $r1"

    hisat2 -x "$hisat2_index" -U "${r1}" -p "$num_t" -S "$bam_dir/${baseprefix}.sam" 
    echo "Alignment done: $baseprefix"

    samtools view -@ "$num_t" -bS "$bam_dir/${baseprefix}.sam" | samtools sort -n -@ "$num_t" -o "$bam_dir/${baseprefix}_nsorted.bam" #indexing not required for name sorted bam files
    echo "sam to bam conversion and name sorting done: $baseprefix"
    rm "$bam_dir/${baseprefix}.sam"

    featureCounts -a "$gtf_dir" -o "$fc_dir/${baseprefix}_counts.txt" -s 0 -C -T "$num_t" "$bam_dir/${baseprefix}_nsorted.bam" #the library is unstranded; for paired-end include -p, -CountReadPairs and -C to exclude the chimeric reads; both multimapping and multi-overlapping reads are excluded
    echo "reads mapped to gene features using featureCounts: $baseprefix"
done