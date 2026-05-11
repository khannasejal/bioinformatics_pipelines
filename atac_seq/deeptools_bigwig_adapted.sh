#/bin/bash
#script is for generating bigwig normalisation files 
bam_dir="/mnt/e/tanmay_lele_19thmay/ATAC/Adapted/bam/b2/removed_bg"
bw_dir="/mnt/e/tanmay_lele_19thmay/ATAC/Adapted/bw/cpm_bg_removed"

num_t=10
for r1 in "$bam_dir"/*.blacklist-filtered_wo_bg.bam; do
    baseprefix=$(basename "$r1" | sed 's/.blacklist-filtered_wo_bg.bam//')
    echo "Processing $r1"
    samtools sort -@ $num_t "$bam_dir/${baseprefix}.blacklist-filtered_wo_bg.bam" -o "$bam_dir/${baseprefix}.blacklist-filtered_wo_bg.sorted.bam"
    samtools index -@ $num_t "$bam_dir/${baseprefix}.blacklist-filtered_wo_bg.sorted.bam" "$bam_dir/${baseprefix}.blacklist-filtered_wo_bg.sorted.bam.bai"
    bamCoverage --bam "$bam_dir/${baseprefix}.blacklist-filtered_wo_bg.sorted.bam" -o "$bw_dir/${baseprefix}.blacklist-filtered.bw" --numberOfProcessors $num_t --binSize 50 --normalizeUsing CPM --smoothLength 200 --effectiveGenomeSize 2913022398
    echo "Big-wig coverage files ready!"
done


