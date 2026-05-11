#!/bin/bash
ulimit -n 10000
num_t=30
fastq_dir="/home/sejal/sejal/tanmay_lele/fastq/done"
bam_bai_dir="/home/sejal/sejal/tanmay_lele/bam_bai"
picard_out_dir="/home/sejal/sejal/tanmay_lele/picard_out"
star_ref="/mnt/RNASeq/hg38_STARIndexed"
for fwd in "$fastq_dir"/*_1_trim.fastq.gz; do
 rev=${fwd/_1_trim/_2_trim}
 filename=$(basename "$fwd")
 bn=$(echo "$filename" | sed 's/_1_trim.*//')
 echo "PROCESSING SAMPLE: $bn"
 STAR --runThreadN $num_t --genomeDir $star_ref --readFilesIn "$fwd" "$rev" --readFilesCommand zcat --outFileNamePrefix "$bam_bai_dir/$bn." --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattrRGline ID:$bn CN:Lab LB:PairedEnd PL:Illumina PU:Unknown SM:$bn --twopassMode Basic
 mv "$bam_bai_dir/$bn.Aligned.sortedByCoord.out.bam" "$bam_bai_dir/$bn.sorted.bam"
 echo "Running Picard MarkDuplicates..."
 picard-java MarkDuplicates I="$bam_bai_dir/$bn.sorted.bam" O="$picard_out_dir/$bn.dupMarked.bam" M="$picard_out_dir/$bn.dup.metrics" CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT
 echo "Completed $bn."
done
