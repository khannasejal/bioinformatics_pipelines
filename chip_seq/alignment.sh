#/bin/bash
mkdir -p bam_bai_folder/
mkdir -p bw_folder_GSE32222_new
search_dir=fastq_folder
num_t=10
for file in "$search_dir"/*;do
 filename=$(basename $file)
 echo $filename
 bowtie2 -t --local -N 1 -k 1 --threads $num_t -x GRCh38_noalt_as/GRCh38_noalt_as -U ${file} -S ${filename}.sam
 echo "Alignment done!"
 samtools view ${filename}.sam -q 30 -bS -@ $num_t > ${filename}.bam
 echo "Conversion to BAM from SAM done!"
 rm ${filename}.sam
 samtools view -b -q 30 -@ $num_t ${filename}.bam > ${filename}.filt.bam
 echo "Filtering BAM file done!"
 rm ${filename}.bam
 samtools sort ${filename}.filt.bam -@ $num_t -o bam_bai_folder/${filename}.sorted.bam
 echo "Conversion to sorted BAM done!"
 rm ${filename}.filt.bam
 samtools index bam_bai_folder/${filename}.sorted.bam bam_bai_folder/${filename}.sorted.bam.bai -@ $num_t
 echo "Indexing of BAM done!"
 bamCoverage --bam bam_bai_folder_GSE32222/bam_bai_folder/${filename}.sorted.bam -o bw_folder_GSE32222_new/${filename}.sorted.CPM.smooth.bw --binSize 50 --normalizeUsing CPM --effectiveGenomeSize 2913022398 --smoothLength 200 -p $num_t
 echo "Coverage plot calculated!"
done    
