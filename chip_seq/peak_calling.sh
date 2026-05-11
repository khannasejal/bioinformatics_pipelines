#/bin/bash
search_dir=/mnt/d/sejal/bam_bai_folder/mcf7
for file in "$search_dir"/*.bam;do
 filename=$(basename $file)
 echo $filename
 macs2 callpeak -t ${file} -c ${search_dir}/input/*.bam -f BAM -g hs -n ${filename}_narrow_peaks -B -q 0.05 --outdir ${search_dir}/narrow
 macs2 callpeak -t ${file} -c ${search_dir}/input/*.bam --broad -f BAM -g hs -n ${filename}_broad_peaks -B -q 0.05 --outdir ${search_dir}/broad
done
  
