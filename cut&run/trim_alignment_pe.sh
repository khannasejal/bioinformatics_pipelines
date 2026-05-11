#!/bin/bash
search_dir="/mnt/e/datasets/19-024-and-19-055"
alignDir="/mnt/e/datasets/19-024-and-19-055/bam_bai"
bw_folder="/mnt/e/datasets/19-024-and-19-055/bw"
genome="/mnt/e/datasets/bowtie2_index/mm10_bowtie2_index/mm10"
spikein="/mnt/e/datasets/bowtie2_index/yeast_bowtie2_index/R64-1-1/R64-1-1"
ecoli="/mnt/e/datasets/bowtie2_index/ecoli_bowtie2_index/genome"

num_t=10
scale=1000

echo -e "Sample\tEcoli_Reads\tSpikeIn_Reads\tMouse_Reads\tScaleFactor" > $alignDir/summary_report.txt

for file in "$search_dir"/*_1.fq.gz; do
    fullfilename=$(basename "$file")
    base=$(echo "$fullfilename" | sed 's/_1\.fq\.gz//')
    
    r1="$search_dir/${base}_1.fq.gz"
    r2="$search_dir/${base}_2.fq.gz"
    
    echo "--------------------------------------"
    echo "Processing Sample: $base"
    echo "--------------------------------------"

    # A. TRIMMING
    fastp -i "$r1" -I "$r2" \
          -o "${base}_1_trim.fq.gz" -O "${base}_2_trim.fq.gz" \
          --html "${base}_fastp.html" -w $num_t

    # B. ALIGN TO MOUSE (Experimental)
    echo "Aligning $base to Mouse..."
    bowtie2 -p $num_t --local --very-sensitive-local --no-unal --no-mixed --no-discordant \
        -I 10 -X 700 -x $genome -1 "${base}_1_trim.fq.gz" -2 "${base}_2_trim.fq.gz" \
        2> ${alignDir}/${base}_mouse_bt2.log | \
        samtools view -hb -f 2 - | samtools sort -o ${alignDir}/${base}_mouse.bam -
    samtools index ${alignDir}/${base}_mouse.bam

    # C. ALIGN TO YEAST (Spike-in)
    echo "Aligning $base to Yeast..."
    bowtie2 -p $num_t --local --very-sensitive-local --no-unal --no-mixed --no-discordant \
        -x $spikein -1 "${base}_1_trim.fq.gz" -2 "${base}_2_trim.fq.gz" \
        2> ${alignDir}/${base}_spike_bt2.log | \
        samtools view -hb -f 2 - > ${alignDir}/${base}_spike.bam

    # D. ALIGN TO ECOLI (Carry-over)
    echo "Aligning $base to Ecoli..."
    bowtie2 -p $num_t --local --very-sensitive-local --no-unal --no-mixed --no-discordant \
        -x $ecoli -1 "${base}_1_trim.fq.gz" -2 "${base}_2_trim.fq.gz" \
        2> ${alignDir}/${base}_ecoli_bt2.log | \
        samtools view -hb -f 2 - > ${alignDir}/${base}_ecoli.bam

    # --- 3. SPIKE-IN CALIBRATION ---
    # Get read counts (using samtools view -c is faster than wc -l on BED)
    n_mouse=$(samtools view -c ${alignDir}/${base}_mouse.bam)
    n_spike=$(samtools view -c ${alignDir}/${base}_spike.bam)
    n_ecoli=$(samtools view -c ${alignDir}/${base}_ecoli.bam)

    # Calculate scaling factor: S = Scale / N_spike
    # Note: If n_spike is 0, we avoid division by zero
    if [ "$n_spike" -gt 0 ]; then
        # Using awk for floating point math
        scale_factor=$(awk -v s=$scale -v n=$n_spike 'BEGIN {print s/n}')
    else
        scale_factor=1
        echo "WARNING: No spike-in reads found for $base. Using factor 1."
    fi

    echo "Scaling Factor for $base: $scale_factor"
    echo -e "$base\t$n_ecoli\t$n_spike\t$n_mouse\t$scale_factor" >> $alignDir/summary_report.txt

    # --- 4. GENERATE BIGWIGS ---
    # 1. Spike-in Normalized BigWig
    bamCoverage -b ${alignDir}/${base}_mouse.bam \
        -o ${bw_folder}/${base}_spikeNorm.bw \
        --binSize 10 --extendReads \
        --scaleFactor $scale_factor

    # 2. Standard CPM Normalized BigWig (for comparison)
    bamCoverage -b ${alignDir}/${base}_mouse.bam \
        -o ${bw_folder}/${base}_cpmNorm.bw \
        --binSize 10 --extendReads \
        --normalizeUsing CPM

    # Cleanup temp trimmed files to save space
    #rm "${base}_1_trim.fq.gz" "${base}_2_trim.fq.gz"

done