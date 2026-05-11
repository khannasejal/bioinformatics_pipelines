#!/bin/bash

INPUT_DIR="/home/sejal/sejal/tanmay_lele/picard_out/undone"            
OUTPUT_DIR="/home/sejal/sejal/tanmay_lele/gatk_preprocessed"
JAVA_CMD="java -Xmx16g -jar /usr/share/java/gatk-bin/GenomeAnalysisTK.jar"

REF="/home/sejal/sejal/tanmay_lele/hg38/Homo_sapiens.GRCh38.dna.primary_assembly_modified.fa"
DBSNP="/home/sejal/sejal/tanmay_lele/known_sites/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
MILLS="/home/sejal/sejal/tanmay_lele/known_sites/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
KGINDELS="/home/sejal/sejal/tanmay_lele/known_sites/1000G_phase1.snps.high_confidence.hg38.vcf.gz"


for bam_file in "$INPUT_DIR"/*.dupMarked.bam; do
 filename=$(basename "$bam_file") 
 BN="${filename%%.dupMarked.bam}"

 echo "========================================"
 echo "Processing Sample: $BN"
 echo "========================================"


 #echo ">>> Step 1: SplitNCigarReads..."
 #$JAVA_CMD SplitNCigarReads -R $REF -I "$bam_file" -O "$OUTPUT_DIR/${BN}.split.bam" 

 echo ">>> Step 2: BaseRecalibrator..."
 $JAVA_CMD BaseRecalibrator -R $REF -I "$OUTPUT_DIR/${BN}.split.bam" --known-sites $KGINDELS --known-sites $MILLS --known-sites $DBSNP -O "$OUTPUT_DIR/${BN}.recal_data.table"

 echo ">>> Step 3: ApplyBQSR..."
 $JAVA_CMD ApplyBQSR -R $REF -I "$OUTPUT_DIR/${BN}.split.bam" --bqsr-recal-file "$OUTPUT_DIR/${BN}.recal_data.table" -O "$OUTPUT_DIR/${BN}.processed.bam"
    
 echo ">>> Finished pre-processing $BN"
done

