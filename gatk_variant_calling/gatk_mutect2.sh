#!/bin/bash

# ==========================================
# 1. CONFIGURATION
# ==========================================
JAVA_CMD="java -Xmx16g -jar /usr/share/java/gatk-bin/GenomeAnalysisTK.jar" 
REF="/home/sejal/sejal/tanmay_lele/hg38/Homo_sapiens.GRCh38.dna.primary_assembly_modified.fa"
GERMLINE="/home/sejal/sejal/tanmay_lele/known_sites/af-only-gnomad.hg38.vcf.gz"
PON="/home/sejal/sejal/tanmay_lele/known_sites/1000g_pon.hg38.vcf.gz"
GATK_OUT="/home/sejal/sejal/tanmay_lele/gatk/st_sel_st_anc"

BAM_DIR="/home/sejal/sejal/tanmay_lele/gatk_preprocessed"
OUTPUT_NAME="${GATK_OUT}/st_sel_st_anc"

# ==========================================
# 2. DEFINE YOUR FILES MANUALLY
# ==========================================
# List your 4 TEST files (Condition A - The "Tumor")
TEST_1="${BAM_DIR}/SRR27970197.processed.bam"
TEST_2="${BAM_DIR}/SRR27970198.processed.bam"
TEST_3="${BAM_DIR}/SRR27970199.processed.bam"
TEST_4="${BAM_DIR}/SRR27970200.processed.bam"

# List your 4 CONTROL files (Condition B - The "Normal")
CTRL_1="${BAM_DIR}/SRR27970205.processed.bam"
CTRL_2="${BAM_DIR}/SRR27970206.processed.bam"
CTRL_3="${BAM_DIR}/SRR27970207.processed.bam"
CTRL_4="${BAM_DIR}/SRR27970208.processed.bam"

# ==========================================
# 3. GET SAMPLE NAMES (CRITICAL)
# ==========================================
# Mutect2 needs the internal "SM" name to know which files are Normal.
# We extract them automatically here.

echo "Extracting Sample Names from Control BAMs..."
N1_NAME=$(samtools view -H $CTRL_1 | grep "^@RG" | sed "s/.*SM:\([^\t]*\).*/\1/" | head -n 1)
N2_NAME=$(samtools view -H $CTRL_2 | grep "^@RG" | sed "s/.*SM:\([^\t]*\).*/\1/" | head -n 1)
N3_NAME=$(samtools view -H $CTRL_3 | grep "^@RG" | sed "s/.*SM:\([^\t]*\).*/\1/" | head -n 1)
N4_NAME=$(samtools view -H $CTRL_4 | grep "^@RG" | sed "s/.*SM:\([^\t]*\).*/\1/" | head -n 1)

echo "Normal Sample Names found: $N1_NAME, $N2_NAME, $N3_NAME, $N4_NAME"

# ==========================================
# 4. RUN MUTECT2 (POOLED NORMAL MODE)
# ==========================================

echo ">>> Starting Multi-Sample Mutect2..."

# We input ALL files (-I), but we specifically flag the control sample names with -normal
$JAVA_CMD Mutect2 -R $REF -I $TEST_1 -I $TEST_2 -I $TEST_3 -I $TEST_4 -I $CTRL_1 -I $CTRL_2 -I $CTRL_3 -I $CTRL_4 -normal "$N1_NAME" -normal "$N2_NAME" -normal "$N3_NAME" -normal "$N4_NAME" --germline-resource $GERMLINE --panel-of-normals $PON --native-pair-hmm-threads 16 --f1r2-tar-gz ${OUTPUT_NAME}.f1r2.tar.gz -O ${OUTPUT_NAME}.unfiltered.vcf.gz --independent-mates --dont-use-soft-clipped-bases

# ==========================================
# 5. POST-PROCESSING (STANDARD STEPS)
# ==========================================

echo ">>> Step 5: LearnReadOrientationModel..."
$JAVA_CMD LearnReadOrientationModel -I ${OUTPUT_NAME}.f1r2.tar.gz -O ${OUTPUT_NAME}.read-orientation-model.tar.gz

echo ">>> Step 6: FilterMutectCalls..."
$JAVA_CMD FilterMutectCalls -R $REF -V ${OUTPUT_NAME}.unfiltered.vcf.gz --ob-priors ${OUTPUT_NAME}.read-orientation-model.tar.gz -O ${OUTPUT_NAME}.filtered.vcf.gz

#echo ">>> Step 7: VariantFiltration..." #this step is not required; it is for haplotype 
#$JAVA_CMD VariantFiltration -R $REF -V ${OUTPUT_NAME}.filtered.vcf.gz -O ${OUTPUT_NAME}.final.vcf.gz --window 35 --cluster 3 --filter-name "FS" --filter-expression "FS > 30.0" --filter-name "QD" --filter-expression "QD < 2.0"

echo ">>> Step 7: VariantsToTable..."
# Note: This will create a VERY wide table with columns for all 8 samples; convert the filtermutect vcf output directly into a tsv file
#$JAVA_CMD VariantsToTable -R $REF -V ${OUTPUT_NAME}.final.vcf.gz -O ${OUTPUT_NAME}.all_variants.tsv -F CHROM -F POS -F REF -F ALT -F FILTER -F TYPE -GF GT -GF AD -GF AF -GF DP

$JAVA_CMD VariantsToTable -R $REF -V ${OUTPUT_NAME}.filtered.vcf.gz -O ${OUTPUT_NAME}.all_variants.tsv -F CHROM -F POS -F REF -F ALT -F FILTER -F TLOD -F NLOD -F DP -GF GT -GF AD -GF AF -GF DP

echo "Done! Check ${OUTPUT_NAME}.all_variants.tsv"
