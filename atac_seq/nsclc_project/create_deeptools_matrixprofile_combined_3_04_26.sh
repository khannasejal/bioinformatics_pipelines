#!/bin/bash

# --- 1. USER VARIABLES ---
# Path to your BED file - Change this variable to run different signatures
BED_FILE="/mnt/e/shreya_nsclc/GSE193258/bw_narrowpeak/bw/osi_dtp/genesets_bed/GRCh38_emt_hallmark_genes.bed" 
# Label for filenames and plot legend (e.g., emt_hallmark, luad, lusc)
BED_LABEL="emt_hallmark"   

# Output directory named by date
OUT_DIR="./results_$(date +%d_%m_%y)"

mkdir -p "$OUT_DIR"
computeMatrix reference-point -p 10 \
        -S GSM5777956_X1a_3562AZ_HCC827_CTRL-ready.bw GSM5777957_X1b_3562AZ_HCC827_CTRL-ready.bw GSM5777958_X1c_3562AZ_HCC827_CTRL-ready.bw GSM5777959_X2a_3562AZ_HCC827_DTP-ready.bw GSM5777960_X2b_3562AZ_HCC827_DTP-ready.bw GSM5777961_X2c_3562AZ_HCC827_DTP-ready.bw \
        -R "$BED_FILE" \
        --referencePoint TSS \
        -b 3000 -a 3000 \
        --skipZeros \
        -o "${OUT_DIR}/HCC827_${BED_LABEL}_matrix.gz"
plotProfile -m "${OUT_DIR}/HCC827_${BED_LABEL}_matrix.gz" -out "${OUT_DIR}/HCC827_${BED_LABEL}_profile.png" --averageType mean --colors blue blue blue red red red

computeMatrix reference-point -p 10 \
        -S GSM5777965_X4a_3562AZ_HCC2935_CTRL-ready.bw GSM5777966_X4b_3562AZ_HCC2935_CTRL-ready.bw GSM5777967_X4c_3562AZ_HCC2935_CTRL-ready.bw GSM5777968_X5a_3562AZ_HCC2935_DTP-ready.bw GSM5777969_X5b_3562AZ_HCC2935_DTP-ready.bw GSM5777970_X5c_3562AZ_HCC2935_DTP-ready.bw \
        -R "$BED_FILE" \
        --referencePoint TSS \
        -b 3000 -a 3000 \
        --skipZeros \
        -o "${OUT_DIR}/HCC2935_${BED_LABEL}_matrix.gz"
plotProfile -m "${OUT_DIR}/HCC2935_${BED_LABEL}_matrix.gz" -out "${OUT_DIR}/HCC2935_${BED_LABEL}_profile.png" --averageType mean --colors blue blue blue red red red

computeMatrix reference-point -p 10 \
        -S GSM5777974_X1a_3519AZ_H1975_DMSO-ready.bw GSM5777975_X1b_3519AZ_H1975_DMSO-ready.bw GSM5777976_X1c_3519AZ_H1975_DMSO-ready.bw GSM5777977_X2a_3519AZ_H1975_9291-ready.bw GSM5777978_X2b_3519AZ_H1975_9291-ready.bw GSM5777979_X2c_3519AZ_H1975_9291-ready.bw \
        -R "$BED_FILE" \
        --referencePoint TSS \
        -b 3000 -a 3000 \
        --skipZeros \
        -o "${OUT_DIR}/H1975_${BED_LABEL}_matrix.gz"
plotProfile -m "${OUT_DIR}/H1975_${BED_LABEL}_matrix.gz" -out "${OUT_DIR}/H1975_${BED_LABEL}_profile.png" --averageType mean --colors blue blue blue red red red

computeMatrix reference-point -p 10 \
        -S GSM5777944_X5a_3519AZ_PC9_DMSO-ready.bw GSM5777945_X5b_3519AZ_PC9_DMSO-ready.bw GSM5777946_X5c_3519AZ_PC9_DMSO-ready.bw GSM5777947_X6a_3519AZ_PC9_9291-ready.bw GSM5777948_X6b_3519AZ_PC9_9291-ready.bw GSM5777949_X6c_3519AZ_PC9_9291-ready.bw \
        -R "$BED_FILE" \
        --referencePoint TSS \
        -b 3000 -a 3000 \
        --skipZeros \
        -o "${OUT_DIR}/PC9_${BED_LABEL}_matrix.gz"
plotProfile -m "${OUT_DIR}/PC9_${BED_LABEL}_matrix.gz" -out "${OUT_DIR}/PC9_${BED_LABEL}_profile.png" --averageType mean --colors blue blue blue red red red
