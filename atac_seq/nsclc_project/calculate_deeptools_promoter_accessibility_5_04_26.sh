#!/bin/bash

OUT_DIR="quantification_results"
#mkdir -p "$OUT_DIR"
BED_DIR="/mnt/e/shreya_nsclc/GSE193258/bw_narrowpeak/bw/osi_dtp/genesets_bed/promoters"

for bedfile in "$BED_DIR"/*GRCh38_All_RefSeq_Promoters.bed; do
    base_name=$(basename "$bedfile" _Promoters.bed)
    echo "---------------------------------------------------"
    echo "Processing: $base_name"
    echo "Input BED: $bedfile"
    multiBigwigSummary BED-file \
      --bwfiles GSM5777956_X1a_3562AZ_HCC827_CTRL-ready.bw GSM5777957_X1b_3562AZ_HCC827_CTRL-ready.bw GSM5777958_X1c_3562AZ_HCC827_CTRL-ready.bw GSM5777959_X2a_3562AZ_HCC827_DTP-ready.bw GSM5777960_X2b_3562AZ_HCC827_DTP-ready.bw GSM5777961_X2c_3562AZ_HCC827_DTP-ready.bw \
      --BED "$bedfile" \
      --labels HCC827_Ctrl1 HCC827_Ctrl2 HCC827_Ctrl3 HCC827_DTP1 HCC827_DTP2 HCC827_DTP3 \
      -o "$OUT_DIR/hcc827_${base_name}.npz" \
      --outRawCounts "$OUT_DIR/hcc827_${base_name}_counts.tab"
    echo "Quantification completed for hcc827_${base_name}.npz and counts saved to hcc827_${base_name}_counts.tab"
    
    multiBigwigSummary BED-file \
      --bwfiles GSM5777965_X4a_3562AZ_HCC2935_CTRL-ready.bw GSM5777966_X4b_3562AZ_HCC2935_CTRL-ready.bw GSM5777967_X4c_3562AZ_HCC2935_CTRL-ready.bw GSM5777968_X5a_3562AZ_HCC2935_DTP-ready.bw GSM5777969_X5b_3562AZ_HCC2935_DTP-ready.bw GSM5777970_X5c_3562AZ_HCC2935_DTP-ready.bw \
      --BED "$bedfile" \
      --labels HCC2935_Ctrl1 HCC2935_Ctrl2 HCC2935_Ctrl3 HCC2935_DTP1 HCC2935_DTP2 HCC2935_DTP3 \
      -o "$OUT_DIR/hcc2935_${base_name}.npz" \
      --outRawCounts "$OUT_DIR/hcc2935_${base_name}_counts.tab"
    echo "Quantification completed for hcc2935_${base_name}.npz and counts saved to hcc2935_${base_name}_counts.tab"

    multiBigwigSummary BED-file \
      --bwfiles GSM5777974_X1a_3519AZ_H1975_DMSO-ready.bw GSM5777975_X1b_3519AZ_H1975_DMSO-ready.bw GSM5777976_X1c_3519AZ_H1975_DMSO-ready.bw GSM5777977_X2a_3519AZ_H1975_9291-ready.bw GSM5777978_X2b_3519AZ_H1975_9291-ready.bw GSM5777979_X2c_3519AZ_H1975_9291-ready.bw \
      --BED "$bedfile" \
      --labels H1975_Ctrl1 H1975_Ctrl2 H1975_Ctrl3 H1975_DTP1 H1975_DTP2 H1975_DTP3 \
      -o "$OUT_DIR/h1975_${base_name}.npz" \
      --outRawCounts "$OUT_DIR/h1975_${base_name}_counts.tab"
    echo "Quantification completed for h1975_${base_name}.npz and counts saved to h1975_${base_name}_counts.tab"

    multiBigwigSummary BED-file \
      --bwfiles GSM5777944_X5a_3519AZ_PC9_DMSO-ready.bw GSM5777945_X5b_3519AZ_PC9_DMSO-ready.bw GSM5777946_X5c_3519AZ_PC9_DMSO-ready.bw GSM5777947_X6a_3519AZ_PC9_9291-ready.bw GSM5777948_X6b_3519AZ_PC9_9291-ready.bw GSM5777949_X6c_3519AZ_PC9_9291-ready.bw \
      --BED "$bedfile" \
      --labels PC9_Ctrl1 PC9_Ctrl2 PC9_Ctrl3 PC9_DTP1 PC9_DTP2 PC9_DTP3 \
      -o "$OUT_DIR/pc9_${base_name}.npz" \
      --outRawCounts "$OUT_DIR/pc9_${base_name}_counts.tab"
    echo "Quantification completed for pc9_${base_name}.npz and counts saved to pc9_${base_name}_counts.tab"
    echo "---------------------------------------------------"

done





