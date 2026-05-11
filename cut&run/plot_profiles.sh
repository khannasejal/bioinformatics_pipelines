#!/bin/bash
BW_DIR="/mnt/e/datasets/19-024-and-19-055/bw/spike_norm/rep2"
MATRIX_OUT="/mnt/e/datasets/cut&run/CUTnRUN_files/matrix_files"
PLOT_OUT="/mnt/e/datasets/cut&run/CUTnRUN_files/plots"

for bw in "$BW_DIR"/*.bw; do
    
    base=$(basename "$bw" .bw)
    echo "Processing sample: ${base}"

    computeMatrix scale-regions -p 10 \
        -S "$bw" \
        -R upregulated_memory_genes.bed upregulated_non_memory_genes.bed upregulated_refractory_genes.bed nonresponsive_genes_random.bed \
        --beforeRegionStartLength 5000 \
        --afterRegionStartLength 5000 \
        -o "${MATRIX_OUT}/${base}_upregulated_matrix.gz"
    echo "Matrix computed: ${base}"
    plotProfile \
        -m "${MATRIX_OUT}/${base}_upregulated_matrix.gz" \
        -out "${PLOT_OUT}/${base}_upregulated_profile.png" \
        --regionsLabel "Memory" "Non-memory" "Refractory" "Non-responsive" \
        --averageType mean \
        --plotType lines \
        --plotTitle "${base}_upregulated_profile" 
    echo "Profile ploted: ${base}"

done
