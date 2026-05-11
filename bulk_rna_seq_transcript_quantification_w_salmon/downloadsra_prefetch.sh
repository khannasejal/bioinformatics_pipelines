#/bin/bash
srr_file="SRR_Acc_List_original.txt"
prefetch_dir="/mnt/e/tanmay_lele_19thmay/RNAseq/mdamb231/sra"
num_t=10
while read -r srr; do
    echo "Processing $srr"
    prefetch -p -O "$prefetch_dir" "$srr"
    echo "Downloaded $srr"
done < "$srr_file"