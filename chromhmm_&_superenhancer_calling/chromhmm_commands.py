#make the cellmarktable
java -mx4000M -jar ChromHMM/ChromHMM.jar BinarizeBam ChromHMM/CHROMSIZES/hg38.txt bam_bai_folder/fasr/k27ac cellmarktable_chromhmm.txt /mnt/d/sejal/chromhmm_analysis/output_binary_f/fasr_2_states
#learn model
java -mx8000M -jar ChromHMM/ChromHMM.jar LearnModel chromhmm_analysis/output_binary_f/tamr_2_states chromhmm_analysis/output_learn_f/tamr/k27ac_2_states 2 hg38
