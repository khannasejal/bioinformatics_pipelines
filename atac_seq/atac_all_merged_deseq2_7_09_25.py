#%%
#use ncbi refseq gene annotations to get the differential accessible regions in the promoter/gene body regions 
#map them with the differentially expressed genes (students' t-test)

#%%
import pandas as pd
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
import pyranges as pr
from scipy.stats import ttest_ind
#%%
all_atac_merged_counts = pd.read_csv("E:\\tanmay_lele_19thmay\\ATAC\\final\\ATAC_raw_counts_reprocessed_unannotated.txt", sep = "\t")
sample_cols = [col for col in all_atac_merged_counts.columns if col not in ['Chromosome', 'Start','End']]

#%%
#read the gene start and end coordinates from the ncbi ref seq 
#try for transcript too
df_TSS = pd.read_csv("grch38_ncbi_refseq_gene_7_09_25.txt", sep="\t")
df_TSS = df_TSS[df_TSS["gene_biotype"] == "protein_coding"]
#based on the strand, consider only the 1kbp upstream from the gene start
# For + strand: (Start - 1000) to Start
window = 1000
def define_region_to_do_deseq2(df,x): #remove x for transcript length quantification
    #based on the strand, consider only the 1kbp upstream from the transcript start
    #For + strand: (Start - 1000) to Start
    df_copy = df.copy()
    #use this code for only upstream region
    '''df_copy.loc[df_copy['Strand'] == "+", 'End']   = df_copy.loc[df_copy['Strand'] == "+", 'Start'] 
    df_copy.loc[df_copy['Strand'] == "+", 'Start'] = df_copy.loc[df_copy['Strand'] == "+", 'Start'] - x

    # For - strand: End to (End + 1000)
    df_copy.loc[df_copy['Strand'] == "-", 'Start'] = df_copy.loc[df_copy['Strand'] == "-", 'End'] 
    df_copy.loc[df_copy['Strand'] == "-", 'End']   = df_copy.loc[df_copy['Strand'] == "-", 'End'] + x''' #comment all these to run for entire transcript length quantification
    
    #use this code for upstream and downstream range
    '''df_copy.loc[df_copy['Strand'] == "+", 'End']   = df_copy.loc[df_copy['Strand'] == "+", 'Start'] + x
    df_copy.loc[df_copy['Strand'] == "+", 'Start'] = df_copy.loc[df_copy['Strand'] == "+", 'Start'] - x

    # For - strand: End to (End + 1000)
    df_copy.loc[df_copy['Strand'] == "-", 'Start'] = df_copy.loc[df_copy['Strand'] == "-", 'End'] - x
    df_copy.loc[df_copy['Strand'] == "-", 'End']   = df_copy.loc[df_copy['Strand'] == "-", 'End'] + x'''
    return df_copy
#do not run any of it for gene body

# %%
df1 = pr.PyRanges(define_region_to_do_deseq2(df_TSS, 1000))
df2 = pr.PyRanges(all_atac_merged_counts)
overlap = df2.join(df1, report_overlap=True) #start_b and end_b coordinates are promoter regions defined from gene start; start and end are peak start and end coordinates; maps all the peaks within the promoter regions
overlap_df = overlap.df.copy()
# %%
grouped = (
    overlap_df
    .groupby("gene_id", as_index=False)
    .agg(
        {**{col: "sum" for col in sample_cols},  
         "Chromosome": "first",                 
         "Start_b": "first",
         "End_b": "first"}
    )
)
final_df = grouped[["Chromosome", "Start_b", "End_b", "gene_id"] + sample_cols]
# %%
deseq_df = final_df[sample_cols].T
conditions_array = ['b1_so_sel']*4+['b1_st_sel']*4+['b2_so_sel']*4+['b2_st_sel']*4+['b1_so_anc','b2_so_anc']*4+['b1_st_anc','b2_st_anc']*4
metadata_df = pd.DataFrame({"condition": conditions_array}, index=sample_cols)

dds = DeseqDataSet(
    counts=deseq_df,
    metadata=metadata_df,
    design_factors="condition"
)

#Perform DESeq2 analysis
dds.deseq2()
print(dds)
# %%
def do_deseq2_for_a_comparison_grp(group1_str, group2_str):
    ds_grp1_grp2 = DeseqStats(dds, contrast=["condition", group1_str, group2_str])
    summary_grp1_grp2 = ds_grp1_grp2.summary()
    results_grp1_grp2 = ds_grp1_grp2.results_df
    results_grp1_grp2['Chromosome'] = final_df['Chromosome'].to_list()
    results_grp1_grp2['Start'] = final_df['Start_b'].to_list()
    results_grp1_grp2['End'] = final_df['End_b'].to_list()
    results_grp1_grp2['gene_id'] = final_df['gene_id'].to_list()
    results_grp1_grp2['-log10pvalue'] = -np.log10(results_grp1_grp2['pvalue'])
    sig_results_df = results_grp1_grp2[results_grp1_grp2['pvalue'] < 0.05]
    #sig_results_df.to_csv(f"{group1_str}_{group2_str}_gene_body_sig_peaks.txt", sep="\t",index=False) #saves the significant/all results for summed up reads across gene body and deseq2 on them
    results_grp1_grp2.to_csv(f"{group1_str}_{group2_str}_gene_body_all_peaks.txt", sep="\t",index=False)
    return sig_results_df
#%%
'''so_sel_so_anc = do_deseq2_for_a_comparison_grp("b2-so-sel", "b2-so-anc")
so_sel_st_anc = do_deseq2_for_a_comparison_grp("b2-so-sel", "b2-st-anc")
st_sel_st_anc = do_deseq2_for_a_comparison_grp("b2-st-sel", "b2-st-anc")
so_anc_st_anc = do_deseq2_for_a_comparison_grp("b2-so-anc", "b2-st-anc")
so_sel_st_sel = do_deseq2_for_a_comparison_grp("b2-so-sel", "b2-st-sel")'''
#%%
#read all the files from the folders and merge them with degs
degs_dir = "E:\\tanmay_lele_19thmay\\ATAC\\final\\ttest_gene_exp"
atac_dir1 = "E:\\tanmay_lele_19thmay\\ATAC\\final\\deseq2_results\\ncbi_refseq\\gene\\1kbp_up_down\\all_peaks"
atac_dir2 = "E:\\tanmay_lele_19thmay\\ATAC\\final\\deseq2_results\\ncbi_refseq\\gene\\gene_body\\all_peaks"

#%%
def merge_atac_gene_exp_dfs(atac_dir, atac_file, degs_file, comparison):
    atac_df = pd.read_csv(os.path.join(atac_dir,atac_file), sep="\t")
    degs_df = pd.read_csv(os.path.join(degs_dir,degs_file), sep="\t")
    merged_df = pd.merge(atac_df,degs_df,left_on="gene_id", right_on="gene_name")
    merged_df.to_csv(f"{comparison}_atac_gene_exp_map.txt", sep="\t", index=False)
    return merged_df

# %%
file_pairs1 = [
    ('b2-so-sel_b2-so-anc_1kbp_up_down_from_start_all_peaks.txt', 'so_sel_mean_so_anc_mean_ttest.txt', 'so_sel_so_anc'),
    ('b2-so-sel_b2-st-anc_1kbp_up_down_from_start_all_peaks.txt', 'so_sel_mean_st_anc_mean_ttest.txt', 'so_sel_st_anc'),
    ('b2-st-sel_b2-st-anc_1kbp_up_down_from_start_all_peaks.txt', 'st_sel_mean_st_anc_mean_ttest.txt', 'st_sel_st_anc'),
    ('b2-so-anc_b2-st-anc_1kbp_up_down_from_start_all_peaks.txt', 'so_anc_mean_st_anc_mean_ttest.txt', 'so_anc_st_anc'),
    ('b2-so-sel_b2-st-sel_1kbp_up_down_from_start_all_peaks.txt', 'so_sel_mean_st_sel_mean_ttest.txt', 'so_sel_st_sel')
]
file_pairs2 = [
    ('b2-so-sel_b2-so-anc_gene_body_all_peaks.txt', 'so_sel_mean_so_anc_mean_ttest.txt', 'so_sel_so_anc'),
    ('b2-so-sel_b2-st-anc_gene_body_all_peaks.txt', 'so_sel_mean_st_anc_mean_ttest.txt', 'so_sel_st_anc'),
    ('b2-st-sel_b2-st-anc_gene_body_all_peaks.txt', 'st_sel_mean_st_anc_mean_ttest.txt', 'st_sel_st_anc'),
    ('b2-so-anc_b2-st-anc_gene_body_all_peaks.txt', 'so_anc_mean_st_anc_mean_ttest.txt', 'so_anc_st_anc'),
    ('b2-so-sel_b2-st-sel_gene_body_all_peaks.txt', 'so_sel_mean_st_sel_mean_ttest.txt', 'so_sel_st_sel')
]
#%%
for peaks_file, expr_file, comparison in file_pairs:
    print(peaks_file)
    merge_atac_gene_exp_dfs(peaks_file, expr_file, comparison)
# %%
modules_dir = "E:\\tanmay_lele_19thmay\\wgcna_23rdmay\\significant_modules_anova"
module_gene_map = pd.read_csv(os.path.join(modules_dir,"module_colors_for_data_sig.txt"), sep='\t')
# %%
def merge_atac_module_genes(atac_dir, atac_file, comparison):
    atac_df = pd.read_csv(os.path.join(atac_dir,atac_file), sep="\t")
    merged_df = pd.merge(atac_df,module_gene_map,left_on="gene_id", right_on="Gene")
    out_file = f"{comparison}_atac_module_map.txt"
    merged_df.to_csv(out_file, sep="\t", index=False)
    return merged_df
# %%
for peaks_file, _, comparison in file_pairs2:
    print(peaks_file)
    merge_atac_module_genes(atac_dir2, peaks_file, comparison)

# %%
#read the deseq2 file
