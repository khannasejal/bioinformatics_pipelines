#%%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import os
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
import pyranges as pr
from scipy.stats import ttest_ind
# %%
deseq2_dir = "E:\\tanmay_lele_19thmay\\ATAC\\final\\deseq2_results"
df_TSS = pd.read_csv("D:/sejal_c/sejal/hic_breast_cancer/TSS_hg38_20_08_25.txt", sep="\t")
df_TSS = df_TSS.dropna(subset=["Gene name"])
df_TSS_filt = df_TSS[df_TSS["Chromosome/scaffold name"].astype(str).isin([str(i) for i in range(1, 23)] + ["X"])]
df_TSS_filt['Chromosome/scaffold name'] = 'chr' + df_TSS['Chromosome/scaffold name'].astype(str)
df_TSS_filt = df_TSS_filt.rename(columns={"Chromosome/scaffold name": "Chromosome"})
#%%
#for the same transcript length, it picks the first instance
df_TSS_maxlen = (
    df_TSS_filt.loc[df_TSS_filt.groupby("Gene name")["Transcript length (including UTRs and CDS)"].idxmax()]
    .reset_index(drop=True)
)
#%%
deseq2_so_sel_so_anc = pd.read_csv(os.path.join(deseq2_dir,"so_sel_so_anc\\b2-so-sel_b2-so-anc_sig_peaks.txt"), sep="\t")
deseq2_st_sel_st_anc = pd.read_csv(os.path.join(deseq2_dir,"st_sel_st_anc\\b2-st-sel_b2-st-anc_sig_peaks.txt"), sep="\t")
#%%
'''deseq2_so_sel_so_anc['peak_length'] = deseq2_so_sel_so_anc['End'] - deseq2_so_sel_so_anc['Start']
plt.hist(deseq2_so_sel_so_anc["peak_length"], bins=30, edgecolor="black")'''

#%%
#check which genes have the same transcript length
'''same_len_genes = []

for gene, group in df_TSS_filt.groupby("Gene name"):
    max_len = group["Transcript length (including UTRs and CDS)"].max()
    max_rows = group[group["Transcript length (including UTRs and CDS)"] == max_len]
    if len(max_rows) > 1:
        same_len_genes.append((gene, max_len, len(max_rows)))'''
#%%
#calculate the overlap percentage wrt to peak length mapped to the gene; pick the peak with the maximum overlap

# %%
def annotate_atac_peaks_for_interval(deseq2_df, comp_group, window, TSS_df):
    TSS_df['Start'] = TSS_df['Transcription start site (TSS)'] - window
    TSS_df['End'] = TSS_df['Transcription start site (TSS)'] + window
    df1 = pr.PyRanges(TSS_df)
    df2 = pr.PyRanges(deseq2_df)
    overlap = df2.join(df1, report_overlap=True)
    overlap_df = overlap.df.copy()
    overlap_df["PeakLength"] = overlap_df["End"] - overlap_df["Start"]
    overlap_df["OverlapPerc"] = (overlap_df["Overlap"] / overlap_df["PeakLength"]) * 100
    #overlap_df_mean = overlap_df.groupby('Gene name', as_index=False)['log2FoldChange'].mean() #for each gene, mapped to peak, take the average 
    #overlap_df_mean['log2FoldChange'] = pd.to_numeric(overlap_df_mean['log2FoldChange'], errors='coerce')
    #overlap_df_mean.to_csv(f"1000bp_annotated_deseq2_{comp_group}_peaks_unique_transcripts.txt", index=False, sep="\t")
    return overlap_df
#%%
'''window = 1000
TSS_df = df_TSS_maxlen
TSS_df['Start'] = TSS_df['Transcription start site (TSS)'] - window
TSS_df['End'] = TSS_df['Transcription start site (TSS)'] + window
df1 = pr.PyRanges(TSS_df)
df2 = pr.PyRanges(deseq2_so_sel_so_anc)
overlap = df2.join(df1, report_overlap=True)
overlap_df = overlap.df.copy()
overlap_df["PeakLength"] = overlap_df["End"] - overlap_df["Start"]
overlap_df["OverlapPerc"] = (overlap_df["Overlap"] / overlap_df["PeakLength"]) * 100

dup_genes = overlap_df[overlap_df["Gene name"].duplicated(keep=False)]
print(dup_genes)
same_overlap = (
    dup_genes.groupby("Gene name")
    .filter(lambda x: x["OverlapPerc"].nunique() < len(x))
)'''

#%%
annotated_so_sel_so_anc = annotate_atac_peaks_for_interval(deseq2_so_sel_so_anc, 'so_sel_so_anc', 1000, df_TSS_filt)
annotated_so_sel_so_anc_unique_transcripts = annotate_atac_peaks_for_interval(deseq2_so_sel_so_anc, 'so_sel_so_anc', 1000, df_TSS_maxlen)
annotated_st_sel_st_anc_unique_transcripts = annotate_atac_peaks_for_interval(deseq2_st_sel_st_anc, 'st_sel_st_anc', 1000, df_TSS_maxlen)

#%%
annotated_so_sel_so_anc_unique_transcripts["PeakCenter"] = (
    (annotated_so_sel_so_anc_unique_transcripts["Start"] +
     annotated_so_sel_so_anc_unique_transcripts["End"]) / 2
)
annotated_so_sel_so_anc_unique_transcripts["DistanceToTSS"] = (
    (annotated_so_sel_so_anc_unique_transcripts["PeakCenter"] -
     annotated_so_sel_so_anc_unique_transcripts["Transcription start site (TSS)"]).abs()
)

closest_peaks = annotated_so_sel_so_anc_unique_transcripts.loc[
    annotated_so_sel_so_anc_unique_transcripts.groupby("Gene name")["DistanceToTSS"].idxmin()
].reset_index(drop=True)
#%%
annotated_st_sel_st_anc_unique_transcripts["PeakCenter"] = (
    (annotated_st_sel_st_anc_unique_transcripts["Start"] +
     annotated_st_sel_st_anc_unique_transcripts["End"]) / 2
)
annotated_st_sel_st_anc_unique_transcripts["DistanceToTSS"] = (
    (annotated_st_sel_st_anc_unique_transcripts["PeakCenter"] -
     annotated_st_sel_st_anc_unique_transcripts["Transcription start site (TSS)"]).abs()
)

closest_peaks_st_sel_st_anc = annotated_st_sel_st_anc_unique_transcripts.loc[
    annotated_st_sel_st_anc_unique_transcripts.groupby("Gene name")["DistanceToTSS"].idxmin()
].reset_index(drop=True)





#%%
closest_peaks_st_sel_st_anc.to_csv(f"1000bp_annotated_deseq2_st_sel_st_anc_peaks_unique_transcripts_closest_peaks.txt", index=False, sep="\t")
#%%












# %%


# %%
annotated_so_sel_so_anc = annotate_atac_peaks_for_interval(deseq2_so_sel_so_anc, 'so_sel_so_anc', 1000, df_TSS_filt)
annotated_st_sel_st_anc = annotate_atac_peaks_for_interval(deseq2_st_sel_st_anc, 'st_sel_st_anc', 1000, df_TSS_filt)
# %%
