#%%
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pyranges as pr
# %%
col_headers = ["Chromosome","Start","End"]
mcf7_tads = pd.read_csv("MCF7_TADs.bed", sep='\t',header = None).drop_duplicates()
fasr_tads = pd.read_csv("FASR_TADs.bed", sep='\t',header = None).drop_duplicates()
tamr_tads = pd.read_csv("TAMR_TADs.bed", sep='\t',header = None).drop_duplicates()
mcf7_tads.columns = col_headers
fasr_tads.columns = col_headers
tamr_tads.columns = col_headers

# %%
df_TSS = pd.read_csv("TSS_hg38.txt", sep="\t")
df_TSS = df_TSS.dropna(subset=["Gene name"])
df_TSS = df_TSS[df_TSS["Chromosome/scaffold name"].isin([f"{i}" for i in range(1, 23)] + ["X"])] #filter
df_TSS["Chromosome/scaffold name"] = df_TSS["Chromosome/scaffold name"].astype(str)
df_TSS["Gene name"] = df_TSS["Gene name"].str.replace('_','-')
df_TSS["gene_transcript"] = df_TSS["Gene name"] + '_' + df_TSS["Transcript stable ID version"]
df_TSS['Chromosome/scaffold name'] = 'chr' + df_TSS['Chromosome/scaffold name'].astype(str)
#%%
print(df_TSS)
# %%
#function to do the overlap of TSS with the start and end of the TAD coordinate 
def overlap_tad_coords_w_tss(tads):
    pr1 = pr.PyRanges(tads)   
    # For df2 (assuming it has Chromosome and TSS)
    df_TSS["Start"] = df_TSS["Transcription start site (TSS)"]   
    df_TSS["End"] = df_TSS["Transcription start site (TSS)"]   
    pr2 = pr.PyRanges(df_TSS.rename(columns={"Chromosome/scaffold name": "Chromosome", "Start": "Start", "End": "End"}))
    overlap = pr1.join(pr2)
    overlap = overlap.df
    return overlap

# %%
'''#check if transcription start site lies exactly within the start and end of the TAD coordinate
count = 0
for index, row in overlap.iterrows():
    if (row["Start"]<=row["Start_b"]) and (row["End"] >= row["Start_b"]):
        count +=1
print(count)'''
# %%
#filter tads which contain atleast 1 ER gene
def filter_tads_containing_atleast_one_ER_gene(overlap):
    er_geneset = set(pd.read_csv("er_corr_genes.txt", header=None)[0])
    updated_er_geneset = {gene.replace("_", "-") for gene in er_geneset}
    filtered_overlap = overlap[overlap["Gene name"].isin(updated_er_geneset)]
    tads_with_atleast_one_ER_gene = filtered_overlap[["Chromosome","Start","End"]].drop_duplicates()
    return tads_with_atleast_one_ER_gene

chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
               'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20','chr21','chr22','chrX']
#%%
def count_fraction_of_tads_w_atleast_one_er_gene_chr_wise(tads):
    cts_er_tads = filter_tads_containing_atleast_one_ER_gene(overlap_tad_coords_w_tss(tads))["Chromosome"].value_counts().reindex(chromosomes).fillna(0).astype(int)
    cts = tads["Chromosome"].value_counts().reindex(chromosomes).fillna(0).astype(int)
    frac = cts_er_tads/cts
    return frac
# %%
#bar plot for fraction of TADs on each chromosome across three conditions
mcf7_frac = count_fraction_of_tads_w_atleast_one_er_gene_chr_wise(mcf7_tads)
fasr_frac = count_fraction_of_tads_w_atleast_one_er_gene_chr_wise(fasr_tads)
tamr_frac = count_fraction_of_tads_w_atleast_one_er_gene_chr_wise(tamr_tads)
bar_width = 0.25
x = np.arange(len(chromosomes))
#plotting
fig, ax = plt.subplots(figsize=(12, 6))
#Plot each group of bars
ax.bar(x - bar_width, mcf7_frac, width=bar_width, label='MCF7', color='navy')
ax.bar(x, fasr_frac, width=bar_width, label='FASR', color='skyblue')
ax.bar(x + bar_width, tamr_frac, width=bar_width, label='TAMR', color='lightgreen')
#Add labels, title, and legend
ax.set_xlabel('Chromosomes')
ax.set_ylabel('Number of TADs containing atleast 1 ER gene / Total number of TADs detected per chromosome')
ax.set_title('TAD Counts by Chromosome for Different Categories')
ax.set_xticks(x)
ax.set_xticklabels(chromosomes, rotation=45)
ax.legend()
#Show the plot
plt.tight_layout()
plt.show()
# %%
