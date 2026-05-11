#%%
#calculate the number of on each chromosome in mcf7/fasr/tamr conditions
import random
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# %%
chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
               'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20','chr21','chr22','chrX']

col_headers = ["Chromosome","Start","End"]
mcf7_tads = pd.read_csv("MCF7_TADs.bed", sep='\t',header = None).drop_duplicates()
fasr_tads = pd.read_csv("FASR_TADs.bed", sep='\t',header = None).drop_duplicates()
tamr_tads = pd.read_csv("TAMR_TADs.bed", sep='\t',header = None).drop_duplicates()
mcf7_tads.columns = col_headers
fasr_tads.columns = col_headers
tamr_tads.columns = col_headers

# %%
mcf7_cts = mcf7_tads["Chromosome"].value_counts().reindex(chromosomes).fillna(0).astype(int)
fasr_cts = fasr_tads["Chromosome"].value_counts().reindex(chromosomes).fillna(0).astype(int)
tamr_cts = tamr_tads["Chromosome"].value_counts().reindex(chromosomes).fillna(0).astype(int)
mcf7_frac = mcf7_cts/len(mcf7_tads)
fasr_frac = fasr_cts/len(fasr_tads)
tamr_frac = tamr_cts/len(tamr_tads)

# %%
#bar plot for fraction of TADs on each chromosome across three conditions
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
ax.set_ylabel('#TADs')
ax.set_title('TAD Counts by Chromosome for Different Categories')
ax.set_xticks(x)
ax.set_xticklabels(chromosomes, rotation=45)
ax.legend()
#Show the plot
plt.tight_layout()
plt.show()


# %%
#plot the distribution of tad sizes across the conditions
mcf7_tads['tad_size'] = mcf7_tads['End'] - mcf7_tads['Start']
mcf7_tads['condition'] = 'MCF7'
fasr_tads['tad_size'] = fasr_tads['End'] - fasr_tads['Start']
fasr_tads['condition'] = 'FASR'
tamr_tads['tad_size'] = tamr_tads['End'] - tamr_tads['Start']
tamr_tads['condition'] = 'TAMR'
df_all = pd.concat([mcf7_tads, fasr_tads, tamr_tads], ignore_index=True)


#%%
# Create the box plot for all three conditions
plt.figure(figsize=(15, 10))
sns.boxplot(data=df_all, x='condition', y='tad_size')
# Customize plot labels and title
plt.title('TAD Size Distribution for Three Conditions')
plt.xlabel('Condition')
plt.ylabel('TAD Size')
# Show plot
plt.tight_layout()
plt.show()


# %%
# Create the KDE plot using sns.displot (distribution plot for size of TADs)
sns.displot(df_all, x='tad_size', hue='condition', kind='kde', fill=True)
# Customize plot labels and title
plt.title('TAD Size Distribution for Three Conditions')
plt.xlabel('TAD Size')
plt.ylabel('Density')
# Show plot
plt.tight_layout()
plt.show()


# %%
#Violin plot for the distribution
sns.violinplot(data=df_all,x = 'condition', y='tad_size', hue='condition', fill=True)
# Customize plot labels and title
plt.title('TAD Size Distribution for Three Conditions')
plt.xlabel('Condition')
plt.ylabel('TAD size')
# Show plot
plt.tight_layout()
plt.show()


# %%
#stripplot plot
sns.stripplot(data=df_all,x='condition', y='tad_size', jitter=True)
# Customize plot labels and title
plt.title('TAD Size Distribution for Three Conditions')
plt.xlabel('Condition')
plt.ylabel('TAD size')
# Show plot
plt.tight_layout()
plt.show()
# %%
