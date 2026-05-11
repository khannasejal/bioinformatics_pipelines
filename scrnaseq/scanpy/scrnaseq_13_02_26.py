#%%
#cleaner version of the previous code 
#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import scipy.sparse
# %%
adata = sc.read_h5ad("GSM6428979_adata.adenonepc.PtRP.final.h5ad")
# %%
def calculate_global_entropy(adata):
    # Use the normalized matrix (X)
    data = adata.X.toarray() if scipy.sparse.issparse(adata.X) else adata.X
    # Normalize each cell to sum to 1
    p = data / data.sum(axis=1, keepdims=True)
    entropy = -np.sum(p * np.log(p + 1e-9), axis=1)
    return entropy
# %%
adata.obs['entropy'] = calculate_global_entropy(adata)
#%%
sample_order = [
    'PRP_8weeks_Intact_R1', 
    'PRP_9weeks_Intact_R1', 
    'PRP_12weeks_Intact_R1', 'PRP_12weeks_Intact_R2', 'PRP_12weeks_Intact_R3', 'PRP_12weeks_Intact_R4',
    'PRP_12weeks_Cas_R1', 'PRP_12weeks_Cas_R2', 'PRP_12weeks_Cas_R3',
    'PRP_12weeks_DHT_R1', 'PRP_12weeks_DHT_R2', 'PRP_12weeks_DHT_R3',
    'PRP_16weeks_Intact_R1'
]

# %%
hallmark_emt_genes = gene_list = [
    "COL3A1","COL5A2","COL5A1","FBN1","COL1A1","FN1","COL6A3","SERPINE1",
    "COL1A2","COL4A1","COL4A2","VCAN","IGFBP3","TGFBI","SPARC","LUM",
    "LAMC1","LOX","LAMC2","CCN2","TAGLN","COL7A1","LOXL2","COL6A2",
    "ITGAV","THBS2","COL16A1","NNMT","TPM1","CDH2","MMP2","COL11A1",
    "THBS1","FAP","BGN","SERPINH1","FSTL1","POSTN","THY1","SPP1","TNC",
    "TFPI2","NID2","ITGB5","MMP3","VIM","LOXL1","FBLN5","COL12A1","ELN",
    "CDH11","COMP","SPOCK1","BMP1","IL32","LAMA3","TIMP1","QSOX1","TIMP3",
    "VCAM1","CCN1","EDIL3","CALD1","MAGEE1","FBLN1","SGCB","ECM1",
    "LAMA2","FSTL3","TPM2","INHBA","DAB2","EMP3","BASP1","ITGA5","MGP",
    "VEGFA","CXCL1","WNT5A","SDC1","PLOD2","PCOLCE","GREM1","ITGB1",
    "COL5A3","RHOB","HTRA1","FGF2","SNTB1","GADD45A","MEST","LRRC15",
    "TNFRSF11B","CD59","ACTA2","EFEMP2","MATN2","PCOLCE2","SERPINE2",
    "GPC1","ABI3BP","FUCA1","SLIT3","LAMA1","PMEPA1","COL8A2","FBN2",
    "IGFBP2","PFN2","SDC4","CD44","GADD45B","CXCL8","GLIPR1","ANPEP",
    "P3H1","VEGFC","MMP14","SGCD","PLOD1","MATN3","MYL9","SLC6A8",
    "CALU","PRRX1","TNFRSF12A","FMOD","ID2","GEM","PLAUR","MYLK",
    "TGFB1","SFRP1","PLOD3","IL6","APLP1","FBLN2","MSX1","PTX3","FZD8",
    "JUN","FERMT2","DKK1","SNAI2","DST","TPM4","DCN","GJA1","PMP22",
    "IGFBP4","COPA","LRP1","ITGA2","FLNA","MFAP5","PTHLH","TGFBR3",
    "SFRP4","LGALS1","RGS4","CDH6","SAT1","NT5E","DPYSL3","PPIB","TGM2",
    "SGCG","ITGB3","PDLIM4","CTHRC1","ECM2","CRLF1","AREG","IL15","MCM7",
    "GAS1","PRSS2","CADM1","OXTR","SCG2","CXCL6","MMP1","TNFAIP3","CAPG",
    "CAP2","MXRA5","FOXC2","NTM","ENO2","FAS","BDNF","ADAM12","PVR",
    "CXCL12","PDGFRB","SLIT2","NOTCH2","COLGALT1","GPX7","WIPF1"
]
available_emt = [g for g in hallmark_emt_genes if g in adata.var_names]
adata.obs['emt_entropy'] = calculate_global_entropy(adata[:, available_emt])
adata.obs['sample_clean'] = adata.obs['sample'].str.replace(r'_R\d+$', '', regex=True)

#%%
def load_signatures_from_gmt(file_path):
    signatures = {}
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            sig_name = parts[0] 
            genes = [g.upper() for g in parts[2:]]
            signatures[sig_name] = genes
    return signatures

path_to_gmt = 'cell_state_signatures.gmt'
all_signatures = load_signatures_from_gmt(path_to_gmt)

for sig_name, gene_list in all_signatures.items():
    # Filter list for genes actually present in your adata
    valid_genes = [g for g in gene_list if g in adata.var_names]
    if len(valid_genes) > 0:
        print(f"Scoring {sig_name} with {len(valid_genes)} genes...")
        sc.tl.score_genes(adata, gene_list=valid_genes, score_name=sig_name)
    else:
        print(f"Warning: No genes found for {sig_name} in adata.var_names.")

new_sig_columns = list(all_signatures.keys())
core_columns = ['pseudotime', 'entropy', 'emt_entropy','JAK_STAT', 'cell_type_final', 'sample_clean','sample']

columns_to_extract = [c for c in core_columns + new_sig_columns if c in adata.obs.columns]
plot_df = adata.obs[columns_to_extract].copy()

print("Plotting DataFrame ready!")
print(plot_df.head())
#%%

x= "Plasticity"
y="entropy"
#plot entropy on one axis and EMT score on the other axis
sns.scatterplot(data=plot_df, x=x, y=y, hue='cell_type_final', s=8, linewidth=0, palette="tab20")

sns.regplot(data=plot_df, x=x, y=y, 
            scatter=False, lowess=True, color='black', line_kws={'linewidth': 3})
plt.legend(title='Samples', bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

plt.tight_layout()
plt.show()
# %%
#plt.figure(figsize=(15, 6))

# We use the 'sample' column to keep replicates separate as requested
# and the 'sample_order' list you created to force the layout
sns.violinplot(
    data=plot_df, 
    x='sample', 
    y='Plasticity', 
    order=sample_order, 
    palette='magma', 
    inner='quartile' # Adds dashed lines for the median and quartiles
)

plt.xticks(rotation=90)
plt.title("Plasticity Signature Distribution across PtRP Samples")
plt.xlabel("Experimental Condition (Time-Course)")
plt.ylabel("Plasticity Score")

plt.tight_layout()
plt.show()
# %%
