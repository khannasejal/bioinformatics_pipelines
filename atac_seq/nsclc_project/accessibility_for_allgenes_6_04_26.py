#%%
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import os
import glob
import numpy as np
# %%
cell_line = "h1975" 
input_folder = cell_line
prom_bed_path = "E:\\shreya_nsclc\\GSE193258\\bw_narrowpeak\\bw\\osi_dtp\\genesets_bed\\promoters"
prom_bed = pd.read_csv(os.path.join(prom_bed_path, "GRCh38_All_RefSeq_Promoters.bed"), sep="\t", header=None, names=["#chr", "start", "end", "gene_id", "score", "strand"])
annot_df = prom_bed[["#chr", "start", "end", "gene_id"]].drop_duplicates()
tab_files = glob.glob(os.path.join(input_folder, "*.tab"))
tab_files_all_genes = glob.glob(os.path.join("E:\\shreya_nsclc\\GSE193258\\bw_narrowpeak\\bw\\osi_dtp\\quantification_results\\all_gene_promoter_quantifications", "*.tab"))

#%%

for file_path in tab_files_all_genes:  #change it to whchever tab files are to be annotated
    print(f"Processing: {os.path.basename(file_path)}")
    df = pd.read_csv(file_path, sep="\t")
    df.columns = [re.sub(r"['\"]", "", c).strip() for c in df.columns]
    
    df = pd.merge(df, annot_df, on=["#chr", "start", "end"], how="left")
    cols = df.columns.tolist()
    if 'gene_id' in cols:
        cols.insert(3, cols.pop(cols.index('gene_id')))
        df = df[cols]

    ctrl_cols = [c for c in df.columns if 'Ctrl' in c]
    dtp_cols = [c for c in df.columns if 'DTP' in c]

    for col in ctrl_cols + dtp_cols:
        df[col] = df[col].astype(str).str.replace("'", "").str.replace('"', "")
        df[col] = pd.to_numeric(df[col], errors='coerce')

    df[ctrl_cols + dtp_cols] = df[ctrl_cols + dtp_cols].fillna(0)
    df['Mean_Control'] = df[ctrl_cols].mean(axis=1)
    df['Mean_DTP'] = df[dtp_cols].mean(axis=1)
    df['Direction'] = np.where(df['Mean_DTP'] > df['Mean_Control'], "Up", "Down")

    def run_ttest(row):
        c_vals = row[ctrl_cols].values.astype(float)
        d_vals = row[dtp_cols].values.astype(float)
        t_stat, p_val = stats.ttest_ind(d_vals, c_vals, equal_var=False)
        return pd.Series({'T_Statistic': t_stat, 'P_Value': p_val})
    
    stats_results = df.apply(run_ttest, axis=1)
    df = pd.concat([df, stats_results], axis=1)
    df.to_csv(file_path, sep="\t", index=False)

print("\nSuccess: All files in the folder are now annotated with gene names and T-test results.")

# %%
