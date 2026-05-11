#%%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import re
# %%
pc9_accessibility = pd.read_csv("pc9_geneset_accessibility_summary.txt", sep="\t")
hcc2935_accessibility = pd.read_csv("hcc2935_geneset_accessibility_summary.txt", sep="\t")
hcc827_accessibility = pd.read_csv("hcc827_geneset_accessibility_summary.txt", sep="\t")
h1975_accessibility = pd.read_csv("h1975_geneset_accessibility_summary.txt", sep="\t")
# %%
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.titlesize'] = 14
dfs_list = [hcc827_accessibility, hcc2935_accessibility, h1975_accessibility, pc9_accessibility]
cell_lines = ["hcc827", "hcc2935", "h1975", "pc9"]


all_data = []

for df, cell in zip(dfs_list, cell_lines):
    temp_df = df.copy()
    
    # Add cell line label
    temp_df['Cell_Line'] = cell.upper()
    
    # Clean the Geneset column: Remove "{cell}_GRCh38_" prefix
    # This ensures "hcc2935_GRCh38_EMT" becomes "EMT"
    prefix_pattern = f"^{cell.lower()}_GRCh38_|^GRCh38_|^"
    temp_df['Clean_Geneset'] = temp_df['Geneset'].apply(
        lambda x: re.sub(r'^[a-z0-9]+_GRCh38_|^GRCh38_', '', x)
    )
    
    all_data.append(temp_df)

master_df = pd.concat(all_data)


# %%

fig, axes = plt.subplots(1, 4, figsize=(16, 8), sharey=True)

# 2. Loop through each cell line's data
for i, cell in enumerate(cell_lines):
    # Get data for this cell line from your master_df
    data = master_df[master_df['Cell_Line'] == cell.upper()].copy()
    
    # Sort by T-statistic so the plot looks organized
    data = data.sort_values('T_Statistic', ascending=False)
    
    # Create a color list: Bold color if p < 0.05, faded if not
    # Red for Up (T > 0), Blue for Down (T < 0)
    colors = []
    for _, row in data.iterrows():
        if row['P_Value'] < 0.05:
            colors.append('#d9534f' if row['T_Statistic'] > 0 else '#8f4adf')
        else:
            colors.append('#f2dede' if row['T_Statistic'] > 0 else "#bb95e5") # Faded versions

    # 3. Plot
    axes[i].barh(data['Clean_Geneset'], data['T_Statistic'], color=colors, edgecolor='none')
    
    # Add a vertical line at 0
    axes[i].axvline(0, color='black', linewidth=0.8)
    
    # Titles and labels
    axes[i].set_title(cell.upper(), fontsize=16, fontweight='bold')
    axes[i].set_xlabel("t-statistic")
    if i == 0:
        axes[i].set_ylabel("Geneset")
    
    sns.despine(ax=axes[i], left=True)

plt.suptitle("Differential Accessibility of Genesets between Controls and DTPs", fontsize=20, y=1.02, fontweight='bold')
plt.tight_layout()
#plt.savefig("Accessibility_Diverging_Bars.png", dpi=600, bbox_inches='tight')
plt.show()
# %%
