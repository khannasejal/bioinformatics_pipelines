#%%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import os
import glob
import textwrap

# %%

cell_line = "hcc2935" 
input_folder = cell_line
#%%
output_folder = f"{cell_line}_plots"
summary_output = f"{cell_line}_geneset_accessibility_summary.txt"

os.makedirs(output_folder, exist_ok=True)

results_list = []
tab_files = glob.glob(os.path.join(input_folder, "*.tab"))
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 14


for file_path in tab_files:
    filename = os.path.basename(file_path)
    # Extract geneset/cell line name from filename
    geneset_name = filename.replace("_counts.tab", "")
    
    # Load data
    df = pd.read_csv(file_path, sep="\t")
    
    # Identify data columns (skip chr, start, end)
    data_cols = df.columns[3:]
    
    # 1. Calculate the mean accessibility for every gene in each sample
    # This gives one value per replicate
    replicate_averages = df[data_cols].mean()
    
    # 2. Organize into a temporary DataFrame to separate groups
    # Note: Handles 'Ctrl' and 'DMSO' as Control, others as DTP
    temp_df = pd.DataFrame({
        'Replicate': replicate_averages.index,
        'Accessibility': replicate_averages.values
    })
    temp_df['Condition'] = temp_df['Replicate'].apply(
        lambda x: 'Control' if any(c in x for c in ['Ctrl', 'DMSO']) else 'DTP'
    )
    
    # 3. Perform T-test
    ctrl_vals = temp_df[temp_df['Condition'] == 'Control']['Accessibility']
    dtp_vals = temp_df[temp_df['Condition'] == 'DTP']['Accessibility']
    
    # Calculate stats
    t_stat, p_val = stats.ttest_ind(dtp_vals, ctrl_vals, equal_var=False)  # Welch's t-test
    
    # Determine direction
    mean_ctrl = ctrl_vals.mean()
    mean_dtp = dtp_vals.mean()
    direction = "Up" if mean_dtp > mean_ctrl else "Down"
    
    # Append to results list
    results_list.append({
        'Geneset': geneset_name,
        'T_Statistic': t_stat,
        'P_Value': p_val,
        'Direction': direction,
        'Mean_Control': mean_ctrl,
        'Mean_DTP': mean_dtp
    })
    
    # 4. Generate and save the Bar Plot
    plt.figure(figsize=(4, 5))
    ax = sns.barplot(
        data=temp_df, 
        x='Condition', 
        y='Accessibility', 
        hue='Condition',      # Assign 'Condition' to hue to satisfy the new requirement
        legend=False,         # Remove the legend since the x-axis already labels the groups
        palette={'Control': "#839099", 'DTP': "#8f4adf"},
        capsize=0.1,
        errorbar='sd',
        edgecolor='grey'
    )
    sns.stripplot(
        data=temp_df, 
        x='Condition', 
        y='Accessibility', 
        color='black', 
        size=8, 
        alpha=0.6
    )

    plt.title(f"{geneset_name}\np = {p_val:.4f}", fontsize=12)
    plt.ylabel("Average CPM Signal")
    plt.xlabel("")
    plt.tight_layout()
    
    plot_name = f"{geneset_name}_barplot.png"
    #plt.savefig(os.path.join(output_folder, plot_name), dpi=600, bbox_inches='tight')
    #plt.show()
    plt.close()


summary_df = pd.DataFrame(results_list)
#summary_df.to_csv(summary_output, sep="\t", index=False)

print(f"Analysis complete.")
print(f"Summary saved to: {summary_output}")
print(f"Plots saved to: {output_folder}/")

#%%
#bar plot for a single file, average accessibility for all genes in the gene set for a particular sample
file_path = "./pc9/pc9_GRCh38_Hallmark_EMT_counts.tab"
df = pd.read_csv(file_path, sep="\t")

data_cols = df.columns[3:]

replicate_averages = df[data_cols].mean()

plot_df = pd.DataFrame({
    'Replicate': replicate_averages.index,
    'Accessibility': replicate_averages.values
})

plot_df['Condition'] = plot_df['Replicate'].apply(
    lambda x: 'Control' if 'Ctrl' in x else 'DTP'
)
plt.rcParams['font.family'] = 'Arial'
plt.figure(figsize=(4, 5))
ax = sns.barplot(
    data=plot_df, 
    x='Condition', 
    y='Accessibility', 
    palette={'Control': '#3498db', 'DTP': '#e74c3c'},
    capsize=0.1,
    errorbar='sd' # Standard Deviation
)
sns.stripplot(
    data=plot_df, 
    x='Condition', 
    y='Accessibility', 
    color='black', 
    size=8, 
    alpha=0.6
)

plt.title(f"Mean Promoter Accessibility\n{file_path.split('/')[-1]}")
plt.ylabel("Average CPM Signal")
plt.xlabel("")
plt.tight_layout()
#plt.savefig("accessibility_barplot.png", dpi=300)
plt.show()
# %%
#box plot for acessibility for a gene

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 14

def plot_gene_across_all_cells(gene_name, cell_lines, base_path="all_gene_promoter_quantifications"):
    """
    Creates subplots showing accessibility for a gene in multiple cell lines.
    """
    # Create the figure with one subplot per cell line
    fig, axes = plt.subplots(1, len(cell_lines), figsize=(4 * len(cell_lines), 5), sharey=False)
    
    # Ensure axes is iterable even if only one cell line is provided
    if len(cell_lines) == 1: axes = [axes]

    for i, cell in enumerate(cell_lines):
        # Adjust path based on your folder structure
        # Assumes: folder_name/folder_name_GRCh38_All_RefSeq_Promoters_counts.tab
        file_path = os.path.join(base_path, f"{cell}_GRCh38_All_RefSeq_counts.tab")
        if not os.path.exists(file_path):
            print(f"Warning: File not found for {cell}: {file_path}")
            axes[i].set_title(f"{cell.upper()}\nMissing File")
            continue
            
        # Load the annotated data
        df = pd.read_csv(file_path, sep="\t")
        
        # Filter for your gene
        gene_data = df[df['gene_id'] == gene_name]
        
        if gene_data.empty:
            print(f"Warning: {gene_name} not found in {cell}")
            axes[i].set_title(f"{cell.upper()}\nGene Not Found")
            continue

        # Identify replicate columns (matching your DMSO/Ctrl/DTP/9291 naming)
        ctrl_cols = [c for c in df.columns if any(w in c for w in ['Ctrl', 'DMSO'])]
        dtp_cols = [c for c in df.columns if any(w in c for w in ['DTP', '9291'])]
        
        # Melt data for Seaborn (Long Format)
        melted = pd.melt(gene_data, id_vars=['gene_id'], value_vars=ctrl_cols + dtp_cols,
                         var_name='Sample', value_name='CPM')
        
        # Define condition labels
        melted['Condition'] = melted['Sample'].apply(
            lambda x: 'Control' if any(c in x for c in ['Ctrl', 'DMSO']) else 'DTP'
        )

        # Draw the Bar Plot
        sns.barplot(
            data=melted, x='Condition', y='CPM', hue='Condition', ax=axes[i],
            palette={'Control': "#839099", 'DTP': "#8f4adf"}, # Grey and Purple
            capsize=0.1, errorbar='sd', legend=False, edgecolor='grey'
        )
        
        # Overlay the individual replicates (dots)
        sns.stripplot(
            data=melted, x='Condition', y='CPM', color='black', ax=axes[i],
            size=7, alpha=0.6, jitter=True
        )

        # Formatting
        p_val = gene_data['P_Value'].values[0]
        axes[i].set_title(f"{cell.upper()}\np = {p_val:.4f}", fontsize=14, fontweight='bold')
        axes[i].set_ylabel("CPM Accessibility" if i == 0 else "")
        axes[i].set_xlabel("")
        sns.despine(ax=axes[i])

    plt.suptitle(f"Promoter Accessibility: {gene_name}", fontsize=18, y=1.05, fontweight='bold')
    plt.tight_layout()
    
    # Save the output
    save_name = f"{gene_name}_cross_cell_line_plot.png"
    #plt.savefig(save_name, dpi=600, bbox_inches='tight')
    print(f"Success! Plot saved as {save_name}")
    plt.show()

# --- RUN THE CODE ---
# List your cell line folders exactly as they appear in your directory


# %%
cell_lines = ["hcc827", "hcc2935", "h1975", "pc9"]
def plot_gene_across_all_cells(gene_name, cell_lines, base_path="all_gene_promoter_quantifications"):
    
    fig, axes = plt.subplots(1, len(cell_lines), figsize=(4 * len(cell_lines), 5), sharey=False)
    
    for i, cell in enumerate(cell_lines):
        file_path = os.path.join(base_path, f"{cell}_GRCh38_All_RefSeq_counts.tab")
        
        df = pd.read_csv(file_path, sep="\t")
        gene_data = df[df['gene_id'] == gene_name]
        
        if gene_data.empty:
            print(f"Warning: {gene_name} not found in {cell}")
            axes[i].set_title(f"{cell.upper()}\nGene Not Found")
            continue

        
        ctrl_cols = [c for c in df.columns if any(w in c for w in ['Ctrl1', 'Ctrl2', 'Ctrl3'])]
        dtp_cols = [c for c in df.columns if any(w in c for w in ['DTP1', 'DTP2', 'DTP3'])]
        
        melted = pd.melt(gene_data, id_vars=['gene_id'], value_vars=ctrl_cols + dtp_cols,
                            var_name='Sample', value_name='CPM')
        
    
        melted['Condition'] = melted['Sample'].apply(
            lambda x: 'Control' if any(c in x for c in ['Ctrl']) else 'DTP'
        )
        

        sns.barplot(
            data=melted, x='Condition', y='CPM', hue='Condition', ax=axes[i],
            palette={'Control': "#839099", 'DTP': "#1fa19b"}, # Grey and Purple
            capsize=0.1, errorbar='sd', legend=False, edgecolor='grey'
        )
        
        sns.stripplot(
            data=melted, x='Condition', y='CPM', color='black', ax=axes[i],
            size=7, alpha=0.6, jitter=True
        )

        p_val = gene_data['P_Value'].values[0]
        axes[i].set_title(f"{cell.upper()}\np = {p_val:.4f}", fontsize=14, fontweight='bold')
        axes[i].set_ylabel("CPM Accessibility" if i == 0 else "")
        axes[i].set_xlabel("")
        sns.despine(ax=axes[i])

    plt.suptitle(f"Promoter Accessibility: {gene_name}", fontsize=18, y=1.05, fontweight='bold')
    plt.tight_layout()
    plt.show()
# %%
cell_lines_to_plot = ["hcc827", "hcc2935", "h1975", "pc9"]
plot_gene_across_all_cells("AXL", cell_lines_to_plot)
# %%
