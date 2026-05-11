#%%
import pyranges as pr
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import statsmodels.api as sm

# %%
#read the segments file from the folder as dataframe and create variables
#filter the enriched regions (for a particular state)
root_dir = "D:\\sejal\\chromhmm_analysis\\chromhmm_superenhancer_analysis_codes\\input_segment_files_for_se\\GSE118716_GSE57498"
col_headers = ['chrom','start','end','state']
#cell_lines = ["76NF2V","AU565","BT549","HCC1937","HCC1954","MB231","MB361","MB436","MCF10A","MDAMB468","SKBR3","SUM149","T47D","UACC812","ZR751"]
cell_lines = ['TAMR','FASR','MCF7']
#get the enriched regions for the given cell line
def get_enriched_regions(cell_line):
    for file_name in os.listdir(root_dir):
        if file_name.startswith(cell_line):
            segments_file_path = os.path.join(root_dir,file_name)
            break
    segments = pd.read_csv(segments_file_path, sep='\t',header = None)
    segments.columns = col_headers
    #mention the state which corresponds to the histone mark enriched regions
    enriched_regions = segments[segments['state'] == 'E2']
    return enriched_regions


df_TSS = pd.read_csv("D:/sejal_c/sejal/hic_breast_cancer/TSS_hg38.txt", sep="\t")
df_TSS = df_TSS.dropna(subset=["Gene name"])
df_TSS = df_TSS[df_TSS["Chromosome/scaffold name"].isin([f"{i}" for i in range(1, 23)] + ["X"])] #filter
df_TSS["Chromosome/scaffold name"] = df_TSS["Chromosome/scaffold name"].astype(str)
df_TSS["Gene name"] = df_TSS["Gene name"].str.replace('_','-')
df_TSS["gene_transcript"] = df_TSS["Gene name"] + '_' + df_TSS["Transcript stable ID version"]
df_TSS['Chromosome/scaffold name'] = 'chr' + df_TSS['Chromosome/scaffold name'].astype(str)


#do the overlap of k27ac enriched regions with the tss regions of all genes (+-2.5kbp from tss) and remove the regions which lie 100% within the promoter region  
# For df1 (assuming it already has Chromosome, Start, and End)
def overlap_function_for_enriched_regions_and_tss_regions(enriched_regions):
    pr1 = pr.PyRanges(enriched_regions.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"}))    
    # For df2 (assuming it has Chromosome and TSS)
    df_TSS["Start"] = df_TSS["Transcription start site (TSS)"] - 2500  # Extend TSS by -2000 bp
    df_TSS["End"] = df_TSS["Transcription start site (TSS)"] + 2500   # Extend TSS by +2000 bp
    pr2 = pr.PyRanges(df_TSS.rename(columns={"Chromosome/scaffold name": "Chromosome", "Start": "Start", "End": "End"}))
    # Now create PyRanges object
    # Find overlaps between the TSS regions (pr2) and active promoter regions (pr1)
    overlap = pr1.join(pr2)
    overlap = overlap.df
    # Calculate overlap percentage with respect to the k27ac enriched regions (percentage of k27ac enriched regions which overlap with the tss/promoter regions)
    overlap['overlap_length'] = overlap[['End', 'End_b']].min(axis=1) - overlap[['Start', 'Start_b']].max(axis=1)
    overlap['k27ac_enriched_regions'] = overlap['End'] - overlap['Start']
    overlap['overlap_percentage'] = (overlap['overlap_length'] / overlap['k27ac_enriched_regions']) * 100   
    #remove the enriched regions which lie completely within the tss region, calculate enhancer size and put rank based on size
    filtered_overlap = overlap[overlap['overlap_percentage'] > 0]
    k27ac_regions_within_tss_region = filtered_overlap[['Chromosome','Start','End','state']]
    #consider only those rows/enriched regions which are not present in tss promoter region
    merged_df = pd.merge(pr1.df, k27ac_regions_within_tss_region, on=list(pr1.columns), how='left', indicator = True)
    merged_df = merged_df[merged_df['_merge'] == 'left_only']
    #calculate enhancer size
    merged_df['size'] = merged_df['End'] - merged_df['Start']
    #sort and assign rank
    merged_df = merged_df.sort_values(by = 'size', ascending=True).reset_index(drop=True)
    merged_df['rank'] = merged_df.index + 1
    #return ranked dataframe for enhancers
    return merged_df


def do_lowess_regression_and_plot_enhancer_size_versus_rank_graph(enhancers):
    lowess = sm.nonparametric.lowess(enhancers['size'],enhancers['rank'],frac =0.1)
    x = lowess[:,0]
    y = lowess[:,1]
    slope = np.gradient(y,x)
    inflection_points = np.where(np.isclose(slope,1,atol = 0.01))[0]
    if inflection_points.size>0:
        inflection_x = x[inflection_points[0]]
        print(inflection_x)
        inflection_y = y[inflection_points[0]]
        print(inflection_y)
    else:
        inflection_x = None
    plt.figure(figsize = (10,10))
    plt.plot(enhancers['rank'],enhancers['size'], marker='o', label = 'Data points')
    plt.plot(x,y,color = 'red', label = 'LOESS fit')
    if inflection_x is not None:
        plt.axvline(x=inflection_x, color = 'green', linestyle = '--' )
        plt.text(inflection_x + 0.5, inflection_y, f'Inflection Point: {inflection_x:.2f}', color = 'green')
    plt.xlabel('Enhancer Rank')
    plt.ylabel('Enhancer size')
    #plt.yscale('log')
    plt.legend(loc='best',fontsize =10)
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    return inflection_x

#superenhancer-gene map
#overlap the se region with the gene region (100kbp from tss)
#df_TSS start and end columns will be replaced with +-100kbp range

def get_se_gene_map(enriched_regions, cell_line):
    enhancers =   overlap_function_for_enriched_regions_and_tss_regions(enriched_regions)
    inflection_x_value = do_lowess_regression_and_plot_enhancer_size_versus_rank_graph(enhancers)
    superenhancers = enhancers[enhancers['rank'] > inflection_x_value]
    #return superenhancers
    superenhancers = superenhancers.reset_index(drop = True)
    superenhancers['name'] = 'enh_' + (superenhancers.index + 1).astype(str)
    se = superenhancers[['Chromosome','Start','End','name']].drop_duplicates()
    output_folder = 'D:/sejal/chromhmm_analysis/chromhmm_superenhancer_analysis_codes/all_se/tamr_fasr_mcf7'
    se.to_csv(f"{output_folder}/{cell_line}_all_se.bed", sep ='\t', index = False)
    return se
    
    pr1 = pr.PyRanges(superenhancers)
    
    #extend the tss region by 100kbp
    df_TSS["Start"] = df_TSS["Transcription start site (TSS)"] - 100000  
    df_TSS["End"] = df_TSS["Transcription start site (TSS)"] + 100000
    pr2 = pr.PyRanges(df_TSS.rename(columns={"Chromosome/scaffold name": "Chromosome", "Start": "Start", "End": "End"}))

    #Now create PyRanges object
    #Find overlaps between the TSS regions (pr2) and active promoter regions (pr1)
    overlap = pr1.join(pr2)
    overlap = overlap.df

    #Calculate overlap percentage with respect to the k27ac enriched regions
    overlap['overlap_length'] = overlap[['End', 'End_b']].min(axis=1) - overlap[['Start', 'Start_b']].max(axis=1)
    overlap['superenhancer_region'] = overlap['End'] - overlap['Start']
    overlap['overlap_percentage'] = (overlap['overlap_length'] / overlap['superenhancer_region']) * 100

    #if 50 percent of se region lies within 100kbp of tss region, assign the gene to superenhancer
    filtered_overlap = overlap[overlap['overlap_percentage'] >= 50]

    #to get se and transcript map, add transcript id version column below
    se_gene_map  = filtered_overlap[['Chromosome','Start','End','name','Gene name']].drop_duplicates()
    #remove the header=False to get the maps wo header
    '''se_gene_map.to_csv(f"{cell_line}_se_gene_map.bed", sep ='\t', index = False)'''
    return se_gene_map
#%%
#number_of_gene_counts = se_gene_map['Gene name'].value_counts()
#number_of_enhancer_counts = se_gene_map['name'].value_counts()
for cell_line in cell_lines:
    enriched_regions = get_enriched_regions(cell_line)
    se = get_se_gene_map(enriched_regions, cell_line)

# %%
#read the se gene map and print the number of enhancers
map_dir = "D:/sejal/chromhmm_analysis/chromhmm_superenhancer_analysis_codes"
def get_the_maps(cell_line):
    for file_name in os.listdir(map_dir):
        if file_name.startswith(cell_line):
            map_file_path = os.path.join(map_dir,file_name)
            break
    map = pd.read_csv(map_file_path, sep='\t')     
    return map

#%%
for cell_line in cell_lines:
    maps = get_the_maps(cell_line)
    superenhancers = maps[['Chromosome','Start','End','name']].drop_duplicates()
    print(cell_line + ':' + str(len(superenhancers)))

# %%
