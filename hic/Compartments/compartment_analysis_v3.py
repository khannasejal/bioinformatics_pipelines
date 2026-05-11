#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import itertools
import random 
import numpy as np
import glob
import os
#%%
def read_comp_file_line_by_line(filename):
    set_of_comps = set()
    with open(filename, 'r') as f:
        for line in f:
            a = line[:-1].split('\t')
            a = '_'.join(a[1:4])
            set_of_comps.add(a)
    return set_of_comps
# %%
a_comp_fasr = read_comp_file_line_by_line('compartmentsA_FASR_copy.txt')
b_comp_fasr = read_comp_file_line_by_line('compartmentsB_FASR_copy.txt')
a_comp_mcf7 = read_comp_file_line_by_line('compartmentsA_MCF7_copy.txt')
b_comp_mcf7 = read_comp_file_line_by_line('compartmentsB_MCF7_copy.txt')
a_comp_tamr = read_comp_file_line_by_line('compartmentsA_TAMR_copy.txt')
b_comp_tamr = read_comp_file_line_by_line('compartmentsB_TAMR_copy.txt')

a_comp_fasr = [f'A_{i}_f' for i in a_comp_fasr]
b_comp_fasr = [f'B_{i}_f' for i in b_comp_fasr]
a_comp_mcf7 = [f'A_{i}_m' for i in a_comp_mcf7]
b_comp_mcf7 = [f'B_{i}_m' for i in b_comp_mcf7]
a_comp_tamr = [f'A_{i}_t' for i in a_comp_tamr]
b_comp_tamr = [f'B_{i}_t' for i in b_comp_tamr]

comp_fasr = a_comp_fasr + b_comp_fasr
comp_mcf7 = a_comp_mcf7 + b_comp_mcf7
comp_tamr = a_comp_tamr + b_comp_tamr

comp_mcf7_fasr = [comp_mcf7,comp_fasr]

df_TSS = pd.read_csv("TSS_hg38.txt", sep="\t") #read the ENSEMBL biomart Hg38 transcription start site information
df_TSS = df_TSS.dropna(subset=["Gene name"])
df_TSS = df_TSS[df_TSS["Chromosome/scaffold name"].isin([f"{i}" for i in range(1, 23)] + ["X"])] #filter
df_TSS["Chromosome/scaffold name"] = df_TSS["Chromosome/scaffold name"].astype(str)
df_TSS["Gene name"] = df_TSS["Gene name"].str.replace('_','-')
# %%
def return_tss_for_genes_within_comp(comp, tss):
    comp = comp.split('_')
    start, end = int(comp[2]), int(comp[3])
    # iterate over all genes that are on the same chromosome as the comp
    #genes = []
    # create a subset of the TSS dataframe for the chromosome of the comp
    tss_for_chr = tss[tss['Chromosome/scaffold name'] == comp[1][3:]]
    # create a TSS dataframe for all genes that are within the TAD
    tss_for_comp = tss_for_chr[(tss_for_chr['Transcription start site (TSS)'] >= start) & (tss_for_chr['Transcription start site (TSS)'] <= end)]
    return tss_for_comp
#format for genes within compartment:: compartmentA/B_chrnum_genename_tss_transcriptid
def return_unique_list_of_genes_within_comp(comp_list, tss, m_f_t):
    comp_num_genes = []
    for comp in comp_list:
        tss_for_comp = return_tss_for_genes_within_comp(comp, tss)
        # print number of unique genes in the comp
        comp_name = list(comp.split('_')[0]) * len(tss_for_comp)
        chr_num = list(tss_for_comp['Chromosome/scaffold name'])
        gene_names = list(tss_for_comp['Gene name'])
        tss_num = list(tss_for_comp['Transcription start site (TSS)'])
        transcript_id = list(tss_for_comp['Transcript stable ID version'])
        gene_chr = [f"{comp_nme}_{chr}_{gene}_{trnscrpt_id}_{transcript_ss}_{m_f_t}" for comp_nme, chr, gene, trnscrpt_id, transcript_ss in zip(comp_name,chr_num,gene_names,transcript_id,tss_num)]
        comp_num_genes.append(list(set(gene_chr)))
    return comp_num_genes  

comp_mcf7_genes = return_unique_list_of_genes_within_comp(comp_mcf7,df_TSS, "m")
comp_fasr_genes = return_unique_list_of_genes_within_comp(comp_fasr,df_TSS, "f")
comp_tamr_genes = return_unique_list_of_genes_within_comp(comp_tamr,df_TSS, "t")
# %%
#filter the comps which do not have any genes
def filter_comps_which_have_genes(list):
    filtered_list = []
    for sublist in list:
        if len(sublist) > 0:
            filtered_list.append(sublist)
    return filtered_list
comp_mcf7_genes = filter_comps_which_have_genes(comp_mcf7_genes)
comp_fasr_genes = filter_comps_which_have_genes(comp_fasr_genes)
comp_tamr_genes = filter_comps_which_have_genes(comp_tamr_genes)

#%%
#group the compartments chr wise
def group_sublists_by_chr(comp_list):
    grouped_sublists = defaultdict(list)
    for sublist in comp_list:
        chr = sublist[0].split('_')[1]
        grouped_sublists[chr].append(sublist)
    grouped_sublists_list = list(grouped_sublists.items())
    grouped_sublists_list.sort(key = lambda x:(x[0]!= 'X', float(x[0]) if x[0].isdigit() else float('inf')))
    return [comp_list for chr, comp_list in grouped_sublists_list]
grouped_comp_mcf7_genes = group_sublists_by_chr(comp_mcf7_genes)
grouped_comp_fasr_genes = group_sublists_by_chr(comp_fasr_genes)
grouped_comp_tamr_genes = group_sublists_by_chr(comp_tamr_genes)
#use random.choice() to replace the gene names and transcripts with other random genes, do the mappings and intersection to get the null value

#%%
#get the mappings of tss in one sample to the other
def return_tss_mappings(comp1,comp2):
    mappings = []
    for group1, group2 in zip(comp1,comp2):
        group2_dict = defaultdict(list)
        for sublist2 in group2:
            for ele2 in sublist2:
                parts = ele2.split('_')
                key = (parts[2], parts[3], parts[4])                
                group2_dict[key].append(ele2)
        group_result =[]
        for sublist1 in group1:
            sublist_result = []
            for ele1 in sublist1:
                parts = ele1.split('_')
                gene_name, tss, transcript_id = parts[2], parts[3], parts[4]
                key = (gene_name, tss, transcript_id)
                matching_elements = group2_dict.get(key, None)
                if matching_elements:
                    sublist_result.append(f"{ele1}:{':'.join(matching_elements)}")
                else:
                    sublist_result.append(f"{ele1}:None")
            group_result.append(sublist_result)
        mappings.append(group_result)
    return mappings

# %%
#the mappings are for all genes within the compartment, the genes which are not mapped to the other compartment, are represented by :None
mcf7_fasr_mappings = return_tss_mappings(grouped_comp_mcf7_genes,grouped_comp_fasr_genes)
mcf7_tamr_mappings = return_tss_mappings(grouped_comp_mcf7_genes,grouped_comp_tamr_genes)
fasr_tamr_mappings = return_tss_mappings(grouped_comp_fasr_genes,grouped_comp_tamr_genes)
tamr_fasr_mappings = return_tss_mappings(grouped_comp_tamr_genes,grouped_comp_fasr_genes)
#%%
def filter_elements_with_prefix(comp_map, str1 , str2):
    C1_C2_switch = []
    for group in comp_map:
        group_result =[]
        for sublist in group:
            filter_list = []
            for element in sublist:
                parts = element.split(':')
                if parts[0].startswith(str1) and parts[1].startswith(str2):
                    filter_list.append(element)
            if len(filter_list) > 0:
                group_result.append(filter_list)
        if len(group_result) > 0:
            C1_C2_switch.append(group_result)
    return C1_C2_switch
# %%
mcf7_fasr_B_B = filter_elements_with_prefix(mcf7_fasr_mappings, "B", "B")
mcf7_fasr_A_A = filter_elements_with_prefix(mcf7_fasr_mappings, "A", "A")
mcf7_fasr_A_B = filter_elements_with_prefix(mcf7_fasr_mappings, "A", "B")
mcf7_fasr_B_A = filter_elements_with_prefix(mcf7_fasr_mappings, "B", "A")

mcf7_tamr_B_B = filter_elements_with_prefix(mcf7_tamr_mappings, "B", "B")
mcf7_tamr_A_A = filter_elements_with_prefix(mcf7_tamr_mappings, "A", "A")
mcf7_tamr_A_B = filter_elements_with_prefix(mcf7_tamr_mappings, "A", "B")
mcf7_tamr_B_A = filter_elements_with_prefix(mcf7_tamr_mappings, "B", "A")

fasr_tamr_B_B = filter_elements_with_prefix(fasr_tamr_mappings, "B", "B")
fasr_tamr_A_A = filter_elements_with_prefix(fasr_tamr_mappings, "A", "A")
fasr_tamr_A_B = filter_elements_with_prefix(fasr_tamr_mappings, "A", "B")
fasr_tamr_B_A = filter_elements_with_prefix(fasr_tamr_mappings, "B", "A")

tamr_fasr_B_B = filter_elements_with_prefix(tamr_fasr_mappings, "B", "B")
tamr_fasr_A_A = filter_elements_with_prefix(tamr_fasr_mappings, "A", "A")
tamr_fasr_A_B = filter_elements_with_prefix(tamr_fasr_mappings, "A", "B")
tamr_fasr_B_A = filter_elements_with_prefix(tamr_fasr_mappings, "B", "A")
#%%
#replace the genes with other random genes
def replace_genes_comp_with_other_genes(mappings, df_TSS):
    random_comp = []
    for group in mappings:
        genes_list = []
        chr_num = group[0][0].split(':')[0].split('_')[1]
        tss_for_chr = df_TSS[df_TSS["Chromosome/scaffold name"].str.startswith(chr_num)]
        genes_in_tss = list(set((tss_for_chr["Gene name"].astype(str)).tolist()))        
        for sublist in group:
            for ele in sublist:
                if ele.split(':')[0].split('_')[2] not in genes_list:
                    genes_list.append(ele.split(':')[0].split('_')[2]) #list of unique genes in the compartment
        for gene in genes_list:
            replacement = random.choice(genes_in_tss)
            genes_in_tss.remove(replacement)
            random_comp.append(f"{replacement}")
    return random_comp

# %%
#get all the genes which change the compartment or are stable
def return_unique_genes_which_switch_comp_and_intersect_w_geneset(comp_list, geneset):
    genes = set()
    for group in comp_list:
        for sublist in group:
            for element in sublist:
                gene_name = element.split(':')[0].split('_')[2]
                genes.add(gene_name)
                intersect = genes.intersection(geneset)
    return intersect
# %%
er_geneset = set(pd.read_csv("er_corr_genes.txt", header=None)[0])
emt_geneset = set(pd.read_csv("emt_corr_genes.txt", header=None)[0])
def read_lines_and_convert_to_set(filename):
    with open(filename,'r') as f:
        first_line = f.readline().strip()
        geneset = set(first_line.split('\t'))
    return geneset
pemt_geneset = read_lines_and_convert_to_set("pemt_geneset.txt")
epithelial_geneset = read_lines_and_convert_to_set("epithelial_geneset.txt")
partial_geneset = read_lines_and_convert_to_set("partial_emt_geneset.txt")
mesenchymal_geneset = read_lines_and_convert_to_set("mesenchymal_geneset.txt")
tnbc_super_enhancer_geneset = read_lines_and_convert_to_set("tnbc_superenhancer_genes.txt")
#%%
folder_path = "C:\\sejal\\hic_breast_cancer\\genesets\\er_emt_mes"
files = glob.glob(os.path.join(folder_path, "*.txt"))
genesets ={}
for file in files:
    base_name = os.path.basename(file).replace(".txt","")
    var_name = base_name
    geneset = set(pd.read_csv(file,header =None)[0])
    
    genesets[base_name] = geneset
conditions = ['B_A','A_B','A_A','B_B']
prefixes = ['fasr_tamr']
#%%
def plot_null_distribution(mappings, df_TSS, geneset, filename):
    true_value = len(return_unique_genes_which_switch_comp_and_intersect_w_geneset(mappings, geneset))
    print(true_value)
    null_distribution = []
    for i in range(1000):
        if i % 10 == 0:
            print(i)
        replaced_gene_list = replace_genes_comp_with_other_genes(mappings, df_TSS)
        null_value = len(set(replaced_gene_list).intersection(geneset))
        null_distribution.append(null_value) 
    mean_val = np.mean(null_distribution)
    std_dev = np.std(null_distribution)  
    fig, ax = plt.subplots()
    plot_label = 'Title: {}\nMean: {:.4f}\nStd Dev: {:.4f}'.format(filename,mean_val,std_dev)
    #plot the histogram
    sns.histplot(null_distribution, kde=True, color='purple', ax=ax, label = plot_label)
    ax.axvline(x= true_value, color='black', label = true_value)
    ax.legend()
    fig.savefig(f'{filename}.jpg', format ='jpg', dpi = 300)
    return fig

# %%
for prefix in prefixes:
    for cond in conditions:
        for geneset_name, geneset in genesets.items():
            var_name = f"{prefix}_{cond.lower()}_{geneset_name.split('_')[0]}_{1000}"
            globals()[var_name] = plot_null_distribution(globals()[f"{prefix}_{cond}"], df_TSS,geneset,var_name)