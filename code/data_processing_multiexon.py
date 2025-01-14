# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 16:45:33 2022

@author: ronnieli
"""

import pandas as pd
import numpy as np
import pyarrow.feather as feather
import os

data_dir = r'F:/lab_data/mtclass_eqtl/multi-exon'

# load data into memory
top_counts = feather.read_feather(os.path.join(data_dir, 'top_exon_TPM.ftr'))
attributes = pd.read_csv(os.path.join(data_dir, 'v8_sample_attributes.txt'), sep='\t', header=0)

# process attributes data
attributes = attributes[['SAMPID','SMTSD']]
attributes['donor'] = attributes['SAMPID'].apply(lambda x: '-'.join(x.split('-')[:2]))
attributes.columns = ['full_id','tissue','donor']
attributes['tissue'] = attributes['tissue'].apply(lambda x: x.replace(' - ', '_').replace(' ','_').replace('(','').replace(')',''))

# filter top counts by mean variance across samples and exons
gene_dict = {'gene':[], 'mean_variance':[]}

for gene in sorted(list(set(top_counts['gene_id']))):
    
    gene_data = top_counts[top_counts['gene_id']==gene]
    gene_counts = gene_data.set_index('gene_id').drop(['Description','Name'], axis=1)
    mean_var = np.mean(np.var(gene_counts, axis=1))
    print((gene, mean_var))
    gene_dict['gene'].append(gene)
    gene_dict['mean_variance'].append(mean_var)

var_df = pd.DataFrame(gene_dict)
var_df.sort_values(by='mean_variance', ascending=False, inplace=True)

# select top half based on variance
sel_df = var_df.iloc[:int(var_df.shape[0]/2), :]
sel_genes = sorted(list(sel_df['gene']))

high_var_counts = top_counts[top_counts['gene_id'].isin(sel_genes)]

gene_info = pd.read_table(r'F:/lab_data/mtclass_eqtl/multi-exon/top_exon_gene_info.txt.gz', sep='\t', header=0)
gene_info = gene_info[gene_info['gene_id'].isin(sel_genes)]
gene_info = gene_info.dropna(how='any', axis=0)

with open(r'F:/lab_data/mtclass_eqtl/multi-exon/high_variance_genes.txt','w') as f:
    for g in sorted(list(set(gene_info['gene_name']))):
        f.write(g + '\n')

#%% Define function to extract genes with n or more exons

def get_top_counts(n_exons, save_file=False):
    
    data_dir = r'F:/lab_data/mtclass_eqtl/multi-exon'
    exon_counts = pd.read_parquet(os.path.join(data_dir, 'GTEx_v8_exon_reads.parquet'))
    top_genes = []
    for gene in sorted(list(set(exon_counts['Description']))):
        
        data = exon_counts[exon_counts['Description']==gene]
        if data.shape[0] >= n_exons:
            top_genes.append(gene)
            
    top_exon_counts = exon_counts[exon_counts['Description'].isin(top_genes)]
    if save_file:
        feather.write_feather(top_exon_counts, os.path.join(data_dir, 'top_exon_TPM.ftr'))
    return top_exon_counts

## extract TPM measures for all GTEx tissues

for tissue in sorted(set(attributes['tissue'])):
    
    results_dir = r'F:/lab_data/mtclass_eqtl/multi-exon/single_tissue_expression'
    
    # get full IDs
    full_ids = attributes[attributes['tissue']==tissue]['full_id']
    tiss_counts = top_counts.set_index('Description').filter(full_ids)
    tiss_counts.columns = list(map(lambda x: '-'.join(x.split('-')[:2]), tiss_counts.columns))
    
    count_dict = {'gene':[], 'count':[]}
    
    for gene in sorted(list(set(tiss_counts.index))):
        
        gene_counts = tiss_counts[tiss_counts.index==gene]
        gene_counts = gene_counts.T
        gene_counts.index.rename('donor', inplace=True)
        gene_counts.columns = ['exon_' + str(i+1) for i in range(len(gene_counts.columns))] # exon_1, exon_2, etc.
        gene_counts = gene_counts.reset_index()
        gene_counts.sort_values(by='donor', inplace=True)
        gene_list = [gene] * gene_counts.shape[0]
    
        count_dict['gene'].extend(gene_list)
        count_dict['count'].append(gene_counts)
        
    final_counts = pd.concat(count_dict['count'], axis=0)
    final_counts['gene'] = count_dict['gene']
    cols = ['gene','donor'] + final_counts.iloc[:,1:-1].columns.tolist()
    final_counts = final_counts.loc[:,cols]
    feather.write_feather(final_counts, os.path.join(results_dir, f'{tissue}.ftr'))
    

'''
Process tissue counts even further. For each gene, determine
whether most exon counts have 0 expression level. If so, exclude that gene.
'''

data_dir = r'F:/lab_data/mtclass_eqtl/multi-exon/single_tissue_expression'
tissue_list = [tissue.split('.ftr')[0] for tissue in os.listdir(data_dir) if '.txt' not in tissue]
brain_tissues = [tissue for tissue in tissue_list if 'Brain_' in tissue]

for tissue in brain_tissues:
    print(f'Now on tissue {tissue}...')
    counts = feather.read_feather(os.path.join(data_dir, f'{tissue}.ftr'))
    
    tissue_dict = {'gene':[], 'prop_zero':[]}
    gene_list = sorted(list(set(counts['gene'])))
    
    for gene in gene_list:
        
        counts_gene = counts[counts['gene']==gene].drop(['gene','donor'], axis=1)
        counts_gene = counts_gene.dropna(axis=1)
        counts_gene = counts_gene.to_numpy()
        
        # count proportion of zeros
        num_donors = counts_gene.shape[0]
        prop_zero = np.count_nonzero(counts_gene==0) / np.size(counts_gene)
        prop_zero = np.around(prop_zero, 4)
        # print(f'Gene {gene} \t Proportion of zeros = {prop_zero}')
        
        # append
        tissue_dict['gene'].append(gene)
        tissue_dict['prop_zero'].append(prop_zero)
    
    tissue_results = pd.DataFrame(tissue_dict)
    tissue_results = tissue_results.sort_values(by='prop_zero', ascending=False)
    print(tissue_results.head())
    
    # Get bottom 25 percentile
    # because we want genes with LOW proportion of zeros
    bottom_percentile = np.percentile(tissue_results['prop_zero'].values, 25)
    top_tissue = tissue_results[tissue_results['prop_zero'] <= bottom_percentile]
    top_tissue = top_tissue.sort_values(by='prop_zero', ascending=True)
    print(top_tissue.head())
    
    num_genes_final = len(top_tissue)
    print('Final number of genes selected: ', num_genes_final)
    
    gene_list = top_tissue['gene'].tolist()
    with open(os.path.join(counts_dir, f'{tissue}_gene_list.txt'), 'w') as f:
        f.write('### Top genes to include in analysis\n')
        f.write('### These are genes with low numbers of zeros in the exon counts\n')
        f.write('gene_name\n')
        for gene in sorted(gene_list):
            f.write(gene)
            f.write('\n')