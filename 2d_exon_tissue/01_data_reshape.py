# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 20:27:48 2023

@author: ronnieli
"""

import pyarrow.feather as feather
import pandas as pd
import numpy as np
import os

def load_expression_level_data(data_dir):
    return feather.read_feather(os.path.join(data_dir, 'gtex_exoneqtl', '9_tissue_exon_expression.ftr'))

def load_attributes(data_dir):
    return pd.read_csv(os.path.join(data_dir,'gtex_exoneqtl','v8_sample_attributes.txt'), sep='\t', header=0)

def load_gene_annotations(data_dir):
    return pd.read_csv(os.path.join(data_dir, 'gtex_exoneqtl', 'top_exon_gene_info.txt.gz'), sep='\t', header=0)

### Define function to reshape expression level data

def prepare_X(gene_id, exp, attr):
    
    gene_exp = exp[exp['gene_id'] == gene_id]
    gene_exp = gene_exp.drop(['Description','Name'], axis=1).set_index('gene_id')
    
    gene_exp = gene_exp.transpose().reset_index()
    gene_exp['donor'] = gene_exp['index'].apply(lambda x: '-'.join(x.split('-')[:2]))
    
    donor_dict = {}
    
    for donor in sorted(list(set(gene_exp['donor']))):
        
        donor_exp = gene_exp[gene_exp['donor'] == donor]
        
        if donor_exp.shape[0] != 9: # if there are not 9 tissues, skip donor
            continue
        
        donor_exp = donor_exp.merge(attr[['SAMPID','SMTSD']], left_on='index', right_on='SAMPID')
        donor_exp = donor_exp.drop(['index','donor','SAMPID'], axis=1).set_index('SMTSD').transpose()
        
        # sort tissues by hierarchical clustering
        sorted_tissues = ['Whole Blood', 'Muscle - Skeletal', 'Skin - Not Sun Exposed (Suprapubic)',
                          'Skin - Sun Exposed (Lower leg)','Lung','Artery - Tibial','Nerve - Tibial',
                          'Adipose - Subcutaneous','Thyroid']
        donor_exp = donor_exp.loc[:,sorted_tissues]
        donor_arr = donor_exp.values
        donor_arr = donor_arr.flatten()
        
        donor_dict[donor] = donor_arr
        
    donor_list = list(donor_dict.keys())
    X = np.vstack(list(donor_dict.values()))
    X_df = pd.DataFrame(X, index=donor_list)
    X_df['gene_id'] = gene_id
    new_cols = [X_df.columns[-1]] + list(X_df.columns[:-1])
    X_df = X_df.loc[:,new_cols]
    X_df = X_df.reset_index()
    X_df = X_df.rename(columns={'index':'donor'})
    
    return X_df

data_dir = r'/projects/compbio/users/yli459'
data_dir = r'E:/lab_data'

exp = load_expression_level_data(data_dir)
attr = load_attributes(data_dir)
annot = load_gene_annotations(data_dir)

all_results = []

for c in range(1,23):
    
    chrom_data = []
    
    exp_genes = sorted(list(set(exp['gene_id'])))
    annot_genes = annot[annot['gene_id'].isin(exp_genes)]
    annot_genes = annot_genes[annot_genes['chr'] == c]
    chrom_genes = sorted(list(annot_genes['gene_id']))
    print(f'There are {len(chrom_genes)} genes in chr{c}')
    
    for i, gene in enumerate(chrom_genes):
        
        if i%100==0:
            print(f'Now on {i}th gene {gene}')
        
        gene_data = prepare_X(gene, exp, attr)
        chrom_data.append(gene_data)
    
    chrom_df = pd.concat(chrom_data)
    chrom_df.to_csv(os.path.join(data_dir, 'gtex_exoneqtl', '2d_expression', 'chr{c}_2d_expression.txt.gz'), sep='\t', index=False)
    print(f'Saved chr{c} 2D expression to file')
        