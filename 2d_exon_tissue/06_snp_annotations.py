# -*- coding: utf-8 -*-
"""
Analysis of top SNPs identified by MTClass
in the 2D exon-tissue study

Created on Sat Apr  8 18:46:25 2023
@author: ronnieli
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pyarrow.feather as feather
import os, sys
import umap
from sklearn.preprocessing import StandardScaler
sys.path.append(r'G:/My Drive/Lab/lab_projects/gtex_exoneqtl/scripts/2d_exon_tissue')
from functions import get_expression_levels_2D, plot_9tissue_exon_expression, plot_umap_2D

#%%
def list_all_eqtls(data_dir=r'E:/lab_data'):
    
    ''' Returns a list of all GTEx eQTLs from all 54 tissues '''
    
    os.chdir(os.path.join(data_dir,'gtex_eqtl','GTEx_Analysis_v8_eQTL'))
    
    eqtls = []
    
    for file in os.listdir():
        if '.signif_variant_gene_pairs' in file:
            f = pd.read_csv(file, sep='\t', header=0)
            eqtls.append(set(f['variant_id']))
    
    all_eqtls = sorted(list(set.union(*eqtls)))
    return all_eqtls

def plot_top_non_eqtls(case):
    
    '''
    Plots top SNPs that are not eQTLs in any tissue.
    Input:
        - case = '9_tissue','brain','2D'
    '''
    
    # Load results
    if case == '9_tissue':
        os.chdir(r'G:/My Drive/Lab/lab_projects/gtex_eqtl/results/9_tissue/')
        mtclass = feather.read_feather('9_tissue_mtclass_ensemble_aggregate123.ftr')
        mphen = feather.read_feather('9_tissue_multiphen_binary.ftr')
        manova = feather.read_feather('9_tissue_manova.ftr')
        f1, mcc = 'f1_macro_median', 'mcc_median'
        threshold = 1.0
        
    elif case == 'brain':
        os.chdir(r'G:/My Drive/Lab/lab_projects/gtex_eqtl/results/brain_tissue/')
        mtclass = feather.read_feather('brain_mtclass_ensemble_aggregate123.ftr')
        mphen = feather.read_feather('brain_multiphen_binary.ftr')
        manova = feather.read_feather('brain_manova.ftr')
        f1, mcc = 'f1_macro_median', 'mcc_median'
        threshold = 0.80
        
    elif case == '2D':
        os.chdir(r'G:/My Drive/Lab/lab_projects/gtex_exoneqtl/results/2d_exon_tissue/')
        mtclass = pd.read_csv('mtclass_2d_NN.txt.gz', sep='\t', header=0)
        mphen = pd.read_csv('multiphen_2d.txt.gz', sep='\t', header=0)
        manova = pd.read_csv('manova_2d.txt.gz', sep='\t', header=0)
        f1, mcc = 'f1_macro', 'mcc'
        threshold = 1.0

    # Select top results from MTClass
    mtclass_top = mtclass[(mtclass[f1] >= threshold) & (mtclass[mcc] >= threshold)]
    top_snps = set(mtclass_top['variant'])
    
    # Use previous function to get all eQTLS from every tissue
    all_eqtls = list_all_eqtls()

    non_eqtls = list(top_snps - all_eqtls)
    print(f'There are {len(non_eqtls)} SNPs that are not eQTLs in any tissue')

    mtclass_non_eqtls = mtclass[mtclass['variant'].isin(non_eqtls)]

    for i, row in mtclass_non_eqtls.iterrows():
    
        gene_id = row['gene_id']
        gene_name = row['gene_name']
        snp = row['variant']

        mphen_pval = mphen[(mphen['gene_id']==gene_id) & (mphen['variant']==snp)]['pval'].values[0]
        manova_pval = manova[(manova['gene_id']==gene_id) & (manova['variant']==snp)]['pval'].values[0]
        
        print(f'{gene_name} \t {snp} \t MultiPhen p-val: {mphen_pval}\t MANOVA p-val: {manova_pval}')
        
        plot_9tissue_exon_expression(gene_id, snp)
        plot_umap_2D(gene_id, snp)

#%%

def reformat_genotypes_exoneqtl(data_dir=r'E:/lab_data', genotype='binary'):
    
    '''
    Reads in genotype data one chromosome at a time,
    then reformats the SNPs to look like GTEx variants,
    then merges with an annotation file to get gene names 
    in addition to Ensembl IDs. Finally, overwrites the
    original genotypes file.

    Input:
        - data_dir = data directory
        - genotype = 'binary' or 'additive'
    '''
    
    def reformat_snp(snp):
        
        snp_list = snp.split(':')
        snp = 'chr'+ '_'.join(snp_list) + '_b38'
        return snp
    
    gene_info = pd.read_csv(os.path.join(data_dir, 'gtex_exoneqtl', 'top_exon_gene_info.txt.gz'), sep='\t', header=0)
    gene_info = gene_info[['gene_id','gene_name']]

    for c in range(1,23):
    
        gt = pd.read_csv(os.path.join(data_dir, 'gtex_exoneqtl', 'genotypes', f'chr{c}_genotypes_10kb_{genotype}.txt.gz'), sep='\t', header=0)
        gt['ID'] = gt['ID'].apply(reformat_snp)
        
        gt_new = gene_info.merge(gt, on='gene_id')
        gt_new.to_csv(os.path.join(data_dir, 'gtex_exoneqtl', 'genotypes', f'chr{c}_genotypes_10kb_{genotype}.txt.gz'), sep='\t', header=True, index=False)
        print(f'Done for chr{c}')
        
def reformat_genotypes_eqtl(genotype):
        
    '''
    Read in genotype data one chromosome at a time,
    then merges with an annotation file to get gene IDs
    in addition to gene names. Finally, writes a new
    .txt.gz file.
    '''
    
    row_metadata = pd.read_pickle(r'E:/lab_data/gtex_eqtl/row_metadata.pkl').reset_index()
    
    for c in range(1,23):
    
        if genotype == 'binary':
            gt = feather.read_feather(rf'E:/lab_data/gtex_eqtl/genotypes/chr{c}_genotypes.ftr')
        elif genotype == 'additive':
            gt = feather.read_feather(rf'E:/lab_data/gtex_eqtl/genotypes/chr{c}_genotypes_additive.ftr')
        
        gt_merge = row_metadata.merge(gt, left_on='Description', right_on='gene')
        gt_merge = gt_merge.drop_duplicates(['gene','ID'])
        
        gt_merge = gt_merge.drop('Description', axis=1)
        
        cols = ['gene_id','gene_name','ID'] + list(gt_merge.columns[3:])
        gt_merge.columns = cols
        gt_merge = gt_merge.sort_values(by=['gene_name','ID'])
        
        if genotype == 'binary':
            gt_merge.to_csv(rf'E:/lab_data/gtex_eqtl/genotypes/chr{c}_genotypes_10kb_binary.txt.gz', index=False, header=True, sep='\t')
        elif genotype == 'additive':
            gt_merge.to_csv(rf'E:/lab_data/gtex_eqtl/genotypes/chr{c}_genotypes_10kb_additive.txt.gz', index=False, header=True, sep='\t')

        print(f'=== Done for chr{c} ===')