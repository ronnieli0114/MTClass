# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 16:16:36 2023

@author: ronnieli
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.stats import hypergeom

def disease_gene_overlap_isoqtl():
    
    os.chdir(r'G:/My Drive/Lab/lab_projects/isoqtl/results')
    mtclass = pd.read_csv('isoqtl_ensemble_aggregate.txt.gz', sep='\t', header=0)
    mphen_binary = pd.read_csv('isoqtl_multiphen_binary.txt.gz', sep='\t', header=0)
    mphen_add = pd.read_csv('isoqtl_multiphen_additive.txt.gz', sep='\t', header=0)
    manova = pd.read_csv('isoqtl_manova.txt.gz', sep='\t', header=0)
        
    mtclass['pair'] = mtclass['gene'] + '-' + mtclass['variant']
    mphen_binary['pair'] = mphen_binary['gene'] + '-' + mphen_binary['variant']
    mphen_add['pair'] = mphen_add['gene'] + '-' + mphen_add['variant']
    manova['pair'] = manova['gene'] + '-' + manova['variant']
    
    all_pairs = [set(data['pair']) for data in [mtclass, mphen_binary, mphen_add, manova]]
    common_pairs = set.intersection(*all_pairs)

    mtclass = mtclass[mtclass['pair'].isin(common_pairs)]
    mphen_binary = mphen_binary[mphen_binary['pair'].isin(common_pairs)]
    mphen_add = mphen_add[mphen_add['pair'].isin(common_pairs)]
    manova = manova[manova['pair'].isin(common_pairs)]

    print('Loaded MTClass, MultiPhen, MANOVA')
    print(f'Using {len(common_pairs)} total gene-SNP pairs for each method')
    
    n_top = 100
    print(f'Selecting the top {n_top} genes')
    gene_dict = {}
    
    for data, name in zip([mtclass, mphen_binary, mphen_add, manova],['MTClass','MultiPhen_binary','MultiPhen_additive','MANOVA']):
        
        top_genes = []
        if name == 'MTClass':
            data_sorted = data.sort_values(by='f1_macro_median', ascending=False)
        else:
            data_sorted = data.sort_values(by='pval', ascending=True)
        
        for gene in data_sorted['gene']:
            if len(set(top_genes)) == n_top:
                break
            else:
                top_genes.append(gene)
        top_genes = sorted(list(set(top_genes)))
        gene_dict[name] = top_genes
        
    data_dir = r'E:/lab_data'
    disease_list = ['AD','ADHD','ALS','autism','PD','schizophrenia']
    print(f'Calculating disease-gene overlap for {", ".join(disease_list)}')
    
    for disease in disease_list:
        
        disease_genes = pd.read_csv(os.path.join(data_dir, f'gtex_exoneqtl/disgenet_lists/{disease}_gene_list.tsv'), sep='\t')
        disease_gene_list = disease_genes['Gene'].tolist()
        background_genes = set(mtclass['gene'])
        disease_gene_list_10_isoforms = set(disease_gene_list) & background_genes
        print(f'There are {len(disease_gene_list_10_isoforms)} genes with at least 10 isoforms in {disease}')
        
        mtclass_num = set(gene_dict['MTClass']) & set(disease_gene_list_10_isoforms)
        
        mtclass_overlap = len(set(gene_dict['MTClass']) & set(disease_gene_list_10_isoforms)) / n_top
        mphen_bin_overlap = len(set(gene_dict['MultiPhen_binary']) & set(disease_gene_list_10_isoforms)) / n_top
        mphen_add_overlap = len(set(gene_dict['MultiPhen_additive']) & set(disease_gene_list_10_isoforms)) / n_top
        manova_overlap = len(set(gene_dict['MANOVA']) & set(disease_gene_list_10_isoforms)) / n_top
        
        # calculate p-value from hypergeometric distribution
        pop_size = len(background_genes)
        num_successes_pop = len(disease_gene_list_10_isoforms)
        sample_size = n_top
        X = len(mtclass_num)
        
        pval = hypergeom.sf(X-1, pop_size, num_successes_pop, sample_size)
        pval = np.round(pval, 3)
        print(f'p-value (MTClass) = {pval}')
        
        # plot
        fig, ax = plt.subplots(figsize=(6,6))
        titles = ['MTClass', 'MultiPhen\n(binary)','MultiPhen\n(additive)','MANOVA']
        overlaps = [mtclass_overlap, mphen_bin_overlap, mphen_add_overlap, manova_overlap]
        colors = ['blue','orange','yellow','green']
        
        b = ax.bar(titles, overlaps, color=colors, alpha=0.7)
        ax.bar_label(b)
        ax.set_xticks(np.arange(len(titles)), titles, fontsize=16)
        ax.set_xlabel('Method', fontsize=16)
        ax.set_ylabel('Percentage overlap with disease genes', fontsize=14)
        ax.set_title(f'{disease}\nisoQTL', fontsize=16)
        ax.annotate(text=f'p-value (MTClass) = {pval}', xy=(-0.15, -0.18), xycoords='axes fraction', fontsize=14)
        plt.tight_layout()
        plt.savefig(rf'C:/Users/ronni/OneDrive/Desktop/{disease}_isoqtl.png', dpi=400)
        plt.show()
        
        
        
        
    