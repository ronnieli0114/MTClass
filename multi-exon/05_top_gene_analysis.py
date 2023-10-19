# -*- coding: utf-8 -*-
"""
Based on the top genes identified by MTClass and top eGenes for single-tissue
approaches, determine which method produces more GWAS-related genes
for neurological and psychiatric related traits, specifically AD and ALS.

Created on Thu Jan  5 22:11:55 2023
@author: ronnieli
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import os
from scipy.stats import hypergeom

def disease_association_single_tissue(tissue, disease):
    
    '''
    Determine whether the top genes identified by MTClass have disease-associations.
    Compare to single-tissue eQTL and sQTL approach.
    '''
    
    metric = 'f1_macro_median'
    n_top_genes = 100
    
    os.chdir(r'G:/My Drive/Lab/lab_projects/gtex_exoneqtl/results')
    data_dir = r'E:/lab_data'
    
    # Analyze results for MTClass
    mtclass = pd.read_csv(f'{tissue}_ensemble_aggregate.txt.gz', sep='\t', header=0)
    
    mtclass_top = mtclass.sort_values(by=metric, ascending=False)
    mtclass_genes = []
    for gene in mtclass_top['gene']:
        if len(set(mtclass_genes)) == n_top_genes:
            break
        else:
            mtclass_genes.append(gene)
    mtclass_genes = list(set(mtclass_genes))
    
    row_metadata = pd.read_pickle(os.path.join(data_dir, 'gtex_eqtl/row_metadata.pkl'))
    genes_10_exons = pd.read_csv(os.path.join(data_dir, 'gtex_exoneqtl/gene_list_10_exons.txt'), header=None)[0].tolist()
    
    # eQTL single tissue
    eqtl = pd.read_csv(os.path.join(data_dir, f'gtex_eqtl/GTEx_Analysis_v8_eQTL/{tissue}.v8.signif_variant_gene_pairs.txt.gz'), sep='\t', header=0)
    eqtl = eqtl.merge(row_metadata, left_on='gene_id', right_index=True)
    eqtl_10_exons = eqtl[eqtl['Description'].isin(genes_10_exons)]
    eqtl_10_exons = eqtl_10_exons.sort_values(by='pval_nominal', ascending=True)
    
    eqtl_top_genes = []
    for gene in eqtl_10_exons['Description']:
        if len(set(eqtl_top_genes)) == n_top_genes:
            break
        else:
            eqtl_top_genes.append(gene)
    eqtl_top_genes = list(set(eqtl_top_genes))
        
    # sQTL single tissue
    sqtl = pd.read_csv(os.path.join(data_dir, f'gtex_eqtl/GTEx_Analysis_v8_sQTL/{tissue}.v8.sqtl_signifpairs.txt.gz'), sep='\t', header=0)
    sqtl['gene_id'] = sqtl['phenotype_id'].apply(lambda x: x.split(':')[-1])
    sqtl = sqtl.merge(row_metadata, left_on='gene_id', right_index=True)
    sqtl_10_exons = sqtl[sqtl['Description'].isin(genes_10_exons)]
    sqtl_10_exons = sqtl_10_exons.sort_values(by='pval_nominal', ascending=True)
    
    sqtl_top_genes = []
    for gene in sqtl_10_exons['Description']:
        if len(set(sqtl_top_genes)) == n_top_genes:
            break
        else:
            sqtl_top_genes.append(gene)
    sqtl_top_genes = list(set(sqtl_top_genes))
    
    # Load disease genes
    disease_genes = pd.read_csv(os.path.join(data_dir, f'gtex_exoneqtl/disgenet_lists/{disease}_gene_list.tsv'), sep='\t')
    disease_gene_list = disease_genes['Gene'].tolist()
    disease_gene_list = list(set(disease_gene_list) & set(genes_10_exons)) # make sure disease genes have 10 exons

    # Determine overlaps
    mtclass_overlap = set(mtclass_genes) & set(disease_gene_list)
    mtclass_percentage = round((len(mtclass_overlap) * 100/ len(mtclass_genes)), 3)
    
    eqtl_overlap = set(eqtl_top_genes) & set(disease_gene_list)
    eqtl_percentage = round((len(eqtl_overlap) * 100 / len(eqtl_top_genes)), 3)

    sqtl_overlap = set(sqtl_top_genes) & set(disease_gene_list)
    sqtl_percentage = round((len(sqtl_overlap) * 100 / len(sqtl_top_genes)), 3)
    print(f'{tissue}\t {disease}\t MTClass: {mtclass_percentage} \t eQTL: {eqtl_percentage} \t sQTL: {sqtl_percentage}')
    
    # calculate p-value from hypergeometric distribution
    pop_size = len(genes_10_exons)
    num_successes_pop = len(disease_gene_list)
    sample_size = n_top_genes
    X = len(mtclass_overlap)
    
    pval = hypergeom.sf(X-1, pop_size, num_successes_pop, sample_size)
    pval = np.round(pval, 3)
    print(f'p-value (MTClass) = {pval}')

    # plot bar graph
    fig, ax = plt.subplots(figsize=(6,6))
    b = ax.bar(['MTClass','eQTL','sQTL'], [mtclass_percentage, eqtl_percentage, sqtl_percentage], 
               color=('blue','pink','yellow'), alpha=0.7)
    ax.bar_label(b)
    ax.set_xticks(np.arange(3), ('MTClass','eQTL','sQTL'), fontsize=20)
    ax.set_xlabel('Method', fontsize=16)
    ax.set_ylabel('Overlap with disease genes', fontsize=16)
    ax.set_title(f'{disease}\n{tissue}\n{n_top_genes} top genes', fontsize=20)
    ax.annotate(text=f'p-value (MTClass) = {pval}', xy=(-0.15, -0.25), xycoords='axes fraction', fontsize=16)
    plt.tight_layout()
    plt.savefig(rf'C:/Users/ronni/OneDrive/Desktop/{disease}_{tissue}_singletissue.png', dpi=400)
    plt.show()
    
    ### Plot three-way Venn Diagrams of MTClass, single-tissue eQTL, single-tissue sQTL approach
    mtclass_set = set(mtclass_genes)
    eqtl_set = set(eqtl_top_genes)
    sqtl_set = set(sqtl_top_genes)
    
    plt.figure(figsize=(5,5))
    venn3([mtclass_set, eqtl_set, sqtl_set], ('MTClass', 'eQTL', 'sQTL'), set_colors=('blue','pink','yellow'), alpha=0.5)
    plt.title(f'{disease}\n{tissue}\n{n_top_genes} top genes', fontsize=20)
    plt.tight_layout()
    plt.savefig(rf'C:/Users/ronni/OneDrive/Desktop/{tissue}_venn3_single_tissue.png', dpi=400)
    plt.show()

def disease_association_multivariate(tissue, disease):
    
    '''
    Determine whether the top genes identified by MTClass have disease-associations.
    Compare to other multivariate approaches (MultiPhen and MANOVA).
    '''
    
    os.chdir(r'G:/My Drive/Lab/lab_projects/gtex_exoneqtl/results')
    data_dir = r'E:/lab_data'

    metric = 'f1_macro_median'
    n_top_genes = 100

    # Make sure they are testing on same gene-SNP pairs
    mtclass = pd.read_csv(f'{tissue}_ensemble_aggregate.txt.gz', sep='\t', header=0)
    mphen = pd.read_csv(f'{tissue}_multiphen_binary.txt.gz', sep='\t', header=0)
    manova = pd.read_csv(f'{tissue}_manova.txt.gz', sep='\t', header=0)

    mtclass['pair'] = mtclass['gene'] + '-' + mtclass['variant']
    mphen['pair'] = mphen['gene'] + '-' + mphen['variant']
    manova['pair'] = manova['gene'] + '-' + manova['variant']
    
    all_pairs = [set(data['pair']) for data in [mtclass, mphen, manova]]
    common_pairs = set.intersection(*all_pairs)
    
    mtclass = mtclass[mtclass['pair'].isin(common_pairs)]
    mphen = mphen[mphen['pair'].isin(common_pairs)]
    manova = manova[manova['pair'].isin(common_pairs)]
    
    print(f'There are {len(common_pairs)} gene-SNP pairs for {tissue}')
    print(f'Selecting {n_top_genes} top genes from each method...')

    # Determine top genes
    mtclass_top = mtclass.sort_values(by=metric, ascending=False)
    mtclass_genes = []
    for gene in mtclass_top['gene']:
        if len(set(mtclass_genes)) == n_top_genes:
            break
        else:
            mtclass_genes.append(gene)
    mtclass_genes = list(set(mtclass_genes))
    
    mphen_top = mphen.sort_values(by='pval', ascending=True)
    mphen_top_genes = []
    for gene in mphen_top['gene']:
        if len(set(mphen_top_genes)) == n_top_genes:
            break
        else:
            mphen_top_genes.append(gene)
    mphen_top_genes = list(set(mphen_top_genes))
    
    manova_top = manova.sort_values(by='pval', ascending=True)
    manova_top_genes = []
    for gene in manova_top['gene']:
        if len(set(manova_top_genes)) == n_top_genes:
            break
        else:
            manova_top_genes.append(gene)
    manova_top_genes = list(set(manova_top_genes))
    
    # Get disease gene list (more than 10 exons)
    print('Loading disease gene list...')
    disease_genes = pd.read_csv(os.path.join(data_dir, f'gtex_exoneqtl/disgenet_lists/{disease}_gene_list.tsv'), sep='\t', header=0)
    disease_gene_list = disease_genes['Gene'].tolist()
    genes_10_exons = pd.read_csv(os.path.join(data_dir, 'gtex_exoneqtl/gene_list_10_exons.txt'), header=None)[0].tolist()
    disease_gene_list = list(set(disease_gene_list) & set(genes_10_exons)) # make sure disease genes have 10 exons
    print('For {}, there are {} associated genes with at least 10 exons'.format(disease, len(disease_gene_list)))
    
    # Determine overlaps with disease genes
    mtclass_overlap = set(mtclass_genes) & set(disease_gene_list)
    mtclass_percentage = round((len(mtclass_overlap) * 100/ len(mtclass_genes)), 3)
    
    mphen_overlap = set(mphen_top_genes) & set(disease_gene_list)
    mphen_percentage = round((len(mphen_overlap) * 100 / len(mphen_top_genes)), 3)
    
    manova_overlap = set(manova_top_genes) & set(disease_gene_list)
    manova_percentage = round((len(manova_overlap) * 100 / len(manova_top_genes)), 3)
    
    # calculate p-value from hypergeometric distribution
    pop_size = len(genes_10_exons)
    num_successes_pop = len(disease_gene_list)
    sample_size = n_top_genes
    X = len(mtclass_overlap)
    
    pval = hypergeom.sf(X-1, pop_size, num_successes_pop, sample_size)
    pval = np.round(pval, 3)
    print(f'p-value (MTClass) = {pval}')
    
    ### Plot bar plot of disease-gene overlap
    fig, ax = plt.subplots(figsize=(6,6))
    b = ax.bar(['MTClass','MultiPhen (binary)','MANOVA'], [mtclass_percentage, mphen_percentage, manova_percentage], 
               color=('blue','orange','green'), alpha=0.7)
    # b = ax.bar(['MTClass','MultiPhen'], [mtclass_percentage, mphen_percentage])
    ax.bar_label(b)
    ax.set_xticks(np.arange(3), ('MTClass','MultiPhen\n(binary)','MANOVA'), fontsize=20)
    ax.set_xlabel('Method', fontsize=16)
    ax.set_ylabel('Overlap with disease genes', fontsize=16)
    ax.set_title(f'{disease}\n{tissue}\n{n_top_genes} top genes', fontsize=20)
    ax.annotate(text=f'p-value (MTClass) = {pval}', xy=(-0.15, -0.35), xycoords='axes fraction', fontsize=16)
    plt.tight_layout()
    plt.savefig(rf'C:/Users/ronni/OneDrive/Desktop/{disease}_{tissue}.png', dpi=400)
    plt.show()

    ### Plot three-way Venn Diagram
    mtclass_set = set(mtclass_genes)
    mphen_set = set(mphen_top_genes)
    manova_set = set(manova_top_genes)
    
    plt.figure(figsize=(5,5))
    venn3([mtclass_set, mphen_set, manova_set], ('MTClass', 'MultiPhen\n(binary)', 'MANOVA'), set_colors=('blue','orange','green'))
    plt.title(f'{disease}\n{tissue}\n{n_top_genes} top genes', fontsize=20)
    plt.tight_layout()
    plt.savefig(rf'C:/Users/ronni/OneDrive/Desktop/{tissue}_venn3_multivariate.png', dpi=400)
    plt.show()