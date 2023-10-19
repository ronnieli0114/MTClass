# -*- coding: utf-8 -*-
"""
Created on Sat Apr  8 21:43:05 2023

@author: ronnieli
"""
import pandas as pd
import numpy as np
import pyarrow.feather as feather
import ensembl_rest
import os
import sys
import subprocess
import matplotlib.pyplot as plt
import umap
from sklearn.preprocessing import StandardScaler

def get_expression_levels_2D(gene_id, snp, data_dir=r'E:/lab_data'):
    
    ''' Gets both multi-exon and multi-tissue expression levels from 9-tissue 
    2D multi-tissue and multi-exon study. '''
    
    # Load expression level data, filter by 9 tissues
    exp = feather.read_feather(os.path.join(data_dir,'gtex_exoneqtl','9_tissue_exon_expression.ftr'))        
    attr = pd.read_csv(os.path.join(data_dir, 'gtex_exoneqtl', 'v8_sample_attributes.txt'), sep='\t', header=0)
    tissue_list = ['Adipose - Subcutaneous','Artery - Tibial','Lung','Muscle - Skeletal','Nerve - Tibial',
                   'Skin - Not Sun Exposed (Suprapubic)','Skin - Sun Exposed (Lower leg)','Thyroid','Whole Blood']
    attr = attr[attr['SMTSD'].isin(tissue_list)][['SMTSD','SAMPID']]
    full_ids = attr['SAMPID'].tolist() 
    
    # reformat expression levels
    X = exp[exp['gene_id']==gene_id]
    X = X.filter(['Name']+full_ids)
    X = X.set_index('Name')
    X = X.stack().reset_index()
    X.columns = ['exon','SAMPID','exp']
    X = X.merge(attr, on='SAMPID')
    X['donor'] = X['SAMPID'].apply(lambda x: '-'.join(x.split('-')[:2]))
    X = X[['donor','SMTSD','exon','exp']].set_index('donor')
    X.columns = ['tissue','exon','exp']
    
    # Load genotypes
    c = snp.split('_')[0][3:]
    gt = pd.read_csv(os.path.join(data_dir,'gtex_exoneqtl','genotypes',f'chr{c}_genotypes_10kb_binary.txt.gz'), sep='\t', header=0)
    y = gt[(gt['ID']==snp) & (gt['gene_id']==gene_id)]
    y = y.iloc[:,3:].T
    y.columns = ['genotype']
    
    # merge together
    Xy = X.merge(y, left_index=True, right_index=True)
    
    # separate genotypes
    exp_0 = Xy[Xy['genotype']==0].iloc[:,:-1]
    exp_1 = Xy[Xy['genotype']==1].iloc[:,:-1]
    return exp_0, exp_1

def plot_9tissue_exon_expression(gene_id, snp, data_dir=r'E:/lab_data'):
    
    exp_0, exp_1 = get_expression_levels_2D(gene_id, snp)
    
    gene_data = ensembl_rest.lookup(gene_id)
    gene_name = gene_data['display_name']
    
    # plot for all 9 tissues
    fig, ax = plt.subplots(3, 3, sharex=True, sharey=False, figsize=(10,10))
    ax = ax.flatten()
    idx = np.arange(len(ax))

    # make separate subplot for each tissue
    
    for i, tissue in zip(idx, sorted(list(set(exp_0['tissue'])))):

        if tissue == 'Skin - Not Sun Exposed (Suprapubic)':
            tissue_name = 'Skin\nNot Sun Exposed\n(Suprapubic)'            
        else:
            tissue_name = tissue

        exp_0_tiss = exp_0[exp_0['tissue']==tissue].reset_index()
        exp_1_tiss = exp_1[exp_1['tissue']==tissue].reset_index()
        
        exp_0_tiss['exon'] = exp_0_tiss['exon'].apply(lambda x: int(x.split('_')[1]))
        exp_1_tiss['exon'] = exp_1_tiss['exon'].apply(lambda x: int(x.split('_')[1]))
        
        exp_0_tiss = exp_0_tiss.pivot(index='index', columns='exon', values='exp')
        exp_1_tiss = exp_1_tiss.pivot(index='index', columns='exon', values='exp')
        
        exp_tiss = pd.concat((exp_0_tiss, exp_1_tiss))
        div0 = exp_0_tiss.shape[0]
        
        # scale data
        scaler = StandardScaler()
        exp_scaled = scaler.fit_transform(exp_tiss)
        exp_scaled = pd.DataFrame(exp_scaled, index=list(exp_tiss.index), columns=list(exp_tiss.columns))
        
        exp_0_tiss_scaled = exp_scaled.iloc[:div0, :]
        exp_1_tiss_scaled = exp_scaled.iloc[div0:, :]
        
        title_list = ['Genotype REF/REF', 'Genotypes REF/ALT & ALT/ALT']
        exp_list = [exp_0_tiss_scaled, exp_1_tiss_scaled]
        shifts = [0, 0.2]
        
        for title, exp, shift in zip(title_list, exp_list, shifts):
            
            x = np.arange(len(exp.columns)) + 1 + shift
            y = np.mean(exp, axis=0)
            yerr = np.std(exp, axis=0)
            
            ax[i].errorbar(x=x, y=y, yerr=yerr, linestyle='', marker='o', label=title)
            ax[i].set_xticks(x, exp.columns, rotation=0)
            ax[i].set_title(f'{tissue_name}', fontsize=20)
            
            # set x and y limits
            ax[i].set_ylim((-3,3))
            if i >= 6:
                ax[i].set_xlabel('Exon number', fontsize=20)
    
    ax[0].legend()
    ax[3].set_ylabel('TPM\n(standardized)', fontsize=20)
    
    # save figure
    plt.suptitle(f'{gene_name}\n{snp}\n2D exon-tissue', fontsize=24, y=1)
    plt.tight_layout()
    save_dir = r'C:/Users/ronni/OneDrive/Desktop/plots'
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)
    plt.savefig(os.path.join(save_dir, f'{gene_name}-{snp}_expression_2D.png'), dpi=500)
    plt.show()

def get_binary_expression_2D(gene_id, snp, data_dir=r'E:/lab_data'):

    ''' Get 2D expression levels based on binary genotypes 
    Inputs:
        - gene_id: gene ID (Ensembl)
        - snp: GTEx formatted SNP
        - data directory: default is E:/lab_data
    '''

    # get chromosome number from SNP
    c = snp.split('_')[0][3:]
    
    # load data (expression and genotypes)
    exp = pd.read_csv(os.path.join(data_dir, 'gtex_exoneqtl', '2d_expression', f'chr{c}_2d_expression.txt.gz'), sep='\t', header=0)
    gt = pd.read_csv(os.path.join(data_dir, 'gtex_exoneqtl', 'genotypes', f'chr{c}_genotypes_10kb_binary.txt.gz'), sep='\t', header=0)
    
    # get expression levels
    X = exp[exp['gene_id']==gene_id]
    X = X.dropna(how='all', axis=1)
    X = X.drop(['gene_id','gene_name'], axis=1).set_index('donor')
    
    # get genotypes
    y = gt[(gt['ID']==snp) & (gt['gene_id']==gene_id)]
    y = y.drop_duplicates(['gene_id','gene_name','ID'])
    y = y.iloc[:,3:].T
    y.columns = ['genotype']
    
    # merge
    Xy = X.merge(y, left_index=True, right_index=True)
    
    exp_0 = Xy[Xy['genotype']==0].iloc[:,:-1]
    exp_1 = Xy[Xy['genotype']==1].iloc[:,:-1]
    return exp_0, exp_1

def plot_umap_2D(gene_id, snp, lab_dir=r'G:/My Drive/Lab', data_dir=r'E:/lab_data', savefig=True):
    
    ''' Plot UMAP for MTClass 2D results '''
    
    # Get gene name
    gene_data = ensembl_rest.lookup(gene_id)
    gene_name = gene_data['display_name']
    
    # get MTClass results
    results_dir = os.path.join(lab_dir,'lab_projects','gtex_exoneqtl','results','2d_exon_tissue')
    os.chdir(results_dir)
    mtclass = pd.read_csv('mtclass_2d_NN.txt.gz', sep='\t', header=0)
    
    f1, mcc = 'f1_macro', 'mcc'
    
    exp_0, exp_1 = get_binary_expression_2D(gene_id, snp)
    exp = pd.concat((exp_0, exp_1))
    
    div1 = exp_0.shape[0]
    
    scaler = StandardScaler()
    exp_scaled = scaler.fit_transform(exp)
    
    # run UMAP
    embed = umap.UMAP(n_neighbors=15, n_components=2).fit_transform(exp_scaled)
    
    fig, ax = plt.subplots(figsize=(5,5))
    ax.scatter(embed[:div1,0], embed[:div1,1], c='blue', label='REF/REF', alpha=0.5)
    ax.scatter(embed[div1:,0], embed[div1:,1], c='orange', label='REF/ALT & ALT/ALT', alpha=0.5)
    
    # annotate with MTClass results  
    f1 = mtclass[(mtclass['gene_id']==gene_id) & (mtclass['variant']==snp)][f1].values[0]
    mcc = mtclass[(mtclass['gene_id']==gene_id) & (mtclass['variant']==snp)][mcc].values[0]
    
    f1 = np.round(f1, 2)
    mcc = np.round(mcc, 2)
    
    x_coord = -50
    y_coord = -50
    
    ax.annotate(f'F1={f1}', xy=(x_coord, y_coord), xycoords='axes points')
    ax.annotate(f'MCC={mcc}', xy=(x_coord, y_coord-14), xycoords='axes points')

    plt.xlabel('UMAP 1', fontsize=14)
    plt.ylabel('UMAP 2', fontsize=14)

    ax.legend()
    ax.set_title(f'{gene_name}\n{snp}\n2D exon-tissue', fontsize=16)
    plt.tight_layout()
    
    if savefig:
        save_dir = r'C:/Users/ronni/OneDrive/Desktop'
        if not os.path.exists(save_dir):
            os.mkdir(save_dir)
        plt.savefig(os.path.join(save_dir, f'{gene_id}_{snp}_umap_2D.png'), dpi=400)
    plt.show()