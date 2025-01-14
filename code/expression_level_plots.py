# -*- coding: utf-8 -*-
"""
Expression level plots
"""
import pandas as pd
import numpy as np
import pyarrow.feather as feather
import os
import sys
import umap
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import TSNE
sys.path.append(r'C:/Users/ronni/OneDrive/Documents/Lab/lab_projects/mtclass_eqtl/multi-tissue')
from variant_functions import get_binary_exp, get_additive_exp

def plot_umap(gene, variant, experiment, savefig=True):
    
    os.chdir(r'C:/Users/ronni/OneDrive/Documents/lab_data/mtclass_eqtl/results')
    if experiment == '9_tissue':
        os.chdir('9_tissue')
        mtclass = pd.read_csv('9_tissue_mtclass.txt.gz', sep='\t', header=0)
        mphen_binary = pd.read_csv('9_tissue_multiphen_binary.txt.gz', sep='\t', header=0)
        manova = pd.read_csv('9_tissue_manova.txt.gz', sep='\t', header=0)
        experiment_title = '9-tissue'
    elif experiment == 'brain':
        os.chdir('brain_tissue')
        mtclass = pd.read_csv('brain_mtclass.txt.gz', sep='\t', header=0)
        mphen_binary = pd.read_csv('brain_multiphen_binary.txt.gz', sep='\t', header=0)
        manova = pd.read_csv('brain_manova.txt.gz', sep='\t', header=0)
        experiment_title = 'Brain tissue'
    elif experiment == '48_tissue':
        os.chdir('48_tissue')
        mtclass = pd.read_csv('48_tissue_mtclass.txt.gz', sep='\t', header=0)
        mphen_binary = pd.read_csv('48_tissue_multiphen_binary.txt.gz', sep='\t', header=0)
        manova = pd.read_csv('48_tissue_manova.txt.gz', sep='\t', header=0)
        experiment_title = '48-tissue'
    else:
        raise ValueError("Experiment name not found")
    print('loaded MTClass, MultiPhen, MANOVA')
    
    f1_key, mcc_key = 'f1_macro_median', 'mcc_median'
    
    print('fetching and standardizing expression levels...')
    counts_0, counts_1 = get_binary_exp(gene, variant, experiment)
    div1 = counts_0.shape[0]
    
    counts = pd.concat((counts_0, counts_1), axis=0)
    scaler = StandardScaler()
    counts_scaled = scaler.fit_transform(counts)
    
    # run UMAP
    print('running UMAP...')
    embed = umap.UMAP(n_neighbors=20, n_components=2).fit_transform(counts_scaled)
   
    fig, ax = plt.subplots(figsize=(6,6))
    ax.scatter(embed[:div1,0], embed[:div1,1], c='blue', label='WT', alpha=0.5)
    ax.scatter(embed[div1:,0], embed[div1:,1], c='orange', label='Mutant', alpha=0.5)
    
    # annotate with MTClass results  
    f1 = mtclass[(mtclass['gene']==gene) & (mtclass['variant']==variant)][f1_key].values[0]
    mcc = mtclass[(mtclass['gene']==gene) & (mtclass['variant']==variant)][mcc_key].values[0]
    
    f1 = np.around(f1, 3)
    mcc = np.around(mcc, 3)
    
    x_coord = -50
    y_coord = -50
    
    ax.annotate(f'F1={f1}', xy=(x_coord, y_coord), xycoords='axes points', fontsize=14)
    ax.annotate(f'MCC={mcc}', xy=(x_coord, y_coord-14), xycoords='axes points', fontsize=14)
        
    # annotate with MultiPhen (binary) results
    pval_mphen_binary = mphen_binary[(mphen_binary['gene']==gene) & (mphen_binary['variant']==variant)]['pval'].values[0]
    if pval_mphen_binary == 'NA':
        ax.annotate('MultiPhen p-val=NA', xy=(x_coord, y_coord-28), xycoords='axes points', fontsize=14)
    else:
        ax.annotate(f'MultiPhen p-val={pval_mphen_binary:.3e}', xy=(x_coord, y_coord-28), xycoords='axes points', fontsize=14)

    # annotate with MANOVA results
    pval_manova = manova[(manova['gene']==gene) & (manova['variant']==variant)]['pval'].values[0]
    if pval_manova == 'NA':
        ax.annotate('MANOVA p-val=NA', xy=(x_coord, y_coord-42), xycoords='axes points', fontsize=14)
    else:
        ax.annotate(f'MANOVA p-val={pval_manova:.3e}', xy=(x_coord, y_coord-42), xycoords='axes points', fontsize=14)
    
    plt.xlabel('UMAP 1', fontsize=14)
    plt.ylabel('UMAP 2', fontsize=14)

    ax.legend()
    ax.set_title(f'{gene}\n{variant}\n{experiment_title}', fontsize=16)
    plt.tight_layout()
    
    plt.savefig(rf'C:/Users/ronni/Desktop/{gene}_{variant}_umap.png', dpi=400)        
    plt.show()


def plot_tsne(gene, variant, experiment, savefig=True):
     
    '''
    Plots t-SNE dimension reduction of expression levels separated by genotype.
    
    Inputs:
        - gene name
        - GTEx formatted variant
        - MTClass, MultiPhen, MANOVA results (for annotation)
        - data directory
        - experiment name ('9_tissue', 'brain', 'all')
        - option to save figure (default True)
    '''
    
    os.chdir(r'C:/Users/ronni/OneDrive/Documents/lab_data/mtclass_eqtl/results')
    if experiment == '9_tissue':
        os.chdir('9_tissue')
        mtclass = pd.read_csv('9_tissue_mtclass_aggregate.txt.gz', sep='\t', header=0)
        mphen_binary = pd.read_csv('9_tissue_multiphen_binary.txt.gz', sep='\t', header=0)
        manova = pd.read_csv('9_tissue_manova.txt.gz', sep='\t', header=0)
        experiment_title = '9-tissue'
    elif experiment == 'brain':
        os.chdir('brain_tissue')
        mtclass = pd.read_csv('brain_mtclass_aggregate.txt.gz', sep='\t', header=0)
        mphen_binary = pd.read_csv('brain_multiphen_binary.txt.gz', sep='\t', header=0)
        manova = pd.read_csv('brain_manova.txt.gz', sep='\t', header=0)
        experiment_title = 'Brain tissue'
    elif experiment == '48_tissue':
        os.chdir('48_tissue')
        mtclass = pd.read_csv('48_tissue_mtclass_aggregate.txt.gz', sep='\t', header=0)
        mphen_binary = pd.read_csv('48_tissue_multiphen_binary.txt.gz', sep='\t', header=0)
        manova = pd.read_csv('48_tissue_manova.txt.gz', sep='\t', header=0)
        experiment_title = '48-tissue'
    else:
        raise ValueError("Experiment name not found")
    print('loaded MTClass, MultiPhen, MANOVA')
    
    f1_key, mcc_key = 'f1_macro_median', 'mcc_median'
    
    print('fetching and standardizing expression levels...')
    counts_0, counts_1 = get_binary_exp(gene, variant, experiment)
    div1 = counts_0.shape[0]
    exp = pd.concat((counts_0, counts_1), axis=0)
    
    scaler = StandardScaler()
    exp_scaled = scaler.fit_transform(exp)
    
    # run TSNE
    embed = TSNE(n_components=2, perplexity=20, random_state=2022, learning_rate='auto').fit_transform(exp_scaled)
    
    fig, ax = plt.subplots(figsize=(6,6))
    ax.scatter(embed[:div1,0], embed[:div1,1], c='blue', label='Genotype 0', alpha=0.5)
    ax.scatter(embed[div1:,0], embed[div1:,1], c='orange', label='Genotype 1', alpha=0.5)

    # annotate with MTClass results 
    f1 = mtclass[(mtclass['gene']==gene) & (mtclass['variant']==variant)][f1_key].values[0]
    mcc = mtclass[(mtclass['gene']==gene) & (mtclass['variant']==variant)][mcc_key].values[0]
    
    f1 = np.round(f1, 3)
    mcc = np.round(mcc, 3)
    
    x_coord = -50
    y_coord = -50
    
    ax.annotate(f'F1={f1}', xy=(x_coord, y_coord), xycoords='axes points')
    ax.annotate(f'MCC={mcc}', xy=(x_coord, y_coord-14), xycoords='axes points')
        
    # annotate with MultiPhen (binary) results
    pval_mphen_binary = mphen_binary[(mphen_binary['gene']==gene) & (mphen_binary['variant']==variant)]['pval'].values[0]
    if pval_mphen_binary == 'NA':
        ax.annotate('MultiPhen p-val=NA', xy=(x_coord, y_coord-28), xycoords='axes points')
    else:
        ax.annotate(f'MultiPhen p-val={pval_mphen_binary:.3e}', xy=(x_coord, y_coord-28), xycoords='axes points')
        
    # annotate with MANOVA results
    pval_manova = manova[(manova['gene']==gene) & (manova['variant']==variant)]['pval'].values[0]
    if pval_manova == 'NA':
        ax.annotate('MANOVA p-val=NA', xy=(x_coord, y_coord-42), xycoords='axes points')
    else:
        ax.annotate(f'MANOVA p-val={pval_manova:.3e}', xy=(x_coord, y_coord-42), xycoords='axes points')

    plt.xlabel('t-SNE 1')
    plt.ylabel('t-SNE 2')

    ax.legend()
    ax.set_title(f'{gene}\n{variant}\n{experiment_title}', fontsize=16)
    plt.tight_layout()
    
    if savefig:
        save_dir = r'C:/Users/ronni/Desktop/plots'
        if not os.path.exists(save_dir):
            os.mkdir(save_dir)
        plt.savefig(os.path.join(save_dir, f'{gene}_{variant}_tsne.png'), dpi=400)
    plt.show()


def plot_heatmap(gene, variant, experiment, savefig=True):
    
    '''
    Plots heatmap of expression levels separated by genotype.
    
    Inputs:
        - gene name
        - GTEx formatted variant
        - experiment name ('9_tissue', 'brain', '48_tissue')
        - option to save figure (default True)
    '''
    
    os.chdir(r'C:/Users/ronni/OneDrive/Documents/lab_data/mtclass_eqtl/results')
    if experiment == '9_tissue':
        os.chdir('9_tissue')
        mtclass = pd.read_csv('9_tissue_mtclass_aggregate.txt.gz', sep='\t', header=0)
        mphen_binary = pd.read_csv('9_tissue_multiphen_binary.txt.gz', sep='\t', header=0)
        manova = pd.read_csv('9_tissue_manova.txt.gz', sep='\t', header=0)
        experiment_title = '9-tissue'
    elif experiment == 'brain':
        os.chdir('brain_tissue')
        mtclass = pd.read_csv('brain_mtclass_aggregate.txt.gz', sep='\t', header=0)
        mphen_binary = pd.read_csv('brain_multiphen_binary.txt.gz', sep='\t', header=0)
        manova = pd.read_csv('brain_manova.txt.gz', sep='\t', header=0)
        experiment_title = 'Brain tissue'
    elif experiment == '48_tissue':
        os.chdir('48_tissue')
        mtclass = pd.read_csv('48_tissue_mtclass_aggregate.txt.gz', sep='\t', header=0)
        mphen_binary = pd.read_csv('48_tissue_multiphen_binary.txt.gz', sep='\t', header=0)
        manova = pd.read_csv('48_tissue_manova.txt.gz', sep='\t', header=0)
        experiment_title = '48-tissue'
    else:
        raise ValueError("Experiment name not found")
    print('loaded MTClass, MultiPhen, MANOVA')
    
    f1_key, mcc_key = 'f1_macro_median', 'mcc_median'
    counts_0, counts_1 = get_binary_exp(gene, variant, experiment)
    div_0 = counts_0.shape[0]
    
    counts = pd.concat((counts_0, counts_1))
    
    scaler = StandardScaler()
    counts_scaled = scaler.fit_transform(counts)
    counts_scaled = pd.DataFrame(counts_scaled, index=list(counts.index), columns=list(counts.columns))
    
    counts_0 = counts_scaled.iloc[:div_0,:]
    counts_1 = counts_scaled.iloc[div_0:,:]
    
    # sort by median of REF/REF expression levels
    means = np.median(counts_0, axis=0)
    inds = np.argsort(means)
    
    counts_0 = counts_0.iloc[:,inds]
    counts_1 = counts_1.iloc[:,inds]
    counts = pd.concat((counts_0, counts_1))

    print('plotting...')
    fig, ax = plt.subplots(1,1, figsize=(10,10))
    ax = sns.heatmap(counts, ax=ax, yticklabels=10, vmin=-1, vmax=1)
    ax.axhline(y=div_0, color='yellow', linewidth=2)
    
    # annotate with MTClass results  
    f1 = mtclass[(mtclass['gene']==gene) & (mtclass['variant']==variant)][f1_key].values[0]
    mcc = mtclass[(mtclass['gene']==gene) & (mtclass['variant']==variant)][mcc_key].values[0]
    
    f1 = np.round(f1, 2)
    mcc = np.round(mcc, 2)
    
    x_coord = -50
    y_coord = -200
    
    ax.annotate(f'F1={f1}', xy=(x_coord, y_coord), xycoords='axes points')
    ax.annotate(f'MCC={mcc}', xy=(x_coord, y_coord-14), xycoords='axes points')
        
    # annotate with MultiPhen (binary) results
    pval_mphen_binary = mphen_binary[(mphen_binary['gene']==gene) & (mphen_binary['variant']==variant)]['pval'].values[0]
    if pval_mphen_binary == 'NA':
        ax.annotate('MultiPhen p-val=NA', xy=(x_coord, y_coord-28), xycoords='axes points')
    else:
        ax.annotate(f'MultiPhen p-val={pval_mphen_binary:.3e}', xy=(x_coord, y_coord-28), xycoords='axes points')

    # annotate with MANOVA results
    pval_manova = manova[(manova['gene']==gene) & (manova['variant']==variant)]['pval'].values[0]
    if pval_manova == 'NA':
        ax.annotate('MANOVA p-val=NA', xy=(x_coord, y_coord-42), xycoords='axes points')
    else:
        ax.annotate(f'MANOVA p-val={pval_manova:.3e}', xy=(x_coord, y_coord-42), xycoords='axes points')
    
    plt.yticks(rotation=0)
    plt.suptitle(f'{gene}\n{variant}\n{experiment_title}', fontsize=20)
    plt.tight_layout()
    
    if savefig:
        save_dir = r'C:/Users/ronni/Desktop/plots'
        if not os.path.exists(save_dir):
            os.mkdir(save_dir)
        plt.savefig(os.path.join(save_dir, f'{gene}_{variant}_heatmap.png'), dpi=400)
    plt.show()


def plot_boxplot(gene, variant, experiment, savefig=True):
    
    '''
    Plots boxplots of expression levels by tissue, separated by genotype.
    
    Inputs:
        - gene name
        - GTEx formatted variant
        - experiment name ('9_tissue', 'brain', '48_tissue')
        - option to save figure (default True)
    '''
 
    os.chdir(r'C:/Users/ronni/OneDrive/Documents/lab_data/mtclass_eqtl/results')
    if experiment == '9_tissue':
        os.chdir('9_tissue')
        mtclass = pd.read_csv('9_tissue_mtclass_aggregate.txt.gz', sep='\t', header=0)
        mphen_binary = pd.read_csv('9_tissue_multiphen_binary.txt.gz', sep='\t', header=0)
        manova = pd.read_csv('9_tissue_manova.txt.gz', sep='\t', header=0)
        experiment_title = '9-tissue'
    elif experiment == 'brain':
        os.chdir('brain_tissue')
        mtclass = pd.read_csv('brain_mtclass_aggregate.txt.gz', sep='\t', header=0)
        mphen_binary = pd.read_csv('brain_multiphen_binary.txt.gz', sep='\t', header=0)
        manova = pd.read_csv('brain_manova.txt.gz', sep='\t', header=0)
        experiment_title = 'Brain tissue'
    elif experiment == '48_tissue':
        os.chdir('48_tissue')
        mtclass = pd.read_csv('48_tissue_mtclass_aggregate.txt.gz', sep='\t', header=0)
        mphen_binary = pd.read_csv('48_tissue_multiphen_binary.txt.gz', sep='\t', header=0)
        manova = pd.read_csv('48_tissue_manova.txt.gz', sep='\t', header=0)
        experiment_title = '48-tissue'
    else:
        raise ValueError("Experiment name not found")
    print('loaded MTClass, MultiPhen, MANOVA')
    
    f1_key, mcc_key = 'f1_macro_median', 'mcc_median'
    
    print('fetching and standardizing expression levels...')
    counts_0, counts_1 = get_binary_exp(gene, variant, experiment)
    div_0 = counts_0.shape[0]
    
    counts = pd.concat((counts_0, counts_1))
    
    scaler = StandardScaler()
    counts_scaled = scaler.fit_transform(counts)
    counts_scaled = pd.DataFrame(counts_scaled, index=list(counts.index), columns=list(counts.columns))
    
    counts_0 = counts_scaled.iloc[:div_0,:]
    counts_1 = counts_scaled.iloc[div_0:,:]

    # sort by median of REF/REF expression levels
    means = np.median(counts_0, axis=0)
    inds = np.argsort(means)
    
    counts_0 = counts_0.iloc[:,inds]
    counts_1 = counts_1.iloc[:,inds]

    print('plotting...')
    # Plot boxplots
    fig, ax = plt.subplots(figsize=(12, 12))
    pos_0 = np.arange(len(counts_scaled.columns))+1
    pos_1 = pos_0 + 0.2
    b0 = ax.boxplot(counts_0, positions=pos_0, widths=0.2, patch_artist=True, sym='')
    b1 = ax.boxplot(counts_1, positions=pos_1, widths=0.2, patch_artist=True, sym='')
    
    # set boxplot colors
    for patch in b0['boxes']:
        patch.set_facecolor('steelblue')
        patch.set_alpha(0.7)
    for patch in b1['boxes']:
        patch.set_facecolor('orange')
        patch.set_alpha(0.7)
    
    # set x and y limits
    ax.set_ylabel('TPM (standardized)', fontsize=20)
    x = np.arange(len(counts_0.columns)) + 1.25
    ax.set_xticks(x, counts_0.columns, rotation=90, fontsize=14)

    # plot legend
    ax.legend(handles=(b0['boxes'][0], b1['boxes'][0]), labels=('Genotype 0','Genotype 1'), fontsize=16)

    # annotate with MTClass results
    f1 = mtclass[(mtclass['gene']==gene) & (mtclass['variant']==variant)][f1_key].values[0]
    mcc = mtclass[(mtclass['gene']==gene) & (mtclass['variant']==variant)][mcc_key].values[0]
    
    f1 = np.round(f1, 3)
    mcc = np.round(mcc, 3)
    
    x_coord = -50
    y_coord = -250
    
    ax.annotate(f'F1={f1}', xy=(x_coord, y_coord), xycoords='axes points')
    ax.annotate(f'MCC={mcc}', xy=(x_coord, y_coord-14), xycoords='axes points')
        
    # annotate with MultiPhen (binary) results
    pval_mphen_binary = mphen_binary[(mphen_binary['gene']==gene) & (mphen_binary['variant']==variant)]['pval'].values[0]
    if pval_mphen_binary == 'NA':
        ax.annotate('MultiPhen p-val=NA', xy=(x_coord, y_coord-28), xycoords='axes points')
    else:
        ax.annotate(f'MultiPhen p-val={pval_mphen_binary:.3e}', xy=(x_coord, y_coord-28), xycoords='axes points')

    # annotate with MANOVA results
    pval_manova = manova[(manova['gene']==gene) & (manova['variant']==variant)]['pval'].values[0]
    if pval_manova == 'NA':
        ax.annotate('MANOVA p-val=NA', xy=(x_coord, y_coord-42), xycoords='axes points')
    else:
        ax.annotate(f'MANOVA p-val={pval_manova:.3e}', xy=(x_coord, y_coord-42), xycoords='axes points')

    plt.suptitle(f'{gene}\n{variant}\n{experiment_title}', fontsize=20)
    
    plt.tight_layout()
    if savefig:
        save_dir = r'C:/Users/ronni/Desktop/plots'
        if not os.path.exists(save_dir):
            os.mkdir(save_dir)
        plt.savefig(os.path.join(save_dir, f'{gene}-{variant}_boxplot.png'), dpi=500)
    plt.show()


def plot_mean_SD(gene, variant, experiment, savefig=True):
    
    '''
    Plots means and variances of expression levels by tissue, grouped by genotype.
    
    Inputs:
        - gene name
        - GTEx formatted variant
        - experiment name ('9_tissue', 'brain', '48_tissue') - to import expression levels
        - option to save figure (default True)
    '''
    
    os.chdir(r'C:/Users/ronni/OneDrive/Documents/lab_data/mtclass_eqtl/results')
    if experiment == '9_tissue':
        os.chdir('9_tissue')
        mtclass = pd.read_csv('9_tissue_mtclass.txt.gz', sep='\t', header=0)
        mphen_binary = pd.read_csv('9_tissue_multiphen_binary.txt.gz', sep='\t', header=0)
        manova = pd.read_csv('9_tissue_manova.txt.gz', sep='\t', header=0)
        experiment_title = '9-tissue'
    elif experiment == 'brain':
        os.chdir('brain_tissue')
        mtclass = pd.read_csv('brain_mtclass.txt.gz', sep='\t', header=0)
        mphen_binary = pd.read_csv('brain_multiphen_binary.txt.gz', sep='\t', header=0)
        manova = pd.read_csv('brain_manova.txt.gz', sep='\t', header=0)
        experiment_title = 'Brain tissue'
    elif experiment == '48_tissue':
        os.chdir('48_tissue')
        mtclass = pd.read_csv('48_tissue_mtclass.txt.gz', sep='\t', header=0)
        mphen_binary = pd.read_csv('48_tissue_multiphen_binary.txt.gz', sep='\t', header=0)
        manova = pd.read_csv('48_tissue_manova.txt.gz', sep='\t', header=0)
        experiment_title = '48-tissue'
    else:
        raise ValueError("Experiment name not found")
    print('loaded MTClass, MultiPhen, MANOVA')
    
    print('fetching and scaling expression levels...')
    counts_0, counts_1 = get_binary_exp(gene, variant, experiment)
    counts_all = pd.concat((counts_0, counts_1))
    
    div_0 = counts_0.shape[0]
    
    scaler = StandardScaler()
    counts_scaled = scaler.fit_transform(counts_all)
    counts_scaled = pd.DataFrame(counts_scaled, index=list(counts_all.index), columns=list(counts_all.columns))
    
    counts_0 = counts_scaled.iloc[:div_0,:]
    counts_1 = counts_scaled.iloc[div_0:,:]
    
    # sort by median of REF/REF expression levels
    means = np.mean(counts_0, axis=0)
    inds = np.argsort(means)
    
    counts_0 = counts_0.iloc[:,inds]
    counts_1 = counts_1.iloc[:,inds]
    
    # plot means and variances as error bars
    print('plotting...')
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    
    title_list = ['WT', 'Mutant']
    counts_list = [counts_0, counts_1]
    shifts = [0, 0.2]
    
    for title, counts, shift in zip(title_list, counts_list, shifts):
        
        counts.loc['mean'] = np.mean(counts, axis=0)
        counts.loc['std'] = np.std(counts, axis=0)
        
        x = np.arange(len(counts.columns)) + 1 + shift
        y = counts.loc['mean'].values
        yerr = counts.loc['std'].values
        
        ax.errorbar(x=x, y=y, yerr=yerr, linestyle='', marker='o', label=title)
        ax.set_xticks(x, counts.columns, rotation=90, fontsize=14)
        
        # set x and y limits
        ax.set_ylim((-3,3))
    
    # # annotate with MTClass results
    # f1_key, mcc_key = 'f1_macro_median', 'mcc_median'

    # f1 = mtclass[(mtclass['gene']==gene) & (mtclass['variant']==variant)][f1_key].values[0]
    # mcc = mtclass[(mtclass['gene']==gene) & (mtclass['variant']==variant)][mcc_key].values[0]
    
    # f1 = np.around(f1, 3)
    # mcc = np.around(mcc, 3)
    
    # x_coord = -50
    # y_coord = -200
    fs = 14
    
    # ax.annotate(f'F1={f1}', xy=(x_coord, y_coord), xycoords='axes points', fontsize=fs)
    # ax.annotate(f'MCC={mcc}', xy=(x_coord, y_coord-14), xycoords='axes points', fontsize=fs)    
    
    # # annotate with MultiPhen (binary) results
    # pval_mphen_binary = mphen_binary[(mphen_binary['gene']==gene) & (mphen_binary['variant']==variant)]['pval'].values[0]
    # if pval_mphen_binary == 'NA':
    #     ax.annotate('MultiPhen p-val=NA', xy=(x_coord, y_coord-28), xycoords='axes points', fontsize=fs)
    # else:
    #     ax.annotate(f'MultiPhen p-val={pval_mphen_binary:.3e}', xy=(x_coord, y_coord-28), xycoords='axes points', fontsize=fs)

    # # annotate with MANOVA results
    # pval_manova = manova[(manova['gene']==gene) & (manova['variant']==variant)]['pval'].values[0]
    # if pval_manova == 'NA':
    #     ax.annotate('MANOVA p-val=NA', xy=(x_coord, y_coord-42), xycoords='axes points', fontsize=fs)
    # else:
    #     ax.annotate(f'MANOVA p-val={pval_manova:.3e}', xy=(x_coord, y_coord-42), xycoords='axes points', fontsize=fs)

    ax.set_xlabel('Tissue', fontsize=20)
    ax.set_ylabel('TPM expression \n (standardized)', fontsize=20)
    plt.suptitle(f'Mean and SD\n{gene}\n{variant}\n{experiment_title}', fontsize=20)
    plt.tight_layout()
    plt.legend(fontsize=fs)
    plt.savefig(rf'C:/Users/ronni/Desktop/{gene}-{variant}_meanSD.png', dpi=400)
    plt.show()
    
#%%

'''
Determine whether HLA genes are enriched among the results
differentially detected by MTClass
'''

df1 = mtclass.merge(mphen_binary, on='pair')
df2 = df1.merge(manova, on='pair', suffixes=('_mphen','_manova'))
df3 = df2.loc[:,['pair','sample_size','mcc_median','f1_macro_median','pval_mphen','pval_manova']]


signif = df3[(df3['pval_mphen'] > 0.05) & (df3['pval_manova'] > 0.05) & (df3['f1_macro_median'] > 0.5) & (df3['mcc_median'] > 0.2)]
not_signif = df3[~df3['pair'].isin(signif['pair'])]

num_hla_signif = len([x for x in signif['pair'] if 'HLA-' in x])
num_hla_nonsignif = len([x for x in not_signif['pair'] if 'HLA-' in x])

spot0 = num_hla_signif
spot1 = len(signif) - num_hla_signif
spot2 = num_hla_nonsignif
spot3 = len(not_signif) - num_hla_nonsignif

df = pd.DataFrame([[spot0,spot1],[spot2,spot3]], columns=['HLA','Not HLA'], index=['Diff_MTClass','Not_Diff_MTClass'])

from scipy.stats import fisher_exact
stat, p = fisher_exact(df)
