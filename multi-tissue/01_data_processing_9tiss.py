# -*- coding: utf-8 -*-
"""
Data processing of GTEx tables to determine
amount of missing data present and
how many tissues and donors to use

Created on Mon Feb 14 13:11:37 2022

@author: ronnieli
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pyarrow.feather as feather
import numpy as np
import time
import os
from cmapPy.pandasGEXpress.parse_gct import parse
pd.set_option('chained_assignment',None)

def parse_GTEx_data(data_dir=r'E:/lab_data', save_files=False):
    
    ''' 
    Parses the GCF file downloaded from GTEx.
    Returns a tuple in this order:
        - Column metadata (donors)
        - Row metadata (gene names and IDs)
        - TPM expression
    '''
    
    # import expression matrix GCF file
    t0 = time.time()
    tpm_data = parse(os.path.join(data_dir, 'gtex_eqtl', 'v8_gene_tpm.gct'))
    t1 = time.time()
    print(f'Took {(t1-t0)/60:.1f} minutes to parse TPM data')
    
    col_metadata = tpm_data.col_metadata_df
    row_metadata = tpm_data.row_metadata_df
    tpm_counts = tpm_data.data_df
    
    # save files to pickle
    if save_files:
        col_metadata.to_pickle(os.path.join(data_dir, 'col_metadata.pkl'))
        row_metadata.to_pickle(os.path.join(data_dir, 'row_metadata.pkl'))
        feather.write_feather(tpm_counts, os.path.join(data_dir,'gene_TPM.ftr'))
    
    return (col_metadata, row_metadata, tpm_counts)

def load_GTEx_data(data_dir=r'E:/lab_data/mtclass_eqtl/multi-tissue/'):
    
    '''
    Returns a tuple in this order:
        - Column metadata (donors)
        - Row metadata (gene names and IDs)
        - TPM expression
        - Donor attributes
    '''
    
    # load GTEx data
    col_metadata = pd.read_pickle(os.path.join(data_dir,'col_metadata.pkl'))
    row_metadata = pd.read_pickle(os.path.join(data_dir,'row_metadata.pkl'))
    tpm = feather.read_feather(os.path.join(data_dir,'gene_TPM.ftr'))
    attributes = pd.read_csv(os.path.join(data_dir,'v8_sample_attributes.txt'), sep='\t', usecols=['SAMPID','SMTSD'])

    return (col_metadata, row_metadata, tpm, attributes)

def sample_availability(data_dir=r'E:/lab_data/mtclass_eqtl/multi-tissue/'):
    
    '''
    Determines the sample availability and determining tissues vs. donors.
    Returns a dataframe of donors x tissues, and the sample ID as the values where they exist.
    '''
    
    col_metadata, row_metadata, _, attributes = load_GTEx_data(data_dir)
    
    ### Make dataframe of sample/tissue availability
    avail_samples = attributes[attributes['SAMPID'].isin(col_metadata.index)]
    avail_samples['donor'] = avail_samples['SAMPID'].apply(lambda x: '-'.join(x.split("-")[:2]))
    avail_samples['count'] = 1
    avail_samples.columns = ['full_id', 'tissue', 'donor', 'count']
    
    samples_pivot = avail_samples.pivot(index='tissue', columns='donor', values='full_id')
    samples_pivot['num_donors'] = samples_pivot.count(axis=1)
    samples_pivot.sort_values(by='num_donors', ascending=False, na_position='last', inplace=True)

    # determine number of donors we'd have for given number of top tissues to include
    for n_tissues in range(1, len(samples_pivot.index)):
        subset = samples_pivot.iloc[:n_tissues, :-1].dropna(axis=1, how='any')
        num_samples = subset.shape[1]
        if num_samples > 0:
            print(f"{n_tissues} tissues:\t {num_samples} donors")

    return samples_pivot
    
def plot_tissue_sample_availability(data_dir=r"E:/lab_data/mtclass_eqtl/multi-tissue"):
    
    samples_pivot = sample_availability(data_dir)
    
    # plot tissues vs. samples
    plot_data = samples_pivot.drop('num_donors', axis=1).isna().astype(int)
    tissues = list(plot_data.index)
    donors = list(plot_data.columns)
    donors = [(i, d) for i, d in enumerate(donors) if i % 30 == 0]
    x = np.array([pos for (pos, donor) in donors])
    d = [d for (pos, d) in donors]
    
    plot_data = plot_data.sort_values(by=tissues[:9], axis=1, ascending=False)
    
    fig, ax = plt.subplots(figsize=(14,14))
    
    sns.heatmap(plot_data, cmap='binary', cbar=False, ax=ax, xticklabels=30)
    
    ax.set_xticks(x+0.5, d, rotation=90, fontsize=14)
    ax.set_yticks(np.arange(len(tissues))+0.5, tissues, fontsize=14)
    
    ax.set_title('Missing data for all GTEx', fontsize=20)
    ax.set_ylabel('Tissue', fontsize=20)
    ax.set_xlabel('Donor', fontsize=20)
    plt.tight_layout()
    plt.savefig(os.path.join(data_dir, 'missing_data_all_GTEx.png'), dpi=400)
    plt.show()

def make_9tissue_features(data_dir=r"E:/lab_data/mtclass_eqtl/multi-tissue/", num_tissues=9):
    
    samples_pivot = sample_availability(data_dir)
    col_metadata, row_metadata, tpm, attributes = load_GTEx_data(data_dir)
    print(f'Selecting top {num_tissues} tissues...')

    avail_samples = attributes[attributes['SAMPID'].isin(col_metadata.index)]
    avail_samples['donor'] = avail_samples['SAMPID'].apply(lambda x: '-'.join(x.split("-")[:2]))

    sel_samples = samples_pivot.iloc[:num_tissues,:-1].dropna(axis=1).values.flatten()

    # get common names for all genes
    tpm = row_metadata.merge(tpm, left_index=True, right_on='rid')

    print("Calculating variances for all genes...")
    
    tpm = tpm.set_index('Description').drop('rid', axis=1)
    
    # calculate top variances and subset using only high variance genes
    var_list = []
    
    # drop genes with TPM of zero in more than half the rows
    for _, row in tpm.iterrows():
        half_row = int(row.shape[0]/2) 
        row_var = np.var(row.values) 
        if np.count_nonzero(row.values==0) <= half_row: 
            var_list.append(row_var)
        else:
            var_list.append('QC_fail')

    tpm.insert(0, "variance", var_list)
    
    # keep genes that pass QC
    tpm_pass = tpm[tpm['variance'] != 'QC_fail']
    tpm_pass.sort_values(by='variance', ascending=False, inplace=True)
    num_high_var_genes = tpm_pass.shape[0]
    print("Number of genes with TPM=0 in less than half the samples:", num_high_var_genes)
    
    # select top half of high variance genes
    top_tpm = tpm_pass.iloc[:int(tpm_pass.shape[0]/2),:]

    # remove duplicate rows based on gene name
    top_tpm = top_tpm.loc[~top_tpm.index.duplicated(), :]
    print('Number of selected genes:', top_tpm.shape[0])

    # create features
    print("Creating feature matrix...")
    feats = top_tpm.filter(sel_samples)
    feats = feats.T.merge(avail_samples, left_index=True, right_on='SAMPID')
    feats = feats.rename(columns={'SMTSD':'tissue','SAMPID':'full_id'})
    feats = feats.set_index(['tissue','donor']).drop('full_id', axis=1)
    feats = feats.T.stack().reset_index().rename(columns={'level_0':'gene'})
    feats = feats.set_index('gene')
    
    feature_matrix = feats.copy()
    print('=== done ===')

    return top_tpm, feature_matrix

def save_features(data_dir=r"E:/lab_data/mtclass_eqtl/multi-tissue"):
    
    top_tpm, features = make_9tissue_features(data_dir)
    print("Saving top TPM and feature matrices...")
    
    feather.write_feather(top_tpm, os.path.join(data_dir, 'top_TPM.ftr'))
    feather.write_feather(features, os.path.join(data_dir, '9_tissue_TPM.ftr'))

def create_TPM_to_impute(data_dir=r"E:/lab_data/mtclass_eqtl/multi-tissue"):
    
    import random
    random.seed(2022)
    
    num_tissues = 9
    
    samples_pivot = sample_availability(data_dir)
    sel_samples = samples_pivot.iloc[:num_tissues,:-1].dropna(axis=1).values.flatten()
    sel_samples = list(sel_samples)
    
    def del_random_elements(input_list, n):
        to_delete = set(random.sample(range(len(input_list)), n))
        return [x for i,x in enumerate(input_list) if i not in to_delete]
    
    random_samples = del_random_elements(sel_samples, int(0.5*len(sel_samples)))
    with open(os.path.join(data_dir, 'samples_impute.txt'), 'w') as f:
        for sample in random_samples:
            f.write(sample)
            f.write('\n')
    
    print("Making feature matrix for imputation testing...")
    top_tpm = feather.read_feather(os.path.join(data_dir, "top_TPM.ftr"))
    
    samples_pivot = sample_availability(data_dir)
    col_metadata, row_metadata, _, attributes = load_GTEx_data(data_dir)
    
    avail_samples = attributes[attributes['SAMPID'].isin(col_metadata.index)]
    avail_samples['donor'] = avail_samples['SAMPID'].apply(lambda x: '-'.join(x.split("-")[:2]))
    
    feats = top_tpm.filter(sel_samples)
    feats = feats.T.merge(avail_samples, left_index=True, right_on='SAMPID', how='left')
    feats = feats.rename(columns={'SMTSD':'tissue','SAMPID':'full_id'})
    feats = feats.set_index(['tissue','donor'])
    feats = feats.mask(~feats['full_id'].isin(random_samples)).drop('full_id', axis=1)
    feats = feats.T.stack(dropna=False).reset_index().rename(columns={'level_0':'gene'})
    feats = feats.set_index('gene')
    
    return feats