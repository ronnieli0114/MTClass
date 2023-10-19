# -*- coding: utf-8 -*-
"""
Fully connected neural network
using 2D data from different exons and tissues

Created on Mon Feb 20 01:25:28 2023
@author: ronnieli
"""

import pyarrow.feather as feather
import pandas as pd
import numpy as np
import os
import gc
import time
import multiprocessing as mp
import matplotlib.pyplot as plt
import tensorflow as tf
from itertools import repeat
from keras.models import Sequential
from keras.layers import Dense, Conv2D, Flatten
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import f1_score, matthews_corrcoef
from scipy.cluster.hierarchy import dendrogram, linkage
tf.config.threading.set_inter_op_parallelism_threads(4)
tf.config.threading.set_intra_op_parallelism_threads(4)
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

c = 1 # chromosome number
data_dir = r'E:/lab_data'
num_cpus = 4
output_filename = f'mtclass_NN_chr{c}'

def hierarchical_clustering(data_dir):
    
    '''
    Performs hierarchical clustering on highest variance genes
    from the 9-tissue case to determine optimal ordering of tissues
    '''
    
    # Determine what genes to use for hierarhical clustering
    top_exp = feather.read_feather(os.path.join(data_dir,'gtex_eqtl','top_tpm_counts.ftr'))
    top_exp = top_exp.sort_values(by='variance', ascending=False)
    top_genes = top_exp[top_exp['variance'] >= np.percentile(top_exp['variance'], 75)]['Description'].tolist() # top 25% of genes by variance
    # 3572 genes

    # read in expression levels
    nine_tiss_exp = feather.read_feather(os.path.join(data_dir, 'gtex_eqtl', 'multi_tissue', 'counts.ftr'))
    nine_tiss_exp = nine_tiss_exp[nine_tiss_exp['gene'].isin(top_genes)]

    data = nine_tiss_exp.drop(['gene','donor'], axis=1).transpose()
    linkage_data = linkage(data, method='ward', metric='euclidean')
    plt.figure(figsize=(8,8))
    dend = dendrogram(linkage_data, labels=list(nine_tiss_exp.columns)[2:], leaf_rotation=90)
    plt.ylabel('Euclidean distance', fontsize=16)
    plt.xlabel('Tissue', fontsize=16)
    plt.title('Dendrogram of 9 tissues', fontsize=20)
    plt.tight_layout()
    plt.savefig(r'C:/Users/ronnieli/Desktop/dendrogram.png', dpi=400)
    plt.show()

def filter_expression_level_data(data_dir):
    
    '''
    Filters top_exon_counts.ftr by donors and tissues from the 9-tissue case 
    (103 donors, 9 somatic tissues with no missing data). Used for a 2D CNN
    analyzing both multi-exon and multi-tissue expression
    '''

    exp = feather.read_feather(os.path.join(data_dir,'gtex_exoneqtl','top_exon_counts.ftr'))
    attr = pd.read_csv(os.path.join(data_dir,'gtex_eqtl','v8_sample_attributes.txt'), sep='\t', header=0)
    donors = pd.read_csv(os.path.join(data_dir, 'gtex_eqtl','multi_tissue','9_tissue_donors.txt'), sep='\t', header=None)[0].tolist()
    gene_annot = pd.read_csv(os.path.join(data_dir, 'gtex_exoneqtl', 'top_exon_gene_info.txt.gz'), sep='\t', header=0)
    exp = exp[exp['gene_id'].isin(gene_annot['gene_id'])]
    print('===== Loaded data =====')

    attr['DONOR'] = attr['SAMPID'].apply(lambda x: '-'.join(x.split('-')[:2]))

    tissues = ['Adipose - Subcutaneous','Artery - Tibial','Lung','Muscle - Skeletal',
               'Nerve - Tibial','Skin - Not Sun Exposed (Suprapubic)',
               'Skin - Sun Exposed (Lower leg)','Thyroid','Whole Blood']

    attr_sub = attr[(attr['DONOR'].isin(donors)) & (attr['SMTSD'].isin(tissues))]
    full_ids = attr_sub['SAMPID'].tolist()
    exp_filter = exp.filter(['gene_id','Description','Name']+full_ids)
    return exp_filter

def load_expression_level_data(data_dir):
    return feather.read_feather(os.path.join(data_dir, 'gtex_exoneqtl', '9_tissue_exon_expression.ftr'))

def load_attributes(data_dir):
    return pd.read_csv(os.path.join(data_dir,'gtex_exoneqtl','v8_sample_attributes.txt'), sep='\t', header=0)

def load_gene_annotations(data_dir):
    return pd.read_csv(os.path.join(data_dir, 'gtex_exoneqtl', 'top_exon_gene_info.txt.gz'), sep='\t', header=0)

### Define function to reshape expression level data

def prepare_dict(gene_id, exp, attr):
    
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

    return donor_dict

### Define function to run CNN

def run_NN(gene_id, gene_annot, exp, attr, data_dir):
    
    # print(f'Now on gene {gene_id}')
    results_dict = {'gene_id':[], 'gene_name':[], 'variant':[], 'f1_micro':[], 'f1_macro':[], 'mcc':[]}
    
    # prepare expression arrays
    donor_dict = prepare_dict(gene_id, exp, attr)
    donor_list = list(donor_dict.keys())
    
    # look up gene
    gene_info = gene_annot[gene_annot['gene_id']==gene_id]
    gene_name = gene_info['gene_name'].values[0]
    chrom = gene_info['chr'].values[0]
        
    # import genotypes
    gt = pd.read_csv(os.path.join(data_dir, 'gtex_exoneqtl', 'genotypes', f'chr{chrom}_genotypes_10kb_binary.txt.gz'), sep='\t', header=0, usecols=['gene_id','ID']+donor_list)
    gene_gt = gt[gt['gene_id'] == gene_id]
    gene_gt = gene_gt.drop_duplicates(subset=['gene_id','ID'], keep='first')

    num_snps = len(gene_gt['ID'])
    print(f'There are {num_snps} SNPs for gene {gene_id}')    

    for snp in gene_gt['ID']:

        tf.keras.backend.clear_session() # clear session to prevent memory leakage
        
        # print(f'Now on pair {gene_id}-{snp}')
        
        snp_gt = gene_gt[gene_gt['ID']==snp]
        snp_gt = snp_gt.filter(donor_list)
        snp_gt = snp_gt.stack().reset_index()
        snp_gt = snp_gt.drop('level_0', axis=1)
        snp_gt.columns = ['donor','genotype']  
        snp_gt = snp_gt.set_index('donor')
        snp_gt = snp_gt.sort_index()
        
        # take care of missing genotypes
        nonmissing_gt = snp_gt[(snp_gt['genotype']==0) | (snp_gt['genotype']==1)]
        nonmissing_donors = list(nonmissing_gt.index)
        
        y = nonmissing_gt['genotype'].values.astype(np.float32)
        y = y.reshape((-1, 1))
        
        if np.count_nonzero(y==0) < 5 or np.count_nonzero(y==1) < 5:
            continue

        # take care of missing expression levels
        nonmissing_array_dict = {donor: array for (donor, array) in donor_dict.items() if donor in nonmissing_donors}
        nonmissing_arrays = list(nonmissing_array_dict.values())

        X = np.vstack(nonmissing_arrays)
        
        num_features = X.shape[1]
            
        f1_micro_list = []
        f1_macro_list = []
        mcc_list = []
            
        # Stratified K-Fold cross validation
        skf = StratifiedKFold(n_splits=5)
        for train_idx, test_idx in skf.split(X, y):

            gc.collect() # garbage collector
            
            X_train, X_test = X[train_idx], X[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]

            # build Keras model
            model = Sequential()
            model.add(Dense(int(num_features/2), activation='relu', input_shape=(num_features,)))
            model.add(Dense(1, activation='sigmoid'))
            
            # compile and fit model
            model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
            model.fit(X_train, y_train, epochs=10, verbose=False)
            
            # make predictions
            preds = model(X_test)
            preds = np.around(preds)
            
            f1_micro = np.around(f1_score(y_test, preds, average='micro'), 3)
            f1_macro = np.around(f1_score(y_test, preds, average='macro'), 3)
            mcc = np.around(matthews_corrcoef(y_test, preds), 3)
            
            f1_micro_list.append(f1_micro)
            f1_macro_list.append(f1_macro)
            mcc_list.append(mcc)
        
        f1_micro_avg = np.around(np.average(f1_micro_list), 3)
        f1_macro_avg = np.around(np.average(f1_macro_list), 3)
        mcc_avg = np.around(np.average(mcc_list), 3)
        print(f'{snp}\t F1 macro={f1_macro_avg} \t MCC={mcc_avg}', flush=True)
        
        results_dict['gene_id'].append(gene_id)
        results_dict['gene_name'].append(gene_name)
        results_dict['variant'].append(snp)
        results_dict['f1_micro'].append(f1_micro_avg)
        results_dict['f1_macro'].append(f1_macro_avg)
        results_dict['mcc'].append(mcc_avg)
        
    return pd.DataFrame(results_dict)

def gene_X(gene_id, gene_annot, exp, attr, data_dir):
    
    data_list = []
    
    # prepare expression arrays
    donor_dict = prepare_dict(gene_id, exp, attr)
    donor_list = list(donor_dict.keys())
    
    # look up gene
    gene_info = gene_annot[gene_annot['gene_id']==gene_id]
    gene_name = gene_info['gene_name'].values[0]
    chrom = gene_info['chr'].values[0]
        
    # import genotypes
    gt = pd.read_csv(os.path.join(data_dir, 'gtex_exoneqtl', 'genotypes', f'chr{chrom}_genotypes_10kb_binary.txt.gz'), sep='\t', header=0, usecols=['gene_id','ID']+donor_list)
    gene_gt = gt[gt['gene_id'] == gene_id]
    gene_gt = gene_gt.drop_duplicates(subset=['gene_id','ID'], keep='first')

    num_snps = len(gene_gt['ID'])
    # print(f'There are {num_snps} SNPs for gene {gene_id}')

    for snp in gene_gt['ID']:

        tf.keras.backend.clear_session() # clear session to prevent memory leakage
        
        # print(f'Now on pair {gene_id}-{snp}')
        
        snp_gt = gene_gt[gene_gt['ID']==snp]
        snp_gt = snp_gt.filter(donor_list)
        snp_gt = snp_gt.stack().reset_index()
        snp_gt = snp_gt.drop('level_0', axis=1)
        snp_gt.columns = ['donor','genotype']  
        snp_gt = snp_gt.set_index('donor')
        snp_gt = snp_gt.sort_index()
        
        # take care of missing genotypes
        nonmissing_gt = snp_gt[(snp_gt['genotype']==0) | (snp_gt['genotype']==1)]
        nonmissing_donors = list(nonmissing_gt.index)
        
        y = nonmissing_gt['genotype'].values.astype(np.float32)
        y = y.reshape((-1, 1))
        
        if np.count_nonzero(y==0) < 5 or np.count_nonzero(y==1) < 5:
            continue

        # take care of missing expression levels
        nonmissing_array_dict = {donor: array for (donor, array) in donor_dict.items() if donor in nonmissing_donors}
        nonmissing_arrays = list(nonmissing_array_dict.values())

        X = np.vstack(nonmissing_arrays)
        
        X_df = pd.DataFrame(X, index=nonmissing_donors)
        X_df['gene_id'] = gene_id
        X_df['gene_name'] = gene_name
        X_df['variant'] = snp
        X_df['chr'] = chrom
        
        new_cols = list(X_df.columns[-4:]) + list(X_df.columns[:-4])
        X_df = X_df.loc[:,new_cols]
        X_df = X_df.reset_index()
        X_df = X_df.rename(columns={'index':'donor'})
        
        data_list.append(X_df)
    
    gene_data = pd.concat(data_list)
    print(f'Done for gene {gene_id}')
    return gene_data

#%% Load data and run CNN

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
    
    for gene in chrom_genes:
        
        gene_data = gene_X(gene, annot, exp, attr, data_dir)
        chrom_data.append(gene_data)
    
    chrom_df = pd.concat(chrom_data)
    chrom_df.to_csv(rf'E:/lab_data/gtex_exoneqtl/2d_expression/chr{c}_2d_expression.txt.gz', sep='\t', index=False)
    print(f'Saved chr{c} 2D expression to file')
        

#%%

if __name__ == '__main__':
    
    all_results = []

    # for gene in chrom_genes:
    #     t0 = time.time()
    #     result = run_NN(gene, annot, exp, attr, data_dir)
    #     t1 = time.time()
    #     print(f'Took {np.around((t1-t0)/60, 1)} min for gene {gene}')
    #     all_results.append(result)

   with mp.Pool(num_cpus) as pool:
       
       items = zip(chrom_genes, repeat(annot), repeat(exp), repeat(attr), repeat(data_dir))
       results = pool.starmap(run_NN, items)
       all_results.append(pd.concat(results))
   
    all_results = pd.concat(all_results)
    print(all_results.head())
    all_results.to_csv(os.path.join(data_dir, 'gtex_exoneqtl', 'results', f'{output_filename}.txt.gz'), sep='\t', index=False)
    print('=== DONE ===')    
  
