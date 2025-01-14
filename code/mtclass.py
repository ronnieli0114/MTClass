# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 14:25:30 2022

@author: ronnieli
"""

import pandas as pd
import numpy as np
import os, sys
import time
import gc
import pyarrow.feather as feather
import multiprocessing as mp
import random
from itertools import repeat
from sklearn.metrics import f1_score, roc_auc_score, matthews_corrcoef, average_precision_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, VotingClassifier

# Define function to classify
def classify_one_variant(g, v, exp_gene, gt_gene, seed=2024, permute=False):
    
    results_dict = {'gene':[], 'variant':[], 'sample_size':[], 'f1_micro':[], 'f1_macro':[], 'f1_weighted':[], 'mcc':[]}
    
    y = gt_gene[gt_gene['ID']==v]
    y = y.drop_duplicates(['gene_name','ID'])
    y = y.iloc[:,3:].T
    y.columns = ['label']
    
    Xy = exp_gene.merge(y, left_index=True, right_index=True)
   
    # drop donors with missing genotype
    Xy['label'] = Xy['label'].replace(9, np.nan)
    Xy.dropna(how='any', subset='label', axis=0, inplace=True)
    sample_size = Xy.shape[0]

    X = Xy.drop('label', axis=1)
    
    X_numpy = X.to_numpy(dtype=float)
    y_numpy = Xy['label'].to_numpy(dtype=int)

    if permute:
        np.random.seed(seed)
        np.random.shuffle(y_numpy)

    if np.count_nonzero(y_numpy==1) < n_splits or np.count_nonzero(y_numpy==0) < n_splits:
        print(f'MAF too small for {gene}-{snp}', flush=True)
        return

    # initialize models
    svc = SVC(kernel=svm_kernel, C=C, probability=True, random_state=seed)
    rf = RandomForestClassifier(n_estimators=n_estimators, max_depth=rf_max_depth, random_state=seed)
    skf = StratifiedKFold(n_splits=n_splits, random_state=seed, shuffle=True)
    scaler = StandardScaler()
    model = VotingClassifier(estimators=[('SVM', svc), ('RF', rf)], voting='soft')

    for train_ind, test_ind in skf.split(X_numpy, y_numpy):

        X_train, y_train = X_numpy[train_ind], y_numpy[train_ind]
        X_test, y_test = X_numpy[test_ind], y_numpy[test_ind]

        # standardize the data
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)

        f1_micro_lst = []
        f1_macro_lst = []
        f1_weighted_lst = []
        mcc_lst = []

        model.fit(X_train, y_train)
    
        y_pred = model.predict(X_test)
        y_proba = model.predict_proba(X_test)
        f1_micro_lst.append(f1_score(y_test, y_pred, average='micro'))
        f1_macro_lst.append(f1_score(y_test, y_pred, average='macro'))
        f1_weighted_lst.append(f1_score(y_test, y_pred, average='weighted'))
        mcc_lst.append(matthews_corrcoef(y_test, y_pred))

    results_dict['gene'].append(g)
    results_dict['variant'].append(v)
    results_dict['sample_size'].append(sample_size)
    results_dict['f1_micro'].append(np.around(np.mean(f1_micro_lst), 3))
    results_dict['f1_macro'].append(np.around(np.mean(f1_macro_lst), 3))
    results_dict['f1_weighted'].append(np.around(np.mean(f1_weighted_lst), 3))
    results_dict['mcc'].append(np.around(np.mean(mcc_lst), 3))

    results_df = pd.DataFrame(results_dict)
    gc.collect()
    return results_df

if __name__ == '__main__':

    # hyperparameters
    C = 1
    n_splits = 4
    n_estimators = 150
    rf_max_depth = 5
    svm_kernel = 'rbf'    
    num_cpus = 12
    permute = False

    data_dir = r'/mnt/d/lab_data/mtclass_eqtl/multi-tissue'
    results_dir = r'/mnt/c/Users/ronni/OneDrive/Documents/mtclass/'

    experiment_list = ['9_tissue','brain_tissue','48_tissue']
    for experiment in experiment_list:
        
        # run MTClass 3 times
        for i in [1,2,3]:

            print(f'Running MTClass on the {experiment} case, run #{i}', flush=True)
            seed = i

            # for each chromosome
            for c in range(1,23):
                
                gt = pd.read_csv(os.path.join(data_dir, f'genotypes/chr{c}_genotypes_10kb_binary.txt.gz'), sep='\t', header=0)
                gt = gt[~pd.isnull(gt['gene_name'])]
                gt_genes = sorted(list(set(gt['gene_name'])))

                if experiment=='9_tissue':
                    exp_path = os.path.join(data_dir, f'{experiment}/TPM_expression.ftr')
                else:
                    exp_path = os.path.join(data_dir, f'{experiment}/PMM_imputed_expression.ftr')
                exp_chrom = feather.read_feather(exp_path)
                if 'gene_name' in exp_chrom.columns:
                    exp_chrom = exp_chrom.rename(columns={'gene_name':'gene'})
                if 'gene_id' in exp_chrom.columns:
                    exp_chrom = exp_chrom.drop('gene_id', axis=1)
                exp_chrom = exp_chrom[exp_chrom['gene'].isin(gt_genes)]
                print(f'loaded data for chr{c}', flush=True)

                out_name = f'{experiment}_{i}_chr{c}'
                genes = sorted(list(set(exp_chrom['gene'])))
                print(f'There are {len(genes)} genes in chr{c}', flush=True)
                
                chr_results_list = []
                
                for j, g in enumerate(genes):
                    
                    exp_gene = exp_chrom[exp_chrom['gene']==g]
                    exp_gene = exp_gene.dropna(how='all', axis=1) # drop columns with all missing values
                    exp_gene = exp_gene.drop('gene', axis=1).set_index('donor')
                    
                    gt_gene = gt[gt['gene_name']==g]
                    variant_list = gt_gene['ID'].tolist()
                    print(f'\t there are {len(variant_list)} variants in #{j} {g}')
                    
                    with mp.Pool(num_cpus) as pool:
                        
                        inputs = zip(repeat(g), variant_list, repeat(exp_gene), repeat(gt_gene), repeat(seed), repeat(permute))
                        result_list = pool.starmap(classify_one_variant, inputs)
                        g_result = pd.concat(result_list)
                        chr_results_list.append(g_result)
                
                chr_results = pd.concat(chr_results_list)
                chr_results.to_csv(os.path.join(results_dir, f'{out_name}.txt'), sep='\t', index=False, header=True)
                print(f'{experiment}\twrote chr{c} to file', flush=True)
                gc.collect()