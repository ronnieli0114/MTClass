# -*- coding: utf-8 -*-
"""
Main script to load the data and run 
the classification model for one tissue

Created on Fri Sep  2 18:54:11 2022
@author: ronnieli
"""

num_iter = 1
chr_range = range(1,23)
model_name = 'RF'
num_cpus = 4
computer = 'pc' # determines location of data
how = 'parallel' # or 'serial', determines processing type

import pandas as pd
import numpy as np
import os
import time
import sys
import pyarrow.feather as feather
import multiprocessing as mp
import random
from itertools import repeat
from sklearn.metrics import f1_score, roc_auc_score, matthews_corrcoef, average_precision_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, VotingClassifier

seed = random.randint(1,1001)
print(f'Randomly generated seed is {seed}')

def load_data(tissue, chrom, computer='pc'):
    
    if computer == 'pc':
        count_dir = r'/mnt/e/lab_data/gtex_exoneqtl/tissue_level_expression'
        gt_dir = r'/mnt/e/lab_data/gtex_exoneqtl/genotypes'
        sys.path.append(r'/mnt/g/My\ Drive/Lab/lab_projects/gtex_exoneqtl/scripts')
    elif computer == 'cluster':
        count_dir = r'/projects/compbio/users/yli459/gtex_exoneqtl/tissue_level_expression'
        gt_dir = r'/projects/compbio/users/yli459/gtex_exoneqtl/genotypes'
        sys.path.append(r'/projects/compbio/users/yli459/gtex_exoneqtl')
     
    counts = feather.read_feather(os.path.join(count_dir, f'{tissue}_counts.ftr'))
    gene_list = pd.read_csv(os.path.join(count_dir, f'{tissue}_gene_list.txt'), header=2)['gene_name'].tolist()
    counts = counts[counts['gene'].isin(gene_list)]

    # get genotypes for chromosome
    gt = pd.read_csv(os.path.join(gt_dir, f'chr{chrom}_genotypes_10kb_binary.txt.gz'), sep='\t', header=0)
    
    # only keep genes in genotypes
    gt = gt[gt['gene_name'].isin(counts['gene'])]
    
    # filter counts by donors available in genotype
    counts = counts[counts['donor'].isin(gt.columns)]
    
    # only keep donors in the counts
    keep_cols = ['gene_name','ID'] + sorted(list(set(counts['donor'])))
    gt = gt.filter(items=keep_cols, axis=1)
    
    num_genes = len(set(gt['gene_name']))
    
    # define counts and labels
    print(f'There are a total of {num_genes} genes in chr{chrom}')
    return counts, gt

def train_model(gene, model_name, counts, gt, seed):

    # hyperparameters
    C = 1
    n_splits = 4
    n_estimators = 150
    rf_max_depth = 5
    svm_kernel = 'rbf'
    
    results_dict = {'model':[], 'variant':[], 'gene':[], 'f1_micro':[], 'f1_macro':[], 'f1_weighted':[], 'mcc':[]}
    
    counts_gene = counts[counts['gene']==gene]
    counts_gene = counts_gene.dropna(how='any', axis=1)   # drop missing data (exons not present for a gene)
    counts_gene = counts_gene.drop('gene', axis=1).set_index('donor')
    num_donors = len(set(counts_gene.index))

    gt_gene = gt[gt['gene_name']==gene]
    num_vars = len(set(gt_gene['ID']))
    
    for snp in gt_gene['ID']:
        
        y = gt_gene[gt_gene['ID']==snp]
        y = y.drop_duplicates(['gene_name','ID'])
        y = y.iloc[:,2:].T
        y.columns = ['label']
        
        Xy = counts_gene.merge(y, left_index=True, right_index=True)
        
        # drop donors with missing genotype
        Xy.replace({9:np.nan}, inplace=True)
        Xy.dropna(how='any', axis=0, inplace=True)
        
        if Xy.shape[0] < 0.9 * num_donors:
            print(f'Sample size too small for {gene}-{snp}: N={Xy.shape[0]}', flush=True)
            continue

        pred_vars = ['label']
        X = Xy.loc[:, [col for col in Xy.columns if col not in pred_vars]]

        X = X.to_numpy().astype(float)
        y = Xy['label'].to_numpy().astype(int)

        if np.count_nonzero(y==1) < n_splits or np.count_nonzero(y==0) < n_splits:
            print(f'MAF too small for {gene}-{snp}', flush=True)
            continue

        # initialize models
        svc = SVC(kernel=svm_kernel, C=C, probability=True, random_state=seed)
        rf = RandomForestClassifier(n_estimators=n_estimators, max_depth=rf_max_depth, random_state=seed)
        skf = StratifiedKFold(n_splits=n_splits, random_state=seed, shuffle=True)
        scaler = StandardScaler()

        if model_name == 'voting_hard':
            model = VotingClassifier(estimators=[('SVM', svc), ('RF', rf)], voting='hard')
        elif model_name == 'voting_soft':
            model = VotingClassifier(estimators=[('SVM', svc), ('RF', rf)], voting='soft')
        elif model_name == 'SVM':
            model = svc
        elif model_name == 'RF':
            model = rf

        for train_ind, test_ind in skf.split(X, y):

            X_train, y_train = X[train_ind], y[train_ind]
            X_test, y_test = X[test_ind], y[test_ind]

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

        results_dict['model'].append(model_name)
        results_dict['gene'].append(gene)
        results_dict['variant'].append(snp)
        results_dict['f1_micro'].append(round(np.mean(f1_micro_lst), 3))
        results_dict['f1_macro'].append(round(np.mean(f1_macro_lst), 3))
        results_dict['f1_weighted'].append(round(np.mean(f1_weighted_lst), 3))
        results_dict['mcc'].append(round(np.mean(mcc_lst), 3))

    results_df = pd.DataFrame(results_dict)
    return results_df

#%% Run the model on all chromosomes

# data_dir = r'/projects/compbio/users/yli459/gtex_exoneqtl/tissue_level_expression'
data_dir = r'/mnt/e/lab_data/gtex_exoneqtl/tissue_level_expression'
tissue_list = [tissue.split('_counts.ftr')[0] for tissue in os.listdir(data_dir) if 'Brain_' in tissue and 'gene_list' not in tissue]

for tissue in tissue_list:
    
    print(f'Now running on {tissue}')
    
    for c in chr_range:
            
        print(f'Now running on chr{c}...')
        t0 = time.time()
       
        counts, gt = load_data(tissue, c, computer='pc') 
        print('=== Loaded data ===')
           
        gene_list = sorted(list(set(gt['gene_name'])))
        num_genes = len(gene_list)
        # results_dir = r'/projects/compbio/users/yli459/gtex_exonqtl/mtclass_results'
        results_dir = r'/mnt/g/My Drive/Lab/lab_projects/gtex_exoneqtl/results'
        results = []
            
        ### serial processing
        if how == 'serial':
            print('Running in serial')
            for idx, gene in enumerate(gene_list):
                print(f'Now on #{idx} gene {gene}...')
                result = train_model(gene, model_name, counts, gt, seed)
                results.append(result)
    
        ### parallel processing
        elif how == 'parallel': 
            print('Running in parallel')
            if __name__ == '__main__':
                
                with mp.Pool(processes=num_cpus) as pool:        
                    
                    items = zip(gene_list, repeat(model_name), repeat(counts), repeat(gt), repeat(seed))
                    result = pool.starmap(train_model, items)
                    results.append(pd.concat(result))
        
        final_results = pd.concat(results, axis=0)
        print(final_results.head())
            
        if not os.path.exists(results_dir):
            os.mkdir(results_dir)
        feather.write_feather(final_results, os.path.join(results_dir, f'{tissue}_chr{c}_{model_name}_{num_iter}.ftr'))
            
        t1 = time.time()
        print(f'Took {(t1-t0)/3600:.1f} hours to run on chr{c}. Wrote results to feather')
        
print('=== DONE ===')
