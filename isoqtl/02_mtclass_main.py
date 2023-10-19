# -*- coding: utf-8 -*-
"""
MTClass main script on 
PsychENCODE data

Created on Mon Jan 30 16:39:27 2023
@author: ronnieli
"""

import pandas as pd
import numpy as np
import pyarrow.feather as feather
import os, sys
import random
import multiprocessing as mp
from itertools import repeat
from sklearn.metrics import f1_score, roc_auc_score, average_precision_score, matthews_corrcoef
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, VotingClassifier

def train_model(gene, model_name, exp, gt, seed):
    
    # hyperparameters
    C = 1
    n_splits = 4
    n_estimators = 150
    rf_max_depth = 5
    svm_kernel = 'rbf'

    results_dict = {'gene':[], 'variant':[], 'f1_micro':[], 'f1_macro':[], 'f1_weighted':[], 'mcc':[]}
    
    exp_gene = exp[exp['gene_name']==gene]
    gt_gene = gt[gt['gene_name']==gene]
    
    num_donors = len(set(exp_gene['genotypingID']))
    num_vars = len(set(gt_gene['ID']))
    print(f'There are {num_vars} variants for gene {gene}')
    
    for snp in gt_gene['ID']:
        
        y = gt_gene[gt_gene['ID']==snp]
        y = y.drop_duplicates('ID')
        y = y.iloc[:,3:].T
        y.columns = ['label']
        y.replace({9:np.nan}, inplace=True)
        
        X = exp_gene.drop(['gene_id','gene_name','individualID','specimenID'], axis=1).set_index('genotypingID')
        X = X.dropna(how='all', axis=1)
        Xy = X.merge(y, left_index=True, right_index=True)
        
        # drop donors with missing genotype
        Xy.dropna(how='any', axis=0, subset='label', inplace=True)
        
        if Xy.shape[0] < 0.8 * num_donors:
            print(f'Sample size too small for {gene}-{snp}, N={Xy.shape[0]}')
        
        X_numpy = Xy.iloc[:,:-1].to_numpy(dtype=float)
        y_numpy = Xy.iloc[:,-1].to_numpy(dtype=int)

        if np.count_nonzero(y_numpy==0) < n_splits or np.count_nonzero(y_numpy==1) < n_splits:
            print(f'MAF too small for {gene}-{snp}')

        # initialize models
        svc = SVC(kernel=svm_kernel, C=C, probability=True, random_state=seed)
        rf = RandomForestClassifier(n_estimators=n_estimators, max_depth=rf_max_depth, random_state=seed)
        
        skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        scaler = StandardScaler()

        if model_name == 'voting_soft':
            model = VotingClassifier(estimators=[('SVM', svc), ('RF', rf)], voting='soft')
        elif model_name == 'SVM':
            model = svc
        elif model_name == 'RF':
            model = rf

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
            
            try:
                model.fit(X_train, y_train)
            except:
                continue

            y_pred = model.predict(X_test)
            f1_micro_lst.append(f1_score(y_test, y_pred, average='micro'))
            f1_macro_lst.append(f1_score(y_test, y_pred, average='macro'))
            f1_weighted_lst.append(f1_score(y_test, y_pred, average='weighted'))
            mcc_lst.append(matthews_corrcoef(y_test, y_pred))

        results_dict['gene'].append(gene)
        results_dict['variant'].append(snp)
        results_dict['f1_micro'].append(round(np.mean(f1_micro_lst), 3))
        results_dict['f1_macro'].append(round(np.mean(f1_macro_lst), 3))
        results_dict['f1_weighted'].append(round(np.mean(f1_weighted_lst), 3))
        results_dict['mcc'].append(round(np.mean(mcc_lst), 3))

        # avg_f1 = round(np.mean(f1_macro_lst),3)
        # avg_mcc = round(np.mean(mcc_lst),3)
        # print(f'SNP {snp} \t macro F1={avg_f1} \t MCC={avg_mcc}')

    results_df = pd.DataFrame(results_dict)
    return results_df


def run_mtclass(computer, out_name, model='ensemble'):
    
    if computer == 'pc':
        data_dir = r'/mnt/e/lab_data/isoqtl/for_mtclass'
        results_dir = r'/mnt/g/My Drive/Lab/lab_projects/isoqtl/results'
        num_cpus = 4
    elif computer == 'cluster':
        data_dir = r'/projects/compbio/users/yli459/isoqtl/data'
        results_dir = r'/projects/compbio/users/yli459/isoqtl/results'
        num_cpus = 8
        
    all_results = []

    for c in range(1,23):
        
        chr_results = []
        
        exp = feather.read_feather(os.path.join(data_dir, 'all_expression.ftr'))
        exp = exp[~exp['gene_name'].isna()]
        gt = pd.read_csv(os.path.join(data_dir, 'genotypes', f'chr{c}_genotypes_binary.txt.gz'), sep='\t', header=0)
        
        exp = exp[exp['gene_name'].isin(gt['gene_name'])]
        
        print(f'=== Loaded data for chr{c} ===')
        
        if __name__ == '__main__':
            
            seed = random.randint(1, 1001)
            print(f'Randomly generated seed is {seed}')
            
            with mp.Pool(num_cpus) as pool:
                
                genes = sorted(list(set(exp['gene_name'])))
                print(f'There are {len(genes)} genes in chr{c}')
                
                if model == 'ensemble':
                    inputs = zip(genes, repeat('voting_soft'), repeat(exp), repeat(gt), repeat(seed))
                elif model == 'RF':
                    inputs = zip(genes, repeat('RF'), repeat(exp), repeat(gt), repeat(seed))
                
                result_list = pool.starmap(train_model, inputs)
                chr_results = pd.concat(result_list)
                all_results.append(chr_results)
    
    all_results = pd.concat(all_results)
    all_results.to_csv(os.path.join(results_dir, f'{out_name}.txt.gz'), sep='\t', index=False, header=True)
    
#%%
run_mtclass('cluster','mtclass_ensemble_1',model='ensemble')
print('=== DONE ===')