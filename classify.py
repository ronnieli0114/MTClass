# -*- coding: utf-8 -*-
"""
Main MTClass script
@author: ronnieli
"""

import pandas as pd
import numpy as np
import os, sys
import random
import multiprocessing as mp
from itertools import repeat
from argparse import ArgumentParser
from sklearn.metrics import f1_score, matthews_corrcoef, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, VotingClassifier

parser = ArgumentParser(prog='MTClass classifier', description='Ensemble classifier for gene-SNP pairs')
parser.add_argument('exp', help='Path to gene expression file', type=str)
parser.add_argument('gt', help='Path to genotypes file', type=str)
parser.add_argument('out_dir', help='Path to output file', type=str)
parser.add_argument('-t', '--num_threads', help='Number of threads to use', default=1, type=int)
parser.add_argument('-v','--verbose', action='store_true', help='Prints verbose output')
args = parser.parse_args()

exp = pd.read_csv(args.exp, sep='\t', header=0)
gt = pd.read_csv(args.gt, sep='\t', header=0)
seed = random.randint(1,1001)
print(f'Randomly generated seed is {seed}')
print(f'Using {args.num_threads} threads')

genes = sorted(list(set(exp['gene']) & set(gt['gene'])))
print(f'There are {len(genes)} genes in common between expression and genotype files')

def classify(gene, exp, gt, seed):
    
    # empty dictionary to store results
    results_dict = {'gene':[], 'variant':[], 'sample_size':[], 'auc':[], 'f1_micro':[], 'f1_macro':[], 'f1_weighted':[], 'mcc':[]}

    # hyperparameters
    C = 1
    n_splits = 4
    n_estimators = 150
    rf_max_depth = 5
    svm_kernel = 'rbf' 

    print(f'Now running on {gene}')
    exp_gene = exp[exp['gene']==gene]
    exp_gene = exp_gene.set_index('donor').drop('gene', axis=1)
    exp_gene = exp_gene.dropna(how='all', axis=1)
    gt_gene = gt[gt['gene']==gene]

    print(f'There are {len(exp_gene.columns)} features')

    for snp in gt_gene['ID']:
        
        y = gt_gene[gt_gene['ID']==snp]
        y = y.iloc[:,2:].T
        y.columns = ['label']
        y.replace({9:np.nan}, inplace=True)
        
        Xy = exp_gene.merge(y, left_index=True, right_index=True)
        
        # drop donors with missing genotype
        Xy.dropna(how='any', axis=0, subset='label', inplace=True)
        sample_size = Xy.shape[0]
        
        pred_vars = ['label']
        X = Xy.loc[:, [x for x in Xy.columns if x not in pred_vars]]
        X = X.to_numpy().astype(float)
        y = Xy['label'].to_numpy().astype(int)

        # require at least 5 samples of each class
        if np.count_nonzero(y==0) < n_splits or np.count_nonzero(y==1) < n_splits:
            print(f'MAF too small for {gene}-{snp}')
            continue

        # initialize models
        svc = SVC(kernel=svm_kernel, C=C, probability=True, random_state=seed)
        rf = RandomForestClassifier(n_estimators=n_estimators, max_depth=rf_max_depth, random_state=seed)
        skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        scaler = StandardScaler()
        model = VotingClassifier([('SVC',svc),('RF',rf)], voting='soft')

        auc_lst = []
        f1_micro_lst = []
        f1_macro_lst = []
        f1_weighted_lst = []
        mcc_lst = []

        for train_ind, test_ind in skf.split(X, y):

            X_train, y_train = X[train_ind], y[train_ind]
            X_test, y_test = X[test_ind], y[test_ind]

            # standardize the data
            X_train = scaler.fit_transform(X_train)
            X_test = scaler.transform(X_test)
        
            model.fit(X_train, y_train)

            y_pred = model.predict(X_test)
            auc_lst.append(roc_auc_score(y_test, y_pred))
            f1_micro_lst.append(f1_score(y_test, y_pred, average='micro'))
            f1_macro_lst.append(f1_score(y_test, y_pred, average='macro'))
            f1_weighted_lst.append(f1_score(y_test, y_pred, average='weighted'))
            mcc_lst.append(matthews_corrcoef(y_test, y_pred))

        if args.verbose:
            f1 = np.around(np.mean(f1_macro_lst), 3)
            print(f'Gene {gene}\tSNP {snp}\tMacro F1 = {f1}')

        results_dict['gene'].append(gene)
        results_dict['variant'].append(snp)
        results_dict['sample_size'].append(sample_size)
        results_dict['auc'].append(np.around(np.mean(auc_lst), 3))
        results_dict['f1_micro'].append(np.around(np.mean(f1_micro_lst), 3))
        results_dict['f1_macro'].append(np.around(np.mean(f1_macro_lst), 3))
        results_dict['f1_weighted'].append(np.around(np.mean(f1_weighted_lst), 3))
        results_dict['mcc'].append(np.around(np.mean(mcc_lst), 3))

    results_df = pd.DataFrame(results_dict)
    results_df.sort_values('f1_macro', ascending=False, inplace=True)
    return results_df

if __name__ == '__main__':
    
    with mp.Pool(args.num_threads) as pool:
    
        inputs = zip(genes, repeat(exp), repeat(gt), repeat(seed))
        result_list = pool.starmap(classify, inputs)
        results = pd.concat(result_list)
        results.to_csv(os.path.join(args.out_dir, 'mtclass_results.txt'), sep='\t', index=False)