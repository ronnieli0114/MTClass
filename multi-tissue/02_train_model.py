# -*- coding: utf-8 -*-
"""
Define function to train model
Voting classifier
"""

import pandas as pd
import numpy as np
from sklearn.metrics import f1_score, roc_auc_score, average_precision_score, matthews_corrcoef
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, VotingClassifier

def classify(gene, model_name, counts, gt, seed):

    # hyperparameters
    C = 1
    n_splits = 4
    n_estimators = 150
    rf_max_depth = 5
    svm_kernel = 'rbf'    

    # empty dictionary to store results
    results_dict = {'model':[], 'gene':[], 'variant':[], 'sample_size':[], 'f1_micro':[], 'f1_macro':[], 'f1_weighted':[], 'mcc':[]}
    
    counts_gene = counts[counts['gene']==gene]
    gt_gene = gt[gt['gene_name']==gene]
    
    num_donors = len(set(counts_gene['donor']))
    num_vars = len(set(gt_gene['ID']))
    
    for snp in gt_gene['ID']:
        
        y = gt_gene[gt_gene['ID']==snp]
        y = y.iloc[:,2:].T
        y.columns = ['label']
        y.replace({9:np.nan}, inplace=True)
        
        Xy = counts_gene.merge(y, left_on='donor', right_index=True)
        
        # drop donors with missing genotype
        Xy.dropna(how='any', axis=0, subset='label', inplace=True)
        sample_size = Xy.shape[0]
        
        pred_vars = ['gene','donor','label']
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
        
        if model_name == 'ensemble':
            model = VotingClassifier(estimators=[('SVM', svc), ('RF', rf)], voting='soft')
        elif model_name == 'SVM':
            model = svc
        elif model_name == 'RF':
            model = rf
        else:
            raise KeyError("Model name must be 'ensemble', 'SVM', or 'RF'!")

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
            f1_micro_lst.append(f1_score(y_test, y_pred, average='micro'))
            f1_macro_lst.append(f1_score(y_test, y_pred, average='macro'))
            f1_weighted_lst.append(f1_score(y_test, y_pred, average='weighted'))
            mcc_lst.append(matthews_corrcoef(y_test, y_pred))

        results_dict['model'].append(model_name)
        results_dict['gene'].append(gene)
        results_dict['variant'].append(snp)
        results_dict['sample_size'].append(sample_size)
        results_dict['f1_micro'].append(round(np.mean(f1_micro_lst), 3))
        results_dict['f1_macro'].append(round(np.mean(f1_macro_lst), 3))
        results_dict['f1_weighted'].append(round(np.mean(f1_weighted_lst), 3))
        results_dict['mcc'].append(round(np.mean(mcc_lst), 3))

    results_df = pd.DataFrame(results_dict)
    return results_df
