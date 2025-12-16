# -*- coding: utf-8 -*-
"""
main MTClass classification algorithm
@author: ronnieli
"""

import pandas as pd
import numpy as np
import gc
import multiprocessing as mp
from itertools import repeat
from sklearn.metrics import f1_score, matthews_corrcoef
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, VotingClassifier
from sklearn.utils import shuffle

# Define function to classify one variant based on gene expression and genotypes
def classify_one_variant(g, v, expression_df, genotype_df, seed=2024, 
                        permutation_test=False, n_permutations=10000):

    '''
    Runs the MTClass ensemble classification algorithm on one gene-variant pair.
    Inputs:
        - g: gene name
        - v: variant ID
        - expression_df: pandas DataFrame of expression values for gene g, in the format of samples x features.
        - genotype_df: pandas DataFrame of genotypes, in the format of variants x samples
        - seed: random seed for initializing weights, default=2024
        - permutation_test: whether to run permutation test for effect size and p-value, default=False
        - n_permutations: number of permutations for permutation test, default=10000
    Returns:
        - a pandas DataFrame with classification results for gene-variant pair g;v. Metrics are F1 scores 
        (macro, micro, weighted), MCC, and optionally effect size and p-value
    '''

    results_dict = {'gene':[], 'variant':[], 'sample_size':[], 'f1_micro':[], 
                    'f1_macro':[], 'f1_weighted':[], 'mcc':[]}
    
    if permutation_test:
        results_dict['effect_size'] = []
        results_dict['p_value'] = []

    # hyperparameters
    C = 1
    n_splits = 4
    n_estimators = 150
    rf_max_depth = 5
    svm_kernel = 'rbf' 

    y = genotype_df.loc[v,:].T.to_frame()
    y.columns = ['label']
    
    Xy = expression_df.merge(y, left_index=True, right_index=True)
   
    # drop donors with missing genotype
    Xy['label'] = Xy['label'].replace(9, np.nan)
    Xy.dropna(how='any', subset='label', axis=0, inplace=True)
    sample_size = Xy.shape[0]

    X = Xy.drop('label', axis=1)
    
    X_numpy = X.to_numpy(dtype=float)
    y_numpy = Xy['label'].to_numpy(dtype=int)

    if np.count_nonzero(y_numpy==1) < n_splits or np.count_nonzero(y_numpy==0) < n_splits:
        print(f'MAF too small for {g}-{v}', flush=True)
        return
        
    # initialize models
    svc = SVC(kernel=svm_kernel, C=C, probability=True, random_state=seed)
    rf = RandomForestClassifier(n_estimators=n_estimators, max_depth=rf_max_depth, random_state=seed)
    skf = StratifiedKFold(n_splits=n_splits, random_state=seed, shuffle=True)
    scaler = StandardScaler()
    model = VotingClassifier(estimators=[('SVM', svc), ('RF', rf)], voting='soft')

    # K-fold cross-validation for classification metrics
    f1_micro_lst = []
    f1_macro_lst = []
    f1_weighted_lst = []
    mcc_lst = []

    for train_ind, test_ind in skf.split(X_numpy, y_numpy):

        X_train, y_train = X_numpy[train_ind], y_numpy[train_ind]
        X_test, y_test = X_numpy[test_ind], y_numpy[test_ind]

        # standardize the data
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)

        model.fit(X_train, y_train)
    
        y_pred = model.predict(X_test)
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

    # Permutation test for effect size and p-value
    if permutation_test:
        # Reset scaler and model for permutation test
        scaler_perm = StandardScaler()
        svc_perm = SVC(kernel=svm_kernel, C=C, probability=True, random_state=seed)
        rf_perm = RandomForestClassifier(n_estimators=n_estimators, max_depth=rf_max_depth, random_state=seed)
        model_perm = VotingClassifier(estimators=[('SVM', svc_perm), ('RF', rf_perm)], voting='soft')
        
        # Single train-test split for permutation test
        X_train, X_test, y_train, y_test = train_test_split(
            X_numpy, y_numpy, test_size=0.20, random_state=seed, stratify=y_numpy
        )
        
        X_train = scaler_perm.fit_transform(X_train)
        X_test = scaler_perm.transform(X_test)
        
        model_perm.fit(X_train, y_train)
        y_pred = model_perm.predict(X_test)
        
        actual_score = matthews_corrcoef(y_test, y_pred)
        
        count = 0
        permuted_scores = []
        
        for _ in range(n_permutations):
            y_test_permuted = shuffle(y_test, random_state=None)
            permuted_score = matthews_corrcoef(y_test_permuted, y_pred)
            permuted_scores.append(permuted_score)
            if permuted_score >= actual_score:
                count += 1
        
        p_value = (count + 1) / (n_permutations + 1)
        effect_size = (actual_score - np.mean(permuted_scores)) / np.max([np.std(permuted_scores), 1/n_permutations])
        
        results_dict['effect_size'].append(np.around(effect_size, 3))
        results_dict['p_value'].append(p_value)

    results_df = pd.DataFrame(results_dict)
    gc.collect()
    return results_df

def run_mtclass(expression_df, genotype_df, iterations=1, num_cores=4, 
                permutation_test=False, n_permutations=10000):

    '''
    Runs the MTClass algorithm using parallel processing (multiprocessing package).
    Inputs:
        - expression_df: pandas DataFrame of gene expression. Can have multiple genes. 
            Needs to have a 'gene' and 'donor' column. The other columns are features (sources of gene expression,
            like multiple tissues or exons).
        - genotype_df: pandas DataFrame of genotypes. Can be made from extract_eqtl_genotypes(). 
            Needs to have a 'gene_name' and 'ID' column. The other columns are the donors, similar to a VCF file.
        - iterations: how many iterations with different random seeds to run MTClass. Default=1.
        - num_cores: how many threads to use in parallel. Default=4.
        - permutation_test: whether to run permutation test for effect size and p-value, default=False
        - n_permutations: number of permutations for permutation test, default=10000
    '''

    # ensure testing on same genes
    common_genes = sorted(list(set(expression_df['gene']) & set(genotype_df['gene_name'])))
    expression_df = expression_df[expression_df['gene'].isin(common_genes)]
    genotype_df = genotype_df[genotype_df['gene_name'].isin(common_genes)]

    # ensure testing on common donors
    common_donors = sorted(list(set(expression_df['donor']) & set(genotype_df.columns)))
    expression_df = expression_df[expression_df['donor'].isin(common_donors)]
    genotype_df = genotype_df.filter(['gene_name','ID']+common_donors)
    print(f'found {len(common_genes)} genes and {len(common_donors)} samples for testing')

    # check number of iterations
    assert iterations >= 1 and num_cores > 0
    
    if permutation_test:
        print(f'Running permutation test with {n_permutations} permutations')
    
    # run MTClass once for each iteration
    result_list = []
    for i in range(iterations):

        print(f'running MTClass iteration {i+1} using {num_cores} cores')
        seed = i+1

        for g in common_genes:
            X = expression_df[expression_df['gene']==g]
            X = X.drop('gene', axis=1).set_index('donor')
            y = genotype_df[genotype_df['gene_name']==g]
            y = y.drop('gene_name', axis=1).set_index('ID')
            y = y[~y.index.duplicated()]
            variants = y.index.tolist()

            with mp.Pool(num_cores) as pool:
                inputs = zip(repeat(g), variants, repeat(X), repeat(y), repeat(seed),
                           repeat(permutation_test), repeat(n_permutations))
                res = pool.starmap(classify_one_variant, inputs)
                res_df = pd.concat(res, axis=0)
                res_df['iteration'] = seed
                result_list.append(res_df)

    final_res = pd.concat(result_list)
    return final_res

def test_mtclass(permutation_test=False):

    '''
    Runs one iteration of MTClass using test data obtained from Ronnie Li's GitHub repository.
    Inputs: 
        - permutation_test: whether to run permutation test, default=False
    Returns: a pandas DataFrame with results on test data
    '''

    expression_df = pd.read_csv('https://github.com/ronnieli0114/MTClass/raw/refs/heads/main/data/test_expression.tsv.gz', sep='\t', header=0)
    genotype_df = pd.read_csv('https://github.com/ronnieli0114/MTClass/raw/refs/heads/main/data/test_genotypes.tsv.gz', sep='\t', header=0)
    res = run_mtclass(expression_df=expression_df, genotype_df=genotype_df, iterations=1, 
                     num_cores=4, permutation_test=permutation_test, n_permutations=1000)
    return res
