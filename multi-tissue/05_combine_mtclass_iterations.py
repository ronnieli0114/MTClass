# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 09:17:01 2022

@author: ronnieli
"""

import pandas as pd
import numpy as np
import os, sys

def combine_iterations_eqtl(experiment):
    
    method = 'ensemble'
    idx_list = ['1','2','3']
    
    if experiment == '9_tissue':
        os.chdir(r'G:/My Drive/Lab/lab_projects/gtex_eqtl/results/9_tissue')
        data_list = [pd.read_csv(f'9_tissue_{method}_{idx}.txt.gz', sep='\t', header=0) for idx in idx_list]
        
    elif experiment == 'brain':
        os.chdir(r'G:/My Drive/Lab/lab_projects/gtex_eqtl/results/brain_tissue')
        data_list = [pd.read_csv(f'brain_{method}_{idx}.txt.gz', sep='\t', header=0) for idx in idx_list]

    elif experiment == '48_tissue':
        os.chdir(r'G:/My Drive/Lab/lab_projects/gtex_eqtl/results/48_tissue')
        data_list = [pd.read_csv(f'48_tissue_{method}_{idx}.txt.gz', sep='\t', header=0) for idx in idx_list]

    pairs = []
    for data in data_list:
        data['pair'] = data['gene'] + '-' + data['variant']
        pairs.append(set(data['pair']))
    common_pairs = set.intersection(*pairs)
    print(f'There are {len(common_pairs)} total gene-SNP pairs')
    
    data_list2 = []
    for data in data_list:
        if 'model' in list(data.columns):
            data.drop('model', axis=1, inplace=True)
        data = data[data['pair'].isin(common_pairs)]
        data = data.drop_duplicates('pair')
        data = data.set_index(['gene','variant','pair','sample_size'])
        data_list2.append(data)
    final = pd.concat(data_list2, axis=1)
    
    final['f1_micro_median'] = np.median(final['f1_micro'].values, axis=1)
    final['f1_micro_mean'] = np.average(final['f1_micro'].values, axis=1)
    final['f1_weighted_median'] = np.median(final['f1_weighted'].values, axis=1)
    final['f1_weighted_mean'] = np.average(final['f1_weighted'].values, axis=1)
    final['f1_macro_median'] = np.median(final['f1_macro'].values, axis=1)
    final['f1_macro_mean'] = np.average(final['f1_macro'].values, axis=1)
    final['mcc_median'] = np.median(final['mcc'].values, axis=1)
    final['mcc_mean'] = np.average(final['mcc'].values, axis=1)
    
    final = final.reset_index()
    cols = ['gene','variant','sample_size','pair','f1_micro_median','f1_micro_mean',
            'f1_macro_median','f1_macro_mean','f1_weighted_median','f1_weighted_mean',
            'mcc_median','mcc_mean']
    final = final[cols]
    print('Aggregated results. Writing to file...')
    final.to_csv(f'{experiment}_{method}_aggregate.txt.gz', sep='\t', index=False)
    
    return final