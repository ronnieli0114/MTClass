# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 09:17:01 2022

@author: ronnieli
"""

import pandas as pd
import numpy as np
import os

def combine_iterations_exoneqtl(tissue):
    
    os.chdir(r'G:/My Drive/Lab/lab_projects/gtex_exoneqtl/results/')
        
    # Read in data from iterations
    method = 'ensemble'
    idx_list = [1,2,3]
    data_list = [pd.read_csv(f'{tissue}_{method}_{idx}.txt.gz', sep='\t', header=0) for idx in idx_list]
    print('Read in data. Determining common set of gene-SNP pairs')
    
    # Determine common gene-SNP pairs
    pairs = []
    for data in data_list:
        data['pair'] = data['gene'] + ';' + data['variant']
        pairs.append(set(data['pair']))
    common_pairs = set.intersection(*pairs)
        
    data_list2 = []
    for data in data_list:
        data = data[data['pair'].isin(common_pairs)]
        data = data.drop_duplicates('pair')
        data = data.set_index(['model','gene','variant','pair','sample_size'])
        data_list2.append(data)
    final = pd.concat(data_list2, axis=1)
    print('Created final dataframe')
    
    final['f1_micro_median'] = np.median(final['f1_micro'].values, axis=1)
    final['f1_micro_mean'] = np.average(final['f1_micro'].values, axis=1)
    final['f1_weighted_median'] = np.median(final['f1_weighted'].values, axis=1)
    final['f1_weighted_mean'] = np.average(final['f1_weighted'].values, axis=1)
    final['f1_macro_median'] = np.median(final['f1_macro'].values, axis=1)
    final['f1_macro_mean'] = np.average(final['f1_macro'].values, axis=1)
    final['mcc_median'] = np.median(final['mcc'].values, axis=1)
    final['mcc_mean'] = np.average(final['mcc'].values, axis=1)
    
    final = final.reset_index()
    cols = ['gene','variant','pair','sample_size','f1_micro_median','f1_micro_mean',
            'f1_macro_median','f1_macro_mean','f1_weighted_median','f1_weighted_mean',
            'mcc_median','mcc_mean']
    final = final.loc[:,cols]
    print('Writing to file...')
    final.to_csv(f'{tissue}_{method}_aggregate.txt.gz', sep='\t', index=False)
    print(f'Done for {tissue}')

    return final