# -*- coding: utf-8 -*-
"""
Anecdotal examples for isoQTL study

Created on Mon Jun 12 08:57:44 2023
@author: ronnieli
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import os
import sys
pd.options.mode.chained_assignment = None
sys.path.append(r"G:/My Drive/Lab/lab_projects/isoqtl/scripts")
from isoqtl_expression_plots import *

def other_methods_NS():
    
    '''
    Returns gene-SNP pairs for which MultiPhen and MANOVA
    both output p-values that are not significant
    '''

    os.chdir(r'G:/My Drive/Lab/lab_projects/isoqtl/results/')
    
    mtclass = pd.read_csv('isoqtl_ensemble_aggregate.txt.gz', sep='\t', header=0)
    mphen_binary = pd.read_csv('isoqtl_multiphen_binary.txt.gz', sep='\t', header=0)
    mphen_add = pd.read_csv('isoqtl_multiphen_additive.txt.gz', sep='\t', header=0)
    manova = pd.read_csv('isoqtl_manova.txt.gz', sep='\t', header=0)

    mtclass['pair'] = mtclass['gene'] + '-' + mtclass['variant']
    mphen_binary['pair'] = mphen_binary['gene'] + '-' + mphen_binary['variant']
    mphen_add['pair'] = mphen_add['gene'] + '-' + mphen_add['variant']
    manova['pair'] = manova['gene'] + '-' + manova['variant']
    
    mtclass.drop_duplicates('pair', inplace=True)
    mphen_binary.drop_duplicates('pair', inplace=True)
    mphen_add.drop_duplicates('pair', inplace=True) 
    manova.drop_duplicates('pair', inplace=True)
    
    all_pairs = [set(data['pair']) for data in [mtclass, mphen_binary, mphen_add, manova]]
    common_pairs = set.intersection(*all_pairs)
    
    mtclass = mtclass[mtclass['pair'].isin(common_pairs)]
    mphen_binary = mphen_binary[mphen_binary['pair'].isin(common_pairs)]
    mphen_add = mphen_add[mphen_add['pair'].isin(common_pairs)]
    manova = manova[manova['pair'].isin(common_pairs)]
    
    print('Loaded MTClass, MultiPhen, MANOVA')
    print(f'There are {len(common_pairs)} total gene-SNP pairs')
    
    metric = 'f1_macro_median'
    threshold = 0.80
    pval_thres = 1e-4
    
    mtclass_top = mtclass[mtclass[metric] >= threshold]
    
    df1 = mtclass_top[['gene','variant','pair',metric]].merge(mphen_binary[['pair','pval']], on='pair')
    df2 = df1.merge(mphen_add[['pair','pval']], on='pair', suffixes=('_binary','_additive'))
    df3 = df2.merge(manova[['pair','pval']], on='pair')
    
    mtclass_only = df3[(df3['pval_binary'] >= pval_thres) & (df3['pval_additive'] >= pval_thres) & (df3['pval'] >= pval_thres)]
    if len(mtclass_only) == 0: # empty dataframe
        mtclass_only = df3[(df3['pval_binary'].isna()) & (df3['pval_additive'].isna()) & (df3['pval'].isna())]
    
    return mtclass_only