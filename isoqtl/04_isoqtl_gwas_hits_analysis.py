# -*- coding: utf-8 -*-
"""
Calculate GWAS hits for isoQTL results

Created on Mon Mar  6 10:10:55 2023
@author: ronnieli
"""

import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
pd.options.mode.chained_assignment = None
sys.path.append(r"G:/My Drive/Lab/lab_projects/gtex_eqtl/scripts")
gwas_functions = __import__('06_gwas_hits_analysis')

def gwas_hits_iterations_isoqtl():
    
    ''' Compare between 3 iterations of MTClass and MultiPhen/MANOVA '''
    
    gwas_table = gwas_functions.load_gwas_table(r'E:/lab_data')
    
    os.chdir(r'G:/My Drive/Lab/lab_projects/isoqtl/results')
    m1 = pd.read_csv('isoqtl_ensemble_1.txt.gz', sep='\t', header=0)
    m2 = pd.read_csv('isoqtl_ensemble_2.txt.gz', sep='\t', header=0)
    m3 = pd.read_csv('isoqtl_ensemble_3.txt.gz', sep='\t', header=0)
    mphen_binary = pd.read_csv('isoqtl_multiphen_binary.txt.gz', sep='\t', header=0)
    mphen_add = pd.read_csv('isoqtl_multiphen_additive.txt.gz', sep='\t', header=0)
    manova = pd.read_csv('isoqtl_manova.txt.gz', sep='\t', header=0)
    title = 'isoQTL'

    m1['pair'] = m1['gene'] + '-' + m1['variant']
    m2['pair'] = m2['gene'] + '-' + m2['variant']
    m3['pair'] = m3['gene'] + '-' + m3['variant']
    mphen_binary['pair'] = mphen_binary['gene'] + '-' + mphen_binary['variant']
    mphen_add['pair'] = mphen_add['gene'] + '-' + mphen_add['variant']
    manova['pair'] = manova['gene'] + '-' + manova['variant']
    
    all_pairs = [set(data['pair']) for data in [m1, m2, m3, mphen_binary, mphen_add, manova]]
    common_pairs = set.intersection(*all_pairs)
    
    m1 = m1[m1['pair'].isin(common_pairs)]
    m2 = m2[m2['pair'].isin(common_pairs)]
    m3 = m3[m3['pair'].isin(common_pairs)]
    mphen_binary = mphen_binary[mphen_binary['pair'].isin(common_pairs)]
    mphen_add = mphen_add[mphen_add['pair'].isin(common_pairs)]
    manova = manova[manova['pair'].isin(common_pairs)]
    
    m1.drop_duplicates('pair', inplace=True)
    m2.drop_duplicates('pair', inplace=True)
    m3.drop_duplicates('pair', inplace=True)
    mphen_binary.drop_duplicates('pair', inplace=True)
    mphen_add.drop_duplicates('pair', inplace=True) 
    manova.drop_duplicates('pair', inplace=True)
    
    num_genes = len(set(m1['gene']))
    
    print('Loaded MTClass, MultiPhen, MANOVA')
    print(f'There are {len(common_pairs)} total gene-SNP pairs in {num_genes} genes')
    print('Calculating GWAS hits')
    
    final_dict = {}
    
    # Calculate GWAS hits
    for metric, name in zip(['f1_macro','mcc'],['Macro F1','MCC']):
        
        result_dict = {'Top variants':[], 'MTClass 1':[], 'MTClass 2':[], 'MTClass 3':[], 
                       'MultiPhen (binary)':[], 'MultiPhen (additive)':[], 'MANOVA':[]}
    
        n_top_list = [100, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000]
        
        for n_top in n_top_list:
            
            m1_count = gwas_functions.gwas_hits(m1, gwas_table, metric=metric, n_top=n_top, bp=10000)
            m2_count = gwas_functions.gwas_hits(m2, gwas_table, metric=metric, n_top=n_top, bp=10000)
            m3_count = gwas_functions.gwas_hits(m3, gwas_table, metric=metric, n_top=n_top, bp=10000)
            mphen_count_add = gwas_functions.gwas_hits(mphen_add, gwas_table, metric='pval', n_top=n_top, bp=10000)
            mphen_count_bin = gwas_functions.gwas_hits(mphen_binary, gwas_table, metric='pval', n_top=n_top, bp=10000)
            manova_count = gwas_functions.gwas_hits(manova, gwas_table, metric='pval', n_top=n_top, bp=10000)
            
            result_dict['Top variants'].append(n_top)
            result_dict['MTClass 1'].append(m1_count)
            result_dict['MTClass 2'].append(m2_count)
            result_dict['MTClass 3'].append(m3_count)
            result_dict['MultiPhen (binary)'].append(mphen_count_bin)
            result_dict['MultiPhen (additive)'].append(mphen_count_add)
            result_dict['MANOVA'].append(manova_count)
            print(pd.DataFrame(result_dict).tail(n=1))
            
        results = pd.DataFrame(result_dict)
        final_dict[metric] = results
        
        # Plot results
        fig, ax = plt.subplots(figsize=(6,6))
        ax.plot(results['Top variants'], results['MTClass 1'], color='blue', linestyle='solid', label='MTClass 1')
        ax.plot(results['Top variants'], results['MTClass 2'], color='blue', linestyle='dotted', label='MTClass 2')
        ax.plot(results['Top variants'], results['MTClass 3'], color='blue', linestyle='dashed', label='MTClass 3')
        ax.plot(results['Top variants'], results['MultiPhen (binary)'], color='orange', linestyle='solid', label='MultiPhen (binary)')
        ax.plot(results['Top variants'], results['MultiPhen (additive)'], color='orange', linestyle='dotted', label='MultiPhen (additive)')
        ax.plot(results['Top variants'], results['MANOVA'], color='green', linestyle='solid', label='MANOVA')
    
        ax.set_title(f'GWAS hits by method\n{name}\n{title}', fontsize=20)
        ax.set_xlabel('Top variants', fontsize=20)
        ax.set_ylabel('GWAS hits', fontsize=20)
        ax.legend(fontsize=12)
        plt.tight_layout()
        plt.savefig(rf'C:/Users/ronni/OneDrive/Desktop/isoqtl_{metric}.png', dpi=400)
        plt.show()
    
    return final_dict

def gwas_hits_aggregate_isoqtl():
    
    ''' Compare GWAS hits between aggregated version of MTClass, MultiPhen, MANOVA '''
    
    gwas_table = gwas_functions.load_gwas_table(r'E:/lab_data')
    
    os.chdir(r'G:/My Drive/Lab/lab_projects/isoqtl/results')
    mtclass = pd.read_csv('isoqtl_ensemble_aggregate.txt.gz', sep='\t', header=0)
    mphen_binary = pd.read_csv('isoqtl_multiphen_binary.txt.gz', sep='\t', header=0)
    mphen_add = pd.read_csv('isoqtl_multiphen_additive.txt.gz', sep='\t', header=0)
    manova = pd.read_csv('isoqtl_manova.txt.gz', sep='\t', header=0)
    title = 'isoQTL'
    
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
    
    num_genes = len(set(mtclass['gene']))
    
    print('Loaded MTClass, MultiPhen, MANOVA')
    print(f'There are {len(common_pairs)} total gene-SNP pairs in {num_genes} genes')
    print('Calculating GWAS hits')
    
    final_dict = {}
    
    for metric, name in zip(['f1_macro_median','mcc_median'],['Median Macro F1','Median MCC']):
        
        results_dict = {'Top variants':[], 'MTClass':[], 'MultiPhen (binary)':[],
                        'MultiPhen (additive)':[], 'MANOVA':[]}
        n_top_list = [100, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000]
        
        for n_top in n_top_list:
            
            mtclass_count = gwas_functions.gwas_hits(mtclass, gwas_table, metric, n_top, bp=10000)
            mphen_binary_count = gwas_functions.gwas_hits(mphen_binary, gwas_table, 'pval', n_top, bp=10000)
            mphen_add_count = gwas_functions.gwas_hits(mphen_add, gwas_table, 'pval', n_top, bp=10000)
            manova_count = gwas_functions.gwas_hits(manova, gwas_table, 'pval', n_top, bp=10000)
            
            results_dict['Top variants'].append(n_top)
            results_dict['MTClass'].append(mtclass_count)
            results_dict['MultiPhen (binary)'].append(mphen_binary_count)
            results_dict['MultiPhen (additive)'].append(mphen_add_count)
            results_dict['MANOVA'].append(manova_count)         
            print(pd.DataFrame(results_dict).tail(n=1))
        
        results = pd.DataFrame(results_dict)
        final_dict[metric] = results
        
        # Plot results
        fig, ax = plt.subplots(figsize=(6,6))
        ax.plot(results['Top variants'], results['MTClass'], color='blue', linestyle='solid', label='MTClass')
        ax.plot(results['Top variants'], results['MultiPhen (binary)'], color='orange', linestyle='solid', label='MultiPhen (binary)')
        ax.plot(results['Top variants'], results['MultiPhen (additive)'], color='orange', linestyle='dotted', label='MultiPhen (additive)')
        ax.plot(results['Top variants'], results['MANOVA'], color='green', linestyle='solid', label='MANOVA')
        
        plt.title(f'GWAS hits by method\n{name}\n{title}', fontsize=20)
        plt.ylabel('GWAS hits', fontsize=20)
        plt.xlabel('Top variants', fontsize=20)
        plt.legend(fontsize=12)
        plt.tight_layout()
        plt.savefig(rf'C:/Users/ronni/OneDrive/Desktop/isoqtl_{metric}.png', dpi=400)
        plt.show()
    
    return final_dict