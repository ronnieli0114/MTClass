# -*- coding: utf-8 -*-
"""
GWAS hits

@author: ronnieli
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pyarrow.feather as feather
import scipy.stats as stats
import os
import sys
pd.options.mode.chained_assignment = None
sys.path.append(r'G:/My Drive/Lab/lab_projects/gtex_eqtl/scripts')
from gwas_functions import gwas_hits, load_gwas_table, gwas_hits_single_tissue_eqtl, gwas_hits_single_tissue_sqtl

#%% Compare MTClass to MultiPhen and MANOVA

def gwas_hits_exoneqtl():
    
    gwas_table = load_gwas_table(data_dir=r'E:/lab_data')
    os.chdir(r'G:/My Drive/Lab/lab_projects/gtex_exoneqtl/results')
    tissues = pd.read_csv(r'E:/lab_data/gtex_exoneqtl/brain_tissues.txt', sep='\t', header=None)[1].tolist()
    
    for tissue in tissues:
        
        mtclass = pd.read_csv(f'{tissue}_ensemble_aggregate.txt.gz', sep='\t', header=0)
        mphen_binary = pd.read_csv(f'{tissue}_multiphen_binary.txt.gz', sep='\t', header=0)
        manova = pd.read_csv(f'{tissue}_manova.txt.gz', sep='\t', header=0)
        
        mtclass['pair'] = mtclass['gene'] + '-' + mtclass['variant']
        mphen_binary['pair'] = mphen_binary['gene'] + '-' + mphen_binary['variant']
        manova['pair'] = manova['gene'] + '-' + manova['variant']
        
        all_pairs = [set(data['pair']) for data in [mtclass, mphen_binary, manova]]
        common_pairs = set.intersection(*all_pairs)
        
        mtclass = mtclass[mtclass['pair'].isin(common_pairs)]
        mphen_binary = mphen_binary[mphen_binary['pair'].isin(common_pairs)]
        manova = manova[manova['pair'].isin(common_pairs)]
        
        print(f'Loaded data for {tissue}. There are {len(common_pairs)} gene-SNP pairs')
        
        for metric, name in zip(['f1_macro_median','mcc_median'],['Median Macro F1','Median MCC']):
            
            result_dict = {'Top variants':[], 'MTClass':[], 'MultiPhen (binary)':[], 'MANOVA':[]}
        
            n_top_list = [100, 500, 1000, 1500, 2000, 2500]
            
            for n_top in n_top_list:
                
                mtclass_count = gwas_hits(mtclass, gwas_table, metric=metric, n_top=n_top, bp=10000)
                mphen_count_bin = gwas_hits(mphen_binary, gwas_table, metric='pval', n_top=n_top, bp=10000)
                manova_count = gwas_hits(manova, gwas_table, metric='pval', n_top=n_top, bp=10000)
                
                result_dict['Top variants'].append(n_top)
                result_dict['MTClass'].append(mtclass_count)
                result_dict['MultiPhen (binary)'].append(mphen_count_bin)
                result_dict['MANOVA'].append(manova_count)
                print(pd.DataFrame(result_dict))
                
            results = pd.DataFrame(result_dict)
            
            # Plot results
            fig, ax = plt.subplots(figsize=(6,6))
            ax.plot(results['Top variants'], results['MTClass'], color='blue', linestyle='solid', label='MTClass')
            ax.plot(results['Top variants'], results['MultiPhen (binary)'], color='orange', linestyle='solid', label='MultiPhen (binary)')
            ax.plot(results['Top variants'], results['MANOVA'], color='green', linestyle='solid', label='MANOVA')
        
            ax.set_title(f'GWAS hits by method\n{name}\n{tissue}', fontsize=20)
            ax.set_xlabel('Top variants', fontsize=20)
            ax.set_ylabel('GWAS hits', fontsize=20)
            ax.legend(fontsize=12)
            plt.tight_layout()
            save_path = r'C:/Users/ronni/OneDrive/Desktop/gwas_hits/'
            if not os.path.exists(save_path):
                os.mkdir(save_path)
            plt.savefig(os.path.join(save_path, f'{tissue}_{metric}.png', dpi=400)
            plt.show()

#%% Compare multi-tissue approach to single-tissue approaches (eQTL and sQTL)

method = 'ensemble' # ensemble or RF

os.chdir(r'C:/Users/ronnieli/OneDrive/Documents/Lab/gtex_exonqtl/results')
tissues = sorted(list(set([tissue.split('_mtclass')[0] for tissue in os.listdir() if 'mtclass' in tissue])))

for tissue in ['Brain_Spinal_cord_cervical_c-1']:
    
    print(f'Now calculating GWAS hits for {tissue}')
    mtclass = feather.read_feather(f'{tissue}_mtclass_{method}.ftr')

    for metric in ['f1_macro']:
        
        results_dict = {'Top variants':[], 'MTClass':[], 'Single-tissue eQTL':[], 'Single-tissue sQTL':[]}

        if metric == 'f1_macro':
            n_top = mtclass[mtclass[metric] >= 0.90].shape[0]
        elif metric == 'mcc':
            n_top = mtclass[mtclass[metric] >= 0.80].shape[0]
        print(f'Using {n_top} top variants')
        
        mtclass_count = gwas_hits(mtclass, gwas_table, 'RF', metric, n_top, get_traits=False)
        eqtl_count = gwas_hits_single_tissue_eqtl(tissue, gwas_table, n_top, 'pc', get_traits=False)
        sqtl_count, sqtl_traits, sqtl_studies, sqtl_snps = gwas_hits_single_tissue_sqtl(tissue, gwas_table, n_top, 'pc', get_traits=True)
        
        results_dict['Top variants'].append(n_top)
        results_dict['MTClass'].append(mtclass_count)
        results_dict['Single-tissue eQTL'].append(eqtl_count)
        results_dict['Single-tissue sQTL'].append(sqtl_count)
        print(pd.DataFrame(results_dict))
        
        results = pd.DataFrame(results_dict)
        title = tissue.split('Brain_')[1]
        
        fig, ax = plt.subplots()
        results.set_index('Top variants').plot.bar(rot=0, figsize=(8,6), ax=ax, fontsize=16)
        for container in ax.containers:
            ax.bar_label(container)
        leg = ax.legend(loc='best', ncol=1, fontsize=12, facecolor='#FFFFFF')
        ax.set_xlabel('Top variants', fontsize=16)
        plt.title(f'GWAS hits by method, MTClass sorted by {metric}\n{tissue}', fontsize=16)
        plt.ylabel('GWAS hits', fontsize=16)
        plt.tight_layout()
        plt.savefig(rf'C:/Users/ronnieli/Desktop/{metric}_{title}.png', dpi=400)
        plt.show()

#%%
os.chdir(r"C:\Users\ronnieli\OneDrive\Documents\Lab\gtex_exonqtl\results")

tissue = 'Brain_Hippocampus'
iter_num = 'aggregate123'
title = tissue.split('Brain_')[1]

mtclass_ensemble = feather.read_feather(f'{tissue}_mtclass_RF.ftr')
mphen = feather.read_feather(f'{tissue}_multiphen.ftr')

gene_list = pd.read_csv(rf'D:/ronnieli/My Drive/lab_data/gtex_exonqtl/by_tissue/{tissue}_gene_list.txt', header=2)['gene_name'].tolist()

names_list = ['MTClass', 'MultiPhen']
data_list = [mtclass_ensemble, mphen]

# Find common gene-SNP pairs
pair_list = []
for data in data_list:
    data = data[data['gene'].isin(gene_list)]
    data['pair'] = data['gene'] + '-' + data['variant']
    pair_list.append(set(data['pair']))
common_pairs = set.intersection(*pair_list)

# Get data ready for analysis
data_dict = {}
for name, data in zip(names_list, data_list):
    data['pair'] = data['gene'] + '-' + data['variant']
    data = data[data['pair'].isin(common_pairs)]
    data.drop_duplicates('pair', inplace=True)
    data_dict[name] = data

print(f'There are {len(set(data_dict["MTClass"]["gene"]))} genes in {tissue}')
print(f'There are {len(set(data_dict["MTClass"]["pair"]))} gene-SNP pairs in {tissue}')

metric_list = ['f1_macro_median','mcc_median','f1_macro_mean','mcc_mean','f1_macro_max','mcc_max']
# metric_list = ['f1_macro','mcc']

for metric in metric_list:
    
    result_dict = {'Top variants':[], 'MTClass':[], 'MultiPhen':[]}
    
    for n_top in range(100, 5101, 500):
        
        result_dict['Top variants'].append(n_top)
    
        for method, data in data_dict.items():   

            if method == 'MultiPhen':
                method_count, _ = gwas_hits(data, gwas_table, None, 'pval', n_top)
            elif method == 'MTClass':
                method_count, _ = gwas_hits(data, gwas_table, 'voting_soft', metric, n_top)
        
            result_dict[method].append(method_count)
            print(f'{metric}, {n_top} top variants \t {method} \t {method_count}')
        
    # result_dict['Top variants'] = sorted(list(set(result_dict['Top variants'])))
    results = pd.DataFrame(result_dict)
    results.plot(x='Top variants', y=list(results.iloc[:,1:].columns), kind='line')
    plt.title(f'GWAS hits by method, MTClass sorted by {metric}\n{tissue}')
    plt.ylabel('GWAS hits')
    plt.tight_layout()
    plt.savefig(rf'C:/Users/ronnieli/Desktop/{metric}_{title}.png', dpi=400)
    plt.show()
    
#%%

