# -*- coding: utf-8 -*-
"""
Impute counts for a given experiment by gene in parallel

Created on Wed Jun 29 17:16:30 2022
@author: ronnieli
"""

import pandas as pd
import os
import sys
import time
import multiprocessing as mp
from itertools import repeat
import pyarrow.feather as feather
pd.options.mode.chained_assignment = None
from statsmodels.imputation import mice
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
from sklearn.ensemble import RandomForestRegressor

computer = 'laptop'
experiment = 'brain'

if experiment == 'brain':
    suffix = '_brain'
elif experiment == 'all':
    suffix = '_all'

if computer == 'pc':
    data_dir = r"D:/ronnieli/My Drive/lab_data/gtex_eqtl"
    results_dir = rf"D:/ronnieli/My Drive/lab_data/gtex_eqtl/multi_tissue{suffix}"
    sys.path.append(r"C:\Users\ronnieli\OneDrive\Documents\Lab\gtex_eqtl_project\scripts")
elif computer == 'cluster':
    data_dir = r'/projects/compbio/users/yli459/gtex_eqtl'
    results_dir = rf'/projects/compbio/users/yli459/gtex_eqtl/multi_tissue{suffix}'
    sys.path.append(r'/projects/compbio/users/yli459/scripts')
elif computer == 'laptop':
    data_dir = r"C:/Users/ronnieli/My Drive/lab_data/gtex_eqtl"
    results_dir = rf"C:/Users/ronnieli/My Drive/lab_data/gtex_eqtl/multi_tissue{suffix}"
    
counts = feather.read_feather(os.path.join(data_dir, f'multi_tissue{suffix}', 'counts.ftr'))

genes = sorted(list(set(counts['gene'])))
num_genes = len(genes)
print(f'There are {num_genes} total genes to impute')

#%% Define function to impute counts for a gene

def impute_counts_iter(gene, counts):
    
    t0 = time.time()
    print(f'Now imputing counts for {gene}')
    counts_gene = counts[counts['gene']==gene]
    donor_list = counts_gene['donor']
    data = counts_gene.drop(['gene','donor'], axis=1)
    cols = list(map(lambda x: x.replace(' - ','_'), data.columns))
    data.columns = cols
    
    imputer = IterativeImputer(estimator = RandomForestRegressor(), max_iter=10, random_state=2022)
    imp_data = imputer.fit_transform(data)
    imp_data = pd.DataFrame(imp_data, columns=cols)
    imp_data['gene'] = gene
    imp_data['donor'] = list(donor_list)
    
    new_cols = ['gene','donor'] + list(imp_data.iloc[:,:-2].columns)
    imp_data = imp_data.loc[:,new_cols]
    
    t1 = time.time()
    print(f'Took {(t1-t0)/60:.1f} min to impute counts for {gene}')
    return imp_data.copy()

def impute_counts_mice(gene, counts):
    
    t0 = time.time()
    
    counts_gene = counts[counts['gene']==gene]
    donor_list = list(counts_gene['donor'])
    data = counts_gene.drop(['gene','donor'], axis=1)
    cols = list(map(lambda x: x.replace('-',''), data.columns))
    data.columns = cols

    imp_sets = []
    imp_data = mice.MICEData(data)
    for _ in range(10):
        imp = imp_data.next_sample()
        imp_sets.append(imp.copy())
    avg_df = sum(imp_sets) / len(imp_sets)
    avg_df['gene'] = gene
    avg_df['donor'] = donor_list
    new_cols = ['gene','donor'] + list(avg_df.iloc[:,:-2].columns)
    avg_df = avg_df.loc[:new_cols]
    
    t1 = time.time()
    print(f'Took {(t1-t0)/60:.1f} min to impute counts for {gene}')
    print(avg_df.head())
    return avg_df

#%% Impute the counts
num_cpus = mp.cpu_count()-4

if __name__ == '__main__':
 
    with mp.Pool(processes=num_cpus) as pool:

        mp_input = zip(genes, repeat(counts))
        result = pool.starmap(impute_counts_mice, mp_input)
        final_results = pd.concat(result, axis=0)
        feather.write_feather(final_results, os.path.join(results_dir, 'MICE_imputed_counts.ftr'))
        print(f'Saved imputed counts to {results_dir}')