# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 20:42:53 2022

@author: ronnieli
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pyarrow.feather as feather
import numpy as np
import time
import os

experiment = 'brain'
data_dir = r"E:/lab_data/gtex_eqtl/"

t0 = time.time()

# load GTEx data
col_metadata = pd.read_pickle(os.path.join(data_dir, 'col_metadata.pkl'))
row_metadata = pd.read_pickle(os.path.join(data_dir, 'row_metadata.pkl'))
tpm_counts = feather.read_feather(os.path.join(data_dir, 'gene_tpm_counts.ftr'))
top_counts = feather.read_feather(os.path.join(data_dir, 'top_tpm_counts.ftr'))
attributes = pd.read_csv(os.path.join(data_dir, 'v8_sample_attributes.txt'), sep='\t', usecols=['SAMPID','SMTSD'])

t1 = time.time()
print(f'Took {t1-t0:.1f} seconds to import all data')

#%% Sample availability and determining tissues vs. donors

### Make dataframe of sample/tissue availability
avail_samples = attributes[attributes['SAMPID'].isin(col_metadata.index)]
avail_samples['donor'] = ['-'.join(x.split("-")[:2]) for x in avail_samples['SAMPID']]
avail_samples['count'] = 1
avail_samples.columns = ['full_id', 'tissue', 'donor', 'count']

samples_pivot = avail_samples.pivot(index='tissue', columns='donor', values='full_id')
samples_pivot['num_donors'] = samples_pivot.count(axis=1)
samples_pivot.sort_values(by='num_donors', ascending=False, na_position='last', inplace=True)

# determine number of donors we'd have for given number of top tissues to include
for n_tissues in range(1, len(samples_pivot.index)):
    subset = samples_pivot.iloc[:n_tissues, :-1].dropna(axis=1, how='any')
    num_samples = subset.shape[1]
    if num_samples > 0:
        print(f"{n_tissues} tissues:\t {num_samples} donors")

if experiment == 'brain':
    samples = avail_samples[avail_samples['tissue'].str.startswith('Brain -')]
    full_ids = samples['full_id']
elif experiment == 'all':
    samples = samples_pivot[samples_pivot['num_donors'] >= 100]
    samples = samples.drop('num_donors', axis=1)
    full_ids = samples.to_numpy().flatten()
    full_ids = pd.Series(full_ids).dropna()

#%% Filter out samples based on full_ids
filtered_counts = top_counts.set_index('Description')
filtered_counts.drop(['rid','variance'], axis=1, inplace=True)

exp_counts = filtered_counts.filter(full_ids).reset_index()
exp_counts = exp_counts.set_index('Description')
exp_counts = exp_counts.T.merge(avail_samples[['full_id','tissue','donor']], left_index=True, right_on='full_id')

df = exp_counts.set_index(['tissue','donor'])
df.drop('full_id', axis=1, inplace=True)
counts = df.T.stack()
counts.reset_index(inplace=True)
counts.rename(columns={'level_0':'gene'}, inplace=True)
counts.columns.rename('', inplace=True)

#%% Visualize missing data for brain tissue

brain_pivot = samples_pivot[samples_pivot.index.str.startswith('Brain')].copy()
brain_pivot = brain_pivot.drop('num_donors', axis=1)

brain_pivot = brain_pivot.isnull().astype(int)
brain_pivot.replace({0:np.nan}, inplace=True)
# brain_pivot = brain_pivot.sort_index(na_position='last')

fig, ax = plt.subplots(figsize=(12,10))
sns.heatmap(brain_pivot.isna(), cmap='YlGnBu', cbar_kws={'label':'Missing Data'})
plt.tight_layout()
plt.savefig(r"C:/Users/ronnieli/Desktop/brain.png", dpi=400)
plt.show()
