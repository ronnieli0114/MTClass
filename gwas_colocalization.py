# -*- coding: utf-8 -*-
"""
GWAS colocalization analysis
@author: ronnieli
"""

import pandas as pd
import numpy as np
import requests
import matplotlib.pyplot as plt
import os
from argparse import ArgumentParser
pd.options.mode.chained_assignment = None

parser = ArgumentParser(prog='MTClass classifier', description='Ensemble classifier for gene-SNP pairs')
parser.add_argument('mtclass', help='Path to MTClass results', type=str)
parser.add_argument('mtclass_metric', help='MTClass classification metric to use', type=str)
parser.add_argument('multiphen', help='Path to MultiPhen results', type=str)
parser.add_argument('manova', help='Path to MANOVA results', type=str)
parser.add_argument('out_dir', help='Path to output folder', type=str)
parser.add_argument('-g','--gwas_catalog', help='Path to downloaded GWAS Catalog', type=str)
parser.add_argument('-d','--download', action='store_true', help='Download GWAS Catalog first')
parser.add_argument('-v','--verbose', action='store_true', help='Prints verbose output')
args = parser.parse_args()

def download_gwas_catalog(data_dir):
    
    ''' Downloads the GWAS Catalog from NHGRI-EBI '''
   
    URL = 'https://www.ebi.ac.uk/gwas/api/search/downloads/full'
    response = requests.get(URL)
    open(os.path.join(data_dir, "gwas_catalog_associations.tsv"), "wb").write(response.content)
    print('Finished downloading GWAS Catalog')

def process_gwas_catalog(gwas_table):

    ''' Processes the GWAS Catalog from NHGRI-EBI '''

    gwas_table = gwas_table[~gwas_table['CHR_POS'].isna()] # remove missing positions
    
    # Get rid of unusually formatted positions
    def check_positions(position): 
        pos_list = position.replace(';',' ').split(' ')
        if len(pos_list) > 1:
            return np.nan
        else:
            return pos_list[0]
        
    gwas_table['position'] = gwas_table['CHR_POS'].apply(check_positions)
    gwas_table = gwas_table[~gwas_table['position'].isna()]
    
    # Define appropriate columns and make GenomicRanges object
    gwas_table['seqnames'] = gwas_table['CHR_ID'].apply(lambda x: "chr"+str(x))
    gwas_table['starts'] = gwas_table['position'].astype(int)
    gwas_table['ends'] = gwas_table['position'].astype(int)
    print('Loaded GWAS Catalog')
    gwas_table = gwas_table.loc[:,['seqnames','starts','ends']]
    return gwas_table

def gwas_hits(results, gwas_table, metric, n_top, bp=10_000):
    
    """ 
    Counts number of GWAS hits in neighborhood window of selected
    variants (both upstream and downstream).
    
    Input: 
        - results = results dataframe (MTClass, MultiPhen, or MANOVA)
        - gwas_table = GWAS Catalog dataframe loaded by load_gwas_table
        - metric = metric to sort by ('pval' if MultiPhen or MANOVA)
        - n_top = number of top SNPs to include
        - bp = +/- base pairs to consider in interval (default 10,000)
    Returns:
        - Total number of GWAS hits nearby
        
    """
    
    if 'pval' in metric: # sort by p-values in ascending order
        results = results.sort_values(by=metric, ascending=True)
        thres = results.iloc[n_top,:][metric]
        thres_list = []
        for metr in results[metric]:
            if metr < thres:
                thres_list.append('include')
            elif metr == thres:
                thres_list.append('dw')
            else:
                thres_list.append('exclude')
                
    else: # assume it's classification metric
        results = results.sort_values(by=metric, ascending=False)
        thres = results.iloc[n_top,:][metric] # value to break ties
        thres_list = [] 
        for metr in results[metric]:
            if metr > thres:
                thres_list.append('include')
            elif metr == thres:
                thres_list.append('dw')
            else:
                thres_list.append('exclude')
    
    results['thres'] = thres_list

    results['seqnames'] = results['variant'].apply(lambda x: "chr"+ str(x.split('_')[0][3:]))
    results['starts'] = results['variant'].apply(lambda x: x.split('_')[1]).astype(int)
    results['ends'] = results['variant'].apply(lambda x: x.split('_')[1]).astype(int)
    
    # Search upstream and downstream 10kb
    results['starts'] = results['starts'].subtract(bp)
    results['ends'] = results['ends'].add(bp)
    
    df_inc = results[results['thres']=='include']
    df_dw = results[results['thres']=='dw']
    
    # count number of GWAS hits to include
    df_inc = df_inc.loc[:,['seqnames','starts','ends']]
    count_inc = 0
    
    for chrom in set(df_inc['seqnames']):
        df_inc_chrom = df_inc[df_inc['seqnames']==chrom]
        gwas_chrom = gwas_table[gwas_table['seqnames']==chrom]
        for _, row in df_inc_chrom.iterrows():
            interval = range(row['starts'], row['ends']+1)
            for gwas_pos in gwas_chrom['starts']:
                if gwas_pos in interval:
                    count_inc += 1
    
    # count number of GWAS hits to down-weight
    df_dw = df_dw.loc[:,['seqnames','starts','ends']]
    count_dw = 0
    
    for chrom in set(df_dw['seqnames']):
        df_dw_chrom = df_dw[df_dw['seqnames']==chrom]
        gwas_chrom = gwas_table[gwas_table['seqnames']==chrom]
        for _, row in df_dw_chrom.iterrows():
            interval = range(row['starts'], row['ends']+1)
            for gwas_pos in gwas_chrom['starts']:
                if gwas_pos in interval:
                    count_dw += 1
    
    # Down-weight the counts with threshold values
    n_have = df_dw.shape[0]
    n_needed = n_top - df_inc.shape[0]
    ratio = n_needed/n_have
    count_dw = count_dw * ratio
    
    # Return total GWAS hits
    total_count = np.around((count_inc + count_dw), 2)
    return total_count


print('Reading results files...')
mtclass = pd.read_csv(args.mtclass, sep='\t', header=0)
multiphen = pd.read_csv(args.multiphen, sep='\t', header=0)
manova = pd.read_csv(args.manova, sep='\t', header=0)

mtclass['pair'] = mtclass['gene'] + ';' + mtclass['variant']
multiphen['pair'] = multiphen['gene'] + ';' + multiphen['variant']
manova['pair'] = manova['gene'] + ';' + manova['variant']

mtclass.drop_duplicates('pair', inplace=True)
multiphen.drop_duplicates('pair', inplace=True)
manova.drop_duplicates('pair', inplace=True)

all_pairs = [set(data['pair']) for data in [mtclass, multiphen, manova]]
common_pairs = set.intersection(*all_pairs)

mtclass = mtclass[mtclass['pair'].isin(common_pairs)]
multiphen = multiphen[multiphen['pair'].isin(common_pairs)]
manova = manova[manova['pair'].isin(common_pairs)]

print('Loaded MTClass, MultiPhen, MANOVA')
print(f'There are {len(common_pairs)} total gene-SNP pairs')

if args.download:
    print('First downloading GWAS Catalog...')
    download_gwas_catalog(args.out_dir)
    gwas_table = pd.read_csv(os.path.join(args.out_dir, 'gwas_catalog_associations.tsv'), sep='\t', low_memory=False)
    gwas_table = process_gwas_catalog(gwas_table)
else:
    print('Reading GWAS Catalog...')
    gwas_table = pd.read_csv(args.gwas_catalog, sep='\t', low_memory=False)
    gwas_table = process_gwas_catalog(gwas_table)
        
if args.mtclass_metric not in list(mtclass.columns):
    raise KeyError('Classification metric not found in column names')
    
results_dict = {'Top variants':[], 'MTClass':[], 'MultiPhen':[], 'MANOVA':[]}
n_top_list = [100, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000]

for n_top in n_top_list:
    
    mtclass_count = gwas_hits(mtclass, gwas_table, args.mtclass_metric, n_top, bp=10000)
    mphen_count = gwas_hits(multiphen, gwas_table, 'pval', n_top, bp=10000)
    manova_count = gwas_hits(manova, gwas_table, 'pval', n_top, bp=10000)
    
    results_dict['Top variants'].append(n_top)
    results_dict['MTClass'].append(mtclass_count)
    results_dict['MultiPhen'].append(mphen_count)
    results_dict['MANOVA'].append(manova_count)

    if args.verbose:       
        print(pd.DataFrame(results_dict).tail(n=1))

results = pd.DataFrame(results_dict)
results.to_csv(os.path.join(args.out_dir, f'gwas_colocalization_{args.mtclass_metric}.txt'), sep='\t', index=False)

# Plot results
fig, ax = plt.subplots(figsize=(6,6))
ax.plot(results['Top variants'], results['MTClass'], color='blue', linestyle='solid', label='MTClass')
ax.plot(results['Top variants'], results['MultiPhen'], color='orange', linestyle='solid', label='MultiPhen')
ax.plot(results['Top variants'], results['MANOVA'], color='green', linestyle='solid', label='MANOVA')

plt.title(f'GWAS hits by method\n{args.mtclass_metric}', fontsize=20)
plt.ylabel('GWAS hits', fontsize=20)
plt.xlabel('Top variants', fontsize=20)
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig(os.path.join(args.out_dir, f'GWAS_colocalization_{args.mtclass_metric}'), dpi=400)
