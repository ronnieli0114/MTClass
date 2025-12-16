# -*- coding: utf-8 -*-
"""
GWAS hits
@author: ronnieli
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
pd.options.mode.chained_assignment = None

def load_gwas_table(data_dir, genome='hg38'):
    
    '''
    Loads and processes the GWAS Catalog
    Inputs:
        - data_dir = data directory
    '''
    gwas_table = pd.read_table(os.path.join(data_dir, 'gwas_catalog_associations_processed.tsv.gz'), header=0, low_memory=False)
    
    # Define appropriate columns and make GenomicRanges object
    gwas_table['chrom'] = gwas_table['CHR_ID'].apply(lambda x: "chr"+str(x))
    gwas_table['start'] = gwas_table[f'pos_{genome}'].astype(int)
    gwas_table['end'] = gwas_table[f'pos_{genome}'].astype(int)
    print('loaded GWAS table')
    return gwas_table.loc[:,['chrom','start','end']]

def gwas_hits(results, gwas_table, metric, n_top, bp=10_000, genome='hg38'):
    
    """ 
    Counts number of GWAS hits in neighborhood window of selected
    variants (both upstream and downstream).
    
    Input: 
        - results = results dataframe (MTClass, MultiPhen, or MANOVA)
        - gwas_table = GWAS Catalog dataframe
        - metric = metric to sort by ('pval' if MultiPhen or MANOVA)
        - n_top = number of top SNPs to include
        - bp = +/- base pairs to consider in interval (default 10,000)
        - genome version that variants are formatted in: default hg38
    Returns:
        - Total number of GWAS hits nearby (adjusted)
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

    if genome=='hg38':
        results['seqnames'] = results['variant'].apply(lambda x: x.split('_')[0])
    elif genome=='hg19':
        results['seqnames'] = results['variant'].apply(lambda x: 'chr'+str(x.split('_')[0]))
    else:
        raise ValueError('genome is not hg38 or hg19. please specify.')
    results['starts'] = results['variant'].apply(lambda x: x.split('_')[1]).astype(int)
    results['ends'] = results['variant'].apply(lambda x: x.split('_')[1]).astype(int)
    
    # Search upstream and downstream 10kb
    results['starts'] = results['starts'].subtract(bp)
    results['ends'] = results['ends'].add(bp)
    
    df_inc = results[results['thres']=='include']
    df_dw = results[results['thres']=='dw']
    
    # use bioframe to count GWAS hits with GWAS table
    bf_inc = df_inc.loc[:,['seqnames','starts','ends']].rename(columns = {'seqnames':'chrom', 'starts':'start', 'ends':'end'})
    bf_dw = df_dw.loc[:,['seqnames','starts','ends']].rename(columns = {'seqnames':'chrom', 'starts':'start', 'ends':'end'})

    count_inc = bf.count_overlaps(bf_inc, gwas_table)['count'].sum()
    count_dw = bf.count_overlaps(bf_dw, gwas_table)['count'].sum()
    
    # Down-weight the counts with threshold values
    n_have = df_dw.shape[0]
    ratio = n_have / n_top
    count_dw = count_dw * ratio
    
    # Return total GWAS hits
    total_count = np.around((count_inc + count_dw), 2)
    return total_count
