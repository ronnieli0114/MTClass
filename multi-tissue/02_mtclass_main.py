# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 14:25:30 2022

@author: ronnieli
"""
import pandas as pd
import multiprocessing as mp
import pyarrow.feather as feather
import os, sys
import time
import random
from itertools import repeat

experiment = '9_tissue'
model_name = 'voting_soft'
permuted_genotypes = False # whether to load permuted genotypes
genotype_type = 'binary' # binary or additive
chr_range = range(1,23) # usually DON'T change this!
idx_range = [1]

computer = 'pc' # determines data directory
cpus = 8
processing_type = 'serial' # whether to use serial or parallel processing

if experiment == '9_tissue':
    suffix = '9_tissue'
elif experiment == 'brain':
    suffix = 'brain_tissue'

if computer == 'pc':
    data_dir = r'E:/lab_data/gtex_eqtl'
    results_dir = rf'G:/My Drive/Lab/lab_projects/gtex_eqtl/results/{suffix}'
    sys.path.append(r"G:/My Drive/Lab/lab_projects/gtex_eqtl/scripts")
elif computer == 'cluster':
    data_dir = "/projects/compbio/users/yli459/gtex_eqtl"
    results_dir = rf'/projects/compbio/users/yli459/gtex_eqtl/{suffix}'
    sys.path.append(r"/projects/compbio/users/yli459/scripts")
    
load_data = __import__('02_load_data')
train_model = __import__('02_train_model')

# Run model
print('Using normal genotypes...')
for idx in idx_range: # loop through model n times
    seed = random.randint(1,1001)
    print(f'Randomly generated seed is {seed}')
    for c in chr_range:
        print(f'Loading data for iteration #{idx}/10, chr{c}')
        counts, gt = load_data.load_data_chrom(c, experiment, computer, genotype_type)
        gene_list = sorted(list(set(gt['gene_name'])))
        num_genes = len(gene_list)
        
        if processing_type == 'serial':
            results = []
            for gene in gene_list:
                gene_result = train_model.train_model(gene, model_name, counts, gt, seed)
                results.append(gene_result)
            chrom_result = pd.concat(results)
            feather.write_feather(chrom_result, os.path.join(results_dir, f'chr{c}_mtclass_{model_name}_iter{idx}.ftr'))
            
        elif processing_type == 'parallel':
            if __name__ == '__main__':
                results = []
                with mp.Pool(processes=cpus) as pool:
                    inputs = zip(gene_list, repeat(model_name), repeat(counts), repeat(gt), repeat(seed))
                    gene_result = pool.starmap(train_model.train_model, inputs)
                    results.append(pd.concat(gene_result))
                chrom_result = pd.concat(results)
                print(chrom_result.head())
                if genotype_type == 'additive':
                    feather.write_feather(chrom_result, os.path.join(results_dir, f'chr{c}_mtclass_{model_name}_iter{idx}_additive.ftr'))
                elif genotype_type == 'binary':
                    feather.write_feather(chrom_result, os.path.join(results_dir, f'chr{c}_mtclass_{model_name}_iter{idx}.ftr'))
print('=== DONE ===')

