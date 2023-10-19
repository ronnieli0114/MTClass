import pandas as pd
import numpy as np
import sys
import os
import time
import pyarrow.feather as feather

# Define function to load data for ONE chromosome
def load_data_chrom(chrom, experiment, computer, genotype_type='binary'):
    
    t0 = time.time()
    
    if experiment == '9_tissue':
        exp_suffix = '9_tissue'
    elif experiment == 'brain':
        exp_suffix = 'brain_tissue'
    
    if computer == 'pc':
        data_dir = r'E:/lab_data/gtex_eqtl'
        sys.path.append(r"G:/My Drive/Lab/lab_projects/gtex_eqtl/scripts")
    elif computer == 'cluster':
        data_dir = r'/projects/compbio/users/yli459/gtex_eqtl'
        sys.path.append(r'/projects/compbio/users/yli459/scripts')
    
    # get GENCODE v26 data
    gencode = feather.read_feather(os.path.join(data_dir, 'gencode_v26.ftr'))
    gencode = gencode[gencode['type']=='gene']
    
    # get counts for experiment
    if experiment == '9_tissue':
        counts = feather.read_feather(os.path.join(data_dir, exp_suffix, 'counts.ftr'))
    elif experiment == 'brain' or experiment == 'all':
        counts = feather.read_feather(os.path.join(data_dir, exp_suffix, 'PMM_imputed_counts.ftr'))
    
    # make sure genes are in GENCODE
    counts = counts[counts['gene'].isin(gencode['gene_name'])]
    
    # get genotypes for chromosome
    if genotype_type == 'binary':
        gt = pd.read_csv(os.path.join(data_dir, 'genotypes', f'chr{chrom}_genotypes_10kb_binary.txt.gz'), sep='\t', header=0)
    elif genotype_type == 'additive':
        gt = pd.read_csv(os.path.join(data_dir, 'genotypes', f'chr{chrom}_genotypes_10kb_additive.txt.gz'), sep='\t', header=0)
    
    # only keep genes in genotypes
    gt = gt[gt['gene_name'].isin(counts['gene'])]
    
    # find list of commons donors between counts and genotypes
    common_donors = sorted(list(set(counts['donor']) & set(gt.iloc[:,2:].columns)))
    
    # only keep common donors
    counts = counts[counts['donor'].isin(common_donors)]
    
    keep_cols = ['gene_name','ID'] + common_donors
    gt = gt.filter(items=keep_cols, axis=1)
    
    # check it's the same set of donors and genes
    assert set(gt.iloc[:,2:].columns) == set(counts['donor'])
    assert set(gt['gene_name']).issubset(set(counts['gene']))
    
    # define counts and labels
    num_genes = len(set(gt['gene_name']))
    num_donors = len(set(counts['donor']))
    
    t1 = time.time()
    print(f"Took {t1-t0:.1f} sec to load data")
    print(f'There are a total of {num_genes} genes and {num_donors} donors in chr{chrom}')
    return counts, gt


# Load data for one chromosome - permuted genotypes

def load_data_chrom_permute(chrom, experiment, computer, index):
    '''
    Loads permuted counts and genotypes.
    
    Inputs:
        - chrom = chromosome number
        - experiment = '9_tissue', 'brain', or 'all'
        - computer = 'pc', 'laptop', or 'cluster'
        - index = permutation number: 1-10
    Returns:
        - counts, genotypes
    '''
    t0 = time.time()
    
    if experiment == '9_tissue':
        suffix = ''
    elif experiment == 'brain':
        suffix = '_brain'
    elif experiment == 'all':
        suffix = '_all'
    
    if computer == 'pc':
        data_dir = r'D:\ronnieli\My Drive\lab_data\gtex_eqtl'
        sys.path.append(r"C:\Users\ronnieli\OneDrive\Documents\Lab\gtex_eqtl\scripts")
    elif computer == 'cluster':
        data_dir = r'/projects/compbio/users/yli459/gtex_eqtl'
        sys.path.append(r'/projects/compbio/users/yli459/scripts')
    elif computer == 'laptop':
        data_dir = r"C:\Users\ronnieli\My Drive\lab_data\gtex_eqtl"
        sys.path.append(r"C:\Users\ronnieli\OneDrive\Documents\Lab\gtex_eqtl\scripts")
    
    # get GENCODE v26 data
    gencode = feather.read_feather(os.path.join(data_dir, 'gencode_v26.ftr'))
    gencode = gencode[gencode['type']=='gene']
    
    # get counts for experiment
    if experiment == '9_tissue':
        counts = feather.read_feather(os.path.join(data_dir, f'multi_tissue{suffix}', 'counts.ftr'))
    elif experiment == 'brain' or experiment == 'all':
        counts = feather.read_feather(os.path.join(data_dir, f'multi_tissue{suffix}', 'MICE_imputed_counts.ftr'))
    
    # make sure genes are in GENCODE
    counts = counts[counts['gene'].isin(gencode['gene_name'])]
    
    # get genotypes - just to identify which SNPs to use per gene
    gt = feather.read_feather(os.path.join(data_dir, 'genotypes', f'chr{chrom}_genotypes.ftr'))
    gt = gt[gt['gene'].isin(counts['gene'])]
    gt = gt[['gene','ID']]
    
    # get permuted genotypes
    gt_permute = feather.read_feather(os.path.join(data_dir, f'multi_tissue{suffix}',
                                                   'genotype_permutations', 
                                                   f'permutation_{index}.ftr'))
    gt_permute = gt_permute[gt_permute['SNP'].str.startswith(f'chr{chrom}_')]
    
    # merge genotypes together
    gt_new = gt.merge(gt_permute, left_on='ID', right_on='SNP')
    gt_new.drop('SNP', axis=1, inplace=True)
    gt_new.sort_values(by='gene', ascending=True, inplace=True)
    
    # only keep common donors
    counts = counts[counts['donor'].isin(list(gt_new.columns[2:]))]
    
    # check it's the same set of donors and genes
    assert set(gt_new.iloc[:,2:].columns) == set(counts['donor'])
    assert set(gt_new['gene']).issubset(set(counts['gene']))
    
    # define counts and labels
    num_genes = len(set(gt_new['gene']))
    num_donors = len(set(counts['donor']))
    
    t1 = time.time()
    print(f"Took {t1-t0:.1f} sec to load data")
    print(f'There are a total of {num_genes} genes and {num_donors} donors in chr{chrom}')
    return counts, gt_new
