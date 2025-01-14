# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 14:35:22 2022
@author: ronnieli
"""

import pandas as pd
import pyarrow.feather as feather
import os

def write_vcf_colon_sep(variant_list, filename):
    
    """ 
    Makes a VCF file out of a list of colon-separated variants.
    
    Input: list of colon-separated variants
        - CHR#:POS#:REF:ALT
    Output: VCF file saved to Desktop
    """
    
    variant_df = pd.DataFrame(variant_list)
    variant_df['chr'] = variant_df[0].apply(lambda x: x.split(':')[0]).astype(int)
    variant_df['pos'] = variant_df[0].apply(lambda x: x.split(':')[1]).astype(int)
    variant_df.sort_values(by=['chr','pos'], ascending=[True, True], inplace=True)
    
    with open(rf"C:/Users/ronni/OneDrive/Desktop/{filename}.vcf", "w") as file:
        
        # write header
        file.write('##fileformat=VCFv4.2')
        file.write('\n')
        file.write('#CHROM\t')
        file.write('POS\t')
        file.write('ID\t')
        file.write('REF\t')
        file.write('ALT\t')
        file.write('QUAL\t')
        file.write('FILTER\t')
        file.write('INFO\t')
        file.write('\n')
        
        # write SNPs
        for variant in variant_df[0].tolist():
            chrom = variant.split(':')[0]
            pos = variant.split(':')[1]
            ref = variant.split(':')[2]
            alt = variant.split(':')[3]
            
            file.write(chrom)
            file.write('\t')
            file.write(pos)
            file.write('\t')
            file.write(variant)
            file.write('\t')
            file.write(ref)
            file.write('\t')
            file.write(alt)
            file.write('\t')
            file.write('.\t')
            file.write('.\t')
            file.write('.\t')
            file.write('\n')


def write_vcf_gtex(variant_list, filename):
    
    """ 
    Makes a VCF file out of a list of GTEx formatted variants.
    
    Input: list of GTEx-formatted variants
    Output: VCF file saved to Desktop
    
    """
    
    variant_df = pd.DataFrame(variant_list)
    variant_df['chr'] = variant_df[0].apply(lambda x: x.split('_')[0][3:]).astype(int)
    variant_df['pos'] = variant_df[0].apply(lambda x: x.split('_')[1]).astype(int)
    variant_df.sort_values(by=['chr','pos'], ascending=[True, True], inplace=True)
    
    with open(rf"C:/Users/ronni/OneDrive/Desktop/{filename}.vcf", "w") as file:
        file.write('##fileformat=VCFv4.2')
        file.write('\n')
        file.write('#CHROM\t')
        file.write('POS\t')
        file.write('ID\t')
        file.write('REF\t')
        file.write('ALT\t')
        file.write('QUAL\t')
        file.write('FILTER\t')
        file.write('INFO\t')
        file.write('\n')
        
        for variant in variant_df[0].tolist():
            chrom = variant.split('_')[0][3:]
            pos = variant.split('_')[1]
            ref = variant.split('_')[2]
            alt = variant.split('_')[3]
            
            file.write(chrom)
            file.write('\t')
            file.write(pos)
            file.write('\t')
            file.write(variant)
            file.write('\t')
            file.write(ref)
            file.write('\t')
            file.write(alt)
            file.write('\t')
            file.write('.\t')
            file.write('.\t')
            file.write('.\t')
            file.write('\n')
    

def get_binary_exp(gene, variant, experiment, data_dir=r'C:/Users/ronni/OneDrive/Documents/lab_data/mtclass_eqtl'):
    
    """
    Retrieves expression levels for a gene and variant, separated by genotype.
    Assumes genotypes are binary, 0 and 1 only (dominant model).
    Multi-tissue study only.
    
    Arguments:
        - gene symbol
        - GTEx formatted variant ID
        - experiment = '9_tissue', 'brain', or '48_tissue' for counts retrieval
        
    Returns:
        - counts_0 (genotype 0, AA)
        - counts_1 (genotype 1, Aa/aa)
    """
        
    if experiment == '9_tissue':
        counts = feather.read_feather(os.path.join(data_dir, 'multi-tissue/expression/9_tissue_TPM.ftr'))
    elif experiment == 'brain':
        counts = feather.read_feather(os.path.join(data_dir, 'multi-tissue/expression/brain_tissue_TPM_PMMimpute.ftr'))
    elif experiment == '48_tissue':
        counts = feather.read_feather(os.path.join(data_dir, 'multi-tissue/expression/48_tissue_TPM_PMMimpute.ftr'))
    else:
        raise ValueError("Experiment name not found")
    
    # get expression levels
    X = counts[counts['gene']==gene]
    X = X.drop('gene', axis=1).set_index('donor')
    
    # get genotypes
    c = variant.split('_')[0][3:] # chromosome number
    gt = pd.read_csv(os.path.join(data_dir, f'multi-tissue/genotypes/chr{c}_genotypes_10kb_binary.txt.gz'), sep='\t', header=0)
    y = gt[(gt['ID']==variant) & (gt['gene_name']==gene)]
    y = y.iloc[:,3:].T
    y.columns = ['genotype']
    
    # merge
    Xy = X.merge(y, left_index=True, right_index=True)
    
    counts_0 = Xy[Xy['genotype']==0].iloc[:,:-1]
    counts_1 = Xy[Xy['genotype']==1].iloc[:,:-1]
    return counts_0, counts_1

def get_additive_exp(gene, variant, experiment, data_dir=r'C:/Users/ronni/OneDrive/Documents/lab_data/mtclass_eqtl'):
    
    """
    Retrieves counts for a gene and variant, separated by genotype.
    Assumes genotypes are 0, 1, and 2 (additive model).
    Multi-tissue study only.
    
    Arguments:
        - gene name
        - variant ID
        - experiment = '9_tissue', 'brain', or '48_tissue'
        
    Returns:
        - counts_0 (genotype 0, AA)
        - counts_1 (genotype 1, Aa)
        - counts_2 (genotype 2, aa)
    """

    # get expression levels
    if experiment == '9_tissue':
        counts = feather.read_feather(os.path.join(data_dir, 'multi-tissue/expression/9_tissue_TPM.ftr'))
    elif experiment == 'brain':
        counts = feather.read_feather(os.path.join(data_dir, 'multi-tissue/expression/brain_tissue_TPM_PMMimpute.ftr'))
    elif experiment == '48_tissue':
        counts = feather.read_feather(os.path.join(data_dir, 'multi-tissue/expression/48_tissue_TPM_PMMimpute.ftr'))
    else:
        raise ValueError("Experiment name not found")
    
    X = counts[counts['gene_name']==gene]
    X = X.drop('gene', axis=1).set_index('donor')
    
    # get genotypes
    c = variant.split('_')[0][3:]
    gt = pd.read_csv(os.path.join(data_dir, f'multi-tissue/genotypes/chr{c}_genotypes_10kb_binary.txt.gz'), sep='\t', header=0)
    y = gt[(gt['ID']==variant) & (gt['gene_name']==gene)]
    y = y.iloc[:,3:].T
    y.columns = ['genotype']
    
    # merge
    Xy = X.merge(y, left_index=True, right_index=True)
    
    counts_0 = Xy[Xy['genotype']==0].iloc[:,:-1]
    counts_1 = Xy[Xy['genotype']==1].iloc[:,:-1]
    counts_2 = Xy[Xy['genotype']==2].iloc[:,:-1]
    return counts_0, counts_1, counts_2

def write_top_snps_to_excel(N=1000):
    
    multi_tissue_dir = r'G:/My Drive/Lab/lab_projects/mtclass_eqtl/multi-tissue/results'
    multi_exon_dir = r'G:/My Drive/Lab/lab_projects/mtclass_eqtl/multi-exon/results'
    isoqtl_dir = r'G:/My Drive/Lab/lab_projects/mtclass_eqtl/isoqtl/results'
    
    out_file = r"C:/Users/ronni/Desktop/MTClass_top_eQTL.xlsx"

    with pd.ExcelWriter(out_file) as writer:
        
        # Write multi-tissue results
        print("Writing multi-tissue results to file...")
        ninetiss = pd.read_csv(os.path.join(multi_tissue_dir, '9_tissue/9_tissue_ensemble_aggregate.txt.gz'), sep='\t', header=0)
        braintiss = pd.read_csv(os.path.join(multi_tissue_dir, 'brain_tissue/brain_ensemble_aggregate.txt.gz'), sep='\t', header=0)
        
        for sheet_name, data in zip(['9_tissue','Brain_tissue'], [ninetiss, braintiss]):
            top_data = data.sort_values('f1_macro_median', ascending=False).iloc[:N,:]
            top_data = top_data[['gene','variant','f1_macro_median','mcc_median']].reset_index(drop=True)
            top_data.to_excel(writer, sheet_name=sheet_name, index=False)
    
        # Write multi-exon results
        print("Writing multi-exon results to file...")
        brain_tissues = [f for f in os.listdir(multi_exon_dir) if 'ensemble_aggregate.txt.gz' in f and 'Brain' in f]
        
        for file in brain_tissues:
            
            tissue_name = file.split('_ensemble_aggregate.txt.gz')[0].replace('Brain_','')
            
            res = pd.read_csv(os.path.join(multi_exon_dir, file), sep='\t', header=0)
            res = res.sort_values('f1_macro_median', ascending=False).iloc[:N,:]
            res = res[['gene','variant','f1_macro_median','mcc_median']].reset_index(drop=True)
            res.to_excel(writer, sheet_name=tissue_name, index=False)
            
        # Write isoQTL results
        print("Writing isoQTL results...")
        isoqtl = pd.read_csv(os.path.join(isoqtl_dir, 'isoqtl_ensemble_aggregate.txt.gz'), sep='\t', header=0)
        
        def rename_snp(snp):
            chrom = snp.split('_')[0][3:]
            pos = snp.split('_')[1]
            ref = snp.split('_')[2]
            alt = snp.split('_')[3]
            return str(chrom)+'_'+str(pos)+'_'+ref+'_'+alt+'_b37'
        
        isoqtl['variant'] = isoqtl['variant'].apply(rename_snp)
        
        res = isoqtl.sort_values('f1_macro_median', ascending=False).iloc[:N,:]
        res = res[['gene','variant','f1_macro_median','mcc_median']].reset_index(drop=True)
        res.to_excel(writer, sheet_name='PFC_isoQTL', index=False)
        
        # Write 2D results
        print("Writing 2D results...")
        exon_tissue = pd.read_csv(os.path.join(multi_exon_dir, '2d_exon_tissue', '2D_exon_tissue_NN.txt.gz'), sep='\t', header=0)
        
        res = exon_tissue.sort_values('f1_macro', ascending=False).iloc[:N,:]
        res = res[['gene_name','variant','f1_macro','mcc']].reset_index(drop=True)
        res.to_excel(writer, sheet_name='2D_Study', index=False)
        
    print('=== done ===')