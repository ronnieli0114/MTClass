# -*- coding: utf-8 -*-
"""
Created on Mon May 22 02:47:44 2023
@author: ronnieli
"""

import pandas as pd
import numpy as np
import pyarrow.feather as feather
import os, sys
import subprocess
import ensembl_rest

def extract_genotypes(gene_list, computer, genotypes='binary'):
    
    '''
    Extracts variants for a list of genes using PLINK
    
    Inputs:
        - gene_list: a list of genes
        - computer name: 'pc' or 'cluster'
        - genotypes: 'binary' or 'additive'
    Returns:
        - A genotype matrix
    '''
    
    all_gts = []
    
    for gene in gene_list:
        
        # define PLINK location and directories
        if computer == 'pc':    
            plink_cmd = 'plink'
            out_dir = r'/mnt/e/lab_data/gtex_exoneqtl'
            out_name = rf'/mnt/e/lab_data/gtex_exoneqtl/{gene}'
        elif computer == 'cluster':
            plink_cmd = "plink"
            out_dir = r"/projects/compbio/users/yli459/gtex_exoneqtl"
            out_name = rf"/projects/compbio/users/yli459/gtex_exoneqtl/{gene}"
        
        # Get Ensembl information about gene
        try:
            data = ensembl_rest.symbol_lookup(species='homo sapiens', symbol=gene)
        except:
            continue
            
        gene_id = data['id']
        chrom = data['seq_region_name']
        start = int(data['start'])-10_000
        end = int(data['end'])-10_000
        
        if computer == 'pc':
            bfile = r'/mnt/e/lab_data/gtex_eqtl/gtex_genotypes'
        elif computer == 'cluster':
            bfile = r'/projects/compbio/users/yli459/gtex_eqtl/gtex_genotypes'
        
        entire_cmd = [plink_cmd, '--bfile', bfile, '--chr', chrom, 
                      '--maf', '0.05', '--mind', '0.1', '--geno', '0.1',
                      '--from-bp', str(start), '--to-bp', str(end), 
                      '--recode', 'vcf', '--out', out_name, 
                      '--memory', '9000000', '--silent']  
        
        process = subprocess.Popen(entire_cmd)
        process.communicate()
        
        if process.returncode == 0:
            print(f'Extracted variants for gene {gene}', flush=True)
    
            # read the file
            gts = pd.read_csv(os.path.join(out_dir, f"{gene}.vcf"), sep='\t', header=6)
            gts['gene_name'] = gene 
            gts['gene_id'] = gene_id
            discard_cols = ['POS','REF','ALT','QUAL','FILTER','INFO','FORMAT']
            gts = gts.loc[:,[col for col in gts.columns if col not in discard_cols]]
            gts.set_index(['#CHROM','gene_name','gene_id','ID'], inplace=True)
            gts.columns = list(map(lambda x: x.split('_')[0], gts.columns))
            gts = gts.reset_index()
            
            # replace genotypes with binary values
            if genotypes == 'binary':
                gts.replace(to_replace={'0/0':0, '0/1':1, '1/1':1, './.':9}, inplace=True)
            elif genotypes == 'additive':
                gts.replace(to_replace={'0/0':0, '0/1':1, '1/1':2, './.':9}, inplace=True)
          
            # remove files generated by PLINK
            os.remove(os.path.join(out_dir, f'{gene}.log'))
            os.remove(os.path.join(out_dir, f'{gene}.nosex'))
            os.remove(os.path.join(out_dir, f'{gene}.vcf'))
            if os.path.exists(os.path.join(out_dir, f'{gene}.irem')):
                os.remove(os.path.join(out_dir, f'{gene}.irem'))
    
            all_gts.append(gts.copy())
        
        else:
            print(f'Extracting variants FAILED for gene {gene}', flush=True)
            print(process.stderr)
            continue
        
    final_df = pd.concat(all_gts)
    return final_df

def list_all_brain_genes(computer='pc'):
    
    if computer == 'pc':
        data_dir = r'E:/lab_data'
    elif computer == 'cluster':
        data_dir = r'/projects/compbio/users/yli459'
    
    files = os.listdir(os.path.join(data_dir, 'gtex_exoneqtl', 'tissue_level_expression'))
    files = [file for file in files if 'Brain_' in file and 'gene_list.txt' not in file]
    brain_tissues = [file.split('_counts.ftr')[0] for file in files]
    
    all_genes = []
    
    for tissue in brain_tissues:
        exp = feather.read_feather(os.path.join(data_dir, 'gtex_exoneqtl', 'tissue_level_expression', f'{tissue}_counts.ftr'))
        all_genes.append(set(exp['gene']))
    
    all_genes = list(set.union(*all_genes))
    return sorted(all_genes)

#%%
comp = 'cluster'
gene_list = list_all_brain_genes(comp)
print('Retrieved gene list', flush=True)
print('Now extracting binary genotypes...', flush=True)
all_gts_binary = extract_genotypes(gene_list, comp, genotypes='binary')
print('Now extracting additive genotypes...', flush=True)
all_gts_additive = extract_genotypes(gene_list, comp, genotypes='additive')

for c, gts in all_gts_binary.groupby('#CHROM'):
    
    gts.drop('#CHROM', axis=1, inplace=True)
    print(gts.head())
    gts.to_csv(rf'/projects/compbio/users/yli459/gtex_exoneqtl/genotypes/chr{c}_genotypes_binary.txt.gz', sep='\t', index=False)
    print(f'Wrote chr{c} binary genotypes to file', flush=True)

for c, gts in all_gts_additive.groupby('#CHROM'):
    
    gts.drop('#CHROM', axis=1, inplace=True)
    print(gts.head())
    gts.to_csv(rf'/projects/compbio/users/yli459/gtex_exoneqtl/genotypes/chr{c}_genotypes_additive.txt.gz', sep='\t', index=False)
    print(f'Wrote chr{c} additive genotypes to file', flush=True)

print('=== DONE ===')