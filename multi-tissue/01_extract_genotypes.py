# -*- coding: utf-8 -*-
"""
Uses PLINK to extract variants based on chromosome
and genomic coordinates, given a gene from the
GENCODE v26 database, the same one used for the GTEx v8 data.

Created on Mon Aug 29 05:57:59 2022
@author: ronnieli 
"""

import pandas as pd
import numpy as np
import pyarrow.feather as feather
import os, sys
import subprocess
import ensembl_rest

def extract_genotypes(gene_list, genotypes='binary'):
    
    '''
    Extracts variants for a list of genes using PLINK
    
    Inputs:
        - gene_list: a list of genes
        - genotypes: 'binary' or 'additive'
    Returns:
        - A genotype matrix
    '''
    
    all_gts = []
    
    for gene in gene_list:
        
        # define PLINK location and directories
        plink_cmd = "plink"
        out_dir = r"/projects/compbio/users/yli459/mtclass_eqtl"
        out_name = rf"/projects/compbio/users/yli459/mtclass_eqtl/{gene}"

        # Get Ensembl information about gene        
        try:
            data = ensembl_rest.symbol_lookup(species='homo sapiens', symbol=gene)
        except:
            continue

        gene_id = data['id']
        chrom = data['seq_region_name']
        start = int(data['start'])-10_000
        end = int(data['end'])-10_000
        
        bfile = r'/projects/compbio/users/yli459/mtclass_eqtl/gtex_genotypes'
        
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

def get_gene_list(data_dir):
    
    top_counts = feather.read_feather(os.path.join(data_dir,'top_TPM.ftr'))
    return sorted(list(set(top_counts['Description'])))

save_dir = r"/peojcts/compbio/users/yli459/mtclass_eqtl/genotypes"
gene_list = get_gene_list("/projects/compbio/users/yli459/mtclass_eqtl/")
print('Retrieved gene list', flush=True)
print('Now extracting binary genotypes...', flush=True)
all_gts_binary = extract_genotypes(gene_list, genotypes='binary')
print('Now extracting additive genotypes...', flush=True)
all_gts_additive = extract_genotypes(gene_list, genotypes='additive')

for c, gts in all_gts_binary.groupby('#CHROM'):
    
    gts.drop('#CHROM', axis=1, inplace=True)
    print(gts.head())
    gts.to_csv(os.path.join(save_dir, f'chr{c}_genotypes_binary.txt.gz'), sep='\t', index=False)
    print(f'Wrote chr{c} binary genotypes to file', flush=True)

for c, gts in all_gts_additive.groupby('#CHROM'):
    
    gts.drop('#CHROM', axis=1, inplace=True)
    print(gts.head())
    gts.to_csv(os.path.join(save_dir, f'chr{c}_genotypes_additive.txt.gz'), sep='\t', index=False)
    print(f'Wrote chr{c} additive genotypes to file', flush=True)

print('=== done ===')