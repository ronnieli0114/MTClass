# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 08:23:47 2023

@author: ronnieli
"""

import pyarrow.feather as feather
import pandas as pd
import os
import ensembl_rest
import subprocess

#%%
def generate_gene_info(data_dir):
    
    '''
    Generate a gene list of genes:
        - from top exon expression levels (at least 10 exons)
        - on autosomes only
        - with genomic coordinates available in Ensembl
    Input:
        - data directory where top_exon_counts.ftr is located
    Returns:
        - pandas DataFrame with gene information
    '''

    exp = feather.read_feather(os.path.join(data_dir, 'top_exon_counts.ftr')) # counts with at least 10 exons
    
    all_genes = sorted(list(set(exp.index)))
    all_genes = [gene.split('.')[0] for gene in all_genes]
    all_genes = sorted(list(set(all_genes)))
    
    print(f'There are {len(all_genes)} genes in total')
    
    gene_dict = {'gene_id':[], 'gene_name':[], 'chr':[], 'gene_start':[], 'gene_end':[]}
    
    for i, gene_id in enumerate(all_genes):
        
        if i % 100 == 0:
            print(f'Now on {i}th gene {gene_id}...')
        
        # lookup Ensembl ID
        try:
            ensembl_dict = ensembl_rest.lookup(gene_id, species='human') # lookup stable ID
        except:
            continue
        chrom = ensembl_dict['seq_region_name']
        start = ensembl_dict['start']
        end = ensembl_dict['end']
        try:
            gene_name = ensembl_dict['display_name']
        except KeyError:
            gene_name = 'NA'
        
        if chrom == 'X' or chrom == 'Y' or chrom == 'MT':
            continue
        
        gene_dict['gene_id'].append(gene_id)
        gene_dict['gene_name'].append(gene_name)
        gene_dict['chr'].append(chrom)
        gene_dict['gene_start'].append(start)
        gene_dict['gene_end'].append(end)
        
    genes = pd.DataFrame(gene_dict)
    genes['chr'] = genes['chr'].astype(int)
    genes['gene_start'] = genes['gene_start'].astype(int)
    genes['gene_end'] = genes['gene_end'].astype(int)
    genes = genes.sort_values(by=['chr','gene_start'], ascending=True)
    return genes

#%%
def extract_genotypes(gene_data, genotypes='binary'):
    
    '''
    Extracts variants for a list of genes using PLINK
    and the genomic coordinates from Ensembl
    
    Inputs:
        - gene_info: a pandas DataFrame of gene information
        - genotypes: 'binary' for dominant model, or 'additive' for dosage
    Returns:
        - A genotype matrix with columns 'gene','ID', and all donors
    '''
    
    gt_dir = '/projects/compbio/users/yli459/gtex_exoneqtl/genotypes/'
    
    for c in range(1,23):
        
        chrom_gts = []
        
        gene_data = gene_data[gene_data['chr']==c]  
    
        for _, row in gene_data.iterrows():
            
            gene_id = row['gene_id']
            
            # get GENCODE information about the gene
            chrom = row['chr']
            start = row['gene_start']-100_000
            end = row['gene_end']+100_000
            if start <= 0:
                start = 1
            
            genofile = r'/projects/compbio/users/yli459/gtex_eqtl/gtex_genotypes'
            
            # define PLINK location and directories
            plink_cmd = 'plink'
            out_dir = r'/projects/compbio/users/yli459/gtex_exoneqtl/data'
            out_name = rf'/projects/compbio/users/yli459/gtex_exoneqtl/data/{gene_id}'
                    
            entire_cmd = [plink_cmd, '--bfile', genofile, '--chr', str(chrom), 
                          '--maf', '0.05', '--mind', '0.1', '--geno', '0.1',
                          '--from-bp', str(start), '--to-bp', str(end), 
                          '--const-fid', '--recode', 'vcf', '--out', out_name, 
                          '--memory', '9000000']  
    
            process = subprocess.Popen(entire_cmd)
            process.communicate()
                    
            if process.returncode == 0:
                print(f'Extracted variants for gene {gene_id}')
    
                # read the file
                gts = pd.read_csv(os.path.join(out_dir, f"{gene_id}.vcf"), sep='\t', header=6)
                gts['gene_id'] = gene_id
                gts['ID'] = gts['#CHROM'].astype(str)+':'+gts['POS'].astype(str)+':'+gts['REF']+':'+gts['ALT']
                discard_cols = ['#CHROM','POS','REF','ALT','QUAL','FILTER','INFO','FORMAT']
                gts = gts.loc[:,[col for col in gts.columns if col not in discard_cols]]
                gts.set_index(['gene_id','ID'], inplace=True)
                
                gts.columns = list(map(lambda x: x.split('_')[0], gts.columns))
                    
                sorted_cols = sorted(gts.columns)
                gts = gts.loc[:,sorted_cols]
                gts = gts.reset_index()
                    
                # replace genotypes with binary values
                if genotypes == 'binary':
                    gts.replace(to_replace={'0/0':0, '0/1':1, '1/1':1, './.':9}, inplace=True)
                elif genotypes == 'additive':
                    gts.replace(to_replace={'0/0':0, '0/1':1, '1/1':2, './.':9}, inplace=True)
              
                # remove files generated by PLINK
                os.remove(os.path.join(out_dir, f'{gene_id}.log'))
                os.remove(os.path.join(out_dir, f'{gene_id}.nosex'))
                os.remove(os.path.join(out_dir, f'{gene_id}.vcf'))
                if os.path.exists(os.path.join(out_dir, f'{gene_id}.irem')):
                    os.remove(os.path.join(out_dir, f'{gene_id}.irem'))
        
                chrom_gts.append(gts.copy())
            
            else:
                print(f'Extracting variants FAILED for gene {gene_id}')
                print(process.stderr)
                continue
            
        chrom_df = pd.concat(chrom_gts)
        chrom_df.to_csv(os.path.join(gt_dir, f'chr{c}_genotypes_{genotypes}.txt.gz'), sep='\t', index=False, header=True)
        print(f'===== WROTE RESULTS FOR CHR{c} =====')
        
#%%
# data_dir = r'E:/lab_data/gtex_exoneqtl'
# gene_info = generate_gene_info(data_dir)
# gene_info.to_csv(os.path.join(data_dir, 'top_exon_gene_info.txt.gz'), sep='\t', header=True, index=False)

data_dir = '/projects/compbio/users/yli459/gtex_exoneqtl'
gene_data = pd.read_csv(os.path.join(data_dir, 'top_exon_gene_info.txt.gz'), sep='\t', header=0)

extract_genotypes(gene_data, genotypes='binary')
extract_genotypes(gene_data, genotypes='additive')