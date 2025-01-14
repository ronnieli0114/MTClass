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

def load_gwas_table(data_dir=r'C:/Users/ronni/OneDrive/Documents/lab_data/mtclass_eqtl', genome='hg38'):
    
    '''
    Loads and processes the GWAS Catalog
    Inputs:
        - data_dir = data directory
    '''
    gwas_table = pd.read_table(os.path.join(data_dir, 'gwas_catalog_associations_processed.tsv'), header=0, low_memory=False)
    
    # Define appropriate columns and make GenomicRanges object
    gwas_table['seqnames'] = gwas_table['CHR_ID'].apply(lambda x: "chr"+str(x))
    gwas_table['starts'] = gwas_table[f'pos_{genome}'].astype(int)
    gwas_table['ends'] = gwas_table[f'pos_{genome}'].astype(int)
    print('loaded GWAS table')
    return gwas_table.loc[:,['seqnames','starts','ends']]

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

def gwas_hits_iterations(experiment):
    
    ''' Compare GWAS hits between 3 iterations of MTClass and MultiPhen/MANOVA '''
    
    gwas_table = load_gwas_table()
    
    if experiment == '9_tissue':
        os.chdir(r'C:/Users/ronni/OneDrive/Documents/Lab/lab_projects/gtex_eqtl/results/9_tissue')
        m1 = pd.read_csv('9_tissue_ensemble_1.txt.gz', sep='\t', header=0)
        m2 = pd.read_csv('9_tissue_ensemble_2.txt.gz', sep='\t', header=0)
        m3 = pd.read_csv('9_tissue_ensemble_3.txt.gz', sep='\t', header=0)
        mphen_binary = pd.read_csv('9_tissue_multiphen_binary.txt.gz', sep='\t', header=0)
        mphen_add = pd.read_csv('9_tissue_multiphen_additive.txt.gz', sep='\t', header=0)
        manova = pd.read_csv('9_tissue_manova.txt.gz', sep='\t', header=0)
        title = '9-tissue'
        sample_size_threshold = 50
        
    elif experiment == 'brain':
        os.chdir(r'C:/Users/ronni/OneDrive/Documents/Lab/lab_projects/gtex_eqtl/results/brain_tissue')
        m1 = pd.read_csv('brain_ensemble_1.txt.gz', sep='\t', header=0)
        m2 = pd.read_csv('brain_ensemble_2.txt.gz', sep='\t', header=0)
        m3 = pd.read_csv('brain_ensemble_3.txt.gz', sep='\t', header=0)
        mphen_binary = pd.read_csv('brain_multiphen_binary.txt.gz', sep='\t', header=0)
        mphen_add = pd.read_csv('brain_multiphen_additive.txt.gz', sep='\t', header=0)
        manova = pd.read_csv('brain_manova.txt.gz', sep='\t', header=0)
        title = 'Brain tissue'
        sample_size_threshold = 150
    
    elif experiment == '48_tissue':
        os.chdir(r'C:/Users/ronni/OneDrive/Documents/Lab/lab_projects/gtex_eqtl/results/48_tissue')
        m1 = pd.read_csv('48_tissue_ensemble_1.txt.gz', sep='\t', header=0)
        m2 = pd.read_csv('48_tissue_ensemble_2.txt.gz', sep='\t', header=0)
        m3 = pd.read_csv('48_tissue_ensemble_3.txt.gz', sep='\t', header=0)
        mphen_binary = pd.read_csv('48_tissue_multiphen_binary.txt.gz', sep='\t', header=0)
        mphen_add = pd.read_csv('48_tissue_multiphen_additive.txt.gz', sep='\t', header=0)
        manova = pd.read_csv('48_tissue_manova.txt.gz', sep='\t', header=0)
        title = '48-tissue'
        sample_size_threshold = 300
    
    else:
        raise ValueError("Experiment name not found!")
    
    # Make sure sample size is sufficient
    print(f'Using sample size threshold of >={sample_size_threshold}')
    m1 = m1[m1['sample_size'] >= sample_size_threshold]
    m2 = m2[m2['sample_size'] >= sample_size_threshold]
    m3 = m3[m3['sample_size'] >= sample_size_threshold]
    
    m1['pair'] = m1['gene'] + '-' + m1['variant']
    m2['pair'] = m2['gene'] + '-' + m2['variant']
    m3['pair'] = m3['gene'] + '-' + m3['variant']
    mphen_binary['pair'] = mphen_binary['gene'] + '-' + mphen_binary['variant']
    mphen_add['pair'] = mphen_add['gene'] + '-' + mphen_add['variant']
    manova['pair'] = manova['gene'] + '-' + manova['variant']
    
    all_pairs = [set(data['pair']) for data in [m1, m2, m3, mphen_binary, mphen_add, manova]]
    common_pairs = set.intersection(*all_pairs)
    
    m1 = m1[m1['pair'].isin(common_pairs)]
    m2 = m2[m2['pair'].isin(common_pairs)]
    m3 = m3[m3['pair'].isin(common_pairs)]
    mphen_binary = mphen_binary[mphen_binary['pair'].isin(common_pairs)]
    mphen_add = mphen_add[mphen_add['pair'].isin(common_pairs)]
    manova = manova[manova['pair'].isin(common_pairs)]
    
    m1.drop_duplicates('pair', inplace=True)
    m2.drop_duplicates('pair', inplace=True)
    m3.drop_duplicates('pair', inplace=True)
    mphen_binary.drop_duplicates('pair', inplace=True)
    mphen_add.drop_duplicates('pair', inplace=True) 
    manova.drop_duplicates('pair', inplace=True)
    
    print('Loaded MTClass, MultiPhen, MANOVA')
    print(f'There are {len(common_pairs)} total gene-SNP pairs')
    
    final_dict = {}
    
    # Calculate GWAS hits
    for metric, name in zip(['f1_macro','mcc'],['Macro F1','MCC']):
        
        result_dict = {'Top variants':[], 'MTClass 1':[], 'MTClass 2':[], 'MTClass 3':[], 
                       'MultiPhen (binary)':[], 'MultiPhen (additive)':[], 'MANOVA':[]}
    
        n_top_list = [100, 500, 1000, 2000, 3000, 4000, 5000]
        
        for n_top in n_top_list:
            
            m1_count = gwas_hits(m1, gwas_table, metric=metric, n_top=n_top, bp=10000)
            m2_count = gwas_hits(m2, gwas_table, metric=metric, n_top=n_top, bp=10000)
            m3_count = gwas_hits(m3, gwas_table, metric=metric, n_top=n_top, bp=10000)
            mphen_count_add = gwas_hits(mphen_add, gwas_table, metric='pval', n_top=n_top, bp=10000)
            mphen_count_bin = gwas_hits(mphen_binary, gwas_table, metric='pval', n_top=n_top, bp=10000)
            manova_count = gwas_hits(manova, gwas_table, metric='pval', n_top=n_top, bp=10000)
            
            result_dict['Top variants'].append(n_top)
            result_dict['MTClass 1'].append(m1_count)
            result_dict['MTClass 2'].append(m2_count)
            result_dict['MTClass 3'].append(m3_count)
            result_dict['MultiPhen (binary)'].append(mphen_count_bin)
            result_dict['MultiPhen (additive)'].append(mphen_count_add)
            result_dict['MANOVA'].append(manova_count)
            print(pd.DataFrame(result_dict).tail(n=1))
            
        results = pd.DataFrame(result_dict)
        final_dict[metric] = results
        
        # Plot results
        fig, ax = plt.subplots(figsize=(6,6))
        ax.plot(results['Top variants'], results['MTClass 1'], color='blue', linestyle='solid', label='MTClass 1')
        ax.plot(results['Top variants'], results['MTClass 2'], color='blue', linestyle='dotted', label='MTClass 2')
        ax.plot(results['Top variants'], results['MTClass 3'], color='blue', linestyle='dashed', label='MTClass 3')
        ax.plot(results['Top variants'], results['MultiPhen (binary)'], color='orange', linestyle='solid', label='MultiPhen (binary)')
        ax.plot(results['Top variants'], results['MultiPhen (additive)'], color='orange', linestyle='dotted', label='MultiPhen (additive)')
        ax.plot(results['Top variants'], results['MANOVA'], color='green', linestyle='solid', label='MANOVA')
    
        ax.set_title(f'GWAS variant colocalization\n{name}', fontsize=20)
        ax.set_xlabel('Top variants', fontsize=20)
        ax.set_ylabel('GWAS hits', fontsize=20)
        ax.legend(fontsize=12)
        plt.tight_layout()
        plt.savefig(rf'C:/Users/ronni/OneDrive/Desktop/{experiment}_{metric}.png', dpi=400)
        plt.show()
        
    return final_dict

def gwas_hits_aggregate(experiment):
    
    ''' Compare GWAS hits between aggregated version of MTClass, MultiPhen, MANOVA '''
    
    gwas_table = load_gwas_table(data_dir=r'C:/Users/ronni/OneDrive/Documents/lab_data/mtclass_eqtl', genome='hg38')
    
    if experiment == '9_tissue':
        os.chdir(r'C:/Users/ronni/OneDrive/Documents/lab_data/mtclass_eqtl/results/9_tissue')
        mtclass = pd.read_csv('9_tissue_mtclass.txt.gz', sep='\t', header=0)
        mphen_binary = pd.read_csv('9_tissue_multiphen_binary.txt.gz', sep='\t', header=0)
        manova = pd.read_csv('9_tissue_manova.txt.gz', sep='\t', header=0)
        sample_size_threshold = 50
        
    elif experiment == 'brain':
        os.chdir(r'C:/Users/ronni/OneDrive/Documents/lab_data/mtclass_eqtl/results/brain_tissue')
        mtclass = pd.read_csv('brain_mtclass.txt.gz', sep='\t', header=0)
        mphen_binary = pd.read_csv('brain_multiphen_binary.txt.gz', sep='\t', header=0)
        manova = pd.read_csv('brain_manova.txt.gz', sep='\t', header=0)
        cca = pd.read_csv('brain_cca.txt.gz', sep='\t', header=0)
        cauchy = pd.read_csv('brain_cauchy.txt.gz', sep='\t', header=0)
        lr = pd.read_csv('brain_logreg.txt.gz', sep='\t', header=0)
        sample_size_threshold = 150
    
    elif experiment == '48_tissue':
        os.chdir(r'C:/Users/ronni/OneDrive/Documents/lab_data/mtclass_eqtl/results/48_tissue')
        mtclass = pd.read_csv('48_tissue_mtclass.txt.gz', sep='\t', header=0)
        mphen_binary = pd.read_csv('48_tissue_multiphen_binary.txt.gz', sep='\t', header=0)
        manova = pd.read_csv('48_tissue_manova.txt.gz', sep='\t', header=0)
        sample_size_threshold = 300
        
    else:
        raise ValueError("experiment name not found!")
        
    # Make sure sample size is sufficient
    # test on same gene-SNP pairs
    # print(f'Using sample size threshold of >={sample_size_threshold}')
    # mtclass = mtclass[mtclass['sample_size'] >= sample_size_threshold]
    # mphen_binary = mphen_binary[mphen_binary['pair'].isin(mtclass['pair'])]
    # manova = manova[manova['pair'].isin(mtclass['pair'])]
    # mtclass.drop_duplicates('pair', inplace=True)
    
    print('loaded MTClass, MultiPhen, MANOVA')
    print([mtclass.shape[0], mphen_binary.shape[0], manova.shape[0]])
    
    final_dict = {}
    
    for metric, name in zip(['f1_macro_median','mcc_median'],['Median Macro F1','Median MCC']):
        
        if experiment == 'brain':
            results_dict = {'Top variants':[], 'MTClass':[], 'MultiPhen':[], 'MANOVA':[], 'CCA':[], 'LR':[], 'Cauchy':[]}
        else:
            results_dict = {'Top variants':[], 'MTClass':[], 'MultiPhen':[], 'MANOVA':[]}
        n_top_list = [100, 500, 1000, 2000, 3000, 4000, 5000]
        
        for n_top in n_top_list:
            
            mtclass_count = gwas_hits(mtclass, gwas_table, metric, n_top, bp=10000, genome='hg38')
            mphen_binary_count = gwas_hits(mphen_binary, gwas_table, 'pval', n_top, bp=10000, genome='hg38')
            manova_count = gwas_hits(manova, gwas_table, 'pval', n_top, bp=10000, genome='hg38')
            if experiment=='brain':
                cca_count = gwas_hits(cca, gwas_table, 'pval', n_top, bp=10000, genome='hg38')
                cauchy_count = gwas_hits(cauchy, gwas_table, 'pval_cauchy', n_top, bp=10000, genome='hg38')
                lr_count = gwas_hits(lr, gwas_table, 'pval', n_top, bp=10000, genome='hg38')
                results_dict['CCA'].append(cca_count)
                results_dict['Cauchy'].append(cauchy_count)
                results_dict['LR'].append(lr_count)
            
            results_dict['Top variants'].append(n_top)
            results_dict['MTClass'].append(mtclass_count)
            results_dict['MultiPhen'].append(mphen_binary_count)
            results_dict['MANOVA'].append(manova_count)         
            print(pd.DataFrame(results_dict).tail(n=1))
        
        results = pd.DataFrame(results_dict)
        final_dict[metric] = results
        
        # Plot results
        fig, ax = plt.subplots(figsize=(6,6))
        ax.plot(results['Top variants'], results['MTClass'], color='blue', linestyle='solid', label='MTClass')
        ax.plot(results['Top variants'], results['MultiPhen'], color='orange', linestyle='solid', label='MultiPhen')
        ax.plot(results['Top variants'], results['MANOVA'], color='green', linestyle='solid', label='MANOVA')
        if experiment == 'brain':
            ax.plot(results['Top variants'], results['CCA'], color='pink', linestyle='solid', label='CCA')
            ax.plot(results['Top variants'], results['Cauchy'], color='red', linestyle='solid', label='Cauchy')
            ax.plot(results['Top variants'], results['LR'], color='brown', linestyle='solid', label='LogReg')
        ax.grid(alpha=0.4)
        
        plt.title(f'GWAS variant colocalization\n{name}', fontsize=16)
        plt.ylabel('GWAS hits', fontsize=16)
        plt.xlabel('Top variants', fontsize=16)
        plt.legend(fontsize=12)
        plt.tight_layout()
        plt.savefig(rf'C:/Users/ronni/Desktop/{experiment}_{metric}.png', dpi=400)
        plt.show()
    return final_dict

def gwas_hits_singletissue(experiment):
    
    ''' Compares GWAS hits between MTClass and top single-tissue eQTLs '''
    
    gwas_table = load_gwas_table()
    
    if experiment == "9_tissue":
        
        sample_size_threshold = 50
        os.chdir(r'C:/Users/ronni/OneDrive/Documents/lab_data/mtclass_eqtl/results/9_tissue')
        mtclass = pd.read_csv('9_tissue_mtclass_aggregate.txt.gz', sep='\t', header=0)
        
        tissue_list = ['Adipose_Subcutaneous','Artery_Tibial','Lung','Muscle_Skeletal',
                       'Nerve_Tibial','Skin_Not_Sun_Exposed_Suprapubic','Skin_Sun_Exposed_Lower_leg',
                       'Thyroid','Whole_Blood']
        title = '9-tissue'
    
    elif experiment == "brain":
        
        sample_size_threshold = 150
        os.chdir(r'C:/Users/ronni/OneDrive/Documents/lab_data/mtclass_eqtl/results/brain_tissue')
        mtclass = pd.read_csv('brain_mtclass_aggregate.txt.gz', sep='\t', header=0)
        
        tissue_list = ['Amygdala','Anterior_cingulate_cortex_BA24','Caudate_basal_ganglia','Cerebellar_Hemisphere',
                       'Cerebellum','Cortex','Frontal_Cortex_BA9','Hippocampus','Hypothalamus',
                       'Nucleus_accumbens_basal_ganglia','Putamen_basal_ganglia','Spinal_cord_cervical_c-1',
                       'Substantia_nigra']
        tissue_list = ['Brain_' + x for x in tissue_list]
        title = "Brain tissue"
        
    elif experiment == '48_tissue':
        
        sample_size_threshold = 300
        os.chdir(r'C:/Users/ronni/OneDrive/Documents/lab_data/mtclass_eqtl/results/48_tissue')
        mtclass = pd.read_csv('48_tissue_mtclass_aggregate.txt.gz', sep='\t', header=0)
        
        tissue_list = [t for t in os.listdir(r'E:/lab_data/mtclass_eqtl/multi-tissue/GTEx_Analysis_v8_eQTL') if 'signif' in t]
        tissue_list = sorted([t.replace('.v8.signif_variant_gene_pairs.txt.gz','') for t in tissue_list])
        tissue_list.remove('Kidney_Cortex')
        title = '48-tissue'

    # Make sure sample size is sufficient
    # print(f'Using sample size threshold of >={sample_size_threshold}')
    # mtclass = mtclass[mtclass['sample_size'] >= sample_size_threshold]
        
    results_dict = {k:[] for k in tissue_list}
    results_dict["MTClass"] = []
    results_dict["Top variants"] = []
    
    n_top_list = [100, 500, 1000, 2000, 3000, 4000, 5000]
    
    for n_top in n_top_list:
        
        mtclass_count = gwas_hits(mtclass, gwas_table, 'f1_macro_median', n_top, bp=10000)
        results_dict['Top variants'].append(n_top)
        results_dict['MTClass'].append(mtclass_count)
        
        for tissue in tissue_list:
            
            tissue_eqtl = pd.read_csv(rf'E:/lab_data/mtclass_eqtl/multi-tissue/GTEx_Analysis_v8_eQTL/{tissue}.v8.signif_variant_gene_pairs.txt.gz', sep='\t', header=0)
            tissue_eqtl = tissue_eqtl.loc[:,['gene_id','variant_id','pval_nominal']]
            tissue_eqtl.columns = ['gene','variant','pval']
            
            tissue_count = gwas_hits(tissue_eqtl, gwas_table, 'pval', n_top, bp=10000)
            results_dict[tissue].append(tissue_count)
        
        print(pd.DataFrame(results_dict).tail(1))
    
    results = pd.DataFrame(results_dict)
    
    if experiment == '48_tissue':
        avg_counts = np.mean(results.drop(['Top variants','MTClass'], axis=1), axis=0)
        avg_counts = avg_counts.sort_values(ascending=False)
        top_tissues = avg_counts.index.tolist()[:10]
        results_filter = results.filter(['Top variants','MTClass'] + top_tissues)
    
    # Plot results
    fig, ax = plt.subplots(figsize=(7,7))
    
    if experiment == '48_tissue':
        ax.plot(results_filter['Top variants'], results_filter['MTClass'], color='blue', linestyle='solid', label='MTClass')
        for tissue in top_tissues:
            ax.plot(results_filter['Top variants'], results_filter[tissue], linestyle='solid', label=tissue)
    else:
        ax.plot(results['Top variants'], results['MTClass'], color='blue', linestyle='solid', label='MTClass')
        for tissue in tissue_list:
            ax.plot(results['Top variants'], results[tissue], linestyle='solid', label=tissue)
    
    plt.title(f'GWAS colocalization\nMedian Macro F1\n{title}', fontsize=14)
    plt.ylabel('GWAS hits', fontsize=14)
    plt.xlabel('Top variants', fontsize=14)
    plt.legend(fontsize=10, loc='upper left')
    plt.tight_layout()
    plt.grid(alpha=0.3)
    plt.savefig(rf'C:/Users/ronni/Desktop/gwas_hits_{experiment}.png', dpi=400)
    plt.show()
    
    return results
    
def gwas_hits_random(experiment):
    
    gwas_table = load_gwas_table(data_dir=r'C:/Users/ronni/OneDrive/Documents/lab_data/mtclass_eqtl')
    
    if experiment == "9_tissue":
        
        sample_size_threshold = 50
        os.chdir(r'C:/Users/ronni/OneDrive/Documents/lab_data/mtclass_eqtl/results/9_tissue')
        mtclass = pd.read_csv('9_tissue_mtclass_aggregate.txt.gz', sep='\t', header=0)
        title = '9-tissue'
        
    elif experiment == "brain":
    
        sample_size_threshold = 150
        os.chdir(r'C:/Users/ronni/OneDrive/Documents/lab_data/mtclass_eqtl/results/brain_tissue')
        mtclass = pd.read_csv('brain_mtclass_aggregate.txt.gz', sep='\t', header=0)
        title = "Brain tissue"
        
    elif experiment == "48_tissue":
    
        sample_size_threshold = 300
        os.chdir(r'C:/Users/ronni/OneDrive/Documents/lab_data/mtclass_eqtl/results/48_tissue')
        mtclass = pd.read_csv('48_tissue_mtclass_aggregate.txt.gz', sep='\t', header=0)
        title = "48-tissue"
    
    # Make sure sample size is sufficient
    print(f'Using sample size threshold of >={sample_size_threshold}')
    mtclass = mtclass[mtclass['sample_size'] >= sample_size_threshold]
    
    n_top_list = [100, 500, 1000, 2000, 3000, 4000, 5000]
    results_dict = {'Top variants':[], 'MTClass':[], 'Random 1':[], 'Random 2':[], 'Random 3':[]}
    
    for n_top in n_top_list:
        
        mtclass_count = gwas_hits(mtclass, gwas_table, 'f1_macro_median', n_top, bp=10000)
        results_dict['Top variants'].append(n_top)
        results_dict['MTClass'].append(mtclass_count)
        
        for seed in [101, 102, 103]:
            
            df_sub = mtclass.sample(n=n_top, replace=False, random_state=seed)
            
            df_sub['seqnames'] = df_sub['variant'].apply(lambda x: "chr"+ str(x.split('_')[0][3:]))
            df_sub['starts'] = df_sub['variant'].apply(lambda x: x.split('_')[1]).astype(int)
            df_sub['ends'] = df_sub['variant'].apply(lambda x: x.split('_')[1]).astype(int)
            
            # Search upstream and downstream 10kb
            df_sub['starts'] = df_sub['starts'].subtract(10000)
            df_sub['ends'] = df_sub['ends'].add(10000)
            
            # count number of GWAS hits to include
            df_sub = df_sub.loc[:,['seqnames','starts','ends']]
            count = 0
            
            for chrom in set(df_sub['seqnames']):
                df_sub_chrom = df_sub[df_sub['seqnames']==chrom]
                gwas_chrom = gwas_table[gwas_table['seqnames']==chrom]
                for _, row in df_sub_chrom.iterrows():
                    interval = range(row['starts'], row['ends']+1)
                    for gwas_pos in gwas_chrom['starts']:
                        if gwas_pos in interval:
                            count += 1
                            
            results_dict[f'Random {seed-100}'].append(count)
        
        print(pd.DataFrame(results_dict).tail(1))
    
    results = pd.DataFrame(results_dict)
    
    # Plot results
    fig, ax = plt.subplots(figsize=(7,7))
    ax.plot(results['Top variants'], results['MTClass'], color='blue', linestyle='solid', label='MTClass')
    for col in ['Random 1','Random 2','Random 3']:
        ax.plot(results['Top variants'], results[col], linestyle='solid', label=col)
    
    plt.title(f'GWAS colocalization\nMedian Macro F1\n{title}', fontsize=14)
    plt.ylabel('GWAS hits', fontsize=14)
    plt.xlabel('Top variants', fontsize=14)
    plt.legend(fontsize=10)
    plt.tight_layout()
    plt.grid(alpha=0.3)
    plt.savefig(rf'C:/Users/ronni/Desktop/gwas_hits_random_{experiment}.png', dpi=400)
    plt.show()
    
    return results