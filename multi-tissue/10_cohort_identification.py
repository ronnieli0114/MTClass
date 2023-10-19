# -*- coding: utf-8 -*-
"""
This script runs MTClass on the SNPs
that were classified perfectly (median macro F1 and MCC = 1.0)
in the 9-tissue case to see whether we can ID
the entire GTEx cohort based solely on the predicted genotypes
from MTClass.

Created on Fri May  5 09:04:18 2023
@author: ronnieli
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import random
import scipy.stats as stats
import pyarrow.feather as feather
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, VotingClassifier
from sklearn.metrics import accuracy_score, f1_score

def write_perfect_snps(results_dir=r'G:/My Drive/Lab/lab_projects'):
    
    os.chdir(os.path.join(results_dir,'gtex_eqtl/results'))
    mtclass = pd.read_csv('9_tissue/9_tissue_ensemble_aggregate.txt.gz', sep='\t', header=0)
    # mtclass_perfect = mtclass[(mtclass['f1_macro_median']==1) & (mtclass['mcc_median']==1)]
    mtclass_perfect = mtclass[(mtclass['f1_macro_median'] >= 0.80) & (mtclass['mcc_median'] >= 0.80)]
    mtclass_perfect = sorted(list(set(mtclass_perfect['variant'])))

    with open(r'C:/Users/ronni/OneDrive/Desktop/perfect_snps.txt','w') as file:
        for snp in mtclass_perfect:
            file.write(snp)
            file.write('\n')
            
def write_top_snps(results_dir=r'G:/My Drive/Lab/lab_projects'):
    
    n_top = 100
    print(f'Writing {n_top} top SNPs to Desktop...')
    
    os.chdir(os.path.join(results_dir,'gtex_eqtl/results/9_tissue'))
    mtclass = pd.read_csv('9_tissue_ensemble_aggregate.txt.gz', sep='\t', header=0)
    
    top_snps = []
    for snp in mtclass.sort_values(by='f1_macro_median', ascending=False)['variant']:
        if len(set(top_snps)) == n_top:
            break
        else:
            top_snps.append(snp)
    top_snps = sorted(list(set(top_snps)))
    
    with open(rf'C:/Users/ronni/OneDrive/Desktop/top_{n_top}_snps.txt','w') as file:
        for snp in top_snps:
            file.write(snp)
            file.write('\n')
    
def load_rename_vcf(data_dir=r'E:/lab_data'):
    
    vcf = pd.read_csv(os.path.join(data_dir,'gtex_eqtl','f1_mcc_80.vcf'), sep='\t', header=0, skiprows=28)
    vcf = vcf.set_index('ID')
    vcf = vcf.drop(['POS','REF','ALT','QUAL','FILTER','INFO','FORMAT'], axis=1)

    vcf.columns = list(map(lambda x: x.split('_')[0], vcf.columns))
    vcf = vcf.replace({'0/0':0, '0/1':1, '1/1':1, './.':9})
    return vcf

def load_top_pairs(results_dir=r'G:/My Drive/Lab/lab_projects'):
    
    mtclass = pd.read_csv(os.path.join(results_dir,'gtex_eqtl/results/9_tissue','9_tissue_ensemble_aggregate.txt.gz'), sep='\t', header=0)
    mtclass_perfect = mtclass[(mtclass['f1_macro_median'] >= 0.80) & (mtclass['mcc_median'] >= 0.80)]
    return mtclass_perfect.loc[:,['gene','variant','f1_macro_median']]

def ID_donors(n_top_snps, random_seed=2023):
    
    random.seed(random_seed)
    
    # perfect_pairs = load_perfect_pairs(results_dir=r'G:/My Drive/Lab/lab_projects')
    top_pairs = load_top_pairs(results_dir=r'G:/My Drive/Lab/lab_projects')
    vcf = load_rename_vcf(data_dir=r'E:/lab_data')
    exp = feather.read_feather(r'E:/lab_data/gtex_eqtl/9_tissue/counts.ftr')
    print('Loaded results, VCF, and expression levels')

    all_preds = {'gene':[], 'variant':[], 'prediction':[]}
 
    true_donors = []
    
    vcf_unique = vcf[~vcf.drop('#CHROM', axis=1).duplicated()]
    top_pairs_unique = top_pairs[top_pairs['variant'].isin(vcf_unique.index)]
    top_pairs_unique = top_pairs_unique.drop_duplicates('variant')

    # Select top N SNPs based on median macro F1
    top_pairs_sort = top_pairs_unique.sort_values('f1_macro_median', ascending=False).reset_index(drop=True)
    # top_pairs_sort = top_pairs_sort.drop_duplicates('gene') # keep top SNP from each gene
    top_n_pairs = top_pairs_sort.iloc[:n_top_snps,:]

    print(f'Now classifying {n_top_snps} gene-SNP pairs')
    for i, row in top_n_pairs.reset_index(drop=True).iterrows():
        
        gene_name = row['gene']
        snp = row['variant']

        gts = vcf[vcf.index==snp].drop('#CHROM', axis=1).T
        gts.columns = ['genotype']
        
        exp_gene = exp[exp['gene']==gene_name]
        exp_gene = exp_gene.set_index('donor').drop('gene', axis=1)
            
        Xy = exp_gene.merge(gts, left_index=True, right_index=True)
        Xy['genotype'] = Xy['genotype'].replace({9:np.nan})
        Xy = Xy.dropna(axis=0, how='any')
        
        X = Xy.iloc[:,:-1]
        y = Xy.iloc[:,-1]
        
        rf = RandomForestClassifier(n_estimators=150, max_depth=5, random_state=random_seed)
        svc = SVC(C=1, kernel='rbf', probability=True, random_state=random_seed)
        clf = VotingClassifier(estimators=[('RF', rf), ('SVM', svc)], voting='soft')
        skf = StratifiedKFold(n_splits=3, shuffle=False)
        
        pred_genotypes = []
        
        for train_index, test_index in skf.split(X, y):
            
            X_train, y_train = X.iloc[train_index,:], y.iloc[train_index]
            X_test, y_test = X.iloc[test_index,:], y.iloc[test_index]
            
            clf.fit(X_train, y_train)
            y_pred = clf.predict(X_test)
            pred_genotypes.extend(list(y_pred))
        
        true_donors.extend(list(Xy.index))
        all_preds['gene'].append(gene_name)
        all_preds['variant'].append(snp)
        all_preds['prediction'].append(pd.Series(pred_genotypes))
        
    true_donors = sorted(list(set(true_donors)))
    num_snps = len(set(all_preds['variant']))
    print(f'Using {num_snps} SNPs to identify {len(true_donors)} donors')
    
    # Make dataframe of all predictions
    preds_df = pd.concat(all_preds['prediction'], axis=1).T
    preds_df.insert(loc=0, column='gene', value=all_preds['gene'])
    preds_df.insert(loc=1, column='variant', value=all_preds['variant'])
    preds_df = preds_df.sort_values(by='variant', ascending=True)
    
    # Use genotypes where there are 103 donors
    all_gts = vcf[vcf.index.isin(all_preds['variant'])]
    all_gts = all_gts.drop('#CHROM', axis=1)
    all_gts = all_gts.sort_index()
    
    assert np.sum(preds_df['variant']==all_gts.index) == num_snps
    
    # for each column of genotypes, check if you can uniquely identify individual
    pred_donors = []
    num_matching_donors = []
    
    for col in preds_df.columns[2:]: # first two columns are gene, variant
        
        pred_gt = preds_df[col].values
        
        # fill in missing values with most common genotype
        pred_gt_no_na = pred_gt[~np.isnan(pred_gt)].astype(int)
        most_common_genotype = np.argmax(np.bincount(pred_gt_no_na))
        pred_gt = np.nan_to_num(pred_gt, nan=most_common_genotype)
        
        overlaps = []
        
        for donor in all_gts.columns:
            
            true_gt = all_gts[donor].values
            overlap, _ = stats.pearsonr(pred_gt, true_gt)
            overlaps.append(overlap)
        
        # Select the donor with highest correlation coefficient that hasn't been called
        donor_inds = np.argsort(overlaps)[::-1]
        for idx in donor_inds:
            pred_donor = np.array(all_gts.columns)[idx]
            if pred_donor not in pred_donors:
                break
        pred_donors.append(pred_donor)
        
        # pred_donor_inds = [i for i, overlap in enumerate(overlaps) if overlap==np.max(overlaps)]
        # all_matching_donors = np.array(all_gts.columns)[pred_donor_inds]
        # num_matching_donors.append(len(all_matching_donors))
        
        # # Require that random called donor has not already been called
        # matching_donors_not_called = [donor for donor in all_matching_donors if donor not in pred_donors]
        # if len(matching_donors_not_called) == 0:
        #     all_donors_not_called = set(all_gts.columns) - set(pred_donors)
        #     random_pred_donor = random.choice(list(all_donors_not_called))
        # else:
        #     random_pred_donor = random.choice(matching_donors_not_called) # randomly select donor from donors with max overlap       
        # pred_donors.append(random_pred_donor)

    # num_matching_donors = np.array(num_matching_donors)
    # num_ones = np.count_nonzero(num_matching_donors==1)

    donors_overlap = set(pred_donors) & set(true_donors)
    num_correct = len(donors_overlap)
    accuracy = np.around(num_correct / len(true_donors), 3)
    # print(f'Number of donors with only 1 match: {num_ones}')
    print(f'Number of donors correct: {num_correct} \t Accuracy: {accuracy}')
    
    return accuracy

#%%
results_dict = {'num_snps':[], 'accuracy':[]}
for n in range(5,501,5):
    acc = ID_donors(n, random_seed=2023)
    results_dict['num_snps'].append(n)
    results_dict['accuracy'].append(acc)
results = pd.DataFrame(results_dict)

fig, ax = plt.subplots(figsize=(5,5))
ax.plot(results['num_snps'], results['accuracy'])
ax.set_ylabel('Accuracy', fontsize=16)
ax.set_xlabel('Top SNPs used', fontsize=16)
ax.set_title('Cohort identification\nMulti-tissue study', fontsize=20)
plt.tight_layout()
plt.savefig(r'C:/Users/ronni/OneDrive/Desktop/multi-tissue_cohort_ID.png', dpi=400)

# max_acc_ind = np.argmax(results['accuracy'])
# best_num_snps = results.iloc[max_acc_ind,:]['num_snps'].astype(int)
# print(f'Best number of SNPs to use: {best_num_snps}')

# seed_results_dict = {'seed':[], 'num_correct':[], 'num_ones':[], 'accuracy':[]}
# for seed in range(0, 201, 20):
#     num_correct, num_ones, acc = ID_donors(best_num_snps, random_seed=seed)
#     seed_results_dict['seed'].append(seed)
#     seed_results_dict['num_correct'].append(num_correct)
#     seed_results_dict['num_ones'].append(num_ones)
#     seed_results_dict['accuracy'].append(acc)
# seed_results = pd.DataFrame(seed_results_dict)

# fig, ax = plt.subplots(figsize=(5,5))
# ax.plot(seed_results['seed'], seed_results['accuracy'])
# ax.set_xlabel('Seed', fontsize=16)
# ax.set_ylabel('Accuracy')
# ax.set_title(f'Cohort identification\nTop {best_num_snps} SNPs', fontsize=20)
# plt.tight_layout()
# plt.savefig(r'C:/Users/ronni/OneDrive/Desktop/cohort_ID_seeds.png', dpi=400)

