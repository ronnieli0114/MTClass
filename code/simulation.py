# -*- coding: utf-8 -*-
"""
Simulation study:
See whether MTClass and MultiPhen differ
in their classification of highly nonlinearly
separable datasets

Created on Fri Jul 29 05:35:41 2022
@author: ronnieli
"""

import pandas as pd
import numpy as np
from math import pi, cos, sin
import os
import sys
import random
import matplotlib.pyplot as plt
import scipy.stats as stats
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, VotingClassifier
from sklearn.metrics import f1_score, roc_auc_score, average_precision_score, matthews_corrcoef, roc_curve
from sklearn.model_selection import StratifiedKFold
sys.path.append(r'G:/My Drive/Lab/lab_projects/gtex_eqtl/scripts')
from expression_level_plots import plot_boxplot, plot_heatmap, plot_tsne, plot_mean_variance

random.seed(2022)
n_samples = 100

#%%
'''
First example of nonlinear data - 
two spheres with the same center (mean)
but different radii (variances)
'''

def point(x0, y0, z0, r):
    theta = random.random() * 2 * pi
    phi = random.random() * pi
    x  = r * cos(theta) * sin(phi) + x0
    y = r * sin(theta) * sin(phi) + y0
    z = r * cos(phi) + z0
    return x, y, z 

xyz0 = [point(20, 20, 20, 20) for _ in range(n_samples)]
xyz1 = [point(20, 20, 20, 5) for _ in range(n_samples)]

x0, y0, z0 = zip(*xyz0)
x1, y1, z1 = zip(*xyz1)

ax = plt.axes(projection='3d')
ax.scatter3D(x0, y0, z0, color='blue', alpha=0.5, label='Genotype AA')
ax.scatter3D(x1, y1, z1, color='orange', alpha=0.5, label='Genotypes Aa/aa')
plt.title('Actual simulated data (3 tissues)', fontsize=20)
plt.tight_layout()
plt.legend()
plt.savefig(r"C:/Users/ronnieli/Desktop/simulated_data.png", dpi=400)
plt.show()

### TSNE on simulated datas
X0 = np.concatenate((x0, y0, z0)).reshape((n_samples, 3))
X1 = np.concatenate((x1, y1, z1)).reshape((n_samples, 3))

y0 = [0] * n_samples
y1 = [1] * n_samples

X = np.vstack((X0, X1))
num_zero = X0.shape[0]
embed = TSNE(n_components=2, perplexity=40, random_state=2022, learning_rate='auto').fit_transform(X)

plt.figure(figsize=(5,5))
plt.scatter(embed[:num_zero,0], embed[:num_zero,1], c='blue', alpha=0.5, label='Genotype AA')
plt.scatter(embed[num_zero:,0], embed[num_zero:,1], c='orange', alpha=0.5, label='Genotypes Aa/aa')
plt.xlabel('Component 1', fontsize=16)
plt.ylabel('Component 2', fontsize=16)
plt.title('TSNE projection of simulated data', fontsize=16)
plt.legend()
plt.tight_layout()
plt.savefig(r"C:/Users/ronnieli/Desktop/TSNE.png", dpi=400)
plt.show()

### CLASSIFY
svc = SVC(kernel='rbf', C=1, probability=True, random_state=2022)
rf = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=2022)

skf = StratifiedKFold(n_splits=5, random_state=2022, shuffle=True)
scaler = StandardScaler()

X = np.vstack((X0, X1))
y = np.concatenate((y0, y1)).reshape(-1)
print(X.shape, y.shape)

### save data for MultiPhen
pheno_data = pd.DataFrame(X)
pheno_data.columns = ['tissue_1','tissue_2','tissue_3']
pheno_data.to_csv(r"C:/Users/ronnieli/Desktop/pheno_data.csv", index=False)

geno_data = pd.DataFrame(y)
geno_data.columns = ['variant_1']
geno_data.to_csv(r"C:/Users/ronnieli/Desktop/geno_data.csv", index=False)

results_dict = {'model':[], 'auc':[], 'f1':[], 'f1_macro':[], 'f1_weighted':[], 'auprc':[], 'mcc':[]}

for train_ind, test_ind in skf.split(X, y):
    
    X_train, y_train = X[train_ind], y[train_ind]
    X_test, y_test = X[test_ind], y[test_ind]
    
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)
    
    f1_lst = []
    f1_macro_lst = []
    f1_weighted_lst = []
    auc_lst = []
    auprc_lst = []
    mcc_lst = []
    
    model = VotingClassifier(estimators=[('SVM', svc), ('RF', rf)], voting='soft')
    model.fit(X_train, y_train)
    
    y_pred = model.predict(X_test)
    y_proba = model.predict_proba(X_test)
    f1_lst.append(f1_score(y_test, y_pred))
    f1_macro_lst.append(f1_score(y_test, y_pred, average='macro'))
    f1_weighted_lst.append(f1_score(y_test, y_pred, average='weighted'))
    auc_lst.append(roc_auc_score(y_test, y_proba[:,1]))
    auprc_lst.append(average_precision_score(y_test, y_proba[:,1]))
    mcc_lst.append(matthews_corrcoef(y_test, y_pred))
    
results_dict['model'].append('SVM')
results_dict['f1'].append(round(np.mean(f1_lst), 3))
results_dict['f1_macro'].append(round(np.mean(f1_macro_lst), 3))
results_dict['f1_weighted'].append(round(np.mean(f1_weighted_lst), 3))
results_dict['mcc'].append(round(np.mean(mcc_lst), 3))
results_dict['auc'].append(round(np.mean(auc_lst),3))
results_dict['auprc'].append(round(np.mean(auprc_lst),3))
results_df = pd.DataFrame(results_dict)
print(results_df)

preds = y_proba[:,1]
fpr, tpr, threshold = roc_curve(y_test, preds)
roc_auc = np.round(np.mean(auc_lst), 3)

plt.figure(figsize=(5,5))
plt.plot(fpr, tpr, 'b', label=f'AUC={roc_auc:.2f}')
plt.legend()
plt.title('ROC curve', fontsize=20)
plt.plot([0,1],'r--')
plt.xlim([-0.05,1.05])
plt.ylim([-0.05,1.05])
plt.ylabel('True Positive Rate', fontsize=16)
plt.xlabel('False Positive Rate', fontsize=16)
plt.tight_layout()
plt.savefig(r'C:/Users/ronnieli/Desktop/roc_curve.png', dpi=400)
plt.show()


#%% Conduct Levene test and independent-samples t-test for each tissue
all_data = pheno_data.merge(geno_data, left_index=True, right_index=True)
stat_dict = {'tissue':[], 'p_ttest':[], 'p_levene':[]}

for col in pheno_data.columns:
    
    counts_0 = all_data[all_data['variant_1']==0][col]
    counts_1 = all_data[all_data['variant_1']==1][col]
    
    # counts_0 = np.log(counts_0 + 1)
    # counts_1 = np.log(counts_1 + 1)
    
    stat_ttest, p_ttest = stats.kruskal(counts_0, counts_1)
    stat_levene, p_levene = stats.levene(counts_0, counts_1)
    
    stat_dict['tissue'].append(col)
    stat_dict['p_ttest'].append(round(p_ttest, 5))
    stat_dict['p_levene'].append(round(p_levene, 5))

stat_df = pd.DataFrame(stat_dict)

#%% Plot means and variances for tissues
fig, ax = plt.subplots(1, 1, figsize=(7,7))

counts_0 = all_data[all_data['variant_1']==0]
counts_1 = all_data[all_data['variant_1']==1]

title_list = ['Genotype REF/REF', 'Genotypes REF/ALT and ALT/ALT']
counts_list = [counts_0, counts_1]
shifts = [0, 0.1]

for title, counts, shift in zip(title_list, counts_list, shifts):
    
    counts = counts.iloc[:,:-1]
    
    counts.loc['mean'] = np.mean(counts, axis=0)
    counts.loc['std'] = np.std(counts, axis=0)
    
    x = np.arange(len(counts.columns)) + 1 + shift
    y = counts.loc['mean']
    yerr = counts.loc['std']
    
    ax.errorbar(x=x, y=y, yerr=yerr, linestyle='', marker='o', label=title)
    ax.set_xticks(x, counts.columns, rotation=90)
    
    ax.set_ylim((0, 35))
    
    f1 = results_df['f1_macro'].values[0]
    mcc = results_df['mcc'].values[0]
    pval = 0.540
    
    x_coord = -50
    y_coord = -100
    
    ax.annotate(f'F1={f1}', xy=(x_coord, y_coord), xycoords='axes points')
    ax.annotate(f'MCC={mcc}', xy=(x_coord, y_coord-14), xycoords='axes points')
    ax.annotate(f'MultiPhen p-val={pval}', xy=(x_coord, y_coord-28), xycoords='axes points')
    ax.legend()

# annotate with p-values for each tissue
x_coord = [1, 1.75, 2.5]
y_coord = [3, 3, 3]
cols = ['tissue_1', 'tissue_2', 'tissue_3']
for x, y, col in zip(x_coord, y_coord, cols):
    
    p_ttest = stat_df[stat_df['tissue']==col]['p_ttest'].values[0]
    p_levene = stat_df[stat_df['tissue']==col]['p_levene'].values[0]
    
    p_ttest = round(p_ttest, 3)
    p_levene = round(p_levene, 3)
    
    ax.text(x, y, f'p (KW test) = {p_ttest}\n p (Levene) = {p_levene}')
    
plt.suptitle('Mean and variance in expression levels per tissue\n Simulated data', fontsize=20)
plt.tight_layout()
save_dir = r'C:/Users/ronnieli/Desktop'
if not os.path.exists(save_dir):
    os.mkdir(save_dir)
plt.savefig(os.path.join(save_dir, 'simulated_results.png'), dpi=400)
plt.show()
    

#%%
'''
Second example of nonlinear data - 
two sine waves intertwined with each other
horizontally shifted by pi
'''

def create_sin(x0, y0, z0, pos=True):
    np.random.seed(2022)
    x = np.random.uniform(2*np.pi, 4*np.pi, n_samples)
    # y = np.sin(x - x0) + y0
    y = np.random.uniform(2*np.pi, 4*np.pi, n_samples)
    if pos:
        z = np.sin(x - x0) + z0
    else:
        z = -np.sin(x - x0) + z0
    noise = np.random.normal(0, 0.05, n_samples)
    return x, y, z + noise

x0, y0, z0 = create_sin(0, 10, 10, pos=True)
x1, y1, z1 = create_sin(np.pi, 10, 10, pos=True)

fig = plt.figure(figsize=(6,6))
ax = plt.axes(projection='3d')
ax.scatter3D(x0, y0, z0, c='blue', alpha=0.5, label='Genotype AA')
ax.scatter3D(x1, y1, z1, c='red', alpha=0.5, label='Genotypes AA/aa')
ax.set_ylabel('y', fontsize=16)
ax.set_xlabel('x', fontsize=16)
ax.set_zlabel('z', fontsize=16)
plt.title('Actual simulated data', fontsize=20)
plt.tight_layout()
plt.legend()
plt.savefig(r"C:/Users/ronnieli/Desktop/data.png", dpi=400)
plt.show()


# labels
label_0 = np.array([0]*n_samples)
label_1 = np.array([1]*n_samples)
y = np.concatenate([label_0, label_1])

# aggregate data
X0 = np.vstack([x0, y0, z0]).transpose()
X1 = np.vstack([x1, y1, z1]).transpose()
X = np.concatenate([X0, X1], axis=0)
num_0 = X0.shape[0]

### TSNE
embed = TSNE(n_components=2, perplexity=40, random_state=2022, learning_rate='auto').fit_transform(X)

plt.figure(figsize=(6,6))
plt.scatter(embed[:num_zero,0], embed[:num_zero,1], c='blue', alpha=0.5, label='Genotype AA')
plt.scatter(embed[num_zero:,0], embed[num_zero:,1], c='red', alpha=0.5, label='Genotypes Aa/aa')
plt.xlabel('Component 1', fontsize=16)
plt.ylabel('Component 2', fontsize=16)
plt.title('TSNE projection of simulated data', fontsize=16)
plt.legend()
plt.tight_layout()
plt.savefig(r"C:/Users/ronnieli/Desktop/TSNE.png", dpi=400)
plt.show()

### CLASSIFY
# svc = SVC(kernel='rbf', C=1, probability=True, random_state=2022)
skf = StratifiedKFold(n_splits=5, random_state=2022, shuffle=True)
scaler = StandardScaler()

print(X.shape, y.shape)

pheno_data = pd.DataFrame(X)
pheno_data.columns = ['tissue1','tissue2','tissue3']
pheno_data.to_csv(r"C:/Users/ronnieli/Desktop/pheno_data.csv", index=False)

geno_data = pd.DataFrame(y)
geno_data.columns = ['variant1']
geno_data.to_csv(r"C:/Users/ronnieli/Desktop/geno_data.csv", index=False)

results_dict = {'model':[], 'auc':[], 'f1':[], 'f1_macro':[], 'f1_weighted':[], 'auprc':[], 'mcc':[]}

for train_ind, test_ind in skf.split(X, y):
    
    X_train, y_train = X[train_ind], y[train_ind]
    X_test, y_test = X[test_ind], y[test_ind]
    
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)
    
    f1_lst = []
    f1_macro_lst = []
    f1_weighted_lst = []
    auc_lst = []
    auprc_lst = []
    mcc_lst = []
    
    # model = VotingClassifier(estimators=[('SVM', svc), ('RF', rf)], voting='soft')
    model = rf
    model.fit(X_train, y_train)
    
    y_pred = model.predict(X_test)
    y_proba = model.predict_proba(X_test)
    f1_lst.append(f1_score(y_test, y_pred))
    f1_macro_lst.append(f1_score(y_test, y_pred, average='macro'))
    f1_weighted_lst.append(f1_score(y_test, y_pred, average='weighted'))
    auc_lst.append(roc_auc_score(y_test, y_proba[:,1]))
    auprc_lst.append(average_precision_score(y_test, y_proba[:,1]))
    mcc_lst.append(matthews_corrcoef(y_test, y_pred))
    
results_dict['model'].append('SVC')
results_dict['f1'].append(round(np.mean(f1_lst), 3))
results_dict['f1_macro'].append(round(np.mean(f1_macro_lst), 3))
results_dict['f1_weighted'].append(round(np.mean(f1_weighted_lst), 3))
results_dict['mcc'].append(round(np.mean(mcc_lst), 3))
results_dict['auc'].append(round(np.mean(auc_lst),3))
results_dict['auprc'].append(round(np.mean(auprc_lst),3))
results_df = pd.DataFrame(results_dict)
print(results_df)

preds = y_proba[:,1]
fpr, tpr, threshold = roc_curve(y_test, preds)
roc_auc = np.round(np.mean(auc_lst), 3)

plt.figure(figsize=(5,5))
plt.plot(fpr, tpr, 'b', label=f'AUC={roc_auc:.2f}')
plt.legend()
plt.title('ROC curve', fontsize=20)
plt.plot([0,1],'r--')
plt.xlim([-0.05,1.05])
plt.ylim([-0.05,1.05])
plt.ylabel('True Positive Rate', fontsize=16)
plt.xlabel('False Positive Rate', fontsize=16)
plt.tight_layout()
plt.savefig(r'C:/Users/ronnieli/Desktop/roc_curve.png', dpi=400)
plt.show()
