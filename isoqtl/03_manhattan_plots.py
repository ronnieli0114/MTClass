# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 22:10:17 2023

@author: ronnieli
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
import pyarrow.feather as feather

results_dir = r'G:/My Drive/Lab/lab_projects/isoqtl/results'
data_dir = r'E:/lab_data/'

def manhattan_plot(results_dir):
    
    '''
    Given a tissue case (9_tissue or brain), plots Manhattan plots
    of aggregated results from 3 iterations of MTClass using both
    median macro F1 and median MCC metrics
    '''
    
    results = pd.read_csv(os.path.join(results_dir, 'mtclass_aggregate123.txt.gz'), sep='\t', header=0)

    results_f1 = results[['gene','variant','f1_macro_median']]
    results_mcc = results[['gene','variant','mcc_median']]
    
    results_f1['CHR'] = results_f1['variant'].apply(lambda x: x.split(':')[0])
    results_f1['POS'] = results_f1['variant'].apply(lambda x: x.split(':')[1])
    results_f1 = results_f1.loc[:,['variant','CHR','POS','f1_macro_median']]
    
    results_mcc['CHR'] = results_mcc['variant'].apply(lambda x: x.split(':')[0])
    results_mcc['POS'] = results_mcc['variant'].apply(lambda x: x.split(':')[1])
    results_mcc = results_mcc.loc[:,['variant','CHR','POS','mcc_median']]

    ### Manhattan plot of F1 values

    results_f1['CHR']= results_f1['CHR'].astype(int)
    results_f1 = results_f1.sort_values('CHR')

    results_f1['ind'] = range(len(results_f1))
    results_grouped = results_f1.groupby(('CHR'))

    fig = plt.figure(figsize=(14, 8)) # Set the figure size
    ax = fig.add_subplot(111)
    colors = ['blue','orange']
    x_labels = []
    x_labels_pos = []
    for num, (name, group) in enumerate(results_grouped):
        group.plot(kind='scatter', x='ind', y='f1_macro_median', color=colors[num % len(colors)], ax=ax, alpha=0.7)
        x_labels.append(name)
        x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)

    # set axis limits
    ax.set_ylim((0, 1.1))

    # x axis label
    ax.set_xlabel('Chromosome', fontsize=16)
    ax.set_ylabel('Median Macro F1', fontsize=16)
    ax.set_title('isoQTL median Macro F1 results', fontsize=20)

    # show the graph
    plt.tight_layout()
    plt.savefig(r'C:/Users/ronnieli/Desktop/isoqtl_f1.png', dpi=400)
    plt.show()

    ### Manhattan plot of MCC values

    results_mcc['CHR']= results_mcc['CHR'].astype(int)
    results_mcc = results_mcc.sort_values('CHR')

    results_mcc['ind'] = range(len(results_mcc))
    results_grouped = results_mcc.groupby(('CHR'))

    # manhattan plot
    fig = plt.figure(figsize=(14, 8)) # Set the figure size
    ax = fig.add_subplot(111)
    colors = ['blue','orange']
    x_labels = []
    x_labels_pos = []
    for num, (name, group) in enumerate(results_grouped):
        group.plot(kind='scatter', x='ind', y='mcc_median', color=colors[num % len(colors)], ax=ax, alpha=0.7)
        x_labels.append(name)
        x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)

    # set axis limits
    ax.set_ylim([-1.2,1.2])

    # x axis label
    ax.set_xlabel('Chromosome', fontsize=16)
    ax.set_ylabel('Median MCC', fontsize=16)
    ax.set_title('isoQTL median MCC results', fontsize=20)

    # show the graph
    plt.tight_layout()
    plt.savefig(rf'C:/Users/ronnieli/Desktop/isoqtl_mcc.png', dpi=400)
    plt.show()
    
manhattan_plot(results_dir)