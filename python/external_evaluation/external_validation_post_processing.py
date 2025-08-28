#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 2023

@author: rintala
"""

import numpy as np
import pandas as pd
import os
import glob

import seaborn as sns

from sklearn import neighbors
from copy import copy
from sklearn.metrics import silhouette_score
from sklearn.metrics.pairwise import pairwise_distances
import importlib
import sys
import re

dsp_spec = importlib.util.spec_from_file_location('dataset_processing', 'datasets/cancer_dataset_processing.py')
dsp = importlib.util.module_from_spec(dsp_spec)
sys.modules['dataset_processing'] = dsp
dsp_spec.loader.exec_module(dsp)

import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

#%% Data
from sys import platform
from socket import gethostname
if gethostname() == 'teemu-pc':
    base_path = '/home/teemu/research_work/superAE_HPO/'
elif platform == 'linux':
    base_path = '/research/work/rintala/superAE_HPO/'
else:
    base_path = '//research/workdir/superAE_HPO/'

#res_path = base_path + '20230821_random_search/brca_test_noclfilter/'
#res_path = base_path + '20231027_random_search/brca_test_noclfilter/'
#res_path = base_path + '20231101_random_search/brca_test_noclfilter/'
#res_path = base_path + '20231103_random_search/brca_test_noclfilter/'
#res_path = base_path + 'random_search_231113/brca_test_noclfilter/'
#res_path = base_path + '231117_random_search/brca_test_noclfilter/'
#res_path = base_path + '20231120_random_search/brca_test_noclfilter/'
#res_path = base_path + '20231121_random_search/brca_test_noclfilter/'
#res_path = base_path + '20231122_random_search/scanb_test_noclfilter/'
#res_path = base_path + '20231124_random_search/pancan_test/'
#res_path = base_path + '20240103_random_search/brca_test_noclfilter/'
#res_path = base_path + '20240108_random_search/brca_test_noclfilter/'
#res_path = base_path + '20240115_random_search/brca_test_noclfilter/'
#res_path = base_path + '20240118_random_search/brca_test_noclfilter/'
#res_path = base_path + '20240125_random_search/brca_test_noclfilter/'
#res_path = base_path + '20240125_random_search/brca_test_noclfilter_nopre/'
#res_path = base_path + '20240206_random_search/brca_test_noclfilter/'
#res_path = base_path + '20240208_random_search/brca_test_noclfilter/'
#res_path = base_path + '20240211_random_search/brca_test_noclfilter/'
#res_path = base_path + '20240213_random_search/brca_test_noclfilter/'
#res_path = base_path + '20240213_random_search/brca_test_noclfilter_pretrain/'
#res_path = base_path + '20240216_random_search/brca_test_noclfilter_pretrain/'
#res_path = base_path + '20240221_random_search/scanb_test_noclfilter/'
#res_path = base_path + '20240301_random_search/scanb_test_noclfilter/'
#res_path = base_path + '20240302_random_search/scanb_test_noclfilter/'
#res_path = base_path + '20240303_random_search/scanb_test_noclfilter/'
#res_path = base_path + '20240328_random_search/scanb_test_noclfilter/'
#res_path = base_path + '20240430_random_search/scanb_test_noclfilter/'
#res_path = base_path + '20240605_random_search/scanb_test_noclfilter/'
#res_path = base_path + '20240610_random_search/scanb_test_noclfilter/'
#res_path = base_path + '20250124_random_search/pancan_test/'
#res_path = base_path + '20250202_random_search/pancan_test/'
#res_path = base_path + '20250205_random_search/scanb_test_noclfilter/'
#res_path = base_path + '20250207_random_search/pancan_test/'
#res_path = base_path + '20250217_random_search/pancan_test/'
#res_path = base_path + '20250219_random_search/pancan_test/'
#res_path = base_path + '20250223_random_search/pancan_test/'
#res_path = base_path + '20250313_random_search/pancan_test/'
#res_path = base_path + '20250402_random_search/pancan_test/'
#res_path = base_path + '20250404_random_search/pancan_test/'
#res_path = base_path + '20250408_random_search/pancan_test/'
res_path = base_path + '20250410_random_search/pancan_test/'
#res_path = base_path + '20250415_random_search/pancan_test/'
#res_path = base_path + '20250506_random_search/pancan_test/' 
#res_path = base_path + '20250520_random_search/pancan_test/'
#res_path = base_path + '20250525_random_search/pancan_test/' # not done
#res_path = base_path + '20250527_random_search/pancan_test/'
#res_path = base_path + '20250530_random_search/pancan_test/'
#res_path = base_path + '20250603_random_search/pancan_test/'

res_path = res_path + 'external_evaluation/'

embedding_files = glob.glob(res_path + '*validation_embeddings.csv.gz')

embeddings_list = []
for f in embedding_files:
    embeddings_list.append(pd.read_csv(f, header = 0))

for i in np.arange(len(embeddings_list)):
    embeddings_list[i].rename(columns = {'Unnamed: 0' : 'id'}, inplace = True)

# Fix dataset column
for i in np.arange(len(embeddings_list)):
    if embeddings_list[i].get('dataset', None) is None:
        if re.search('internal_survival', embedding_files[i]) is not None:
            embeddings_list[i]['dataset'] = 'tcga'
        if re.search('external_survival', embedding_files[i]) is not None:
            embeddings_list[i]['dataset'] = 'scanb'

embeddings_df = pd.concat(embeddings_list, axis = 0)
embeddings_df.value_counts('dataset')

'''
embeddings_dict = dict(zip(
    [embeddings_list[i]['dataset'][0] for i in np.arange(len(embeddings_list))],
    embeddings_list))
'''

prediction_files = glob.glob(res_path + '*validation_predictions.csv.gz')

prediction_list = []
for f in prediction_files:
    prediction_list.append(pd.read_csv(f, header = 0))

for i in np.arange(len(prediction_list)):
    prediction_list[i].rename(columns = {'Unnamed: 0' : 'id'}, inplace = True)

# Fix dataset column
import re
for i in np.arange(len(prediction_list)):
    if prediction_list[i].get('dataset', None) is None:
        if re.search('internal_survival', prediction_files[i]) is not None:
            prediction_list[i]['dataset'] = 'tcga'
        if re.search('external_survival', prediction_files[i]) is not None:
            prediction_list[i]['dataset'] = 'scanb'

'''
prediction_dict = dict(zip(
    [prediction_list[i]['task'][0] for i in np.arange(len(prediction_list))],
    prediction_list))
'''

predictions_df = pd.concat(prediction_list, axis = 0)
predictions_df.value_counts('dataset')

#%% TCGA and CCLE data
if gethostname() == 'teemu-pc':
    patient_expression_root_dir = '/home/teemu/research_work/tcga/pan_cancer/data_full/'
    cell_line_expression_root_dir = '/home/teemu/research_work/ccle/'
    cell_line_drug_response_root_dir = '/home/teemu/research_work/ctrp/'
    bruna_path = '/home/teemu/research_work/breast_pdtc_bruna/'
    drugbank_path = '/home/teemu/research_fortino_group/Common_databases/DrugBank/'
    gao_path = "~/research_work/pdtc_gao/"
    home = True
elif platform == 'linux':
    patient_expression_root_dir = '/research/work/rintala/tcga/pan_cancer/data_full/'
    cell_line_expression_root_dir = '/research/users/rintala/ccle/'
    cell_line_drug_response_root_dir = '/research/work/rintala/ctrp/'
    bruna_path = '/research/work/rintala/breast_pdtc_bruna/'
    drugbank_path = '/research/groups/fortino/Common_databases/DrugBank/'
    gao_path = "/research/work/rintala/pdtc_gao/"
    home = False
else:
    patient_expression_root_dir = '//research/workdir/tcga/pan_cancer/data_full/'
    cell_line_expression_root_dir = '//research/rintala/ccle/'
    cell_line_drug_response_root_dir = '//research/workdir/rintala/ctrp/'
    bruna_path = '//research/workdir/breast_pdtc_bruna/'
    drugbank_path = '/research/groups/fortino/Common_databases/DrugBank/'
    gao_path = "//research/workdir/pdtc_gao/"
    home = False


#%% CTRP or CCLE?
prediction_cols = predictions_df.filter(regex = 'dr_pred_[0-9]+', axis = 1).shape[1]
if prediction_cols == 24:
    response_dataset = 'ccle'
    raise NotImplementedError('CCLE based validation is not implemented yet.')
elif prediction_cols == 545:
    response_dataset = 'ctrp'
    with open(f"{cell_line_drug_response_root_dir}drug_names.txt", 'r') as f:
        train_drugs = f.read().splitlines()
    train_drugs = map(str.lower, train_drugs)
    train_drugs = [re.sub('\\(.*?\\)', '', i) for i in train_drugs]
    train_drugs = [re.sub('\\s', '', i) for i in train_drugs]
    train_drugs = [re.sub(':', '+', i) for i in train_drugs]
    train_drugs = np.array(['tipifarnib' if re.search('tipifarnib', i) else i for i in train_drugs])
elif prediction_cols == 544:
    response_dataset = 'ctrp' # Xia data
    with open(f"{base_path}../drug_response_dataset/processed_drug_names.txt", 'r') as f:
        train_drugs = f.read().splitlines()
    train_drugs = map(str.lower, train_drugs)
    train_drugs = [re.sub('\\(.*?\\)', '', i) for i in train_drugs]
    train_drugs = [re.sub('\\s', '', i) for i in train_drugs]
    train_drugs = [re.sub(':', '+', i) for i in train_drugs]
    train_drugs = np.array(['azd-2281' if re.search('olaparib', i) else i for i in train_drugs])
    train_drugs = np.array(['5-fu' if re.search('fluorouracil', i) else i for i in train_drugs])
    train_drugs = np.array(['plx-4032' if re.search('vemurafenib', i) else i for i in train_drugs])
    train_drugs = np.array(['platin' if re.search('cisplatin', i) else i for i in train_drugs])
else:
    raise ValueError('Drug response prediction column number does not match known datasets.')


#%% Bruna drug response data

pdtc_drug_response = pd.read_csv(
    f"{bruna_path}DrugResponsesAUCModels.txt", 
    header = 0, sep = "\t")

pdtx_drug_response = pd.read_csv(
    f"{bruna_path}DrugResponsesAUCSamples.txt", 
    header = 0, sep = "\t")

#pdtx_drug_response.columns
pdtx_drug_response['drug'] = pdtx_drug_response['Drug'].map(str.lower)
pdtx_drug_response['drug'] = [re.sub('\\(.*?\\)', '', i) for i in pdtx_drug_response['drug']]
pdtx_drug_response['drug'] = [re.sub('\\s', '', i) for i in pdtx_drug_response['drug']]
pdtc_drug_response['drug'] = pdtc_drug_response['Drug'].map(str.lower)
pdtc_drug_response['drug'] = [re.sub('\\(.*?\\)', '', i) for i in pdtc_drug_response['drug']]
pdtc_drug_response['drug'] = [re.sub('\\s', '', i) for i in pdtc_drug_response['drug']]

if not np.isin(pdtc_drug_response['drug'].unique(), pdtx_drug_response['drug'].unique()).all():
    raise(ValueError())
if not np.isin(pdtx_drug_response['drug'].unique(), pdtc_drug_response['drug'].unique()).all():
    raise(ValueError())

#%% Gao response data

gao_drug_response = pd.read_excel(
    f"{gao_path}41591_2015_BFnm3954_MOESM10_ESM.xlsx", 
    sheet_name = "PCT curve metrics"
)
gao_drug_response['drug'] = gao_drug_response['Treatment'].map(str.lower)
gao_drug_response['drug'] = [re.sub('\\s', '', i) for i in gao_drug_response['drug']]
gao_drug_response['drug'].unique()

#%% PDTX name matching

def root_match(str1, str2):
    a = ''
    for i,c in enumerate(str1):
        if len(str2) > i and str2[i] == c:
            a += c
        else:
            break
    return a

expression_pdtx_names = [re.sub('\\.', '-', i) for i in predictions_df['id'].loc[predictions_df['dataset'].isin(['bruna_pdtx'])]]

root_strings = []
for s1 in np.unique(pdtx_drug_response['ID']):
    temp_str = []
    for s2 in expression_pdtx_names:
        temp_str.append(root_match(s1, s2))
    m = np.argmax(np.char.str_len(temp_str))
    
    root = temp_str[m]
    s2 = expression_pdtx_names[m]
    print(f"{s1} ~ {s2}, root: {root}")
    
    
    root_strings.append(temp_str)

best_match = np.argmax(np.char.str_len(np.array(root_strings)), axis = 0)
for i, s1, m in zip(np.arange(len(root_strings)), np.unique(pdtx_drug_response['ID']), best_match):
    root = root_strings[i][m]
    s2 = expression_pdtx_names[m]
    print(f"{s1} ~ {s2}, root: {root}")
#%% Evaluation
from scipy.stats import spearmanr, pearsonr
from sklearn.metrics import r2_score

val_set_list = [
    ('bruna_pdtc', pdtc_drug_response, 'AUC', 'Model'),
    ('gao_pdtx', gao_drug_response, 'BestAvgResponse', 'Model')
]
dr_perf_ext_list = []
for dataset, dr_true_df, target_col, id_col in val_set_list:
    dr_pred = predictions_df.loc[predictions_df['dataset'] == dataset]
    target_drugs = dr_true_df['drug'].unique()
    shared_drugs = target_drugs[np.isin(target_drugs, train_drugs)]
    for di in shared_drugs:
        di_ind = np.argwhere(train_drugs == di)
        di_dr_pred = dr_pred[f"dr_pred_{di_ind[0,0]}"]
        di_dr_pred_id = [re.sub('\\.', '-', i) for i in dr_pred['id']]
        di_dr_pred_id_map = dict(zip(di_dr_pred_id, np.arange(dr_pred.shape[0])))
        
        di_dr_true = dr_true_df.loc[dr_true_df['drug'] == di, ]
        
        di_dr_pred_ind = np.array([di_dr_pred_id_map.get(i, np.nan) for i in di_dr_true[id_col]])
        
        nnan_ind = np.logical_not(np.isnan(di_dr_pred_ind))
        
        di_dr_pred_val = di_dr_pred.iloc[di_dr_pred_ind[nnan_ind]]
        di_dr_true_val = di_dr_true[target_col].to_numpy()[nnan_ind]
        
        di_dr_r2 = r2_score(di_dr_true_val, di_dr_pred_val)
        di_dr_pearson, ppval = pearsonr(di_dr_true_val, di_dr_pred_val)
        di_dr_spearman, spval = spearmanr(di_dr_true_val, di_dr_pred_val)
        
        dr_perf_ext_list.append(pd.DataFrame({
            'dataset' : dataset, 
            'drug' : di, 
            'R2' : di_dr_r2, 
            'PearsonR' : di_dr_pearson, 
            'SpearmanR' : di_dr_spearman}, 
            pd.RangeIndex(0,1)))

dr_perf_ext_df = pd.concat(dr_perf_ext_list)
dr_perf_ext_df.to_csv(f"{res_path}external_drug_response_validation_performance.csv")

