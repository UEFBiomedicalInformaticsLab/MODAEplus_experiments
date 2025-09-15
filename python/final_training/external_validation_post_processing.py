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

import importlib
import sys
import re

script_path = os.environ.get('MODAE_SCRIPT_PATH', default = None)
if script_path is None:
    raise ValueError('Please define MODAE_SCRIPT_PATH')

dsp_spec = importlib.util.spec_from_file_location(
    'dataset_processing', 
    os.path.join(script_path, 'python/utilities/cancer_dataset_processing.py')
)
dsp = importlib.util.module_from_spec(dsp_spec)
sys.modules['dataset_processing'] = dsp
dsp_spec.loader.exec_module(dsp)

import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

#%% Data
output_path = os.environ.get('MODAE_OUTPUT_PATH', default = None)
if output_path is None:
    raise ValueError('Please define MODAE_OUTPUT_PATH')
data_root = os.environ.get('MODAE_DATA_PATH', default = None)
if data_root is None:
    raise ValueError('Please define MODAE_DATA_PATH')

res_path = os.path.join(output_path, '20250410_random_search/pancan_test/')

res_path = os.path.join(res_path, 'external_evaluation/')

embedding_files = glob.glob(os.path.join(res_path, '*validation_embeddings.csv.gz'))

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

prediction_files = glob.glob(os.path.join(res_path, '*validation_predictions.csv.gz'))

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

predictions_df = pd.concat(prediction_list, axis = 0)
predictions_df.value_counts('dataset')

#%% TCGA and CCLE data

patient_expression_root_dir = os.path.join(data_root, 'tcga/')
cell_line_expression_root_dir = os.path.join(data_root, 'ccle/')
cell_line_drug_response_root_dir = os.path.join(data_root, 'ctrp/')
xia_path = os.path.join(data_root, 'drug_sensitivity_xia/')
bruna_path = os.path.join(data_root, 'breast_pdtc_bruna/')
gao_path = os.path.join(data_root, 'pan_cancer_pdx_gao/')

#%% CTRP or CCLE?
'''
Get appropriate drug names for predictions. 
'''
prediction_cols = predictions_df.filter(regex = 'dr_pred_[0-9]+', axis = 1).shape[1]
if prediction_cols == 24:
    response_dataset = 'ccle'
    raise NotImplementedError('CCLE based validation is not implemented yet.')
elif prediction_cols == 545:
    response_dataset = 'ctrp'
    with open(os.path.join(cell_line_drug_response_root_dir, 'drug_names.txt'), 'r') as f:
        train_drugs = f.read().splitlines()
    train_drugs = map(str.lower, train_drugs)
    train_drugs = [re.sub('\\(.*?\\)', '', i) for i in train_drugs]
    train_drugs = [re.sub('\\s', '', i) for i in train_drugs]
    train_drugs = [re.sub(':', '+', i) for i in train_drugs]
    train_drugs = np.array(['tipifarnib' if re.search('tipifarnib', i) else i for i in train_drugs])
elif prediction_cols == 544:
    response_dataset = 'ctrp' # Xia data
    with open(os.path.join(xia_path, 'processed_drug_names.txt'), 'r') as f:
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
    os.path.join(bruna_path, 'DrugResponsesAUCModels.txt'), 
    header = 0, sep = "\t")

pdtx_drug_response = pd.read_csv(
    os.path.join(bruna_path, 'DrugResponsesAUCSamples.txt'), 
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
    os.path.join(gao_path, '41591_2015_BFnm3954_MOESM10_ESM.xlsx'),  
    sheet_name = "PCT curve metrics"
)
gao_drug_response['drug'] = gao_drug_response['Treatment'].map(str.lower)
gao_drug_response['drug'] = [re.sub('\\s', '', i) for i in gao_drug_response['drug']]
gao_drug_response['drug'].unique()

#%% PDTX name matching
'''
Left here as a loose end. Bruna has PDX response data, but gene-expression 
only for PDTC. Matching IDs does not seem to be trivial. 
'''

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
dr_perf_ext_df.to_csv(os.path.join(res_path, 'external_drug_response_validation_performance.csv'))

