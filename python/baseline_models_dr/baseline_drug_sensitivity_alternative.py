#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 14:36:40 2023

@author: teemu
"""


import numpy as np
import pandas as pd
import time
import os
import warnings
import re
import matplotlib.pyplot as plt

from sys import platform
from copy import copy

#%% load pre-computed results
from sys import platform
from socket import gethostname
if gethostname() == 'teemu-pc':
    root_dir = '/home/teemu/research_work/drug_response_dataset/'
    home = True
elif platform == 'linux':
    root_dir = '/research/work/rintala/drug_response_dataset/'
    home = False
else:
    root_dir = '//research/workdir/drug_response_dataset/'
    home = False

eln_aac = pd.read_csv(root_dir + 'elasticnet_aac1.csv', header = 0, index_col = 0)
eln_aac_old = pd.read_csv(root_dir + 'elasticnet_aac_old.csv', header = 0, index_col = 0)
eln_aac_old_fixed = pd.read_csv(root_dir + 'elasticnet_aac_old_fixed.csv', header = 0, index_col = 0)
fs_eln_aac_old_fixed = pd.read_csv(root_dir + 'fs_elasticnet_aac_old_fixed.csv', header = 0, index_col = 0)
pca10_eln_aac_old_fixed = pd.read_csv(root_dir + 'pca10_elasticnet_aac_old.csv', header = 0, index_col = 0)

#%%
import seaborn as sns
sns.histplot(eln_aac['train_r2'])
sns.histplot(eln_aac['test_r2'])

sns.histplot(eln_aac_old['train_r2'])
sns.histplot(eln_aac_old['test_r2'])

sns.histplot(eln_aac_old_fixed['train_r2'], bins = 30)
sns.histplot(eln_aac_old_fixed['test_r2'], bins = 30)

fs_eln_aac_old_fixed_filtered = fs_eln_aac_old_fixed.loc[fs_eln_aac_old_fixed['test_r2'] > -1, :]
sns.histplot(fs_eln_aac_old_fixed_filtered['train_r2'], bins = 30)
sns.histplot(fs_eln_aac_old_fixed_filtered['test_r2'], bins = 30)

# fix folds
eln_aac_old_fixed['fold'] = eln_aac_old_fixed['fold'] // 545
eln_aac_old_fixed.groupby(['fold']).apply(lambda x: x['drug']).duplicated() # correct
fs_eln_aac_old_fixed['fold'] = fs_eln_aac_old_fixed['fold'] // 545
fs_eln_aac_old_fixed.groupby(['fold']).apply(lambda x: x['drug']).duplicated() # correct
pca10_eln_aac_old_fixed['fold'] = pca10_eln_aac_old_fixed['fold'] // 545
pca10_eln_aac_old_fixed.groupby(['fold']).apply(lambda x: x['drug']).duplicated() # correct

nofs_train = eln_aac_old_fixed.groupby(['fold'])['train_r2'].quantile(0.9).mean()
nofs_test = eln_aac_old_fixed.groupby(['fold'])['test_r2'].quantile(0.9).mean()
fs_train = fs_eln_aac_old_fixed.groupby(['fold'])['train_r2'].quantile(0.9).mean()
fs_test = fs_eln_aac_old_fixed.groupby(['fold'])['test_r2'].quantile(0.9).mean()
pca10_train = pca10_eln_aac_old_fixed.groupby(['fold'])['train_r2'].quantile(0.9).mean()
pca10_test = pca10_eln_aac_old_fixed.groupby(['fold'])['test_r2'].quantile(0.9).mean()

print(f"No FS - train R2: {nofs_train}, test R2: {nofs_test}\nFS - train R2: {fs_train}, test R2: {fs_test}\nPCA10 - train R2: {pca10_train}, test R2: {pca10_test}")

#%% brca drug performnace
import json

eln_res_mean_df = eln_aac

best_drug = eln_res_mean_df.index[eln_res_mean_df['test_r2'].argmax()]
print(eln_res_mean_df.loc[best_drug])
print(drug_info.loc[best_drug])

if platform == 'linux':
    patient_expression_root_dir = '/research/work/rintala/tcga/pan_cancer/data_full/'
    cell_line_expression_root_dir = '/research/users/rintala/ccle/'
    cell_line_drug_response_root_dir = '/research/users/rintala/ctrp/'
else:
    patient_expression_root_dir = '//research/workdir/tcga/pan_cancer/data_full/'
    cell_line_expression_root_dir = '//research/rintala/ccle/'
    cell_line_drug_response_root_dir = '//research/rintala/ctrp/'
if False:
    patient_expression_root_dir = '/home/teemu/research_work/tcga/pan_cancer/data_full/'
    cell_line_expression_root_dir = '/home/teemu/research/ccle/'
    cell_line_drug_response_root_dir = '/home/teemu/research/ctrp/'
with open(cell_line_drug_response_root_dir + 'brca_drugs.json', 'r') as f:
    brca_drugs_json = f.readline()
brca_drugs = json.loads(brca_drugs_json)

brca_drug_match = []
for bd in brca_drugs:
    if type(bd) is str:
        #print('str')
        ind = drug_info['NAME'].str.contains(bd, case = False)
        if ind.any():
            brca_drug_match.append(drug_info.index[ind].to_list())
        else:
            brca_drug_match.append(None)
    if type(bd) is list:
        #print('list')
        hits = []
        for bdi in bd:
            ind = drug_info['NAME'].str.contains(bdi, case = False)
            if ind.any():
                hits.append(drug_info.index[ind].to_list())
        if len(hits):
            brca_drug_match.append(hits)
        else:
            brca_drug_match.append(None)


ind = drug_info['NAME'].str.contains('Erlotinib', case = False)
drug_info.loc[ind]

eln_res_mean_df.loc[np.intersect1d(X_key['DRUG'].unique(), drug_info.loc[ind].index)]

#%% load data

# RNA-Seq file seems to be a .tsv
with open(root_dir + 'combined_rnaseq_data', 'r') as f:
    lines = [f.readline() for i in np.arange(2)]

cl_metadata = pd.read_csv(root_dir + 'combined_cl_metadata', sep = '\t', index_col = 0)
#cl_metadata.set_index('sample_name', inplace = True)
cl_rnaseq = pd.read_csv(root_dir + 'combined_rnaseq_data', sep = '\t', index_col = 0)
#cl_rnaseq.set_index('Sample', inplace = True)
#check = [np.isin(cl_metadata.index[i], cl_rnaseq.index) for i in np.arange(cl_metadata.shape[0])]
#np.unique(check, return_counts = True)
cl_rnaseq = cl_rnaseq.loc[cl_metadata.index]
cl_mapping = pd.read_csv(
    root_dir + 'cl_mapping', 
    sep = '\t', 
    header = None, 
    names = ['cl1', 'cl2'])
cl_response = pd.read_csv(
    root_dir + 'combined_single_response_agg', 
    sep = '\t', 
    dtype = {
        'SOURCE' : 'str', 
        'CELL' : 'str', 
        'DRUG' : 'str', 
        'STUDY' : 'str', 
        'AUC' : 'float64', 
        'IC50' : 'float64', 
        'EC50' : 'float64', 
        'EC50se' : 'float64', 
        'R2fit' : 'float64', 
        'Einf' : 'float64', 
        'HS' : 'float64', 
        'AAC1' : 'float64', 
        'AUC1' : 'float64', 
        'DSS1' : 'float64'
    })
cl_response.to_csv(
    f"{root_dir}/combined_single_response_agg.csv.gz", 
    header = True, 
    index = False, 
    compression = 'gzip'
)

drug_dragon_desc = pd.read_csv(
    root_dir + 'Combined_PubChem_dragon7_descriptors.tsv', 
    sep = '\t', 
    na_values = ['na'], 
    index_col = 0)
drug_dragon_ecfp = pd.read_csv(
    root_dir + 'Combined_PubChem_dragon7_ECFP.tsv', 
    sep = '\s+', skiprows = 1, header = None, index_col = 0)
drug_dragon_pfp = pd.read_csv(
    root_dir + 'Combined_PubChem_dragon7_PFP.tsv', 
    sep = '\s+', skiprows = 1, header = None, index_col = 0)

drug_mordred = pd.read_csv(
    root_dir + 'extended_combined_mordred_descriptors', 
    sep = '\t', index_col = 0)
    #sep = '\s+', skiprows = 1, header = None, index_col = 0)
test = drug_mordred.astype('str')
def check_isnumeric(x):
    return pd.to_numeric(x, errors = 'coerce').isna()
test_num = test.apply(check_isnumeric, axis = 0)
test_sums = test_num.apply(pd.Series.value_counts)
test_mixed = test_sums.apply(pd.Series.isna).any()
mixed_string_ind = test_mixed.index[np.logical_not(test_mixed)]
for i in mixed_string_ind:
    #drug_mordred.loc[test_num.loc[:,i],i]
    drug_mordred.loc[:,i] = pd.to_numeric(drug_mordred.loc[:,i], errors = 'coerce')

drug_info = pd.read_csv(
    root_dir + 'drug_info', 
    sep = '\t', index_col = 0)
if False:
    if False:
        with open(f"{root_dir}processed_drug_columns.txt", 'w') as f:
            f.writelines(np.char.add(data_dict['dr_table_cols'], '\n'))
    with open(f"{root_dir}processed_drug_columns.txt", 'r') as f:
        processed_drug_cols = f.readlines()
    processed_drug_cols = [re.sub('\n', '', i) for i in processed_drug_cols]
    if False:
        with open(f"{root_dir}processed_drug_names.txt", 'r') as f:
            processed_drug_names = f.readlines()
    with open(f"{root_dir}processed_drug_names.txt", 'w') as f:
        f.writelines(np.char.add(drug_info.loc[processed_drug_cols,'NAME'].to_list(), '\n'))
#%% response exploration
import seaborn as sns

cl_response['SOURCE'].value_counts()
cl_response['STUDY'].value_counts()

sns.histplot(cl_response['AUC'])
#sns.histplot(cl_response['IC50'])
print(cl_response['IC50'].min())
print(cl_response['IC50'].max())
sns.histplot(cl_response['AAC1'])
sns.histplot(cl_response['AUC1'])
sns.histplot(cl_response['DSS1'])

#%% Build training matrices

# Response data
X_key = cl_response.loc[cl_response['SOURCE'] == 'CTRP', ['CELL', 'DRUG', 'AAC1']]
X_key.duplicated(subset = ['CELL', 'DRUG']).value_counts()

# Filter based on availability
X_key_nonmissing_exp = X_key['CELL'].isin(cl_rnaseq.index)
X_key = X_key.loc[X_key_nonmissing_exp]

X_key['DRUG'].isin(drug_info.index).value_counts()

X_key['drug_pubchem'] = drug_info.loc[X_key['DRUG'], 'PUBCHEM'].reset_index(drop = True)
X_key['drug_pubchem'].isna().value_counts() # > 30% missing
X_key.dropna(how = 'any', inplace = True)

X_key['drug_pubchem'] = [f"PubChem.CID.{i}" for i in X_key['drug_pubchem'].astype('int64')]

X_key_nonmissing_drug_ecfp = X_key['drug_pubchem'].isin(drug_dragon_ecfp.index)
X_key_nonmissing_drug_ecfp.value_counts() # 0 missing
X_key = X_key.loc[X_key_nonmissing_drug_ecfp]

X_key_nonmissing_drug_pfp = X_key['drug_pubchem'].isin(drug_dragon_pfp.index)
X_key_nonmissing_drug_pfp.value_counts() # 0 missing
X_key = X_key.loc[X_key_nonmissing_drug_pfp]

X_key_nonmissing_drug_mordred = X_key['DRUG'].isin(drug_mordred.index)
X_key_nonmissing_drug_mordred.value_counts() # > 10% missing
X_key = X_key.loc[X_key_nonmissing_drug_mordred]

# Gather features
#X_cl_exp = cl_rnaseq.loc[X_key['CELL']]
X_drug_ecfp = drug_dragon_ecfp.loc[X_key['drug_pubchem']].dropna(axis = 1).to_numpy()
X_drug_pfp = drug_dragon_pfp.loc[X_key['drug_pubchem']].dropna(axis = 1).to_numpy()
X_drug_mordred = drug_mordred.loc[X_key['DRUG']].dropna(axis = 1).to_numpy()

#%% Expression exploration
import seaborn as sns
sns.histplot(cl_rnaseq.loc[X_key['CELL'].unique()].to_numpy().flatten())

#%% Evaluation function for arbitrary models with cross-validation excluding cell-lines
from sklearn.model_selection import GroupKFold, ParameterGrid
from sklearn.decomposition import PCA
from threadpoolctl import threadpool_limits, threadpool_info
from sklearn.utils._testing import ignore_warnings
from sklearn.exceptions import ConvergenceWarning

#cv_outer = GroupKFold(n_splits = 5)
#cv_inner = GroupKFold(n_splits = 5)

@ignore_warnings(category=ConvergenceWarning)
def nested_train_eval_drugwise(
        model, 
        param_grid, 
        key_df, 
        y_col, 
        drug_col, 
        cell_col, 
        cl_exp_df, 
        pca_instance = None, 
        cv_outer = GroupKFold(n_splits = 5), 
        cv_inner = GroupKFold(n_splits = 5), 
        group_col = 'cell', 
        thread_limit = 1):
    if group_col == 'cell':
        cv_outer_args = {'groups' : key_df[cell_col].to_numpy()}
    else:
        cv_outer_args = {}
    cv_results = []
    for train_ind, test_ind in cv_outer.split(key_df, **cv_outer_args):
        cl_train_names = key_df[cell_col].iloc[train_ind].unique()
        cl_train_map = pd.Series(np.arange(len(cl_train_names)), index = cl_train_names)
        cl_test_names = key_df[cell_col].iloc[test_ind].unique()
        cl_test_map = pd.Series(np.arange(len(cl_test_names)), index = cl_test_names)
        
        if pca_instance is not None:
            cl_pca = copy(pca_instance)
            cl_pcs_train = cl_pca.fit_transform(cl_exp_df.loc[cl_train_names].to_numpy())
            cl_pcs_test = cl_pca.transform(cl_exp_df.loc[cl_test_names].to_numpy())
        
        key_df_train = key_df.iloc[train_ind]
        key_df_test = key_df.iloc[test_ind]
        
        drug_names = key_df[drug_col].unique()
        for dn in drug_names:
            key_df_train_di = key_df_train.loc[key_df_train[drug_col] == dn]
            key_df_test_di = key_df_test.loc[key_df_test[drug_col] == dn]
            
            if pca_instance is not None:
                X_cl_train = cl_pcs_train[cl_train_map.loc[key_df_train_di[cell_col]].to_numpy(), :]
                X_cl_test = cl_pcs_test[cl_test_map.loc[key_df_test_di[cell_col]].to_numpy(), :]
            else:
                X_cl_train = cl_exp_df.loc[key_df_train_di[cell_col], :].to_numpy()
                X_cl_test = cl_exp_df.loc[key_df_test_di[cell_col], :].to_numpy()
            
            X_train = X_cl_train
            X_test = X_cl_test
            
            y_train = key_df_train_di[y_col]
            y_test = key_df_test_di[y_col]
            
            grid_search_scores = []
            if group_col == 'cell':
                cv_inner_args = {'groups' : key_df_train_di[cell_col].to_numpy()}
            else:
                cv_inner_args = {}
            for par_dict in ParameterGrid(param_grid):
                res = []
                for train_ind_inner, test_ind_inner in cv_outer.split(key_df_train_di, **cv_inner_args):
                    X_train_inner = X_train[train_ind_inner, :]
                    X_test_inner = X_train[test_ind_inner, :]
                    y_train_inner = key_df_train_di[y_col].iloc[train_ind_inner].to_numpy()
                    y_test_inner = key_df_train_di[y_col].iloc[test_ind_inner].to_numpy()
                    
                    model_inner = copy(model)
                    model_inner = model_inner.set_params(**par_dict)
                    with threadpool_limits(limits = thread_limit, user_api = 'blas'):
                        model_inner.fit(X_train_inner, y_train_inner)
                    y_r2_train_inner = model_inner.score(X_train_inner, y_train_inner)
                    y_r2_test_inner = model_inner.score(X_test_inner, y_test_inner)
                    res.append(pd.DataFrame({
                        'train_r2' : y_r2_train_inner, 
                        'test_r2' : y_r2_test_inner}, 
                        index = pd.RangeIndex(0,1)))
                res_df = pd.concat(res, axis = 0)
                grid_search_scores.append(pd.DataFrame({
                    **par_dict, 
                    'train_r2' : res_df['train_r2'].mean(), 
                    'test_r2' : res_df['test_r2'].mean()}, 
                    index = pd.RangeIndex(0,1)))
            
            grid_search_scores_df = pd.concat(grid_search_scores, axis = 0)
            best_params = grid_search_scores_df.astype('object').iloc[grid_search_scores_df['test_r2'].argmax()]
            best_par_dict = best_params[param_grid.keys()].to_dict()
            
            model_best = copy(model)
            model_best = model_best.set_params(**best_par_dict)
            with threadpool_limits(limits = thread_limit, user_api = 'blas'):
                model_best.fit(X_train, y_train.to_numpy())
            y_r2_train = model_best.score(X_train, y_train.to_numpy())
            y_r2_test = model_best.score(X_test, y_test.to_numpy())
            
            cv_results.append(pd.DataFrame({
                'drug' : dn, 
                'fold' : len(cv_results), 
                **par_dict, 
                'train_r2' : y_r2_train, 
                'test_r2' : y_r2_test, 
                'inner_train_r2' : best_params['train_r2'], 
                'inner_test_r2' : best_params['test_r2']}, 
                index = pd.RangeIndex(0,1)))
    
    cv_results_df = pd.concat(cv_results, axis = 0)
    
    return cv_results_df

#%% Elastic-net
from sklearn.linear_model import ElasticNet
eln_model = ElasticNet(max_iter = 100, random_state = 2)

eln_param_grid = {
    'alpha' : 10**np.linspace(-3, 0, 11), 
    'l1_ratio' : np.power(np.linspace(0.1**4, 0.99**4, 11), 1./4.)}

eln_res = []
for target in ['AAC1']:
    res = nested_train_eval_drugwise(
        model = eln_model, 
        param_grid = eln_param_grid, 
        key_df = X_key, 
        y_col = target, 
        drug_col = 'DRUG', 
        cell_col = 'CELL', 
        cl_exp_df = cl_rnaseq, 
        #pca_instance = PCA(n_components = 100, svd_solver = 'arpack'), 
        pca_instance = None, 
        cv_outer = GroupKFold(n_splits = 5), 
        cv_inner = GroupKFold(n_splits = 5), 
        group_col = 'cell', 
        thread_limit = 8)
    res['target'] = target
    eln_res.append(res)

#%% Elastic-net with feature selection
from sklearn.linear_model import ElasticNet
from sklearn.feature_selection import SelectKBest, f_regression
from sklearn.pipeline import Pipeline
eln_model = ElasticNet(max_iter = 100, random_state = 2)

fs_eln_model = Pipeline(steps = [
    ('fs', SelectKBest(score_func = f_regression)), 
    ('eln', eln_model)])

pipeline_param_grid = {
    'fs__k' : [2,10], #np.linspace(10, 100, 11), 
    'eln__alpha' : 10**np.linspace(-3, 0, 4), 
    'eln__l1_ratio' : np.power(np.linspace(0.1**4, 0.99**4, 3), 1./4.)}

eln_res = []
for target in ['AAC1']:
    res = nested_train_eval_drugwise(
        model = fs_eln_model, 
        param_grid = pipeline_param_grid, 
        key_df = X_key, 
        y_col = target, 
        drug_col = 'DRUG', 
        cell_col = 'CELL', 
        cl_exp_df = cl_rnaseq, 
        #pca_instance = PCA(n_components = 100, svd_solver = 'arpack'), 
        pca_instance = None, 
        cv_outer = GroupKFold(n_splits = 5), 
        cv_inner = GroupKFold(n_splits = 5), 
        group_col = 'cell', 
        thread_limit = 8)
    res['target'] = target
    eln_res.append(res)

#%% R2 plots
import seaborn as sns
res_cols = ['train_r2', 'test_r2', 'inner_train_r2', 'inner_test_r2']

eln_res_mean_df = eln_res[0].loc[:,['drug'] + res_cols].groupby(['drug']).mean()

eln_res_mean_df_stacked = eln_res_mean_df.stack().reset_index()
eln_res_mean_df_stacked.rename({'level_1' : 'data_set', 0 : 'R2'}, axis = 1, inplace = True)

sns.histplot(eln_res_mean_df_stacked, x = 'R2', hue = 'data_set')

#%% Evaluation function for arbitrary models with cross-validation excluding cell-lines
from sklearn.model_selection import GroupKFold, ParameterGrid
from sklearn.decomposition import PCA
from threadpoolctl import threadpool_limits, threadpool_info

#cv_outer = GroupKFold(n_splits = 5)
#cv_inner = GroupKFold(n_splits = 5)

@ignore_warnings(category=ConvergenceWarning)
def nested_train_eval(model, param_grid, key_df, y_col, cell_col, cl_exp_df, drug_features, 
                      cv_outer, cv_inner, group_col = 'cell', npcs = 100, thread_limit = 1):
    if group_col == 'cell':
        cv_outer_args = {'groups' : key_df[cell_col].to_numpy()}
    else:
        cv_outer_args = {}
    cv_results = []
    for train_ind, test_ind in cv_outer.split(key_df, **cv_outer_args):
        cl_pca = PCA(n_components = npcs, svd_solver = 'arpack')
        cl_train_names = key_df[cell_col].iloc[train_ind].unique()
        cl_pcs_train = cl_pca.fit_transform(cl_exp_df.loc[cl_train_names].to_numpy())
        cl_train_map = pd.Series(np.arange(len(cl_train_names)), index = cl_train_names)
        X_cl_train = cl_pcs_train[cl_train_map.loc[key_df[cell_col].iloc[train_ind]].to_numpy(), :]
        
        cl_test_names = key_df[cell_col].iloc[test_ind].unique()
        cl_pcs_test = cl_pca.transform(cl_exp_df.loc[cl_test_names].to_numpy())
        cl_test_map = pd.Series(np.arange(len(cl_test_names)), index = cl_test_names)
        X_cl_test = cl_pcs_test[cl_test_map.loc[key_df[cell_col].iloc[test_ind]].to_numpy(), :]
        
        X_train = np.concatenate((X_cl_train, drug_features[train_ind, :]), axis = 1)
        X_test = np.concatenate((X_cl_test, drug_features[test_ind, :]), axis = 1)
        
        key_df_train = key_df.iloc[train_ind]
        key_df_test = key_df.iloc[test_ind]
        
        y_train = key_df_train[y_col]
        y_test = key_df_test[y_col]
        
        grid_search_scores = []
        if group_col == 'cell':
            cv_inner_args = {'groups' : key_df_train[cell_col].to_numpy()}
        else:
            cv_inner_args = {}
        for par_dict in ParameterGrid(param_grid):
            res = []
            for train_ind_inner, test_ind_inner in cv_outer.split(key_df_train, **cv_inner_args):
                X_train_inner = X_train[train_ind_inner, :]
                X_test_inner = X_train[test_ind_inner, :]
                y_train_inner = key_df_train[y_col].iloc[train_ind_inner].to_numpy()
                y_test_inner = key_df_train[y_col].iloc[test_ind_inner].to_numpy()
                
                model_inner = copy(model)
                model_inner = model_inner.set_params(**par_dict)
                with threadpool_limits(limits = thread_limit, user_api = 'blas'):
                    model_inner.fit(X_train_inner, y_train_inner)
                y_r2_train_inner = model_inner.score(X_train_inner, y_train_inner)
                y_r2_test_inner = model_inner.score(X_test_inner, y_test_inner)
                res.append(pd.DataFrame({
                    'train_r2' : y_r2_train_inner, 
                    'test_r2' : y_r2_test_inner}, 
                    index = pd.RangeIndex(0,1)))
            res_df = pd.concat(res, axis = 0)
            grid_search_scores.append(pd.DataFrame({
                **par_dict, 
                'train_r2' : res_df['train_r2'].mean(), 
                'test_r2' : res_df['test_r2'].mean()}, 
                index = pd.RangeIndex(0,1)))
        
        grid_search_scores_df = pd.concat(grid_search_scores, axis = 0)
        best_params = grid_search_scores_df.iloc[grid_search_scores_df['test_r2'].argmax()]
        best_par_dict = best_params[param_grid.keys()].to_dict()
        
        model_best = copy(model)
        model_best = model_best.set_params(**best_par_dict)
        with threadpool_limits(limits = thread_limit, user_api = 'blas'):
            model_best.fit(X_train, y_train.to_numpy())
        y_r2_train = model_best.score(X_train, y_train.to_numpy())
        y_r2_test = model_best.score(X_test, y_test.to_numpy())
        
        cv_results.append(pd.DataFrame({
            'fold' : len(cv_results), 
            **par_dict, 
            'train_r2' : y_r2_train, 
            'test_r2' : y_r2_test, 
            'inner_train_r2' : best_params['train_r2'], 
            'inner_test_r2' : best_params['test_r2']}, 
            index = pd.RangeIndex(0,1)))
    
    cv_results_df = pd.concat(cv_results, axis = 0)
    
    return cv_results_df

#%% Elastic-net
from sklearn.linear_model import ElasticNet
eln_model = ElasticNet(max_iter = 100, random_state = 2)

eln_param_grid = {
    'alpha' : 10**np.linspace(-3, 0, 11), 
    'l1_ratio' : np.power(np.linspace(0.1**4, 0.99**4, 11), 1./4.)}

eln_res = []
for dfi in [('dragon_pfp', X_drug_pfp)]:
    for target in ['AAC1']:
        res = nested_train_eval(
            model = eln_model, 
            param_grid = eln_param_grid, 
            key_df = X_key, 
            y_col = target, 
            cell_col = 'CELL', 
            cl_exp_df = cl_rnaseq, 
            drug_features = dfi[1], 
            cv_outer = GroupKFold(n_splits = 5), 
            cv_inner = GroupKFold(n_splits = 5), 
            group_col = 'cell', 
            npcs = 100, 
            thread_limit = 20)
        res['target'] = target
        res['features'] = dfi[0]
        eln_res.append(res)

#%% GBT
from sklearn.ensemble import GradientBoostingRegressor

gbr_param_grid = {
    'max_depth' : np.arange(3, 10), 
    'ccp_alpha' : 10**np.linspace(-3, 0, 4),
    'learning_rate' : 10**np.linspace(-3, -1, 5)}
gbr_model = GradientBoostingRegressor(
    n_estimators = 200, 
    loss = 'squared_error',
    random_state = 2)

gbr_res = []
for dfi in [('dragon_pfp', X_drug_pfp)]:
    for target in ['AAC1']:
        res = nested_train_eval(
            model = gbr_model, 
            param_grid = gbr_param_grid, 
            key_df = X_key, 
            y_col = target, 
            cell_col = 'CELL', 
            cl_exp_df = cl_rnaseq, 
            drug_features = dfi[1], 
            cv_outer = GroupKFold(n_splits = 5), 
            cv_inner = GroupKFold(n_splits = 5), 
            group_col = 'cell', 
            npcs = 100, 
            thread_limit = 20)
        res['target'] = target
        res['features'] = dfi[0]
        gbr_res.append(res)