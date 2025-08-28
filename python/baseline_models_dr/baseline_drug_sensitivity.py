#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 14:33:45 2023

@author: rintala
"""

import tensorflow as tf
import numpy as np
import pandas as pd
import time
import os
import warnings
import re

import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC
import sklearn.metrics as metrics
from training_utilities2 import train_aecl_with_pretraining, get_model_losses
from data_utilities import rdata_loader
from sys import platform
from data_utilities import complete_data_loader

#%% new data loader test
if platform == 'linux':
    patient_expression_root_dir = '/research/work/rintala/tcga/pan_cancer/data_full/'
    cell_line_expression_root_dir = '/research/users/rintala/ccle/'
    cell_line_drug_response_root_dir = '/research/users/rintala/ctrp/'
else:
    patient_expression_root_dir = '//research/workdir/tcga/pan_cancer/data_full/'
    cell_line_expression_root_dir = '//research/rintala/ccle/'
    cell_line_drug_response_root_dir = '//research/rintala/ctrp/'
if False:
    patient_expression_root_dir = '~/research_work/tcga/pan_cancer/data_full/'
    cell_line_expression_root_dir = '~/research/ccle/'
    cell_line_drug_response_root_dir = '~/research/ctrp/'
patient_expression_cancer_list = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 
                                  'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 
                                  'KIRC', 'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD',
                                  'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 
                                  'PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 
                                  'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 
                                  'UVM']
patient_expression_file_name = 'mrna.csv.gz'
patient_expression_log2 = True
patient_survival_file = 'survival.csv.gz'
cell_line_drug_response_file = 'screen_aac.csv.gz'
cell_line_drug_response_row_info_file = 'screen_rowinfo.csv.gz'
cell_line_drug_response_row_map_file = 'ctrp_ccle_clid_map.csv.gz'
cell_line_expression_file = 'CCLE_expression.csv'
cell_line_gene_mapping_file = 'CCLE_gene_mapping.csv'
patient_gene_mapping_file = 'TCGA_gene_mapping.csv'
gene_harmonization_union = False
separate_gene_standardization = True

if True:
    patient_expression_cancer_list = ['BRCA']
    if True:
        cell_line_oncotree_mapping_file = None
        cell_line_oncotree_mapping_column = None
        cell_line_oncotree_cancer_list = []
    else:
        cell_line_oncotree_mapping_file = 'Model_oncotree.csv'
        cell_line_oncotree_mapping_column = 'level_2'
        cell_line_oncotree_cancer_list = ['BRCA']
    
    if False:
        cl_data = pd.read_csv(cell_line_expression_root_dir + cell_line_expression_file, header = 0, index_col = 0)
        cl_data_ind = cl_data.index.to_numpy()
        cl_oncotree = pd.read_csv(cell_line_expression_root_dir + cell_line_oncotree_mapping_file, header = 0, index_col = 1)
        cl_oncotree_code = cl_oncotree.loc[cl_data_ind, cell_line_oncotree_mapping_column]
        cl_oncotree_filter = np.logical_not(cl_oncotree_code.isna().to_numpy())
        cl_oncotree_filter[cl_oncotree_filter] = np.isin(cl_oncotree_code[cl_oncotree_filter].to_numpy(), cell_line_oncotree_cancer_list)
        cl_data_ind = cl_data_ind[cl_oncotree_filter]
        cl_data = cl_data.loc[cl_oncotree_filter,:]
elif False:
    patient_expression_cancer_list = ['LUAD', 'LUSC']
    cell_line_oncotree_mapping_file = 'Model_oncotree.csv'
    cell_line_oncotree_mapping_column = 'level_3'
    cell_line_oncotree_cancer_list = ['LUAD', 'LUSC']
else:
    cell_line_oncotree_mapping_file = None
    cell_line_oncotree_mapping_column = None
    cell_line_oncotree_cancer_list = []

#%% load data

data_dict = complete_data_loader(
    patient_expression_root_dir = patient_expression_root_dir,
    patient_expression_cancer_list = patient_expression_cancer_list,
    patient_expression_file_name = patient_expression_file_name,
    patient_gene_mapping_file = patient_gene_mapping_file, 
    patient_expression_log2 = True, 
    patient_survival_file = patient_survival_file, 
    patient_survival_event_col = 'OS', 
    patient_survival_time_col = 'OS.time', 
    patient_survival_covar_cols = ['gender', 'age_at_initial_pathologic_diagnosis', 'type'],
    patient_survival_covar_onehot = [True, False, True], 
    cell_line_expression_root_dir = cell_line_expression_root_dir,
    cell_line_expression_file = cell_line_expression_file,
    cell_line_gene_mapping_file = cell_line_gene_mapping_file, 
    gene_harmonization_union = False, 
    #separate_gene_standardization = True, 
    cell_line_drug_response_root_dir = cell_line_drug_response_root_dir,
    cell_line_drug_response_file = cell_line_drug_response_file,
    cell_line_drug_response_row_info_file = cell_line_drug_response_row_info_file, 
    #cell_line_drug_response_colid_file = None, #
    cell_line_drug_response_row_map_file = cell_line_drug_response_row_map_file,
    cell_line_drug_response_target_column = 'aac_recomputed',
    cell_line_oncotree_mapping_file = cell_line_oncotree_mapping_file, 
    cell_line_oncotree_mapping_column = cell_line_oncotree_mapping_column, 
    cell_line_oncotree_cancer_list = cell_line_oncotree_cancer_list)

#%% Manual run of complete_data_loader
if False:
    cell_line_drug_response_target_column = 'aac_recomputed'
    
    cl_exp_flag = (cell_line_expression_root_dir is not None) and (cell_line_expression_file is not None)
    cl_dr_flag = cl_exp_flag and (cell_line_drug_response_root_dir is not None) and (cell_line_drug_response_file is not None)
    p_exp_flag = (patient_expression_root_dir is not None) and (patient_expression_file_name is not None) and (len(patient_expression_cancer_list)>0)
    p_surv_flag = p_exp_flag and (patient_survival_file is not None)
    
    p_exp_flag = p_surv_flag = False
    
    if cl_exp_flag:
        cl_data = pd.read_csv(cell_line_expression_root_dir + cell_line_expression_file, header = 0, index_col = 0)
        cl_data_ind = cl_data.index.to_numpy()
        if len(cell_line_oncotree_cancer_list) > 0:
            cl_oncotree = pd.read_csv(cell_line_expression_root_dir + cell_line_oncotree_mapping_file, header = 0, index_col = 1)
            cl_oncotree_code = cl_oncotree.loc[cl_data_ind, cell_line_oncotree_mapping_column]
            cl_oncotree_filter = np.logical_not(cl_oncotree_code.isna().to_numpy())
            cl_oncotree_filter[cl_oncotree_filter] = np.isin(cl_oncotree_code[cl_oncotree_filter].to_numpy(), cell_line_oncotree_cancer_list)
            if np.sum(cl_oncotree_filter) == 0:
                raise ValueError('Oncotree filter returned 0 cell-lines.')
            cl_data_ind = cl_data_ind[cl_oncotree_filter]
            cl_data = cl_data.loc[cl_oncotree_filter,:]
    
    if cl_dr_flag:
        dr_data = pd.read_csv(cell_line_drug_response_root_dir + cell_line_drug_response_file, header = 0, index_col = 0)
        dr_row_info = pd.read_csv(cell_line_drug_response_root_dir + cell_line_drug_response_row_info_file, header = 0, index_col = 0)
        dr_screen_cl_map = pd.read_csv(cell_line_drug_response_root_dir + cell_line_drug_response_row_map_file, header = 0, index_col = 0)
        
        nan_cl_ind = np.any(dr_row_info.isna().to_numpy(), axis = 1)
        dr_cl_id_mapped = dr_screen_cl_map.loc[dr_row_info['sampleid'], 'CCLE_model_id']
        dr_cl_id_mapped_not_nan_ind = np.logical_not(dr_cl_id_mapped.isna().to_numpy())
        no_expr_cl = np.setdiff1d(dr_cl_id_mapped.to_numpy()[dr_cl_id_mapped_not_nan_ind], cl_data_ind)
        no_expr_cl_ind = np.isin(dr_cl_id_mapped.to_numpy(), no_expr_cl)
        
        dr_filter = np.logical_and(np.logical_not(nan_cl_ind), dr_cl_id_mapped_not_nan_ind) # exclude dr with mismatched cl ids
        dr_filter = np.logical_and(dr_filter, np.logical_not(no_expr_cl_ind)) # exclude dr without cl exp
        dr_missing_aac = dr_data[cell_line_drug_response_target_column].isna().to_numpy()
        dr_filter = np.logical_and(dr_filter, np.logical_not(dr_missing_aac)) # exclude dr with missing AAC
        dr_cl_vec = dr_cl_id_mapped[dr_filter].to_numpy(dtype = 'str')
        
        # Summarize mediums
        dr_treat_vec = dr_row_info.loc[dr_filter, 'treatmentid_fixed'].to_numpy(dtype = 'str')
        dr_cl_unique, dr_cl_inverse = np.unique(dr_cl_vec, return_inverse = True)
        dr_treat_unique, dr_treat_inverse = np.unique(dr_treat_vec, return_inverse = True) 
        dr_table = np.full((dr_cl_unique.shape[0], dr_treat_unique.shape[0]), fill_value = 0.)
        dr_table_count = np.full(dr_table.shape, fill_value = 0.)
        dr_vec = dr_data[cell_line_drug_response_target_column].to_numpy()[dr_filter]
        if np.any(np.isnan(dr_vec)):
            raise ValueError('Unhandled missing values in drug response data.')
        
        # Compute average of drug response matrix
        for dr, t, cl in zip(dr_vec, dr_treat_inverse, dr_cl_inverse):
            dr_table[cl,t] += dr
            dr_table_count[cl,t] += 1.
        dr_table_count[dr_table_count == 0] = np.nan
        dr_table = dr_table / dr_table_count
        
        # Encode missing values as negative
        dr_table_missing = np.isnan(dr_table_count)
        dr_table[dr_table_missing] = -1
        dr_table_mask = np.logical_not(dr_table_missing)
        
        # Cell-line expression to dr index
        cl_name_map = dict(zip(dr_cl_unique, np.arange(dr_cl_unique.shape[0], dtype = 'int64')))
        cl_exp_ind = np.array([cl_name_map.get(i, -1) for i in cl_data_ind])
        
        if np.any(np.unique(cl_exp_ind[cl_exp_ind > -1], return_counts = True)[1] > 1):
            raise ValueError('Cell line expression ids match multiple drug response profiles.')

#%% other tests

if False:
    # sigmoid output in NN?
    # What is the distribution of logit?
    import seaborn as sns
    def logit(x):
        return np.log(x) - np.log(1.-x)
    dr_vec = data_dict['dr_table'][data_dict['dr_table_mask']]
    sns.histplot(logit(dr_vec), bins = 30) # Decent

#%% baseline elastic net regression model

from sklearn.model_selection import RepeatedKFold, GridSearchCV, cross_val_score
from sklearn.linear_model import ElasticNet, Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline
from sklearn.metrics import r2_score
from copy import copy
from threadpoolctl import threadpool_limits, threadpool_info

import seaborn as sns

piped_eln_param_grid = {
    'elnet__alpha' : 10**np.linspace(-3, 0, 11), 
    'elnet__l1_ratio' : np.power(np.linspace(0.1**4, 0.99**4, 11), 1./4.)
    }
piped_pca_param_grid = {
    'pca__n_components' : np.arange(2, 20)
    }
rkf_inner = RepeatedKFold(n_splits=5, n_repeats=1, random_state=1)
rkf_outer = RepeatedKFold(n_splits=5, n_repeats=1, random_state=2)
scale_pipe = Pipeline(steps = [
    #('pca', PCA(svd_solver = 'randomized', whiten = True)), # random is default solver for this data dimension
    #('pca', PCA(svd_solver = 'arpack', whiten = True, n_components = 10)), 
    #('scale', StandardScaler(with_mean = True, with_std = True)), 
    ('elnet', ElasticNet(max_iter = 100)),
    #('elnet', Ridge()),
    ])
gs = GridSearchCV(
    scale_pipe, 
    #param_grid = {**piped_eln_param_grid, **piped_pca_param_grid}, 
    param_grid = piped_eln_param_grid, 
    cv = rkf_inner, 
    scoring = 'r2',
    n_jobs = 1)

dr_true = data_dict['dr_table']
dr_true[np.logical_not(data_dict['dr_table_mask'])] = np.nan
sns.histplot(np.nanmean(dr_true, axis = 0))
X = data_dict['cl_exp']

from sklearn.utils._testing import ignore_warnings
from sklearn.exceptions import ConvergenceWarning

@ignore_warnings(category=ConvergenceWarning)
def eval_fun(gs, data_dict = data_dict, X = X, dr_true = dr_true, rkf_outer = rkf_outer):
    r2_df_list_eln = []
    a = 0
    for outer_train_ind, outer_test_ind in rkf_outer.split(X):
        exp_pca = PCA(svd_solver = 'arpack', whiten = True, n_components = 100)
        X_pca_train = exp_pca.fit_transform(X[outer_train_ind, :])
        X_pca_test = exp_pca.transform(X[outer_test_ind, :])
        for coli, col in enumerate(data_dict['dr_table_cols']):
            if coli + 1 > len(r2_df_list_eln):
                dri = dr_true[:, coli]
                dri_train = dri[outer_train_ind]
                dri_test = dri[outer_test_ind]
                nna_ind_train = np.logical_not(np.isnan(dri_train))
                nna_ind_test = np.logical_not(np.isnan(dri_test))
                            
                dri_train_nna = np.expand_dims(dri_train[nna_ind_train], axis = -1)
                dri_test_nna = np.expand_dims(dri_test[nna_ind_test], axis = -1)
                
                Xi_pca_train = X_pca_train[nna_ind_train, :]
                Xi_pca_test = X_pca_test[nna_ind_test, :]
                
                with threadpool_limits(limits = 8, user_api = 'blas'):
                    modeli = copy(gs)
                    modeli.fit(Xi_pca_train, dri_train_nna)
                    predi_train = modeli.predict(Xi_pca_train)
                    predi_test = modeli.predict(Xi_pca_test)
                    r2i_train = r2_score(dri_train_nna, predi_train)
                    r2i_test = r2_score(dri_test_nna, predi_test)
                    
                    #sns.histplot(dri_train_nna)
                    #sns.histplot(dri_test_nna)
                    #sns.histplot(predi_train)
                    #sns.histplot(predi_test)
                    
                    r2_col_df = pd.DataFrame({
                        'run' : a % rkf_outer.get_n_splits(), 
                        'fold' : a // rkf_outer.get_n_splits(), 
                        'drug' : col, 
                        'dr_r2_train' : r2i_train, 
                        'dr_r2_test' : r2i_test},
                        index = pd.RangeIndex(0,1))
                    r2_df_list_eln.append(r2_col_df)
        a += 1
    
    r2_df_eln = pd.concat(r2_df_list_eln, axis = 0)
    return r2_df_eln

r2_df_eln = eval_fun(gs)

#%% plots

r2_df_eln_whitened = copy(r2_df_eln)
sns.histplot(r2_df_eln_whitened['dr_r2_train'])
sns.histplot(r2_df_eln_whitened['dr_r2_test'])

sns.histplot(r2_df_eln['dr_r2_train'])
sns.histplot(r2_df_eln['dr_r2_test'])

#%% baseline gradient boosting regressor

from sklearn.model_selection import RepeatedKFold, GridSearchCV, cross_val_score
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline
from sklearn.metrics import r2_score
from copy import copy
from threadpoolctl import threadpool_limits, threadpool_info

piped_gbr_param_grid = {
    'gbr__max_depth' : np.arange(3, 10), 
    'gbr__ccp_alpha' : 10**np.linspace(-3, 0, 4),
    'gbr__learning_rate' : 10**np.linspace(-3, -1, 5)}
rkf_inner = RepeatedKFold(n_splits=5, n_repeats=1, random_state=1)
rkf_outer = RepeatedKFold(n_splits=5, n_repeats=1, random_state=2)
scale_pipe = Pipeline(steps = [
    #('scale', StandardScaler(with_mean = True, with_std = True)), 
    ('gbr', GradientBoostingRegressor(
        n_estimators = 200, 
        loss = 'squared_error',
        random_state = 0))])
gs = GridSearchCV(
    scale_pipe, 
    param_grid = piped_gbr_param_grid, 
    cv = rkf_inner, 
    n_jobs = 8)

r2_df_list_gbr = []
a = 0
for outer_train_ind, outer_test_ind in rkf_outer.split(X):
    exp_pca = PCA(svd_solver = 'arpack', whiten = True, n_components = 100)
    X_pca_train = exp_pca.fit_transform(X[outer_train_ind, :])
    X_pca_test = exp_pca.transform(X[outer_test_ind, :])
    for coli, col in enumerate(data_dict['dr_table_cols']):
        if coli + 1 > len(r2_df_list_gbr):
            dri = dr_true[:, coli]
            dri_train = dri[outer_train_ind]
            dri_test = dri[outer_test_ind]
            nna_ind_train = np.logical_not(np.isnan(dri_train))
            nna_ind_test = np.logical_not(np.isnan(dri_test))
                        
            dri_train_nna = dri_train[nna_ind_train]
            dri_test_nna = dri_test[nna_ind_test]
            
            Xi_pca_train = X_pca_train[nna_ind_train, :]
            Xi_pca_test = X_pca_test[nna_ind_test, :]
            
            with threadpool_limits(limits = 1, user_api = 'blas'):
                modeli = copy(gs)
                modeli.fit(Xi_pca_train, dri_train_nna)
                predi_train = modeli.predict(Xi_pca_train)
                predi_test = modeli.predict(Xi_pca_test)
                r2i_train = r2_score(dri_train_nna, predi_train)
                r2i_test = r2_score(dri_test_nna, predi_test)
                
                #sns.histplot(dri_train_nna)
                #sns.histplot(dri_test_nna)
                #sns.histplot(predi_train)
                #sns.histplot(predi_test)
                
                r2_col_df = pd.DataFrame({
                    'run' : a % rkf_outer.get_n_splits(), 
                    'fold' : a // rkf_outer.get_n_splits(), 
                    'drug' : col, 
                    'dr_r2_train' : r2i_train, 
                    'dr_r2_test' : r2i_test},
                    index = pd.RangeIndex(0,1))
                r2_df_list_gbr.append(r2_col_df)
    a += 1

r2_df_gbt = pd.concat(r2_df_list_gbr, axis = 0)
