#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 15:39:10 2023

@author: teemu
"""

import argparse
import os
import sys
import re

import pandas as pd
import numpy as np
import importlib

dsc_spec = importlib.util.spec_from_file_location('dataset_collections', 'datasets/cancer_dataset_collections.py')
dsc = importlib.util.module_from_spec(dsc_spec)
sys.modules['dataset_collections'] = dsc
dsc_spec.loader.exec_module(dsc)

drd_spec = importlib.util.spec_from_file_location('baseline_dr_data', 'baseline_models_dr/baseline_drug_sensitivity_data.py')
drd = importlib.util.module_from_spec(drd_spec)
sys.modules['baseline_dr_data'] = drd
drd_spec.loader.exec_module(drd)

drm_spec = importlib.util.spec_from_file_location('baseline_dr_models', 'baseline_models_dr/baseline_drug_sensitivity_models.py')
drm = importlib.util.module_from_spec(drm_spec)
sys.modules['baseline_dr_models'] = drm
drm_spec.loader.exec_module(drm)

#%%

from dataset_collections import get_tcga_brca_ctrp_ccle_full
from baseline_dr_data import load_sensitivity_data, load_old_sensitivity_data
from baseline_dr_models import nested_train_eval_drugwise

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.linear_model import ElasticNet
from sklearn.ensemble import GradientBoostingRegressor

from sklearn.feature_selection import SelectKBest, f_regression
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GroupKFold

from threadpoolctl import threadpool_limits, threadpool_info

#thread_limit = 20

#%%

def eval_genes(X, y, gene_lists, dr_model):
    drij_r2 = []
    for genes in gene_lists:
        with threadpool_limits(limits = 1, user_api = 'blas'):
            m = dr_model.fit(X[:, genes], y)
            drij_r2.append(m.score(X[:, genes], y))
    return drij_r2

#%% Main
def main(
        root_dir,
        cl_exp_root = None, 
        cl_dr_root = None, 
        identical_data = False,
        old_data = False, 
        gene_preselection = False, 
        drug_features = '', 
        response_sources = 'CTRP', 
        feature_selection_only = False,
        multivariate_fs = False, 
        fs_n_features = 5000, 
        threads = 1, 
        target = 'AAC1', 
        result_file = '',
        standardize = False, 
        use_pca = False, 
        pca_npc = 100, 
        model = 'elasticnet'
):
    os.makedirs(root_dir, exist_ok=True)
    if identical_data:
        data_dict_combined = get_tcga_brca_ctrp_ccle_full(
            home = False, 
            gene_preselection = gene_preselection
        )
        ind = np.argwhere(data_dict_combined['dr_table_mask'])
        ccle_ids = data_dict_combined['cl_exp_rows']
        drug_ids = data_dict_combined['dr_table_cols']
        cl_exp_df = pd.DataFrame(
            data_dict_combined['cl_exp'], 
            index = ccle_ids, 
            columns = data_dict_combined['gene_ids'])
        key_df = pd.DataFrame({
            'CCLE_ID' : ccle_ids[ind[:,0]], 
            'CTRP_TREAT_ID' : drug_ids[ind[:,1]], 
            'AAC' : data_dict_combined['dr_table'][ind[:,0], ind[:,1]]
        })
        data_dict = {'X_key' : key_df, 'cl_rnaseq' : cl_exp_df}
    elif old_data:
        data_dict = load_old_sensitivity_data(
            cell_line_expression_root_dir = cl_exp_root, 
            cell_line_drug_response_root_dir = cl_dr_root
        )
    else:
        data_dict = load_sensitivity_data(
            root_dir = root_dir, 
            drug_features = [i.strip() for i in drug_features.split(',')], 
            response_sources = [i.strip() for i in response_sources.split(',')])
    
    if feature_selection_only:
        fs_dict = {}
        key_df = data_dict['X_key']
        exp_mat = data_dict['cl_rnaseq']
        drug_col = 'CTRP_TREAT_ID'
        clid_col = 'CCLE_ID'
        res_col = 'AAC'
        drug_names = key_df[drug_col].unique()
        
        if multivariate_fs:
            # Manual forward selection
            from sklearn.linear_model import ElasticNet
            from sklearn.preprocessing import StandardScaler
            from multiprocessing import Pool
            from copy import copy
            
            gene_names = exp_mat.columns.to_numpy()
            in_gene_list = np.arange(exp_mat.shape[1]).tolist()
            out_gene_list = []
            dr_r2_df_list = []
            dr_model = ElasticNet(alpha = 1e-3, l1_ratio = 0.9, max_iter = 100)
            
            X_scaler = StandardScaler()
            X_all = X_scaler.fit_transform(exp_mat.to_numpy())
            X_list = []
            y_list = []
            for dn in drug_names:
                res_df = key_df.loc[key_df[drug_col] == dn,:]
                cl_ind = [exp_mat.index.get_loc(i) for i in res_df[clid_col]]
                X_list.append(X_all[cl_ind, :])
                y_list.append(res_df[res_col].to_numpy())
            for i in np.arange(fs_n_features):
                dri_r2 = []
                gene_lists = [out_gene_list + [j] for j in in_gene_list]
                with Pool(processes = threads) as mp:
                    map_iterator = [(X, y, gene_lists, copy(dr_model)) for X, y in zip(X_list, y_list)]
                    for res in mp.starmap(eval_genes, map_iterator):
                        dri_r2.append(res)
                dri_r2 = np.array(dri_r2)
                dri_r2_mean = dri_r2.mean(axis = 0)
                j_max = np.argmax(dri_r2_mean)
                j_max_gene = in_gene_list.pop(j_max)
                out_gene_list.append(j_max_gene)
                res_df = pd.DataFrame({
                    'gene' : [gene_names[j_max_gene]], 
                    'R2_mean' : [dri_r2_mean[j_max]], 
                })
                fn = f"{root_dir}{result_file}"
                res_df.to_csv(fn, mode = 'a')
        else:
            for dn in drug_names:
                res_df = key_df.loc[key_df[drug_col] == dn,:]
                with threadpool_limits(limits = threads, user_api = 'blas'):
                    f,p = f_regression(
                        exp_mat.loc[res_df[clid_col]].to_numpy(), 
                        res_df[res_col].to_numpy())
                #fs_dict[dn] = exp_mat.columns[p < 0.01]
                fs_dict[dn] = exp_mat.columns[np.argsort(p)[:10]]
            fs_all = np.concatenate([i for i in fs_dict.values()])
            fs_all_uq = np.unique(fs_all)
            fs_symbols = [re.sub(' \\([0-9]+\\)', '', i) for i in fs_all_uq]
            with open(f"{root_dir}top10_univariate_genes_combined.txt", 'w') as f:
                f.writelines([f"{i}\n" for i in fs_symbols])
    
    if standardize:
        scaler_instance = StandardScaler()
    else:
        scaler_instance = None
    if use_pca:
        pca_instance = PCA(n_components = pca_npc, whiten = True)
    else:
        pca_instance = None
    
    if model == 'elasticnet':
        model = ElasticNet(max_iter = 100, random_state = 2)
        param_grid = {
            'alpha' : 10**np.linspace(-3, 0, 11), 
            'l1_ratio' : np.power(np.linspace(0.1**4, 0.99**4, 11), 1./4.)}
    elif model == 'fs_elasticnet':
        eln_model = ElasticNet(max_iter = 100, random_state = 2)
        model = Pipeline(steps = [
            ('fs', SelectKBest(score_func = f_regression)), 
            ('eln', eln_model)])
        param_grid = {
            'fs__k' : np.linspace(100, 1000, 3).astype('int64'), 
            'eln__alpha' : 10**np.linspace(-3, -1, 3), 
            'eln__l1_ratio' : np.power(np.linspace(0.1**4, 0.99**4, 3), 1./4.)}
    elif model == 'gbt':
        param_grid = {
            'max_depth' : np.arange(3, 10), 
            'ccp_alpha' : 10**np.linspace(-3, 0, 4),
            'learning_rate' : 10**np.linspace(-3, -1, 5)}
        model = GradientBoostingRegressor(
            n_estimators = 200, 
            loss = 'squared_error',
            random_state = 2)
    elif model == 'fs_gbt':
        gbr_model = GradientBoostingRegressor(
            n_estimators = 200, 
            loss = 'squared_error',
            random_state = 2)
        model = Pipeline(steps = [
            ('fs', SelectKBest(score_func = f_regression)), 
            ('gbt', gbr_model)])
        param_grid = {
            'fs__k' : np.linspace(100, 1000, 3).astype('int64'), 
            'gbt__max_depth' : np.arange(2, 9, 2), 
            'gbt__ccp_alpha' : 10**np.linspace(-3, 0, 4), 
            'gbt__learning_rate' : np.power(np.linspace(0.1**4, 0.99**4, 3), 1./4.)}
    
    if old_data or identical_data:
        res = nested_train_eval_drugwise(
            model = model, 
            param_grid = param_grid, 
            key_df = data_dict['X_key'], 
            y_col = 'AAC', 
            drug_col = 'CTRP_TREAT_ID', 
            cell_col = 'CCLE_ID', 
            cl_exp_df = data_dict['cl_rnaseq'], 
            scaler_instance = scaler_instance, 
            pca_instance = pca_instance, 
            cv_outer = GroupKFold(n_splits = 5), 
            cv_inner = GroupKFold(n_splits = 5), 
            group_col = 'cell', 
            thread_limit = threads)
        res['target'] = 'aac_old'
    else:
        res = nested_train_eval_drugwise(
            model = model, 
            param_grid = param_grid, 
            key_df = data_dict['X_key'], 
            y_col = target, 
            drug_col = 'DRUG', 
            cell_col = 'CELL', 
            cl_exp_df = data_dict['cl_rnaseq'], 
            scaler_instance = scaler_instance, 
            pca_instance = pca_instance, 
            cv_outer = GroupKFold(n_splits = 5), 
            cv_inner = GroupKFold(n_splits = 5), 
            group_col = 'cell', 
            thread_limit = threads)
        res['target'] = target
    
    res.to_csv(f"{root_dir}{result_file}")

if __name__ == '__main__':
    desc_str = 'Command line tool for evaluating several drug sensitivity models. \
    \
    Uses data from Xia et al. 2021 study and models from scikit-learn. '
    #%% Parse command line arguments
    parser = argparse.ArgumentParser(
        prog = 'Baseline drug sensitivity model evaluation', 
        description = desc_str
    )
    parser.add_argument('--root_dir', type = str, default = '')
    parser.add_argument('--result_file', type = str, default = '')
    parser.add_argument('--old_data', action = 'store_true')
    parser.add_argument('--identical_data', action = 'store_true')
    parser.add_argument('--gene_preselection', action = 'store_true')
    parser.add_argument('--multivariate_fs', action = 'store_true')
    parser.add_argument('--response_sources', type = str, default = 'CTRP')
    parser.add_argument('--standardize', action = 'store_true')
    parser.add_argument('--use_pca', action = 'store_true')
    parser.add_argument('--pca_npc', type = int, default = 100)
    parser.add_argument('--model', type = str, default = 'elasticnet')
    parser.add_argument('--target', type = str, default = 'AAC1')
    parser.add_argument('--threads', type = int, default = 8)
    parser.add_argument('--feature_selection_only', action = 'store_true')
    
    kwargs = parser.parse_args()
    
    main(**vars(kwargs))