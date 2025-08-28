#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 15:40:12 2023

@author: teemu
"""

import numpy as np
import pandas as pd

from copy import copy

from sklearn.model_selection import GroupKFold, ParameterGrid
from sklearn.metrics import mean_squared_error
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
        scaler_instance = None, 
        cv_outer = GroupKFold(n_splits = 5), 
        cv_inner = GroupKFold(n_splits = 5), 
        group_col = 'cell', 
        thread_limit = 1, 
        return_predictions = False):
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
        
        X_train = cl_exp_df.loc[cl_train_names].to_numpy()
        X_test = cl_exp_df.loc[cl_test_names].to_numpy()
        
        if scaler_instance is not None:
            scaler = copy(scaler_instance)
            X_train = scaler.fit_transform(X_train)
            X_test = scaler.transform(X_test)
        
        if pca_instance is not None:
            cl_pca = copy(pca_instance)
            cl_pcs_train = cl_pca.fit_transform(X_train)
            cl_pcs_test = cl_pca.transform(X_test)
        
        key_df_train = key_df.iloc[train_ind]
        key_df_test = key_df.iloc[test_ind]
        
        drug_names = key_df[drug_col].unique()
        for dn in drug_names:
            key_df_train_di = key_df_train.loc[key_df_train[drug_col] == dn]
            key_df_test_di = key_df_test.loc[key_df_test[drug_col] == dn]
            
            cl_train_ind = cl_train_map.loc[key_df_train_di[cell_col]].to_numpy()
            cl_test_ind = cl_test_map.loc[key_df_test_di[cell_col]].to_numpy()
            
            if pca_instance is not None:
                X_cl_train = cl_pcs_train[cl_train_ind, :]
                X_cl_test = cl_pcs_test[cl_test_ind, :]
            else:
                X_cl_train = X_train[cl_train_ind, :]
                X_cl_test = X_test[cl_test_ind, :]
            
            y_cl_train = key_df_train_di[y_col].to_numpy()
            y_cl_test = key_df_test_di[y_col].to_numpy()
            
            grid_search_scores = []
            if group_col == 'cell':
                cv_inner_args = {'groups' : key_df_train_di[cell_col].to_numpy()}
            else:
                cv_inner_args = {}
            for par_dict in ParameterGrid(param_grid):
                res = []
                for train_ind_inner, test_ind_inner in cv_inner.split(key_df_train_di, **cv_inner_args):
                    X_train_inner = X_cl_train[train_ind_inner, :]
                    X_test_inner = X_cl_train[test_ind_inner, :]
                    y_train_inner = y_cl_train[train_ind_inner]
                    y_test_inner = y_cl_train[test_ind_inner]
                    
                    model_inner = copy(model)
                    model_inner = model_inner.set_params(**par_dict)
                    with threadpool_limits(limits = thread_limit, user_api = 'blas'):
                        model_inner.fit(X_train_inner, y_train_inner)
                    y_r2_train_inner = model_inner.score(X_train_inner, y_train_inner)
                    y_r2_test_inner = model_inner.score(X_test_inner, y_test_inner)
                    y_mse_train_inner = mean_squared_error(model_inner.predict(X_train_inner), y_train_inner)
                    y_mse_test_inner = mean_squared_error(model_inner.predict(X_test_inner), y_test_inner)
                    
                    res.append(pd.DataFrame({
                        'train_r2' : y_r2_train_inner, 
                        'test_r2' : y_r2_test_inner, 
                        'train_mse' : y_mse_train_inner, 
                        'test_mse' : y_mse_test_inner}, 
                        index = pd.RangeIndex(0,1)))
                res_df = pd.concat(res, axis = 0)
                grid_search_scores.append(pd.DataFrame({
                    **par_dict, 
                    'train_r2' : res_df['train_r2'].mean(), 
                    'test_r2' : res_df['test_r2'].mean(), 
                    'train_mse' : res_df['train_mse'].mean(), 
                    'test_mse' : res_df['test_mse'].mean()}, 
                    index = pd.RangeIndex(0,1)))
            
            grid_search_scores_df = pd.concat(grid_search_scores, axis = 0)
            best_params = grid_search_scores_df.astype('object').iloc[grid_search_scores_df['test_r2'].argmax()]
            best_par_dict = best_params[param_grid.keys()].to_dict()
            
            model_best = copy(model)
            model_best = model_best.set_params(**best_par_dict)
            with threadpool_limits(limits = thread_limit, user_api = 'blas'):
                model_best.fit(X_cl_train, y_cl_train)
            y_r2_train = model_best.score(X_cl_train, y_cl_train)
            y_r2_test = model_best.score(X_cl_test, y_cl_test)
            y_mse_train = mean_squared_error(model_best.predict(X_cl_train),  y_cl_train)
            y_mse_test = mean_squared_error(model_best.predict(X_cl_test), y_cl_test)
            
            cv_results.append(pd.DataFrame({
                'drug' : dn, 
                'fold' : len(cv_results), 
                **best_par_dict, 
                'train_r2' : y_r2_train, 
                'test_r2' : y_r2_test, 
                'train_mse' : y_mse_train, 
                'test_mse' : y_mse_test, 
                'inner_train_r2' : best_params['train_r2'], 
                'inner_test_r2' : best_params['test_r2'], 
                'inner_train_mse' : best_params['train_mse'], 
                'inner_test_mse' : best_params['test_mse']}, 
                index = pd.RangeIndex(0,1)))
    
    if return_predictions:
        pass
    
    cv_results_df = pd.concat(cv_results, axis = 0)
    
    return cv_results_df