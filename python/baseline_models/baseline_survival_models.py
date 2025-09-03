# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 14:46:30 2023

@author: rintala
"""

import numpy as np
import pandas as pd

from copy import copy

from sklearn.model_selection import RepeatedKFold, ParameterGrid
from threadpoolctl import threadpool_limits, threadpool_info
from sklearn.utils._testing import ignore_warnings
from sklearn.exceptions import ConvergenceWarning

#cv_outer = GroupKFold(n_splits = 5)
#cv_inner = GroupKFold(n_splits = 5)

@ignore_warnings(category=ConvergenceWarning)
def nested_survival_eval(
        model, 
        param_grid, 
        X, 
        y, 
        y_mask,
        X_covar = None, 
        strat_var = None, 
        scaler_instance = None,
        covar_scaler_instance = None,
        pca_instance = None, 
        cv_outer = RepeatedKFold(n_splits = 5, n_repeats = 2), 
        cv_inner = RepeatedKFold(n_splits = 5, n_repeats = 1), 
        thread_limit = 1):
    cv_results = []
    if strat_var is not None:
        cv_outer_args = {'groups' : strat_var}
    else:
        cv_outer_args = {}
    for train_ind, test_ind in cv_outer.split(X, y, **cv_outer_args):
        if scaler_instance is not None:
            scaler = copy(scaler_instance)
            X_train = scaler.fit_transform(X[train_ind,:])
            X_test = scaler.transform(X[test_ind,:])
        else:
            X_train = X[train_ind,:]
            X_test = X[test_ind,:]
        if pca_instance is not None:
            pca = copy(pca_instance)
            with threadpool_limits(limits = thread_limit, user_api = 'blas'):
                pcs_train = pca.fit_transform(X_train)
                pcs_test = pca.transform(X_test)
            X_train = pcs_train[y_mask[train_ind], :]
            X_test = pcs_test[y_mask[test_ind], :]
        else:
            X_train = X_train[y_mask[train_ind], :]
            X_test = X_test[y_mask[test_ind], :]
        
        if X_covar is not None:
            if covar_scaler_instance is not None:
                covar_scaler = copy(covar_scaler_instance)
                X_covar_train = covar_scaler.fit_transform(X_covar[train_ind[y_mask[train_ind]],:])
                X_covar_test = covar_scaler.transform(X_covar[test_ind[y_mask[test_ind]],:])
            else:
                X_covar_train = X_covar[train_ind[y_mask[train_ind]],:]
                X_covar_test = X_covar[test_ind[y_mask[test_ind]],:]
            X_train = np.concatenate((X_train, X_covar_train), axis = 1)
            X_test = np.concatenate((X_test, X_covar_test), axis = 1)
        
        y_train = y[train_ind][y_mask[train_ind]]
        y_test = y[test_ind][y_mask[test_ind]]
        
        if strat_var is not None:
            cv_inner_args = {'groups' : strat_var[train_ind][y_mask[train_ind]]}
        else:
            cv_inner_args = {}
        
        grid_search_score_d = grid_search(
            model = model, 
            X = X_train, 
            y = y_train, 
            param_grid = param_grid, 
            cv = cv_inner, 
            cv_args = cv_inner_args, 
            thread_limit = thread_limit)
        best_i = np.argmax(grid_search_score_d['test_score'])
        best_par_dict = dict([(k, grid_search_score_d[k][best_i]) for k in grid_search_score_d.keys()])
        best_par_train_score = best_par_dict.pop('train_score')
        best_par_test_score = best_par_dict.pop('test_score')
        
        model_best = copy(model)
        model_best = model_best.set_params(**best_par_dict)
        with threadpool_limits(limits = thread_limit, user_api = 'blas'):
            model_best.fit(X_train, y_train)
        y_c_train = model_best.score(X_train, y_train)
        y_c_test = model_best.score(X_test, y_test)
        
        cv_results.append(pd.DataFrame({
            'fold' : len(cv_results), 
            **best_par_dict, 
            'train_c' : y_c_train, 
            'test_c' : y_c_test, 
            'inner_train_c' : best_par_train_score, 
            'inner_test_c' : best_par_test_score}, 
            index = pd.RangeIndex(0,1)))
    
    cv_results_df = pd.concat(cv_results, axis = 0)
    
    return cv_results_df

def grid_search(model, X, y, param_grid, cv, cv_args = {}, thread_limit = 1):
    grid_search_scores = []
    for par_dict in ParameterGrid(param_grid):
        res_train = []
        res_test = []
        for train_ind, test_ind in cv.split(X, y, **cv_args):
            X_train = X[train_ind, :]
            X_test = X[test_ind, :]
            y_train = y[train_ind]
            y_test = y[test_ind]
            
            model_i = copy(model)
            model_i = model_i.set_params(**par_dict)
            try:
                with threadpool_limits(limits = thread_limit, user_api = 'blas'):
                    model_i.fit(X_train, y_train)
                score_train = model_i.score(X_train, y_train)
                score_test = model_i.score(X_test, y_test)
            except:
                score_train = np.nan
                score_test = np.nan
            res_train.append(score_train)
            res_test.append(score_test)
        grid_search_scores.append({
            **par_dict, 
            'train_score' : np.mean(res_train), 
            'test_score' : np.mean(res_test)})
    
    grid_search_score_d = {}
    for key in grid_search_scores[0].keys():
        grid_search_score_d[key] = [i[key] for i in grid_search_scores]
    return grid_search_score_d
