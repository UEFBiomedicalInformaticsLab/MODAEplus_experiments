#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 09:57:47 2023

@author: rintala
"""

import numpy as np
import pandas as pd
import os
import glob
import re

import seaborn as sns

from sklearn import neighbors
from copy import copy
from sklearn.metrics import silhouette_score
from sklearn.metrics.pairwise import pairwise_distances

import sys
import importlib.util

script_path = os.environ.get('MODAE_SCRIPT_PATH', default = None)
if script_path is None:
    raise ValueError('Please define MODAE_SCRIPT_PATH')

et_spec = importlib.util.spec_from_file_location(
    'evaluation_tools', 
    f"{script_path}/python/utilities/evaluation_tools.py"
)
et = importlib.util.module_from_spec(et_spec)
sys.modules['evaluation_tools'] = et
et_spec.loader.exec_module(et)

import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

n_workers = 100

#%% Data
output_path = os.environ.get('MODAE_OUTPUT_PATH', default = None)
if output_path is None:
    raise ValueError('Please define MODAE_OUTPUT_PATH')
data_root = os.environ.get('MODAE_DATA_PATH', default = None)
if data_root is None:
    raise ValueError('Please define MODAE_DATA_PATH')

res_path = f"{output_path}20250410_random_search/pancan_test/"

embedding_files = glob.glob(res_path + '*test_cv_*embeddings*.csv.gz')
ps_embedding_files = glob.glob(res_path + '*ps_cv_*embeddings*.csv.gz')

from evaluation_tools import load_results

embeddings_dict = load_results(embedding_files)
ps_embeddings_dict = load_results(ps_embedding_files)

#%% Parameter dicts

parameter_files = glob.glob(res_path + '*parameters*.csv')

if len(parameter_files):
    from evaluation_tools import get_kwargs
    parameters_list = [get_kwargs(i) for i  in parameter_files]
    parameters_dict = dict([(i['task_id'], i) for i in parameters_list])
else:
    parameter_files = glob.glob(res_path + '*parameters*.json')
    from modae.parsing_utilities import modae_args_json_decoder
    import json
    parameters_list = []
    for fn in parameter_files:    
        param_handle = open(fn)
        search_kwargs = json.load(param_handle)
        param_handle.close()
        search_kwargs = modae_args_json_decoder(search_kwargs)
        parameters_list.append(search_kwargs)
    parameters_dict = dict([(i['task_id'], i) for i in parameters_list])

#dict([(i, ) for i in embeddings_dict.keys()])

#%%
from sklearn.svm import LinearSVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV, cross_validate, KFold
from sklearn import metrics
from sklearn import preprocessing
from multiprocessing import Pool
from modae.evaluation import DSC

def mapper(e_df, par_d):
    i = e_df['task'].iloc[0]
    encoder_layers = par_d['model_args']['encoder_layers']
    encoder_depth = encoder_layers[len(encoder_layers)-1]
    shared_embeddings = np.array([f"z{j}" for j in np.arange(encoder_depth)+1])
    batch_df_list = []
    for fold in np.unique(e_df['fold']):
        for run in np.unique(e_df['run']):
            ind = np.logical_and(e_df['fold'] == fold, e_df['run'] == run)
            X = e_df.loc[ind, shared_embeddings].to_numpy()
            y_raw = e_df.loc[ind, 'dataset'].to_numpy()
            y_map = {
                'cl_train' : 1, 
                'cl_test' : 1, 
                'patient_train' : 0, 
                'patient_test' : 0}
            y = np.array([y_map[j] for j in y_raw])
            #y_oh = preprocessing.OneHotEncoder().fit_transform(y)
            
            silh = metrics.silhouette_score(X, y, metric = 'euclidean')
            dsc = DSC(X, y).numpy()
            
            eval_cv = KFold(n_splits = 5, shuffle = True, random_state = 0)
            
            svm = LinearSVC(
                penalty = 'l2', 
                loss = 'squared_hinge', 
                dual = False, 
                tol = 1e-8, 
                C = 1., 
                fit_intercept = True)
            svm = svm.fit(X, y)
            metrics_list = [
                'average_precision', 
                'f1', 
                'balanced_accuracy', 
                'roc_auc', 
                'accuracy'
            ]
            svm_cv = cross_validate(svm, X, y, scoring = metrics_list, cv = eval_cv)
            
            
            rfc = RandomForestClassifier(n_estimators = 100)
            grid_cv = KFold(n_splits = 5, shuffle = True, random_state = 0)
            rfc_grid = GridSearchCV(
                rfc, 
                param_grid = {'ccp_alpha' : np.logspace(-4, -1, num = 10)}, 
                cv = grid_cv, 
                scoring = 'f1')
            rfc_cv = cross_validate(rfc_grid, X, y, scoring = metrics_list, cv = eval_cv)
            
            batch_dfi = pd.DataFrame({
                'task' : i, 
                'run' : run, 
                'fold' : fold, 
                'silh' : silh, 
                'dsc' : dsc, 
                'svm_apr' : np.mean(svm_cv['test_average_precision']), 
                'svm_f1' : np.mean(svm_cv['test_f1']), 
                'svm_bacc' : np.mean(svm_cv['test_balanced_accuracy']), 
                'svm_acc' : np.mean(svm_cv['test_accuracy']), 
                'svm_auc' : np.mean(svm_cv['test_roc_auc']), 
                'rfc_apr' : np.mean(rfc_cv['test_average_precision']), 
                'rfc_f1' : np.mean(rfc_cv['test_f1']), 
                'rfc_bacc' : np.mean(rfc_cv['test_balanced_accuracy']), 
                'rfc_acc' : np.mean(rfc_cv['test_accuracy']), 
                'rfc_auc' : np.mean(rfc_cv['test_roc_auc'])}, 
                index = pd.RangeIndex(0,1))
            
            batch_df_list.append(batch_dfi)
    
    return pd.concat(batch_df_list, axis = 0)

batch_df_list = []
with Pool(processes = n_workers, maxtasksperchild = 1) as mp:
    for res in mp.starmap(
        mapper, 
        [(embeddings_dict[i], parameters_dict[i]) for i in embeddings_dict.keys()]
    ):
        batch_df_list.append(res)
batch_df = pd.concat(batch_df_list, axis = 0)
batch_df.to_csv(res_path + 'batch_predictions.csv')

#%%

ps_batch_df_list = []
with Pool(processes = n_workers, maxtasksperchild = 1) as mp:
    for res in mp.starmap(
        mapper, 
        [(ps_embeddings_dict[i], parameters_dict[i]) for i in ps_embeddings_dict.keys()]
    ):
        ps_batch_df_list.append(res)
ps_batch_df = pd.concat(ps_batch_df_list, axis = 0)
ps_batch_df.to_csv(res_path + 'ps_batch_predictions.csv')
