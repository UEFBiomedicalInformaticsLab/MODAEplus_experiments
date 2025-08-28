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

et_spec = importlib.util.spec_from_file_location('evaluation_tools', 'utilities/evaluation_tools.py')
et = importlib.util.module_from_spec(et_spec)
sys.modules['evaluation_tools'] = et
et_spec.loader.exec_module(et)

import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

n_workers = 100

#%% Data
from sys import platform
from socket import gethostname
if gethostname() == 'teemu-pc':
    base_path = '/home/teemu/research_work/'
elif platform == 'linux':
    base_path = '/research/work/rintala/'
else:
    base_path = '//research/workdir/'

#res_path = base_path + 'superAE_HPO/20230821_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231027_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231101_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231103_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/random_search_231113/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/231117_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231120_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231121_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231122_random_search/scanb_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231124_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20231128_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231129_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231130_random_search_nostandard/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231201_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231201_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20231218_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231218_random_search_relu/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231219_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231219_random_search_relu/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231219_random_search_variational/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231219_random_search_variational_relu/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231220_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231220_random_search_no_joint/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231221_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231221_random_search_no_joint/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240103_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240108_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240115_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240118_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240119_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240123_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240124_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240124_random_search/brca_test_noclfilter_alternative/'
#res_path = base_path + 'superAE_HPO/20240125_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240125_random_search/brca_test_noclfilter_nopre/'
#res_path = base_path + 'superAE_HPO/20240206_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240208_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240211_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240213_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240213_random_search/brca_test_noclfilter_pretrain/'
#res_path = base_path + 'superAE_HPO/20240215_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20240216_random_search/brca_test_noclfilter_pretrain/'
#res_path = base_path + 'superAE_HPO/20240221_random_search/scanb_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240301_random_search/scanb_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240302_random_search/scanb_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240303_random_search/scanb_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240305_random_search/scanb_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240305_2_random_search/scanb_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240306_random_search/scanb_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240310_random_search/scanb_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240315_random_search/scanb_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240319_random_search/scanb_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240328_random_search/scanb_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240430_random_search/scanb_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240605_random_search/scanb_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240610_random_search/scanb_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240614_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20240729_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20240816_random_search/scanb_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20241210_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20241218_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20241231_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20250108_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20250110_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20250117_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20250120_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20250121_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20250123_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20250124_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20250127_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20250128_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20250202_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20250205_random_search/scanb_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20250207_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20250216_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20250217_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20250219_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20250223_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20250313_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20250402_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20250404_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20250408_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20250410_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20250415_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20250506_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20250520_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20250525_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20250527_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20250530_random_search/pancan_test/'
res_path = base_path + 'superAE_HPO/20250603_random_search/pancan_test/'

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
    from sae.parsing_utilities import modae_args_json_decoder
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
from sae.evaluation import DSC

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
