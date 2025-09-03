#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 09:57:47 2023

@author: rintala
"""

import numpy as np
import pandas as pd
import scipy
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

dsc_spec = importlib.util.spec_from_file_location('dataset_collections', 'datasets/cancer_dataset_collections.py')
dsc = importlib.util.module_from_spec(dsc_spec)
sys.modules['dataset_collections'] = dsc
dsc_spec.loader.exec_module(dsc)

et_spec = importlib.util.spec_from_file_location('evaluation_tools', 'utilities/evaluation_tools.py')
et = importlib.util.module_from_spec(et_spec)
sys.modules['evaluation_tools'] = et
et_spec.loader.exec_module(et)

import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
warnings.simplefilter(action='ignore', category=scipy.stats.ConstantInputWarning)

#%% Data
output_path = os.environ.get('MODAE_OUTPUT_PATH', default = None)
if output_path is None:
    raise ValueError('Please define MODAE_OUTPUT_PATH')
data_root = os.environ.get('MODAE_DATA_PATH', default = None)
if data_root is None:
    raise ValueError('Please define MODAE_DATA_PATH')

res_path = f"{output_path}20250410_random_search/pancan_test/"

prediction_files = glob.glob(res_path + '*test_cv_*predictions*.csv.gz')
ps_prediction_files = glob.glob(res_path + '*ps_cv_*predictions*.csv.gz')

from evaluation_tools import load_results, result_file_dict

prediction_dict = result_file_dict(prediction_files)
ps_prediction_dict = result_file_dict(ps_prediction_files)

incomplete_runs = {}
for k in prediction_dict.keys():
    n_files = np.unique(prediction_dict[k]).shape[0]
    if n_files > 1:
        incomplete_runs[k] = n_files
for k in incomplete_runs.keys():
    _ = prediction_dict.pop(k)

p_df = next(iter(load_results(prediction_dict['1']).values()))
n_classes = p_df.filter(regex = 'class_pred_[0-9]+').shape[1]

#%% Cell-line meta-stasis info 
fn = f"{data_root}ccle/Model_augmented.csv"
cl_info = pd.read_csv(fn, index_col = 'ModelID')

cl_info['PrimaryOrMetastasis'].value_counts()

cl_metas_ptr = cl_info.index[cl_info['PrimaryOrMetastasis'] == 'Metastatic']
cl_prims_ptr = cl_info.index[cl_info['PrimaryOrMetastasis'] == 'Primary']

#%% Define comparison function

from sklearn.metrics import adjusted_rand_score
from scipy.special import softmax

def stability_evaluator(d1, d2):
    preds1 = np.argmax(d1, axis = 1)
    preds2 = np.argmax(d2, axis = 1)
    
    ari = adjusted_rand_score(preds1, preds2)
    
    d1_sm = softmax(d1, axis = 1)
    d2_sm = softmax(d2, axis = 1)
    dot_multi = np.mean(np.sum(d1_sm * d2_sm, axis = 1))
    
    return ari, dot_multi


#%% Within run stability (identical initialization, different data)



within_class_stability_df_list = []

for i in prediction_dict.keys():
    p_df = next(iter(load_results(prediction_dict[i]).values()))
    n_run = np.unique(p_df['run']).shape[0]
    n_fold = np.unique(p_df['fold']).shape[0]
    
    for run1 in np.arange(n_run):
        for fold1 in np.arange(n_fold - 1):
            ind1 = np.logical_and(p_df['fold'] == fold1, p_df['run'] == run1)
            class_preds1 = p_df.loc[ind1]
            class_preds1.set_index(p_df.loc[ind1, 'id'], inplace = True)
            class_preds1_grouped = class_preds1.groupby(
                by = 'dataset', 
                as_index = False, 
                sort = False
            )
            for fold2 in np.arange(fold1+1, n_fold):
                ind2 = np.logical_and(p_df['fold'] == fold2, p_df['run'] == run1)
                class_preds2 = p_df.loc[ind2]
                class_preds2.set_index(p_df.loc[ind2, 'id'], inplace = True)
                class_preds2_grouped = class_preds2.groupby(
                    by = 'dataset', 
                    as_index = False, 
                    sort = False
                )
                for (data_name1, dataset1), (data_name2, dataset2) in zip(class_preds1_grouped, class_preds2_grouped):
                    if data_name1 != data_name2:
                        raise ValueError('Pandas groupby order mismatch.')
                    
                    common_idx = dataset1.index.intersection(dataset2.index)
                    
                    if common_idx.shape[0] > 0:
                        dataset1 = dataset1.loc[common_idx, [f"class_pred_{i}" for i in np.arange(n_classes, dtype = 'int')]]
                        dataset2 = dataset2.loc[common_idx, [f"class_pred_{i}" for i in np.arange(n_classes, dtype = 'int')]]
                        
                        ari, dot_multi = stability_evaluator(
                            dataset1.to_numpy(), 
                            dataset2.to_numpy()
                        )
                        
                        within_class_stability_df_list.append(pd.DataFrame({
                            'task' : [i],
                            'run' : [run1], 
                            'fold1' : [fold1], 
                            'fold2' : [fold2],
                            'dataset' : [data_name1],
                            'ari' : [ari], 
                            'soft_max_mean_dot' : [dot_multi]
                        }))
                        
                        if data_name1 == 'cl_train':
                            d1_metas_ptr = dataset1.index.intersection(cl_metas_ptr)
                            d2_metas_ptr = dataset2.index.intersection(cl_metas_ptr)
                            if np.any(d1_metas_ptr != d2_metas_ptr):
                                raise ValueError('Index mismatch')
                            meta_ari, meta_dot_multi = stability_evaluator(
                                dataset1.loc[d1_metas_ptr].to_numpy(), 
                                dataset2.loc[d2_metas_ptr].to_numpy()
                            )
                            within_class_stability_df_list.append(pd.DataFrame({
                                'task' : [i],
                                'run' : [run1], 
                                'fold1' : [fold1], 
                                'fold2' : [fold2],
                                'dataset' : [f"{data_name1}_metastatic"],
                                'ari' : [meta_ari], 
                                'soft_max_mean_dot' : [meta_dot_multi]
                            }))
                            
                            d1_prims_ptr = dataset1.index.intersection(cl_prims_ptr)
                            d2_prims_ptr = dataset2.index.intersection(cl_prims_ptr)
                            if np.any(d1_prims_ptr != d2_prims_ptr):
                                raise ValueError('Index mismatch')
                            primary_ari, primary_dot_multi = stability_evaluator(
                                dataset1.loc[d1_prims_ptr].to_numpy(), 
                                dataset2.loc[d2_prims_ptr].to_numpy()
                            )
                            within_class_stability_df_list.append(pd.DataFrame({
                                'task' : [i],
                                'run' : [run1], 
                                'fold1' : [fold1], 
                                'fold2' : [fold2],
                                'dataset' : [f"{data_name1}_primary"],
                                'ari' : [primary_ari], 
                                'soft_max_mean_dot' : [primary_dot_multi]
                            }))
                        
within_class_stability_df = pd.concat(within_class_stability_df_list, axis = 0)

within_class_stability_df.to_csv(f"{res_path}tissue_classifier_stability_within_run.csv.gz")

#%% Between run stability (different initialization, different_data)

if False:
    raise NotImplementedError('This has not been implemented fully')
    
    between_class_stability_df_list = []
    if n_run > 1:
        for i in prediction_dict.keys():
            p_df = next(iter(load_results(prediction_dict[i]).values()))
            n_run = np.unique(p_df['run']).shape[0]
            n_fold = np.unique(p_df['fold']).shape[0]
            for run1 in np.arange(n_run):
                for run2 in np.arange(n_run):
                    # Combine test folds for more precise evaluation
                    for data_name in ['cl_test', 'patient_test']:
                        next
                    # Loop through train folds
                    for data_name in ['cl_train', 'patient_train']:
                        for fold1 in np.arange(n_fold):
                            ind1 = np.logical_and(p_df['fold'] == fold1, p_df['run'] == run1)
                            class_preds1 = p_df.loc[ind1]
                            class_preds1.set_index(p_df.loc[ind1, 'id'], inplace = True)
                            class_preds1_grouped = class_preds1.groupby(
                                by = 'dataset', 
                                as_index = False, 
                                sort = False
                            )
                            ind2 = np.logical_and(p_df['fold'] == fold2, p_df['run'] == run2)
                            class_preds2 = p_df.loc[ind2]
                            class_preds2.set_index(p_df.loc[ind2, 'id'], inplace = True)
                            class_preds2_grouped = class_preds2.groupby(
                                by = 'dataset', 
                                as_index = False, 
                                sort = False
                            )
                            for (data_name1, dataset1), (data_name2, dataset2) in zip(class_preds1_grouped, class_preds2_grouped):
                                if data_name1 != data_name2:
                                    raise ValueError('Pandas groupby order mismatch.')
                                common_idx = dataset1.index.intersection(dataset2.index)
                                
                                if common_idx.shape[0] > 0:
                                    dataset1 = dataset1.loc[common_idx, [f"class_pred_{i}" for i in np.arange(n_classes, dtype = 'int')]]
                                    dataset2 = dataset2.loc[common_idx, [f"class_pred_{i}" for i in np.arange(n_classes, dtype = 'int')]]
                                    
                                    preds1 = np.argmax(dataset1.to_numpy(), axis = 1)
                                    preds2 = np.argmax(dataset2.to_numpy(), axis = 1)
                                    
                                    ari = adjusted_rand_score(preds1, preds2)
                                    
                                    d1_sm = scipy.special.softmax(dataset1.to_numpy(), axis = 1)
                                    d2_sm = scipy.special.softmax(dataset2.to_numpy(), axis = 1)
                                    dot_multi = np.mean(np.sum(d1_sm * d2_sm, axis = 1))
                                    
                                    between_class_stability_df_list.append(pd.DataFrame({
                                        'task' : [i],
                                        'run' : [run1], 
                                        'fold1' : [fold1], 
                                        'fold2' : [fold2],
                                        'dataset' : [data_name1],
                                        'ari' : [ari], 
                                        'soft_max_mean_dot' : [dot_multi]
                                    }))
                                    
                                    if data_name1 == 'cl_train':
                                        d1_metas_ptr = dataset1.index.intersection(cl_metas_ptr)
                                        d2_metas_ptr = dataset2.index.intersection(cl_metas_ptr)
                                        if np.any(d1_metas_ptr != d2_metas_ptr):
                                            raise ValueError('Index mismatch')
                                        meta_ari, meta_dot_multi = stability_evaluator(
                                            dataset1.loc[d1_metas_ptr].to_numpy(), 
                                            dataset2.loc[d2_metas_ptr].to_numpy()
                                        )
                                        between_class_stability_df_list.append(pd.DataFrame({
                                            'task' : [i],
                                            'run' : [run1], 
                                            'fold1' : [fold1], 
                                            'fold2' : [fold2],
                                            'dataset' : [f"{data_name1}_metastatic"],
                                            'ari' : [meta_ari], 
                                            'soft_max_mean_dot' : [meta_dot_multi]
                                        }))
                                        
                                        d1_prims_ptr = dataset1.index.intersection(cl_prims_ptr)
                                        d2_prims_ptr = dataset2.index.intersection(cl_prims_ptr)
                                        if np.any(d1_prims_ptr != d2_prims_ptr):
                                            raise ValueError('Index mismatch')
                                        primary_ari, primary_dot_multi = stability_evaluator(
                                            dataset1.loc[d1_prims_ptr].to_numpy(), 
                                            dataset2.loc[d2_prims_ptr].to_numpy()
                                        )
                                        between_class_stability_df_list.append(pd.DataFrame({
                                            'task' : [i],
                                            'run' : [run1], 
                                            'fold1' : [fold1], 
                                            'fold2' : [fold2],
                                            'dataset' : [f"{data_name1}_primary"],
                                            'ari' : [primary_ari], 
                                            'soft_max_mean_dot' : [primary_dot_multi]
                                        }))
                                            
        between_class_stability_df = pd.concat(between_class_stability_df_list, axis = 0)
        between_class_stability_df.to_csv(f"{res_path}tissue_classifier_stability_between_run.csv.gz")

#%% Between task stability
if False:
    class_df_list = []
    
    tasks = list(prediction_dict.keys())
    
    for i in np.arange(len(tasks)-1):
        p1_df = next(iter(load_results(prediction_dict[tasks[i]]).values()))
        n_run = np.unique(p1_df['run']).shape[0]
        n_fold = np.unique(p1_df['fold']).shape[0]
        
        for j in np.arange(i+1, len(tasks)):
            p2_df = next(iter(load_results(prediction_dict[tasks[j]]).values()))
            
            for run in np.arange(n_run):
                for fold in np.arange(n_fold):
                    ind1 = np.logical_and(p1_df['fold'] == fold, p1_df['run'] == run)
                    class_preds1 = p1_df.loc[ind1]
                    class_preds1.set_index(p1_df.loc[ind1, 'id'], inplace = True)
                    class_preds1_grouped = class_preds1.groupby(
                        by = 'dataset', 
                        as_index = False, 
                        sort = False
                    )
                    ind2 = np.logical_and(p2_df['fold'] == fold, p2_df['run'] == run)
                    class_preds2 = p2_df.loc[ind2]
                    class_preds2.set_index(p2_df.loc[ind2, 'id'], inplace = True)
                    class_preds2_grouped = class_preds2.groupby(
                        by = 'dataset', 
                        as_index = False, 
                        sort = False
                    )
                    for (data_name1, dataset1), (data_name2, dataset2) in zip(class_preds1_grouped, class_preds2_grouped):
                        if data_name1 != data_name2:
                            raise ValueError('Pandas groupby order mismatch.')
                        common_idx = dataset1.index.intersection(dataset2.index)
                        
                        if common_idx.shape[0] > 0:
                            dataset1 = dataset1.loc[common_idx, [f"class_pred_{i}" for i in np.arange(n_classes, dtype = 'int')]]
                            dataset2 = dataset2.loc[common_idx, [f"class_pred_{i}" for i in np.arange(n_classes, dtype = 'int')]]
                            
                            preds1 = np.argmax(dataset1.to_numpy(), axis = 1)
                            preds2 = np.argmax(dataset2.to_numpy(), axis = 1)
                            
                            ari = adjusted_rand_score(preds1, preds2)
                            
                            d1_sm = scipy.special.softmax(dataset1.to_numpy(), axis = 1)
                            d2_sm = scipy.special.softmax(dataset2.to_numpy(), axis = 1)
                            dot_multi = np.mean(np.sum(d1_sm * d2_sm, axis = 1))
                            
                            class_df_list.append(pd.DataFrame({
                                'task1' : [tasks[i]],
                                'task2' : [tasks[j]],
                                'run' : [run], 
                                'fold' : [fold], 
                                'dataset' : [data_name1],
                                'ari' : [ari], 
                                'soft_max_mean_dot' : [dot_multi]
                            }))
                            
                            if data_name1 == 'cl_train':
                                d1_metas_ptr = dataset1.index.intersection(cl_metas_ptr)
                                d2_metas_ptr = dataset2.index.intersection(cl_metas_ptr)
                                if np.any(d1_metas_ptr != d2_metas_ptr):
                                    raise ValueError('Index mismatch')
                                meta_ari, meta_dot_multi = stability_evaluator(
                                    dataset1.loc[d1_metas_ptr].to_numpy(), 
                                    dataset2.loc[d2_metas_ptr].to_numpy()
                                )
                                class_df_list.append(pd.DataFrame({
                                    'task1' : [tasks[i]],
                                    'task2' : [tasks[j]],
                                    'run' : [run], 
                                    'fold' : [fold], 
                                    'dataset' : [f"{data_name1}_metastatic"],
                                    'ari' : [meta_ari], 
                                    'soft_max_mean_dot' : [meta_dot_multi]
                                }))
                                
                                d1_prims_ptr = dataset1.index.intersection(cl_prims_ptr)
                                d2_prims_ptr = dataset2.index.intersection(cl_prims_ptr)
                                if np.any(d1_prims_ptr != d2_prims_ptr):
                                    raise ValueError('Index mismatch')
                                primary_ari, primary_dot_multi = stability_evaluator(
                                    dataset1.loc[d1_prims_ptr].to_numpy(), 
                                    dataset2.loc[d2_prims_ptr].to_numpy()
                                )
                                class_df_list.append(pd.DataFrame({
                                    'task1' : [tasks[i]],
                                    'task2' : [tasks[j]],
                                    'run' : [run], 
                                    'fold' : [fold], 
                                    'dataset' : [f"{data_name1}_primary"],
                                    'ari' : [primary_ari], 
                                    'soft_max_mean_dot' : [primary_dot_multi]
                                }))
    
    class_df = pd.concat(class_df_list, axis = 0)
    class_df.to_csv(f"{res_path}tissue_classifier_stability_between_task.csv.gz")