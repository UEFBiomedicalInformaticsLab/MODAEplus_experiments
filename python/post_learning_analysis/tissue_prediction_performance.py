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
from sys import platform
from socket import gethostname
if gethostname() == 'teemu-pc':
    base_path = '/home/teemu/research_work/'
    home = True
elif platform == 'linux':
    base_path = '/research/work/rintala/'
    home = False
else:
    base_path = '//research/workdir/'
    home = False

n_workers = 10

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

run_date = re.sub('[a-z_]*/[a-z_]*/$', '', res_path)
run_date = re.sub('^.*/', '', run_date)
dss_sensitivity = int(run_date) >= 20250313
include_metas_labels = int(run_date) >= 20250404 and int(run_date) <= 20250410
include_metas_labels = include_metas_labels or (int(run_date) >= 20250506 and int(run_date) <= 20250520) or (int(run_date) >= 20250603 and int(run_date) <= 20250603)
exclude_metas_data = int(run_date) >= 20250527 and int(run_date) <= 20250527

#embedding_files = glob.glob(res_path + '*test_cv_*embeddings*.csv.gz')
prediction_files = glob.glob(res_path + '*test_cv_*predictions*.csv.gz')
ps_prediction_files = glob.glob(res_path + '*ps_cv_*predictions*.csv.gz')

if False:
    prediction_files = prediction_files[:10]
    ps_prediction_files = ps_prediction_files[:10]

if False:
    complete_preds = [i for i in prediction_files if 'fold' not in i]
    complete_ps_preds = [i for i in ps_prediction_files if 'fold' not in i]
    cpred_tasks = [re.sub('.*?task', '', i) for i in complete_preds]
    cpred_tasks = [re.sub('\\.csv\\.gz', '', i) for i in cpred_tasks]
    
    cpspred_tasks = [re.sub('.*?task', '', i) for i in complete_ps_preds]
    cpspred_tasks = [re.sub('\\.csv\\.gz', '', i) for i in cpspred_tasks]
    
    full_preds = set(cpspred_tasks).intersection(cpred_tasks)
    incomplete_preds = set([str(i) for i in range(1,101)]).difference(full_preds)

from evaluation_tools import load_results, result_file_dict

#embeddings_dict = result_file_dict(embedding_files)
prediction_dict = result_file_dict(prediction_files)
ps_prediction_dict = result_file_dict(ps_prediction_files)

#%% CTRP or CCLE?
prediction_test = next(iter(load_results(next(iter(prediction_dict.values()))).values())).filter(regex = 'dr_pred_[0-9]+')
if prediction_test.shape[1] == 24:
    response_dataset = 'ccle'
elif prediction_test.shape[1] == 545:
    response_dataset = 'ctrp'
elif prediction_test.shape[1] == 544:
    response_dataset = 'ctrp' # Xia data
else:
    raise ValueError('Drug response prediction column number does not match known datasets.')
#%% Measured data
date = int(re.search('[0-9]+', res_path).group(0))
drug_response_scaling = (date == 20231219)

if re.search('brca', res_path) is not None:
    if response_dataset == 'ccle':
        from dataset_collections import get_tcga_brca_ccle_full
        data_dict = get_tcga_brca_ccle_full(
            home = home, 
            drug_response_maxscale = drug_response_scaling, 
            tissue_classifier = True, 
            include_metas_labels = include_metas_labels
        )
    elif response_dataset == 'ctrp':
        from dataset_collections import get_tcga_brca_ctrp_ccle_full
        data_dict = get_tcga_brca_ctrp_ccle_full(
            home = home, 
            drug_response_maxscale = drug_response_scaling, 
            tissue_classifier = True, 
            include_metas_labels = include_metas_labels
        )
elif re.search('pancan', res_path) is not None:
    if response_dataset == 'ccle':
        from dataset_collections import get_tcga_pancan_ccle_solid
        data_dict = get_tcga_pancan_ccle_solid(
            home = home, 
            drug_response_maxscale = drug_response_scaling, 
            tissue_classifier = True, 
            include_metas_labels = include_metas_labels
        )
    elif response_dataset == 'ctrp':
        if dss_sensitivity:
            from dataset_collections import get_xia_ctrp_data
            data_dict = get_xia_ctrp_data(
                home = home, 
                tissue_classifier = True, 
                include_metas_labels = include_metas_labels, 
                exclude_metas_data = exclude_metas_data
            )
        else:
            from dataset_collections import get_tcga_pancan_ctrp_ccle_solid
            data_dict = get_tcga_pancan_ctrp_ccle_solid(
                home = home, 
                drug_response_maxscale = drug_response_scaling, 
                tissue_classifier = True, 
                include_metas_labels = include_metas_labels
            )
elif re.search('scanb', res_path) is not None:
    if response_dataset == 'ccle':
        from dataset_collections import get_scanb_ccle_full
        data_dict = get_scanb_ccle_full(
            home = home, 
            drug_response_maxscale = drug_response_scaling, 
            tissue_classifier = True, 
            include_metas_labels = include_metas_labels
        )
    elif response_dataset == 'ctrp':
        from dataset_collections import get_scanb_ctrp_ccle_full
        data_dict = get_scanb_ctrp_ccle_full(
            home = home, 
            drug_response_maxscale = drug_response_scaling, 
            tissue_classifier = True, 
            include_metas_labels = include_metas_labels
        )

patient_class_table = pd.DataFrame(
    copy(data_dict['patient_class'])[data_dict['patient_class_mask']], 
    columns = ['class'],
    index = data_dict['patient_rows'][data_dict['patient_class_mask']])
cl_class_table = pd.DataFrame(
    copy(data_dict['cl_class'])[data_dict['cl_class_mask']], 
    columns = ['class'],
    index = data_dict['cl_exp_rows'][data_dict['cl_class_mask']])

class_table = pd.concat((patient_class_table, cl_class_table), axis = 0)
n_classes = class_table['class'].max()+1

#%% Classification accuracy

from sklearn.metrics import (
    average_precision_score,
    f1_score, 
    balanced_accuracy_score, 
    accuracy_score, 
    roc_auc_score
)
from sklearn.preprocessing import OneHotEncoder

from multiprocessing import Pool
from sae.evaluation import DSC
from scipy.special import softmax

oh_enc = OneHotEncoder(
    categories = [np.arange(n_classes, dtype = 'int')], 
    sparse_output = False
)

def mapper(p_f, class_true = class_table):
    p_df = next(iter(load_results(p_f).values()))
    i = p_df['task'].iloc[0]
    class_df_list = []
    for fold in np.unique(p_df['fold']):
        for run in np.unique(p_df['run']):
            ind = np.logical_and(p_df['fold'] == fold, p_df['run'] == run)
            class_preds = p_df.loc[ind]
            #class_preds = p_df.filter(regex = 'class_pred_[0-9]+').loc[ind]
            #class_preds['dataset'] = p_df.loc[ind, 'dataset'].to_numpy()
            class_preds.set_index(p_df.loc[ind, 'id'], inplace = True)
            class_preds = class_preds.loc[np.intersect1d(class_preds.index.to_numpy(), class_table.index.to_numpy())]
            class_preds_grouped = class_preds.groupby(
                by = 'dataset', 
                as_index = False, 
                sort = False
            )
            for data_name, dataset in class_preds_grouped:
                #dataset.drop('dataset', axis = 1, inplace = True)
                '''Ensure column order and drop other columns'''
                dataset = dataset[[f"class_pred_{i}" for i in np.arange(n_classes, dtype = 'int')]]
                preds = np.argmax(dataset.to_numpy(), axis = 1)[:,np.newaxis]
                #probs = softmax(dataset.to_numpy(), axis = 1)
                '''
                apr = average_precision_score(
                    y_true = class_true.loc[dataset.index].to_numpy(), 
                    y_score = dataset.to_numpy(), 
                    average = 'macro'
                )
                '''
                f1 = f1_score(
                    y_true = class_true.loc[dataset.index].to_numpy(), 
                    y_pred = preds, 
                    average = 'macro'
                )
                bacc = balanced_accuracy_score(
                    y_true = class_true.loc[dataset.index].to_numpy(), 
                    y_pred = preds
                )
                acc = accuracy_score(
                    y_true = class_true.loc[dataset.index].to_numpy(), 
                    y_pred = preds
                )
                '''
                auroc = roc_auc_score(
                    y_true = oh_enc.fit_transform(class_true.loc[dataset.index].to_numpy()), 
                    y_score = probs, 
                    average = 'macro', 
                    multi_class = 'ovr', 
                    labels = np.arange(n_classes, dtype = 'int')
                )
                '''
                class_dfi = pd.DataFrame({
                    'dataset' : data_name, 
                    'task' : i, 
                    'run' : run, 
                    'fold' : fold, 
                    #'apr' : apr, 
                    'f1' : f1, 
                    'bacc' : bacc, 
                    'acc' : acc, 
                    #'auroc' : np.mean(svm_cv['test_roc_auc'])
                    }, 
                    index = pd.RangeIndex(0,1))
                class_df_list.append(class_dfi)
    
    return pd.concat(class_df_list, axis = 0)

from multiprocessing import Pool
class_df_list = []
with Pool(processes = n_workers) as mp:
    for res in mp.map(
            mapper, 
            [prediction_dict[i] for i in prediction_dict.keys()]):
        class_df_list.append(res)
class_df = pd.concat(class_df_list, axis = 0)

class_df.to_csv(res_path + 'tissue_classifier_performance.csv.gz')

class_df_list = []
with Pool(processes = n_workers) as mp:
    for res in mp.map(
            mapper, 
            [ps_prediction_dict[i] for i in ps_prediction_dict.keys()]):
        class_df_list.append(res)
class_df = pd.concat(class_df_list, axis = 0)

class_df.to_csv(res_path + 'ps_tissue_classifier_performance.csv.gz')
