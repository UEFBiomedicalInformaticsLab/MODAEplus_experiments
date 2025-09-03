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

dsc_spec = importlib.util.spec_from_file_location(
    'dataset_collections', 
    'datasets/cancer_dataset_collections.py'
)
dsc = importlib.util.module_from_spec(dsc_spec)
sys.modules['dataset_collections'] = dsc
dsc_spec.loader.exec_module(dsc)

et_spec = importlib.util.spec_from_file_location(
    'evaluation_tools', 
    'utilities/evaluation_tools.py'
)
et = importlib.util.module_from_spec(et_spec)
sys.modules['evaluation_tools'] = et
et_spec.loader.exec_module(et)

import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
warnings.simplefilter(action='ignore', category=scipy.stats.ConstantInputWarning)

n_workers = 100 # match with SLURM config

#%% Data
output_path = os.environ.get('MODAE_OUTPUT_PATH', default = None)
if output_path is None:
    raise ValueError('Please define MODAE_OUTPUT_PATH')
data_root = os.environ.get('MODAE_DATA_PATH', default = None)
if data_root is None:
    raise ValueError('Please define MODAE_DATA_PATH')

res_path = f"{output_path}20250410_random_search/pancan_test/"

run_date = re.sub('[a-z_]*/[a-z_]*/$', '', res_path)
run_date = re.sub('^.*/', '', run_date)
#dss_sensitivity = int(run_date) >= 20250313
exclude_metas_data = int(run_date) >= 20250527 and int(run_date) <= 20250527

#embedding_files = glob.glob(res_path + '*test_cv_*embeddings*.csv.gz')
prediction_files = glob.glob(res_path + '*test_cv_*predictions*.csv.gz')
ps_prediction_files = glob.glob(res_path + '*ps_cv_*predictions*.csv.gz')

if False:
    complete_preds = [i for i in prediction_files if 'fold' not in i]
    complete_ps_preds = [i for i in ps_prediction_files if 'fold' not in i]
    cpred_tasks = [re.sub('.*?task', '', i) for i in complete_preds]
    cpred_tasks = [re.sub('\\.csv\\.gz', '', i) for i in cpred_tasks]
    
    cpspred_tasks = [re.sub('.*?task', '', i) for i in complete_ps_preds]
    cpspred_tasks = [re.sub('\\.csv\\.gz', '', i) for i in cpspred_tasks]
    
    full_preds = set(cpspred_tasks).intersection(cpred_tasks)
    incomplete_preds = set([str(i) for i in range(1,101)]).difference(full_preds)

from evaluation_tools import load_results

#embeddings_dict = load_results(embedding_files)
prediction_dict = load_results(prediction_files)
ps_prediction_dict = load_results(ps_prediction_files)

#%% CTRP or CCLE?
prediction_test = prediction_dict[1].filter(regex = 'dr_pred_[0-9]+')
if prediction_test.shape[1] == 24:
    response_dataset = 'ccle'
elif prediction_test.shape[1] == 545:
    response_dataset = 'ctrp'
    dss_sensitivity = False
elif prediction_test.shape[1] == 544:
    response_dataset = 'ctrp' # Xia data
    dss_sensitivity = True
else:
    raise ValueError('Drug response prediction column number does not match known datasets.')
#%% Measured data
date = int(re.search('[0-9]+', res_path).group(0))
drug_response_scaling = (date == 20231219)

if re.search('brca', res_path) is not None:
    if response_dataset == 'ccle':
        from dataset_collections import get_tcga_brca_ccle_full
        data_dict = get_tcga_brca_ccle_full(
            data_root = data_root, 
            drug_response_maxscale = drug_response_scaling
        )
    elif response_dataset == 'ctrp':
        from dataset_collections import get_tcga_brca_ctrp_ccle_full
        data_dict = get_tcga_brca_ctrp_ccle_full(
            data_root = data_root, 
            drug_response_maxscale = drug_response_scaling
        )
elif re.search('pancan', res_path) is not None:
    if response_dataset == 'ccle':
        from dataset_collections import get_tcga_pancan_ccle_solid
        data_dict = get_tcga_pancan_ccle_solid(
            data_root = data_root, 
            drug_response_maxscale = drug_response_scaling
        )
    elif response_dataset == 'ctrp':
        if dss_sensitivity:
            from dataset_collections import get_xia_ctrp_data
            data_dict = get_xia_ctrp_data(
                data_root = data_root, 
                exclude_metas_data = exclude_metas_data
            )
        else:
            from dataset_collections import get_tcga_pancan_ctrp_ccle_solid
            data_dict = get_tcga_pancan_ctrp_ccle_solid(
                data_root = data_root, 
                drug_response_maxscale = drug_response_scaling
            )
elif re.search('scanb', res_path) is not None:
    if response_dataset == 'ccle':
        from dataset_collections import get_scanb_ccle_full
        data_dict = get_scanb_ccle_full(
            data_root = data_root, 
            drug_response_maxscale = drug_response_scaling
        )
    elif response_dataset == 'ctrp':
        from dataset_collections import get_scanb_ctrp_ccle_full
        data_dict = get_scanb_ctrp_ccle_full(
            data_root = data_root, 
            drug_response_maxscale = drug_response_scaling
        )

dr_temp = copy(data_dict['dr_table'])
dr_temp[np.logical_not(data_dict['dr_table_mask'])] = np.nan
cl_variable_table = pd.DataFrame(
    dr_temp, 
    columns = data_dict['dr_table_cols'], 
    index = data_dict['cl_exp_rows'])

#%% R^2 scores

from sklearn.metrics import r2_score, mean_squared_error

dr_true = copy(cl_variable_table)

def mapper(p_df, dr_true = dr_true):
    i = p_df['task'].iloc[0]
    r2_df_list = []
    mse_df_list = []
    for fold in np.unique(p_df['fold']):
        for run in np.unique(p_df['run']):
            p_df_train_subset = p_df.loc[
                (p_df['fold'] == fold) & 
                (p_df['run'] == run) & 
                (p_df['dataset'] == 'cl_train')]
            p_df_test_subset = p_df.loc[
                (p_df['fold'] == fold) & 
                (p_df['run'] == run) & 
                (p_df['dataset'] == 'cl_test')]
            dr_pred_train = p_df_train_subset.filter(regex = 'dr_pred_[0-9]+')
            dr_pred_test = p_df_test_subset.filter(regex = 'dr_pred_[0-9]+')
            
            cl_id_train = p_df_train_subset['id']
            cl_id_test = p_df_test_subset['id']
            
            dr_true_train = dr_true.loc[cl_id_train]
            dr_true_test = dr_true.loc[cl_id_test]
            
            for coli, col in enumerate(dr_true_train.columns):
                dri = dr_true_train[col]
                train_nna_ind = np.logical_not(dri.isna().to_numpy())
                train_nna_indi = np.logical_not(dr_pred_train.iloc[:, coli].isna().to_numpy())
                train_nna_ind = np.logical_and(train_nna_ind, train_nna_indi)
                dr_true_train_nna = dri.loc[train_nna_ind]
                dr_pred_train_nna = dr_pred_train.iloc[:, coli].loc[train_nna_ind]
                drii = dr_true_test[col]
                test_nna_ind = np.logical_not(drii.isna().to_numpy())
                test_nna_indi = np.logical_not(dr_pred_test.iloc[:, coli].isna().to_numpy())
                test_nna_ind = np.logical_and(test_nna_ind, test_nna_indi)
                dr_true_test_nna = drii.loc[test_nna_ind]
                dr_pred_test_nna = dr_pred_test.iloc[:, coli].loc[test_nna_ind]
                if train_nna_ind.any():
                    r2_train = r2_score(
                        dr_true_train_nna.to_numpy(), 
                        dr_pred_train_nna.to_numpy())
                    mse_train = mean_squared_error(
                        dr_true_train_nna.to_numpy(), 
                        dr_pred_train_nna.to_numpy())
                else:
                    r2_train = np.nan
                    mse_train = np.nan
                if test_nna_ind.any():
                    r2_test = r2_score(
                        dr_true_test_nna.to_numpy(), 
                        dr_pred_test_nna.to_numpy())
                    mse_test = mean_squared_error(
                        dr_true_test_nna.to_numpy(), 
                        dr_pred_test_nna.to_numpy())
                else:
                    r2_test = np.nan
                    mse_test = np.nan
                
                r2_col_df = pd.DataFrame({
                    'task' : i, 
                    'run' : run, 
                    'fold' : fold, 
                    'drug' : col, 
                    'dataset' : ['cl_train', 'cl_test'], 
                    'dr_r2' : [r2_train, r2_test]},
                    index = pd.RangeIndex(0,2))
                r2_df_list.append(r2_col_df)
                
                mse_col_df = pd.DataFrame({
                    'task' : i, 
                    'run' : run, 
                    'fold' : fold, 
                    'drug' : col, 
                    'dataset' : ['cl_train', 'cl_test'], 
                    'dr_mse' : [mse_train, mse_test]},
                    index = pd.RangeIndex(0,2))
                mse_df_list.append(mse_col_df)
    return (pd.concat(r2_df_list, axis = 0), 
            pd.concat(mse_df_list, axis = 0))

from multiprocessing import Pool
r2_df_list = []
mse_df_list = []
with Pool(processes = n_workers) as mp:
    for res_r2, res_mse in mp.map(
            mapper, 
            [prediction_dict[i] for i in prediction_dict.keys()]):
        r2_df_list.append(res_r2)
        mse_df_list.append(res_mse)
r2_df = pd.concat(r2_df_list, axis = 0)
mse_df = pd.concat(mse_df_list, axis = 0)

#%% save r2

r2_df.to_csv(res_path + 'drug_response_r2.csv.gz')
mse_df.to_csv(res_path + 'drug_response_mse.csv.gz')

#%% parameter search R^2 scores

from multiprocessing import Pool
r2_df_list = []
mse_df_list = []
with Pool(processes = n_workers) as mp:
    for res_r2, res_mse in mp.map(
            mapper, 
            [ps_prediction_dict[i] for i in ps_prediction_dict.keys()]):
        r2_df_list.append(res_r2)
        mse_df_list.append(res_mse)
r2_df = pd.concat(r2_df_list, axis = 0)
mse_df = pd.concat(mse_df_list, axis = 0)

#%% save parameter search r2

r2_df.to_csv(res_path + 'ps_drug_response_r2.csv.gz')
mse_df.to_csv(res_path + 'ps_drug_response_mse.csv.gz')

#%% Mean MSE

dr_true_mean = np.nanmean(dr_true, axis = 0)
dr_mean_mse = np.mean((dr_true - np.expand_dims(dr_true_mean, axis = 0))**2)

#%% overall R^2 scores

if False:
    import tensorflow as tf
    
    dr_true_tf = tf.constant(data_dict['dr_table'], dtype = 'float32')
    dr_mask_tf = tf.constant(data_dict['dr_table_mask'], dtype = 'bool')
    dr_true_masked = tf.ragged.boolean_mask(dr_true_tf, dr_mask_tf)
    cl_index = dict(zip(data_dict['cl_exp_rows'], np.arange(dr_true_tf.shape[0])))
    
    def mapper(p_df, 
               dr_mask_tf = dr_mask_tf, 
               dr_true_masked = dr_true_masked, 
               cl_index = cl_index):
        i = p_df['task'].iloc[0]
        r2_df_list = []
        for fold in np.unique(p_df['fold']):
            for run in np.unique(p_df['run']):
                p_df_train_subset = p_df.loc[
                    (p_df['fold'] == fold) & 
                    (p_df['run'] == run) & 
                    (p_df['dataset'] == 'cl_train')]
                p_df_test_subset = p_df.loc[
                    (p_df['fold'] == fold) & 
                    (p_df['run'] == run) & 
                    (p_df['dataset'] == 'cl_test')]
                dr_pred_train = p_df_train_subset.filter(regex = 'dr_pred_[0-9]+')
                dr_pred_test = p_df_test_subset.filter(regex = 'dr_pred_[0-9]+')
                
                cl_id_train = p_df_train_subset['id']
                cl_id_test = p_df_test_subset['id']
                
                dr_pred_train_cl_ind = np.array([cl_index[i] for i in cl_id_train])
                dr_pred_test_cl_ind = np.array([cl_index[i] for i in cl_id_test])
                
                min_pred_train = dr_pred_train.min(axis = None)
                max_pred_train = dr_pred_train.max(axis = None)
                min_pred_test = dr_pred_test.min(axis = None)
                max_pred_test = dr_pred_test.max(axis = None)
                
                dr_pred_train_tf = tf.constant(dr_pred_train.to_numpy(), dtype = 'float32')
                dr_pred_test_tf = tf.constant(dr_pred_test.to_numpy(), dtype = 'float32')
                dr_pred_train_tf = tf.clip_by_value(dr_pred_train_tf, 0., 1.)
                dr_pred_test_tf = tf.clip_by_value(dr_pred_test_tf, 0., 1.)
                
                dr_pred_train_mask = tf.gather(dr_mask_tf, dr_pred_train_cl_ind, axis = 0)
                dr_pred_test_mask = tf.gather(dr_mask_tf, dr_pred_test_cl_ind, axis = 0)
                
                dr_pred_train_tf_masked = tf.ragged.boolean_mask(
                    dr_pred_train_tf, dr_pred_train_mask)
                dr_pred_test_tf_masked = tf.ragged.boolean_mask(
                    dr_pred_test_tf, dr_pred_test_mask)
                
                dr_true_train_tf = tf.gather(dr_true_masked, dr_pred_train_cl_ind, axis = 0)
                dr_true_test_tf = tf.gather(dr_true_masked, dr_pred_test_cl_ind, axis = 0)
                
                r2_train = tf.keras.metrics.R2Score(class_aggregation='uniform_average')
                r2_test = tf.keras.metrics.R2Score(class_aggregation='uniform_average')
                
                r2_train.update_state(
                    tf.expand_dims(dr_true_train_tf.merge_dims(0,-1), -1), 
                    tf.expand_dims(dr_pred_train_tf_masked.merge_dims(0,-1), -1))
                r2_test.update_state(
                    tf.expand_dims(dr_true_test_tf.merge_dims(0,-1), -1), 
                    tf.expand_dims(dr_pred_test_tf_masked.merge_dims(0,-1), -1))
                
                if False:
                    import seaborn as sns
                    sns.histplot(dr_true_train_tf.merge_dims(0,-1))
                    sns.histplot(dr_true_test_tf.merge_dims(0,-1))
                
                r2_df = pd.DataFrame({
                    'task' : i, 
                    'run' : run, 
                    'fold' : fold, 
                    'dataset' : ['cl_train', 'cl_test'], 
                    'dr_r2' : [r2_train.result().numpy(), 
                               r2_test.result().numpy()],
                    'min_pred' : [min_pred_train, min_pred_test],
                    'max_pred' : [max_pred_train, max_pred_test]},
                    index = pd.RangeIndex(0,2))
                r2_df_list.append(r2_df)
        return pd.concat(r2_df_list, axis = 0)
    
    from multiprocessing import Pool
    r2_df_list = []
    if False:
        with Pool(processes = n_workers) as mp:
            for res in mp.map(
                    mapper, 
                    [prediction_dict[i] for i in prediction_dict.keys()]):
                r2_df_list.append(res)
    else:
        for i in prediction_dict.keys():
            r2_df_list.append(mapper(prediction_dict[i]))
    r2_df = pd.concat(r2_df_list, axis = 0)
    r2_df.to_csv(res_path + 'adjusted_overall_drug_response_r2.csv.gz')

#%% drug-wise R^2 scores
if False:
    import tensorflow as tf
    
    dr_true_tf = tf.constant(data_dict['dr_table'], dtype = 'float32')
    dr_mask_tf = tf.constant(data_dict['dr_table_mask'], dtype = 'bool')
    dr_true_masked = tf.ragged.boolean_mask(dr_true_tf, dr_mask_tf)
    cl_index = dict(zip(data_dict['cl_exp_rows'], np.arange(dr_true_tf.shape[0])))
    
    def mapper(p_df, 
               dr_mask_tf = dr_mask_tf, 
               dr_true_tf = dr_true_tf, 
               cl_index = cl_index):
        i = p_df['task'].iloc[0]
        r2_df_list = []
        for fold in np.unique(p_df['fold']):
            for run in np.unique(p_df['run']):
                p_df_train_subset = p_df.loc[
                    (p_df['fold'] == fold) & 
                    (p_df['run'] == run) & 
                    (p_df['dataset'] == 'cl_train')]
                p_df_test_subset = p_df.loc[
                    (p_df['fold'] == fold) & 
                    (p_df['run'] == run) & 
                    (p_df['dataset'] == 'cl_test')]
                dr_pred_train = p_df_train_subset.filter(regex = 'dr_pred_[0-9]+')
                dr_pred_test = p_df_test_subset.filter(regex = 'dr_pred_[0-9]+')
                
                cl_id_train = p_df_train_subset['id']
                cl_id_test = p_df_test_subset['id']
                
                dr_pred_train_cl_ind = np.array([cl_index[i] for i in cl_id_train])
                dr_pred_test_cl_ind = np.array([cl_index[i] for i in cl_id_test])
                
                min_pred_train = dr_pred_train.min(axis = None)
                max_pred_train = dr_pred_train.max(axis = None)
                min_pred_test = dr_pred_test.min(axis = None)
                max_pred_test = dr_pred_test.max(axis = None)
                
                dr_pred_train_tf = tf.constant(dr_pred_train.to_numpy(), dtype = 'float32')
                dr_pred_test_tf = tf.constant(dr_pred_test.to_numpy(), dtype = 'float32')
                dr_pred_train_tf = tf.clip_by_value(dr_pred_train_tf, 0., 1.)
                dr_pred_test_tf = tf.clip_by_value(dr_pred_test_tf, 0., 1.)
                
                dr_pred_train_mask = tf.gather(dr_mask_tf, dr_pred_train_cl_ind, axis = 0)
                dr_pred_test_mask = tf.gather(dr_mask_tf, dr_pred_test_cl_ind, axis = 0)
                
                dr_pred_train_tf_masked = tf.ragged.boolean_mask(
                    tf.transpose(dr_pred_train_tf), tf.transpose(dr_pred_train_mask))
                dr_pred_test_tf_masked = tf.ragged.boolean_mask(
                    tf.transpose(dr_pred_test_tf), tf.transpose(dr_pred_test_mask))
                
                dr_true_train_tf = tf.gather(dr_true_tf, dr_pred_train_cl_ind, axis = 0)
                dr_true_test_tf = tf.gather(dr_true_tf, dr_pred_test_cl_ind, axis = 0)
                
                dr_true_train_tf_masked = tf.ragged.boolean_mask(
                    tf.transpose(dr_true_train_tf), tf.transpose(dr_pred_train_mask))
                dr_true_test_tf_masked = tf.ragged.boolean_mask(
                    tf.transpose(dr_true_test_tf), tf.transpose(dr_pred_test_mask))
                
                for coli, col in enumerate(data_dict['dr_table_cols']):
                    r2_train = tf.keras.metrics.R2Score(class_aggregation=None)
                    r2_test = tf.keras.metrics.R2Score(class_aggregation=None)
                    
                    r2_train.update_state(
                        tf.expand_dims(dr_true_train_tf_masked[coli], -1), 
                        tf.expand_dims(dr_pred_train_tf_masked[coli], -1))
                    r2_test.update_state(
                        tf.expand_dims(dr_true_test_tf_masked[coli], -1), 
                        tf.expand_dims(dr_pred_test_tf_masked[coli], -1))
                    
                    r2_df = pd.DataFrame({
                        'task' : i, 
                        'run' : run, 
                        'fold' : fold, 
                        'drug' : col, 
                        'dataset' : ['cl_train', 'cl_test'], 
                        'dr_r2' : [r2_train.result().numpy(), 
                                   r2_test.result().numpy()],
                        'min_pred' : [min_pred_train, min_pred_test],
                        'max_pred' : [max_pred_train, max_pred_test]},
                        index = pd.RangeIndex(0,2))
                    r2_df_list.append(r2_df)
        return pd.concat(r2_df_list, axis = 0)
    
    from multiprocessing import Pool
    r2_df_list = []
    if False:
        with Pool(processes = n_workers) as mp:
            for res in mp.map(
                    mapper, 
                    [prediction_dict[i] for i in prediction_dict.keys()]):
                r2_df_list.append(res)
    else:
        for i in prediction_dict.keys():
            r2_df_list.append(mapper(prediction_dict[i]))
    r2_df = pd.concat(r2_df_list, axis = 0)
    
    #%% save r2
    
    r2_df.to_csv(res_path + 'adjusted_drug_response_r2.csv.gz')
    
#%% prediction correlations

#from sklearn.metrics import r2_score, mean_squared_error
from scipy.stats import spearmanr
#spearmanr(np.random.rand(1e2), np.random.rand(1e2))

def mapper(p_df):
    i = p_df['task'].iloc[0]
    cor_df_list = []
    for fold in np.unique(p_df['fold']):
        for run in np.unique(p_df['run']):
            p_df_train_subset = p_df.loc[
                (p_df['fold'] == fold) & 
                (p_df['run'] == run) & 
                (p_df['dataset'] == 'cl_train')]
            p_df_test_subset = p_df.loc[
                (p_df['fold'] == fold) & 
                (p_df['run'] == run) & 
                (p_df['dataset'] == 'cl_test')]
            dr_pred_train = p_df_train_subset.filter(regex = 'dr_pred_[0-9]+')
            dr_pred_test = p_df_test_subset.filter(regex = 'dr_pred_[0-9]+')
            
            train_cor = []
            test_cor = []
            di = []
            dj = []
            for coli, colnamei in enumerate(dr_pred_train.columns[:-1]):
                for colnamej in dr_pred_train.columns[(coli+1):]:
                    di.append(colnamei)
                    dj.append(colnamej)
                    train_cor_ij = spearmanr(
                        dr_pred_train[colnamei], 
                        dr_pred_train[colnamej])
                    test_cor_ij = spearmanr(
                        dr_pred_test[colnamei], 
                        dr_pred_test[colnamej])
                    train_cor.append(train_cor_ij.statistic)
                    test_cor.append(test_cor_ij.statistic)
                    if False:
                        cor_ij_df = pd.DataFrame({
                            'task' : i, 
                            'run' : run, 
                            'fold' : fold, 
                            'drug_i' : colnamei, 
                            'drug_j' : colnamej, 
                            'dataset' : ['cl_train', 'cl_test'], 
                            'cor' : [
                                train_cor_ij.statistic, 
                                test_cor_ij.statistic
                            ]},
                            index = pd.RangeIndex(0,2))
                        cor_df_list.append(cor_ij_df)
            cor_ij_df = pd.DataFrame({
                'task' : i, 
                'run' : run, 
                'fold' : fold, 
                'drug_i' : di, 
                'drug_j' : dj, 
                'dr_cor_train' : train_cor,
                'dr_cor_test' : test_cor})
            cor_df_list.append(cor_ij_df)
    return pd.concat(cor_df_list, axis = 0)

from multiprocessing import Pool
cor_df_list = []
with Pool(processes = n_workers, maxtasksperchild = 1) as mp:
    for res in mp.map(
        mapper, 
        [prediction_dict[i] for i in prediction_dict.keys()]
    ):
        cor_df_list.append(res)
cor_df = pd.concat(cor_df_list, axis = 0)

#%% save r2

cor_df.to_csv(res_path + 'drug_response_cor.csv.gz')
