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

dsc_spec = importlib.util.spec_from_file_location(
    'dataset_collections', 
    f"{script_path}/python/utilities/cancer_dataset_collections.py"
)
dsc = importlib.util.module_from_spec(dsc_spec)
sys.modules['dataset_collections'] = dsc
dsc_spec.loader.exec_module(dsc)

et_spec = importlib.util.spec_from_file_location(
    'evaluation_tools', 
    f"{script_path}/python/utilities/evaluation_tools.py"
)
et = importlib.util.module_from_spec(et_spec)
sys.modules['evaluation_tools'] = et
et_spec.loader.exec_module(et)

import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

#%% Data
output_path = os.environ.get('MODAE_OUTPUT_PATH', default = None)
if output_path is None:
    raise ValueError('Please define MODAE_OUTPUT_PATH')
data_root = os.environ.get('MODAE_DATA_PATH', default = None)
if data_root is None:
    raise ValueError('Please define MODAE_DATA_PATH')

res_path = f"{output_path}20250410_random_search/pancan_test/"

embedding_files = glob.glob(res_path + '*test_cv_*embeddings*.csv.gz')
prediction_files = glob.glob(res_path + '*test_cv_*predictions*.csv.gz')

from evaluation_tools import load_results

embeddings_dict = load_results(embedding_files)
prediction_dict = load_results(prediction_files)

#%% 

def neighborhood_preds_stnr(preds, neighbor_mat_a, neighbor_mat_b):
    a_mean = np.mean(preds[neighbor_mat_a,:], axis = 1)
    b_mean = np.mean(preds[neighbor_mat_b,:], axis = 1)
    a_std = np.std(preds[neighbor_mat_a,:], axis = 1)
    b_std = np.std(preds[neighbor_mat_a,:], axis = 1)
    
    stnr = np.abs(a_mean - b_mean) / (a_std + b_std)
    
    return stnr

#%% Measured data
import re
date = int(re.search('[0-9]+', res_path).group(0))
drug_response_scaling =  date == 20231219

if re.search('brca', res_path) is not None:
    from dataset_collections import get_tcga_brca_ctrp_ccle_full
    data_dict = get_tcga_brca_ctrp_ccle_full(
        data_root = data_root, 
        drug_response_maxscale = drug_response_scaling
    )
elif re.search('pancan', res_path) is not None:
    from dataset_collections import get_tcga_pancan_ctrp_ccle_solid
    data_dict = get_tcga_pancan_ctrp_ccle_solid(
        data_root = data_root, 
        drug_response_maxscale = drug_response_scaling
    )

#%% Preparation

patient_variable_table = pd.DataFrame(
    {
     'event' : data_dict['survival_event'],
     'time' : data_dict['survival_time']
    }, 
    index = data_dict['patient_rows'])

patient_variable_table['event_time'] = patient_variable_table['time']
patient_variable_table.loc[patient_variable_table['event'] == 0, 'event_time'] = np.nan

def bin_variable(df, col, n_quants = 10):
    quants = np.nanquantile(
        df[col], 
        np.arange(1, n_quants) / n_quants)
    df[col + '_binned'] = np.digitize(
        df[col], 
        bins = quants, 
        right = True)
    df.loc[df[col].isna(), col + '_binned'] = np.nan
    df.drop(columns = col, inplace = True)
    return df

patient_variable_table_binned = copy(patient_variable_table)
bin_variable(patient_variable_table_binned, 'time', 2)
bin_variable(patient_variable_table_binned, 'event_time', 2)

dr_temp = copy(data_dict['dr_table'])
dr_temp[np.logical_not(data_dict['dr_table_mask'])] = np.nan
cl_variable_table = pd.DataFrame(
    dr_temp, 
    columns = data_dict['dr_table_cols'], 
    index = data_dict['cl_exp_rows'])

cl_variable_table_binned = copy(cl_variable_table)
for col in data_dict['dr_table_cols']:
    bin_variable(cl_variable_table_binned, col, 2)

def neighborhood_preds_stnr(preds, neighbor_mat_a, neighbor_mat_b):
    a_mean = np.mean(preds[neighbor_mat_a,:], axis = 1)
    b_mean = np.mean(preds[neighbor_mat_b,:], axis = 1)
    a_std = np.std(preds[neighbor_mat_a,:], axis = 1)
    b_std = np.std(preds[neighbor_mat_a,:], axis = 1)
    
    stnr = np.abs(a_mean - b_mean) / (a_std + b_std)
    
    return stnr

def neighborhood_mean_diff(x, neighbor_mat_a, neighbor_mat_b):
    a_mean = np.mean(x[neighbor_mat_a,:], axis = 1)
    b_mean = np.mean(x[neighbor_mat_b,:], axis = 1)
    
    mdiff = np.abs(a_mean - b_mean)
    
    return mdiff

#%% Silhouette

from multiprocessing import Pool

def mapper(
        e_df, 
        cl_variable_table = cl_variable_table, 
        cl_variable_table_binned = cl_variable_table_binned, 
        patient_variable_table = patient_variable_table, 
        patient_variable_table_binned = patient_variable_table_binned
        ):
    i = e_df['task'].iloc[0]
    cl_silhouette_list = []
    patient_silhouette_list = []
    for fold in np.unique(e_df['fold']):
        for run in np.unique(e_df['run']):
            for dataset in [
                    ('cl_test', 'cl_train', cl_variable_table, cl_variable_table_binned), 
                    ('patient_test', 'patient_train', patient_variable_table, patient_variable_table_binned)]:
                e_df_subset = e_df.loc[
                    (e_df['fold'] == fold) & 
                    (e_df['run'] == run) & 
                    (e_df['dataset'].isin(dataset))]
                X = e_df_subset.filter(regex = 'z[0-9]+')
                X_dist = pairwise_distances(X, metric = 'euclidean')
                np.fill_diagonal(X_dist, 0)
                sample_type = dataset[0].replace('_test', '')
                for col, col_binned in zip(dataset[2].columns, dataset[3].columns):
                    nnan_ind = np.logical_not(dataset[2][col].isna().to_numpy())
                    X_dist_nnan = X_dist[nnan_ind,:][:,nnan_ind]
                    
                    if np.unique(dataset[3][col_binned].loc[nnan_ind].astype('int64')).shape[0] > 1:
                        silh = silhouette_score(
                            X_dist_nnan, 
                            dataset[3][col_binned].loc[nnan_ind].astype('int64').to_numpy(),
                            metric = 'precomputed')
                        
                        silh_df = pd.DataFrame({
                            'task' : i, 
                            'run' : run, 
                            'fold' : fold, 
                            'sample_type' : sample_type, 
                            'col' : col, 
                            'silhouette' : silh},
                            index = pd.RangeIndex(0,1))
                        if sample_type == 'cl':
                            cl_silhouette_list.append(silh_df)
                        if sample_type == 'patient':
                            patient_silhouette_list.append(silh_df)
    return {'cl' : pd.concat(cl_silhouette_list, axis = 0), 
            'patient' : pd.concat(patient_silhouette_list, axis = 0)}

cl_silh_list = []
patient_silh_list = []

if True:
    with Pool(processes = 20) as mp:
        for res in mp.map(
                mapper, 
                [embeddings_dict[i] for i in embeddings_dict.keys()]):
            cl_silh_list.append(res['cl'])
            patient_silh_list.append(res['patient'])
else:
    for i in embeddings_dict.keys():
        res = mapper(embeddings_dict[i])
        cl_silh_list.append(res['cl'])
        patient_silh_list.append(res['patient'])

cl_silh_df = pd.concat(cl_silh_list, axis = 0)
patient_silh_df = pd.concat(patient_silh_list, axis = 0)

#%% Single-threaded silhouette calculation

if False:
    cl_silhouette_list = []
    patient_silhouette_list = []
    
    for i in embeddings_dict.keys():
        e_df = embeddings_dict[i]
        for fold in np.unique(embeddings_dict[i]['fold']):
            for run in np.unique(embeddings_dict[i]['run']):
                for dataset in [('cl_test', 'cl_train', cl_variable_table, cl_variable_table_binned), 
                                ('patient_test', 'patient_train', patient_variable_table, patient_variable_table_binned)]:
                    e_df_subset = e_df.loc[
                        (e_df['fold'] == fold) & 
                        (e_df['run'] == run) & 
                        (e_df['dataset'].isin(dataset))]
                    X = e_df_subset.filter(regex = 'z[0-9]+')
                    X_dist = pairwise_distances(X, metric = 'euclidean')
                    np.fill_diagonal(X_dist, 0)
                    for col, col_binned in zip(dataset[2].columns, dataset[3].columns):
                        sample_type = dataset[0].replace('_test', '')
                        nnan_ind = np.logical_not(dataset[2][col].isna().to_numpy())
                        X_dist_nnan = X_dist[nnan_ind,:][:,nnan_ind]
                        
                        if np.unique(dataset[3][col_binned].loc[nnan_ind].astype('int64')).shape[0] > 1:
                            silh = silhouette_score(
                                X_dist_nnan, 
                                dataset[3][col_binned].loc[nnan_ind].astype('int64').to_numpy(),
                                metric = 'precomputed')
                            
                            silh_df = pd.DataFrame({
                                'task' : i, 
                                'run' : run, 
                                'fold' : fold, 
                                'sample_type' : sample_type, 
                                'col' : col, 
                                'silhouette' : silh},
                                index = pd.RangeIndex(0,1))
                            if sample_type == 'cl':
                                cl_silhouette_list.append(silh_df)
                            if sample_type == 'patient':
                                patient_silhouette_list.append(silh_df)
    
    cl_silh_df = pd.concat(cl_silhouette_list)
    patient_silh_df = pd.concat(patient_silhouette_list)

#%% Save Silhouette

cl_silh_df.to_csv(res_path + 'cl_silhouette.csv.gz')
patient_silh_df.to_csv(res_path + 'patient_silhouette.csv.gz')


#%% STNR

dr_stnr_list = []
sr_stnr_list = []

for i in embeddings_dict.keys():
    e_df = embeddings_dict[i]
    for fold in np.unique(embeddings_dict[i]['fold']):
        for run in np.unique(embeddings_dict[i]['run']):
            for dataset in [('cl_test', 'cl_train', cl_variable_table, cl_variable_table_binned), 
                            ('patient_test', 'patient_train', patient_variable_table, patient_variable_table_binned)]:
                e_df_subset = e_df.loc[
                    (e_df['fold'] == fold) & 
                    (e_df['run'] == run) & 
                    (e_df['dataset'].isin(dataset))]
                X = e_df_subset.filter(regex = 'z[0-9]+')
                X_dist = pairwise_distances(X, metric = 'euclidean')
                np.fill_diagonal(X_dist, 0)
                for col, col_binned in zip(dataset[2].columns, dataset[3].columns):
                    sample_type = dataset[0].replace('_test', '')
                    nnan_ind = np.logical_not(dataset[2][col].isna().to_numpy())
                    X_dist_nnan = X_dist[nnan_ind,:][:,nnan_ind]

                    if np.unique(dataset[3][col_binned].loc[nnan_ind].astype('int64')).shape[0] > 1:
                        X_dist_rank = np.argsort(X_dist_nnan, axis = 1)

                        col_stnr = neighborhood_preds_stnr(
                            np.expand_dims(dataset[2][col].loc[nnan_ind].to_numpy(), axis = 1), 
                            X_dist_rank[:,:11], 
                            X_dist_rank[:,:101])
                        
                        stnr_df = pd.DataFrame({
                            'task' : i, 
                            'run' : run, 
                            'fold' : fold, 
                            'sample_type' : sample_type, 
                            'col' : col, 
                            'stnr' : np.mean(col_stnr)},
                            index = pd.RangeIndex(0,1))
                        if sample_type == 'cl':
                            dr_stnr_list.append(stnr_df)
                        if sample_type == 'patient':
                            sr_stnr_list.append(stnr_df)

dr_stnr_df = pd.concat(dr_stnr_list)
sr_stnr_df = pd.concat(sr_stnr_list)

#%% Save stnr

dr_stnr_df.to_csv(res_path + 'drug_response_neighborhood_signal_to_noise.csv.gz')
sr_stnr_df.to_csv(res_path + 'survival_event_neighborhood_signal_to_noise.csv.gz')

#%% neighborhood mean difference

from multiprocessing import Pool

def mapper(
        e_df, 
        cl_variable_table = cl_variable_table, 
        cl_variable_table_binned = cl_variable_table_binned, 
        patient_variable_table = patient_variable_table, 
        patient_variable_table_binned = patient_variable_table_binned
        ):
    i = e_df['task'].iloc[0]
    dr_mdiff_list = []
    sr_mdiff_list = []
    for fold in np.unique(e_df['fold']):
        for run in np.unique(e_df['run']):
            for dataset in [
                    ('cl_test', 'cl_train', cl_variable_table, cl_variable_table_binned), 
                    ('patient_test', 'patient_train', patient_variable_table, patient_variable_table_binned)]:
                e_df_subset = e_df.loc[
                    (e_df['fold'] == fold) & 
                    (e_df['run'] == run) & 
                    (e_df['dataset'].isin(dataset))]
                X = e_df_subset.filter(regex = 'z[0-9]+')
                # It is better to calculate the full distance matrix and subset 
                # for non-missing values. 
                X_dist = pairwise_distances(X, metric = 'euclidean')
                np.fill_diagonal(X_dist, 0)
                y_vals = dataset[2]
                sample_type = dataset[0].replace('_test', '')
                for col in y_vals.columns:
                    nnan_ind = np.logical_not(y_vals[col].isna().to_numpy())
                    X_dist_nnan = X_dist[nnan_ind,:][:,nnan_ind]
                    X_dist_rank = np.argsort(X_dist_nnan, axis = 1)
                    col_mdiff = neighborhood_mean_diff(
                        np.expand_dims(y_vals[col].loc[nnan_ind].to_numpy(), axis = 1), 
                        X_dist_rank[:,:11], 
                        X_dist_rank[:,:101])
                    mdiff_df = pd.DataFrame({
                        'task' : i, 
                        'run' : run, 
                        'fold' : fold, 
                        'sample_type' : sample_type, 
                        'col' : col, 
                        'mdiff' : np.mean(col_mdiff)},
                        index = pd.RangeIndex(0,1))
                    if sample_type == 'cl':
                        dr_mdiff_list.append(mdiff_df)
                    if sample_type == 'patient':
                        sr_mdiff_list.append(mdiff_df)
        
    return {'dr' : dr_mdiff_list, 'sr' : sr_mdiff_list}

dr_mdiff_list = []
sr_mdiff_list = []

if True:
    with Pool(processes = 20) as mp:
        for res in mp.map(
                mapper, 
                [embeddings_dict[i] for i in embeddings_dict.keys()]):
            dr_mdiff_list.append(res['dr'])
            sr_mdiff_list.append(res['sr'])
else:
    for i in embeddings_dict.keys():
        res = mapper(embeddings_dict[i])
        dr_mdiff_list.append(res['dr'])
        sr_mdiff_list.append(res['sr'])

dr_mdiff = []
sr_mdiff = []
for i in np.arange(len(dr_mdiff_list)):
    dr_mdiff.append(pd.concat(dr_mdiff_list[i], axis = 0))
for i in np.arange(len(sr_mdiff_list)):
    sr_mdiff.append(pd.concat(sr_mdiff_list[i], axis = 0))

dr_mdiff_df = pd.concat(dr_mdiff, axis = 0)
sr_mdiff_df = pd.concat(sr_mdiff, axis = 0)

#%% Save real mdiff

dr_mdiff_df.to_csv(res_path + 'drug_response_neighborhood_mean_diff.csv.gz')
sr_mdiff_df.to_csv(res_path + 'survival_event_neighborhood_mean_diff.csv.gz')

#%% Prediction signal to noise ratios between local and global neighborhood
dr_stnr_list = []
sr_stnr_list = []

for i in embeddings_dict.keys():
    e_df = embeddings_dict[i]
    p_df = prediction_dict[i]
    for fold in np.unique(embeddings_dict[i]['fold']):
        for run in np.unique(embeddings_dict[i]['run']):
            for dataset in [('cl_test', 'cl_train'), ('patient_test', 'patient_train')]:
                e_df_subset = e_df.loc[
                    (e_df['fold'] == fold) & 
                    (e_df['run'] == run) & 
                    (e_df['dataset'].isin(dataset))]
                X = e_df_subset.filter(regex = 'z[0-9]+')
                p_df_subset = p_df.loc[
                    (p_df['fold'] == fold) & 
                    (p_df['run'] == run) & 
                    (p_df['dataset'].isin(dataset))]
                
                
                # Ensure correct row mapping
                prediction_map = dict(zip(p_df_subset['id'], np.arange(p_df_subset.shape[0])))
                prediction_ind = np.array([prediction_map[ind] for ind in e_df_subset['id']])
                p_df_subset = p_df_subset.iloc[prediction_ind,:]
                
                # https://jakevdp.github.io/blog/2013/04/29/benchmarking-nearest-neighbor-searches-in-python/
                kdtree = neighbors.KDTree(X.to_numpy(), leaf_size = 15) 
                
                neighbor_inds_10 = kdtree.query(X.to_numpy(), k = 10+1, return_distance = False)
                neighbor_inds_100 = kdtree.query(X.to_numpy(), k = 100+1, return_distance = False)
                
                y_dr = p_df_subset.filter(regex = 'dr_pred_[0-9]+')
                y_dr_mean_10 = np.mean(y_dr.to_numpy()[neighbor_inds_10,:], axis = 1)
                
                dr_stnr = neighborhood_preds_stnr(y_dr.to_numpy(), neighbor_inds_10, neighbor_inds_100)
                
                random_uniform_stnr = neighborhood_preds_stnr(np.random.uniform(size = (X.shape[0], 1000)), neighbor_inds_10, neighbor_inds_100)
                random_std_stnr = neighborhood_preds_stnr(np.random.standard_normal(size = (X.shape[0], 1000)), neighbor_inds_10, neighbor_inds_100)
                
                dr_stnr_df = pd.DataFrame({
                    'task' : i, 
                    'run' : run, 
                    'fold' : fold, 
                    'sample_type' : dataset[0].replace('_test', ''), 
                    **dict(zip(y_dr.columns.to_numpy(), np.mean(dr_stnr, axis = 0)))}, 
                    index = pd.RangeIndex(0,1))
                dr_stnr_list.append(dr_stnr_df)
                
                # survival
                y_sr = p_df_subset.filter(regex = 'survival_risk_[0-9]+')
                y_sr_nna_ind = np.logical_not(np.any(y_sr.isna().to_numpy(), axis = 1))
                if np.any(y_sr_nna_ind):
                    y_sr = y_sr.loc[y_sr_nna_ind,:]
                    
                    X = X.loc[y_sr_nna_ind,:]
                    kdtree = neighbors.KDTree(X.to_numpy(), leaf_size = 15) 
                    
                    neighbor_inds_10 = kdtree.query(X.to_numpy(), k = 10+1, return_distance = False)
                    neighbor_inds_100 = kdtree.query(X.to_numpy(), k = 100+1, return_distance = False)
                    
                    y_sr_mean_10 = np.mean(y_sr.to_numpy()[neighbor_inds_10,:], axis = 1)
                    sr_stnr = neighborhood_preds_stnr(y_sr.to_numpy(), neighbor_inds_10, neighbor_inds_100)
                    
                    sr_stnr_df = pd.DataFrame({
                        'task' : i, 
                        'run' : run, 
                        'fold' : fold, 
                        'sample_type' : dataset[0].replace('_test', ''), 
                        **dict(zip(y_sr.columns.to_numpy(), np.mean(sr_stnr, axis = 0)))}, 
                        index = pd.RangeIndex(0,1))
                    sr_stnr_list.append(sr_stnr_df)
                
                #sns.histplot(np.mean(dr_stnr, axis = 0))
                #sns.histplot(np.mean(sr_stnr, axis = 0))
                #sns.histplot(np.mean(random_uniform_stnr, axis = 0))
                #sns.histplot(np.mean(random_std_stnr, axis = 0))
                
dr_stnr_df = pd.concat(dr_stnr_list)
sr_stnr_df = pd.concat(sr_stnr_list)

#%% Save STNR

dr_stnr_df.to_csv(res_path + 'predicted_drug_response_neighborhood_signal_to_noise.csv.gz')
sr_stnr_df.to_csv(res_path + 'predicted_survival_risk_neighborhood_signal_to_noise.csv.gz')

#%% Prediction based mean difference
from multiprocessing import Pool

def mapper(
        e_df, 
        p_df
        ):
    i = e_df['task'].iloc[0]
    dr_mdiff_list = []
    sr_mdiff_list = []
    for fold in np.unique(e_df['fold']):
        for run in np.unique(e_df['run']):
            for dataset in [
                    ('cl_test', 'cl_train'), 
                    ('patient_test', 'patient_train')]:
                e_df_subset = e_df.loc[
                    (e_df['fold'] == fold) & 
                    (e_df['run'] == run) & 
                    (e_df['dataset'].isin(dataset[:2]))]
                X = e_df_subset.filter(regex = 'z[0-9]+')
                p_df_subset = p_df.loc[
                    (p_df['fold'] == fold) & 
                    (p_df['run'] == run) & 
                    (p_df['dataset'].isin(dataset[:2]))]
                
                # Ensure correct row mapping
                prediction_map = dict(zip(
                    p_df_subset['id'], 
                    np.arange(p_df_subset.shape[0])))
                prediction_ind = np.array([prediction_map[ind] for ind in e_df_subset['id']])
                p_df_subset = p_df_subset.iloc[prediction_ind,:]
                
                # https://jakevdp.github.io/blog/2013/04/29/benchmarking-nearest-neighbor-searches-in-python/
                kdtree = neighbors.KDTree(X.to_numpy(), leaf_size = 15) 
                
                neighbor_inds_10 = kdtree.query(X.to_numpy(), k = 10+1, return_distance = False)
                neighbor_inds_100 = kdtree.query(X.to_numpy(), k = 100+1, return_distance = False)
                
                y_dr = p_df_subset.filter(regex = 'dr_pred_[0-9]+')
                
                dr_mdiff = neighborhood_mean_diff(y_dr.to_numpy(), neighbor_inds_10, neighbor_inds_100)
                
                dr_mdiff_df = pd.DataFrame({
                    'task' : i, 
                    'run' : run, 
                    'fold' : fold, 
                    'sample_type' : dataset[0].replace('_test', ''), 
                    **dict(zip(y_dr.columns.to_numpy(), np.mean(dr_mdiff, axis = 0)))}, 
                    index = pd.RangeIndex(0,1))
                dr_mdiff_list.append(dr_mdiff_df)
                
                # survival
                y_sr = p_df_subset.filter(regex = 'survival_risk_[0-9]+')
                y_sr_nna_ind = np.logical_not(np.any(y_sr.isna().to_numpy(), axis = 1))
                if np.any(y_sr_nna_ind):
                    y_sr = y_sr.loc[y_sr_nna_ind,:]
                    
                    X = X.loc[y_sr_nna_ind,:]
                    kdtree = neighbors.KDTree(X.to_numpy(), leaf_size = 15) 
                    
                    neighbor_inds_10 = kdtree.query(X.to_numpy(), k = 10+1, return_distance = False)
                    neighbor_inds_100 = kdtree.query(X.to_numpy(), k = 100+1, return_distance = False)
                    
                    #y_sr_mean_10 = np.mean(y_sr.to_numpy()[neighbor_inds_10,:], axis = 1)
                    sr_mdiff = neighborhood_mean_diff(y_sr.to_numpy(), neighbor_inds_10, neighbor_inds_100)
                    
                    sr_mdiff_df = pd.DataFrame({
                        'task' : i, 
                        'run' : run, 
                        'fold' : fold, 
                        'sample_type' : dataset[0].replace('_test', ''), 
                        **dict(zip(y_sr.columns.to_numpy(), np.mean(sr_mdiff, axis = 0)))}, 
                        index = pd.RangeIndex(0,1))
                    sr_mdiff_list.append(sr_mdiff_df)
    return {'dr' : pd.concat(dr_mdiff_list, axis = 0), 
            'sr' : pd.concat(sr_mdiff_list, axis = 0)}

dr_mdiff_list = []
sr_mdiff_list = []

if True:
    with Pool(processes = 20) as mp:
        for res in mp.starmap(
                mapper, 
                [(embeddings_dict[i], prediction_dict[i]) for i in embeddings_dict.keys()]):
            dr_mdiff_list.append(res['dr'])
            sr_mdiff_list.append(res['sr'])

dr_mdiff_df = pd.concat(dr_mdiff_list)
sr_mdiff_df = pd.concat(sr_mdiff_list)

#%% save predicted mdiff

dr_mdiff_df.to_csv(res_path + 'predicted_drug_response_neighborhood_mean_diff.csv.gz')
sr_mdiff_df.to_csv(res_path + 'predicted_survival_event_neighborhood_mean_diff.csv.gz')

#%% Prediction-based Silhouette

from multiprocessing import Pool

def mapper(
        e_df, 
        p_df
        ):
    i = e_df['task'].iloc[0]
    dr_silh_list = []
    sr_silh_list = []
    for fold in np.unique(e_df['fold']):
        for run in np.unique(e_df['run']):
            for dataset in [
                    ('cl_test', 'cl_train'), 
                    ('patient_test', 'patient_train')]:
                e_df_subset = e_df.loc[
                    (e_df['fold'] == fold) & 
                    (e_df['run'] == run) & 
                    (e_df['dataset'].isin(dataset[:2]))]
                X = e_df_subset.filter(regex = 'z[0-9]+')
                p_df_subset = p_df.loc[
                    (p_df['fold'] == fold) & 
                    (p_df['run'] == run) & 
                    (p_df['dataset'].isin(dataset[:2]))]
                
                # Ensure correct row mapping
                prediction_map = dict(zip(
                    p_df_subset['id'], 
                    np.arange(p_df_subset.shape[0])))
                prediction_ind = np.array([prediction_map[ind] for ind in e_df_subset['id']])
                p_df_subset = p_df_subset.iloc[prediction_ind,:]
                
                
                
                y_dr = p_df_subset.filter(regex = 'dr_pred_[0-9]+')
                y_dr_binned = copy(y_dr)
                for col in y_dr_binned.columns:
                    bin_variable(y_dr_binned, col, 2)
                
                X_dist = pairwise_distances(X, metric = 'euclidean')
                np.fill_diagonal(X_dist, 0)
                for col in y_dr_binned.columns:
                    if np.unique(y_dr_binned[col].astype('int64')).shape[0] > 1:
                        silh = silhouette_score(
                            X_dist, 
                            y_dr_binned[col].astype('int64').to_numpy(),
                            metric = 'precomputed')
                        silh_df = pd.DataFrame({
                            'task' : i, 
                            'run' : run, 
                            'fold' : fold, 
                            'sample_type' : dataset[0].replace('_test', ''), 
                            'col' : col, 
                            'silhouette' : silh},
                            index = pd.RangeIndex(0,1))
                        dr_silh_list.append(silh_df)
                
                # survival
                y_sr = p_df_subset.filter(regex = 'survival_risk_[0-9]+')
                for col in y_sr.columns:
                    y_sr_nna_ind = np.logical_not(y_sr[col].isna().to_numpy())
                    if np.any(y_sr_nna_ind):
                        y_sr_i = y_sr.loc[y_sr_nna_ind, [col]]
                        y_sr_i_binned = copy(y_sr_i)
                        y_sr_i_binned = bin_variable(y_sr_i_binned, col, 2)
                        if np.unique(y_sr_i_binned.astype('int64')).shape[0] > 1:
                            X_dist_nnan = X_dist[y_sr_nna_ind,:][:,y_sr_nna_ind]
                            silh = silhouette_score(
                                X_dist_nnan, 
                                y_sr_i_binned[col + '_binned'].astype('int64').to_numpy(),
                                metric = 'precomputed')
                            silh_df = pd.DataFrame({
                                'task' : i, 
                                'run' : run, 
                                'fold' : fold, 
                                'sample_type' : dataset[0].replace('_test', ''), 
                                'col' : col + '_binned', 
                                'silhouette' : silh},
                                index = pd.RangeIndex(0,1))
                            sr_silh_list.append(silh_df)
    
    out = {}
    if len(dr_silh_list):
        out['dr'] = pd.concat(dr_silh_list, axis = 0)
    if len(sr_silh_list):
        out['sr'] = pd.concat(sr_silh_list, axis = 0)
    return out

dr_silh_list = []
sr_silh_list = []

with Pool(processes = 20) as mp:
    for res in mp.starmap(
            mapper, 
            [(embeddings_dict[i], prediction_dict[i]) for i in embeddings_dict.keys()]):
        dr_silh_list.append(res.get('dr', []))
        sr_silh_list.append(res.get('sr', []))

dr_silh_list = [i for i in dr_silh_list if len(i) > 0]
sr_silh_list = [i for i in sr_silh_list if len(i) > 0]

dr_silh_df = pd.concat(dr_silh_list)
sr_silh_df = pd.concat(sr_silh_list)

#%% Save predicted silhouette

dr_silh_df.to_csv(res_path + 'predicted_drug_response_silhouette.csv.gz')
sr_silh_df.to_csv(res_path + 'predicted_survival_event_silhouette.csv.gz')
