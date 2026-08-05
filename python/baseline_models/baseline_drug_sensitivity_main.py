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
import json

import pandas as pd
import numpy as np
import importlib

script_path = os.environ.get('MODAE_SCRIPT_PATH', default = None)
if script_path is None:
    raise ValueError('Please define MODAE_SCRIPT_PATH')

dsc_spec = importlib.util.spec_from_file_location(
    'dataset_collections', 
    os.path.join(script_path, 'python/utilities/cancer_dataset_collections.py')
)
dsc = importlib.util.module_from_spec(dsc_spec)
sys.modules['dataset_collections'] = dsc
dsc_spec.loader.exec_module(dsc)

drd_spec = importlib.util.spec_from_file_location(
    'baseline_dr_data', 
    os.path.join(script_path, 'python/baseline_models/baseline_drug_sensitivity_data.py')
)
drd = importlib.util.module_from_spec(drd_spec)
sys.modules['baseline_dr_data'] = drd
drd_spec.loader.exec_module(drd)

drm_spec = importlib.util.spec_from_file_location(
    'baseline_dr_models', 
    os.path.join(script_path, 'python/baseline_models/baseline_drug_sensitivity_models.py')
)
drm = importlib.util.module_from_spec(drm_spec)
sys.modules['baseline_dr_models'] = drm
drm_spec.loader.exec_module(drm)

#%% environment
output_path = os.environ.get('MODAE_OUTPUT_PATH', default = None)
if output_path is None:
    raise ValueError('Please define MODAE_OUTPUT_PATH')
data_root = os.environ.get('MODAE_DATA_PATH', default = None)
if data_root is None:
    raise ValueError('Please define MODAE_DATA_PATH')

#%% imports

from dataset_collections import get_tcga_brca_ctrp_ccle_full
from baseline_dr_data import load_sensitivity_data, load_old_sensitivity_data
from baseline_dr_models import nested_train_eval_drugwise

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.linear_model import ElasticNet
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor

from sklearn.feature_selection import SelectKBest, f_regression
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GroupKFold, ParameterGrid
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.base import clone

from threadpoolctl import threadpool_limits, threadpool_info

from modae.evaluation_utilities import cross_validation_index
from modae.data_utilities import data_serialization

#%% Main
def main(
        data_root_dir,
        result_file = '',
        data_standardize = True, 
        use_pca = False, 
        pca_npc = 100, 
        model_name = 'elasticnet'
):
    os.makedirs(output_path, exist_ok=True)
    # Load and process data
    from dataset_collections import get_xia_ctrp_data_combat
    data_dict = get_xia_ctrp_data_combat(
        data_root = data_root, 
        tissue_classifier = True, # Matching processing
        include_metas_labels = True, 
        exclude_metas_data = False
    )
    # Use data serialization function to replicate MODAE splits
    cv_data = data_serialization(
        data_dict,
        nruns = 1,
        nfolds = 5,
        stratify_survival = True,
        stratify_subtype = True,
        stratify_dr = True,
        stratify_time_quantiles = 10,
        cv_seed = 0, # first run + 0 starting seed
        model_args = {'drug_response_model': True},
        data_standardize = True,
        shared_scaler = False,
        omics_layer = 'mrna',
        ps_validation_sets = True,
        ps_test_sets = True,
        return_data = True
    )
    # Define model
    if use_pca:
        pca = PCA(n_components = pca_npc, whiten = True)
    if model_name == 'elasticnet':
        model = ElasticNet(max_iter = 100, random_state = 2)
        param_grid = {
            'alpha' : 10**np.linspace(-3, 0, 11), 
            'l1_ratio' : np.power(np.linspace(0.1**4, 0.99**4, 11), 1./4.)
        }
    elif model_name == 'fs_elasticnet':
        eln_model = ElasticNet(max_iter = 100, random_state = 2)
        model = Pipeline(steps = [
            ('fs', SelectKBest(score_func = f_regression)), 
            ('eln', eln_model)
        ])
        param_grid = {
            'fs__k' : np.linspace(100, 1000, 3).astype('int64').tolist(), 
            'eln__alpha' : 10**np.linspace(-3, -1, 3), 
            'eln__l1_ratio' : np.power(np.linspace(0.1**4, 0.99**4, 3), 1./4.)
        }
    elif model_name == 'gbt':
        param_grid = {
            'max_depth' : np.arange(3, 10).tolist(), 
            'ccp_alpha' : 10**np.linspace(-3, 0, 4),
            'learning_rate' : 10**np.linspace(-3, -1, 5)
        }
        model = GradientBoostingRegressor(
            n_estimators = 200, 
            loss = 'squared_error',
            random_state = 2
        )
    elif model_name == 'fs_gbt':
        gbr_model = GradientBoostingRegressor(
            n_estimators = 200, 
            loss = 'squared_error',
            random_state = 2
        )
        model = Pipeline(steps = [
            ('fs', SelectKBest(score_func = f_regression)), 
            ('gbt', gbr_model)
        ])
        param_grid = {
            'fs__k' : np.linspace(100, 1000, 3).astype('int64').tolist(), 
            'gbt__max_depth' : np.arange(2, 9, 2), 
            'gbt__ccp_alpha' : 10**np.linspace(-3, 0, 4), 
            'gbt__learning_rate' : np.power(np.linspace(0.1**4, 0.99**4, 3), 1./4.)
        }
    # Run CV with nested train test split
    result_list = []
    for par_dict in ParameterGrid(param_grid):
        # Outer CV
        for run in cv_data.keys():
            for fold in cv_data[run].keys():
                for (trainset_name, testset_name), setkey in [(('train', 'validation'), 'ps_cv_'), (('retrain', 'test'), 'test_cv_')]:
                    sliced_data = cv_data[run][fold][setkey]['cell_line']
                    X = sliced_data['train']['expression']
                    X_test = sliced_data['test']['expression']
                    if use_pca:
                        pca_inst = clone(pca)
                        X = pca_inst.fit_transform(X)
                        X_test = pca_inst.transform(X_test)
                    n_drugs = sliced_data['train']['drug_response'].shape[1]
                    for drugi in np.arange(n_drugs):
                        model_inst = clone(model)
                        model_inst = model_inst.set_params(**par_dict)
                        maski = sliced_data['train']['drug_response_mask'][:,drugi]
                        maski_test = sliced_data['test']['drug_response_mask'][:,drugi]
                        Xi = X[maski,:]
                        Xi_test = X_test[maski_test,:]
                        yi = sliced_data['train']['drug_response'][maski, drugi]
                        yi_test = sliced_data['test']['drug_response'][maski_test, drugi]
                        model_inst = model_inst.fit(Xi, yi)
                        yi_pred = model_inst.predict(Xi)
                        yi_test_pred = model_inst.predict(Xi_test)
                        mse = mean_squared_error(yi, yi_pred)
                        mse_test = mean_squared_error(yi_test, yi_test_pred)
                        r2 = mean_squared_error(yi, yi_pred)
                        r2_test = mean_squared_error(yi_test, yi_test_pred)
                        predi_df = pd.DataFrame({
                            'run': run,
                            'fold': fold,
                            'set': [trainset_name, testset_name],
                            'params': json.dumps(par_dict),
                            'drugi': drugi,
                            'mse': [mse, mse_test],
                            'r2': [r2, r2_test]
                        })
                        result_list.append(predi_df)
    model_str = f'pca{pca_npc}_' if use_pca else ''
    model_str += model_name
    result_file = os.path.join(output_path, f'baseline_results/drug_response/ctrp/{model_str}_metrics.csv.gz')
    result_df = pd.concat(result_list)
    result_df.to_csv(result_file)
    prediction_file = os.path.join(output_path, f'baseline_results/drug_response/ctrp/{model_str}_predictions.csv.gz')
    if prediction_file:
        # Select best model for final training
        test_df = result_df.loc[result_df['set'] == 'test']
        r2_q90 = test_df.groupby(['params', 'run', 'fold']).agg(
            func = {
                'r2': lambda x: np.quantile(x, 0.9)
            },
            axis = 0
        )
        r2_q90_mean = r2_q90.groupby(['params']).agg(func =  {'r2': 'mean'}, axis = 0)
        best_params_str = r2_q90_mean['r2'].idxmax(axis = 0)
        best_params = json.loads(best_params_str)
        # Final training
        prediction_list = []
        X = data_dict['cl_exp']
        scaler = StandardScaler().fit(X)
        X = scaler.transform(X)
        X_target = data_dict['patient_exp']
        X_target = scaler.transform(X_target)
        if use_pca:
            pca_inst = clone(pca)
            X = pca_inst.fit_transform(X)
            X_target = pca_inst.transform(X_target)
        n_drugs = data_dict['dr_table'].shape[1]
        preds = np.full((X.shape[0], n_drugs), np.nan)
        preds_target = np.full((X_target.shape[0], n_drugs), np.nan)
        for drugi in np.arange(n_drugs):
            model_inst = clone(model)
            model_inst = model_inst.set_params(**best_params)
            maski = data_dict['dr_table_mask'][:,drugi]
            Xi = X[maski,:]
            yi = data_dict['dr_table'][maski, drugi]
            model_inst = model_inst.fit(Xi, yi)
            yi_pred = model_inst.predict(X)
            yi_target_pred = model_inst.predict(X_target)
            preds[:,drugi] = yi_pred
            preds_target[:,drugi] = yi_target_pred
        pred_df = pd.DataFrame(
            preds,
            index = data_dict['cl_exp_rows'],
            columns = data_dict['dr_table_cols']
        )
        pred_target_df = pd.DataFrame(
            preds_target,
            index = data_dict['patient_rows'],
            columns = data_dict['dr_table_cols']
        )
        pred_df = pd.concat([pred_df, pred_target_df], axis = 0)
        pred_df.to_csv(prediction_file)
    #res.to_csv(os.path.join(data_root_dir, result_file))

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