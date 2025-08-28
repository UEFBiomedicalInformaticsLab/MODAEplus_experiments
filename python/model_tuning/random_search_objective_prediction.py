#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 11:02:27 2023

@author: teemu
"""

import numpy as np
import pandas as pd
import time
import os
import warnings
import re

import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV
import sklearn.metrics as metrics
from sys import platform
from copy import copy

#%%
from sys import platform
from socket import gethostname
if gethostname() == 'teemu-pc':
    base_path = '/home/teemu/research_work/superAE_HPO/'
    home = True
elif platform == 'linux':
    base_path = '/research/work/rintala/superAE_HPO/'
    home = False
else:
    base_path = '//research/workdir/superAE_HPO/'
    home = False

#res_dir = base_path + '20231027_random_search/brca_test_noclfilter/'
#res_dir = base_path + '20231101_random_search/brca_test_noclfilter/'
#res_dir = base_path + '20231103_random_search/brca_test_noclfilter/'
#res_dir = base_path + 'random_search_231113/brca_test_noclfilter/'
#res_dir = base_path + '231117_random_search/brca_test_noclfilter/'
#res_dir = base_path + '20231120_random_search/brca_test_noclfilter/'
#res_dir = base_path + '20231121_random_search/brca_test_noclfilter/'
#res_dir = base_path + '20231122_random_search/scanb_test_noclfilter/'
#res_dir = base_path + '20231128_random_search/brca_test_noclfilter/'
#res_dir = base_path + '20231201_random_search/pancan_test/'
#res_path = base_path + '20240124_random_search/brca_test_noclfilter_alternative/'
#res_path = base_path + '20240125_random_search/brca_test_noclfilter/'
#res_path = base_path + '20240125_random_search/brca_test_noclfilter_nopre/'
#res_path = base_path + '20240206_random_search/brca_test_noclfilter/'
#res_path = base_path + '20240208_random_search/brca_test_noclfilter/'
#res_path = base_path + '20240211_random_search/brca_test_noclfilter/'
#res_path = base_path + '20250110_random_search/pancan_test/'
#res_path = base_path + '20250117_random_search/pancan_test/'
#res_path = base_path + '20250121_random_search/pancan_test/'
#res_path = base_path + '20250124_random_search/pancan_test/'
#res_path = base_path + '20250202_random_search/pancan_test/'
#res_path = base_path + '20250205_random_search/scanb_test_noclfilter/'
#res_path = base_path + '20250207_random_search/pancan_test/'
#res_path = base_path + '20250216_random_search/pancan_test/'
#res_path = base_path + '20250217_random_search/pancan_test/'
#res_path = base_path + '20250219_random_search/pancan_test/'
#res_path = base_path + '20250223_random_search/pancan_test/'
#res_path = base_path + '20250506_random_search/pancan_test/'
#res_path = base_path + '20250520_random_search/pancan_test/'
res_path = base_path + '20250603_random_search/pancan_test/'

#%%

X = pd.read_csv(res_path + 'objective_x.csv', index_col = 0, header = 0)
y = pd.read_csv(res_path + 'objective_y.csv', index_col = 0, header = 0)
y_ps = pd.read_csv(res_path + 'ps_objective_y.csv', index_col = 0, header = 0)

X.dropna(axis = 1, inplace = True)
#%%

from sklearn.model_selection import RepeatedKFold, GridSearchCV, cross_val_score
from sklearn.linear_model import ElasticNet
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.pipeline import Pipeline
from sklearn.metrics import r2_score
from copy import copy
from threadpoolctl import threadpool_limits, threadpool_info

import seaborn as sns

if not np.all(np.isfinite(X['log10_alpha'])):
    X.drop(columns = 'log10_alpha', inplace = True)

piped_eln_param_grid = {
    'elnet__alpha' : 10**np.linspace(-3, 0, 11), 
    'elnet__l1_ratio' : np.power(np.linspace(0.1**4, 0.99**4, 11), 1./4.)
    }
rkf_inner = RepeatedKFold(n_splits=5, n_repeats=1, random_state=1)
scale_pipe = Pipeline(steps = [
    ('scale', StandardScaler(with_mean = True, with_std = True)), 
    ('polynomials', PolynomialFeatures(degree = 1, include_bias = False)), 
    ('elnet', ElasticNet(max_iter = 1000)),
    ])
gs = GridSearchCV(
    scale_pipe, 
    param_grid = piped_eln_param_grid, 
    cv = rkf_inner, 
    scoring = 'r2',
    n_jobs = 1)

eln_res = {}
eln_res_coef = {}
ps_eln_res = {}
ps_eln_res_coef = {}
for i, j in zip(y_ps.columns, y.columns):
    ps_gsi = copy(gs)
    ps_y_filter = np.isfinite(y_ps[i]).to_numpy()
    ps_y_filter = np.logical_and(ps_y_filter, np.logical_not(y_ps[i].isna().to_numpy()))
    ps_gsi = ps_gsi.fit(X.loc[ps_y_filter,:].to_numpy(), np.expand_dims(y_ps.loc[ps_y_filter,i].to_numpy(), -1))
    gsi = copy(gs)
    gsi = gsi.fit(X.to_numpy(), np.expand_dims(y[j].to_numpy(), -1))
    
    eln_res[j] = gsi.best_score_
    ps_eln_res[i] = ps_gsi.best_score_
    eln_res_coef[j] = pd.DataFrame({
        'hyper_parameter' : gsi.best_estimator_.named_steps['polynomials'].get_feature_names_out(X.columns.to_numpy()),
        'target' : j,  
        'coefficient' : gsi.best_estimator_.named_steps['elnet'].coef_})
    ps_eln_res_coef[i] = pd.DataFrame({
        'hyper_parameter' : ps_gsi.best_estimator_.named_steps['polynomials'].get_feature_names_out(X.columns.to_numpy()),
        'target' : i,  
        'coefficient' : ps_gsi.best_estimator_.named_steps['elnet'].coef_})

eln_res_coef = pd.concat(eln_res_coef, axis = 0)
ps_eln_res_coef = pd.concat(ps_eln_res_coef, axis = 0)

#%% test
fig, ax = plt.subplots(figsize=(10, 12))
sns.set_style("whitegrid")
#sns.barplot(eln_res_coef_filtered, y = 'hyper_parameter', x = 'coefficient', hue = 'target')
sns.barplot(eln_res_coef, x = 'coefficient', y = 'hyper_parameter', hue = 'target')
#ax.set_xticklabels(ax.get_xticklabels(), rotation=90, horizontalalignment='center')
plt.grid()
#plt.show()
plt.tight_layout()
plt.savefig(f"{res_path}../plots/model_tuning_test_coefficients.png", 
            dpi = 300)

#%% ps
fig, ax = plt.subplots(figsize=(10, 12))
sns.set_style("whitegrid")
#sns.barplot(eln_res_coef_filtered, y = 'hyper_parameter', x = 'coefficient', hue = 'target')
sns.barplot(ps_eln_res_coef, x = 'coefficient', y = 'hyper_parameter', hue = 'target')
#ax.set_xticklabels(ax.get_xticklabels(), rotation=90, horizontalalignment='center')
plt.grid()
#plt.show()
plt.tight_layout()
plt.savefig(f"{res_path}../plots/model_tuning_ps_coefficients.png", 
            dpi = 300)

#%%
#fig, ax = plt.subplots(figsize=(10, 12))
g = sns.FacetGrid(eln_res_coef, row = 'hyper_parameter', height = 0.5, aspect = 2)
g.map_dataframe(sns.barplot, x = 'coefficient', hue = 'target')
g.add_legend()
#ax.set_xticklabels(ax.get_xticklabels(), rotation=90, horizontalalignment='center')
#plt.show()

#%%
eln_res_coef_filtered =  eln_res_coef.loc[eln_res_coef['coefficient'].abs() > 0.01]
