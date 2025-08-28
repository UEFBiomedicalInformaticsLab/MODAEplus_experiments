#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 12:10:18 2024

@author: rintala
"""

import numpy as np
import pandas as pd
import os
import glob
from copy import copy
import importlib.util
import sys

dsc_spec = importlib.util.spec_from_file_location('dataset_collections', 'datasets/cancer_dataset_collections.py')
dsc = importlib.util.module_from_spec(dsc_spec)
sys.modules['dataset_collections'] = dsc
dsc_spec.loader.exec_module(dsc)

dsp_spec = importlib.util.spec_from_file_location('dataset_processing', 'datasets/cancer_dataset_processing.py')
dsp = importlib.util.module_from_spec(dsp_spec)
sys.modules['dataset_processing'] = dsp
dsp_spec.loader.exec_module(dsp)

et_spec = importlib.util.spec_from_file_location('evaluation_tools', 'utilities/evaluation_tools.py')
et = importlib.util.module_from_spec(et_spec)
sys.modules['evaluation_tools'] = et
et_spec.loader.exec_module(et)

import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

pd.set_option('display.max_columns', None)

#%% Dataset matching a given MODAE result
from sys import platform
from socket import gethostname
if gethostname() == 'teemu-pc':
    base_path = '/home/teemu/research_work/superAE_HPO/'
    home = True
elif platform == 'linux':
    base_path = '/research/work/rintala/superAE_HPO/'
    home = False
else:
    base_path = '//research.uefad.uef.fi/workdir/superAE_HPO/'
    home = False

#res_path = base_path + '20230821_random_search/brca_test_noclfilter/'
#res_path = base_path + '20231027_random_search/brca_test_noclfilter/'
#res_path = base_path + '20231101_random_search/brca_test_noclfilter/'
#res_path = base_path + '20231103_random_search/brca_test_noclfilter/'
#res_path = base_path + '20231120_random_search/brca_test_noclfilter/'
#res_path = base_path + '20231121_random_search/brca_test_noclfilter/'
#res_path = base_path + '20231124_random_search/pancan_test/'
#res_path = base_path + '20240103_random_search/brca_test_noclfilter/'
#res_path = base_path + '20240108_random_search/brca_test_noclfilter/'
#res_path = base_path + '20240115_random_search/brca_test_noclfilter/'
#res_path = base_path + '20240118_random_search/brca_test_noclfilter/'
#res_path = base_path + '20240119_random_search/brca_test_noclfilter/'
#res_path = base_path + '20240123_random_search/brca_test_noclfilter/'
#res_path = base_path + '20240124_random_search/brca_test_noclfilter/'
#res_path = base_path + '20240124_random_search/brca_test_noclfilter_alternative/'
#res_path = base_path + '20240125_random_search/brca_test_noclfilter/'
#res_path = base_path + '20240125_random_search/brca_test_noclfilter_nopre/'
#res_path = base_path + '20240206_random_search/brca_test_noclfilter/'
#res_path = base_path + '20240208_random_search/brca_test_noclfilter/'
#res_path = base_path + '20240211_random_search/brca_test_noclfilter/'
#res_path = base_path + '20240213_random_search/brca_test_noclfilter/'
#res_path = base_path + '20240213_random_search/brca_test_noclfilter_pretrain/'
#res_path = base_path + '20240215_random_search/pancan_test/'
#res_path = base_path + '20240216_random_search/brca_test_noclfilter_pretrain/'
#res_path = base_path + '20240221_random_search/scanb_test_noclfilter/'
res_path = base_path + '20240301_random_search/scanb_test_noclfilter/'
#res_path = base_path + '20240302_random_search/scanb_test_noclfilter/'
#res_path = base_path + '20240303_random_search/scanb_test_noclfilter/'

gene_preselection = False
if 'pancan' in res_path:
    from dataset_collections import get_tcga_pancan_ctrp_ccle_solid
    data_dict = get_tcga_pancan_ctrp_ccle_solid(home = home, gene_preselection = gene_preselection)
    internal_dataset = 'tcga'
    external_dataset = None
elif 'brca' in res_path:
    from dataset_collections import get_tcga_brca_ctrp_ccle_full
    data_dict = get_tcga_brca_ctrp_ccle_full(home = home, gene_preselection = gene_preselection)
    internal_dataset = 'tcga_brca'
    external_dataset = 'scanb'
elif 'scanb' in res_path:
    from dataset_collections import get_scanb_ctrp_ccle_full
    data_dict = get_scanb_ctrp_ccle_full(home = home, gene_preselection = gene_preselection)
    internal_dataset = 'scanb'
    external_dataset = 'tcga_brca'
elif 'pancan_ccle' in res_path:
    from dataset_collections import get_tcga_pancan_ccle_solid
    data_dict = get_tcga_pancan_ccle_solid(home = home, gene_preselection = gene_preselection)
    internal_dataset = 'tcga'
    external_dataset = None
elif 'brca_ccle' in res_path:
    from dataset_collections import get_tcga_brca_ccle_full
    data_dict = get_tcga_brca_ccle_full(home = home, gene_preselection = gene_preselection)
    internal_dataset = 'tcga_brca'
    external_dataset = 'scanb'
elif 'scanb_ccle' in res_path:
    from dataset_collections import get_scanb_ccle_full
    data_dict = get_scanb_ccle_full(home = home, gene_preselection = gene_preselection)
    internal_dataset = 'scanb'
    external_dataset = 'tcga_brca'
else:
    raise ValueError('res_path format does not match known dataset')

#%%

ind = np.argwhere(data_dict['dr_table_mask'])
ccle_ids = data_dict['cl_exp_rows']
drug_ids = data_dict['dr_table_cols']
cl_exp_df = pd.DataFrame(
    data_dict['cl_exp'], 
    index = ccle_ids, 
    columns = data_dict['gene_ids'])
key_df = pd.DataFrame({
    'CCLE_ID' : ccle_ids[ind[:,0]], 
    'CTRP_TREAT_ID' : drug_ids[ind[:,1]], 
    'AAC' : data_dict['dr_table'][ind[:,0], ind[:,1]]
})
data_ds = {'X_key' : key_df, 'cl_rnaseq' : cl_exp_df}

#%% Best params for 10 PC EN

from threadpoolctl import threadpool_limits
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.linear_model import ElasticNet

pc10_en_fn = f"{base_path}../drug_response_dataset/pca10_elasticnet_aac_old.csv"
result_df = pd.read_csv(pc10_en_fn, header = 0)

scaler_instance = StandardScaler()
pca_instance = PCA(n_components = 10, whiten = True)

x_z = scaler_instance.fit_transform(data_dict['cl_exp'])
#x_z_df = pd.DataFrame(x_z, index = data_dict['cl_exp_rows'])
x_pc = pca_instance.fit_transform(x_z)

x_pc_df = pd.DataFrame(x_pc, index = data_dict['cl_exp_rows'])

preds = []
drug_names = result_df['drug'].unique()
for drugi in drug_names:
    di_df = result_df.loc[result_df['drug'] == drugi]
    alpha = np.exp(np.mean(np.log(di_df['alpha'])))
    l1_ratio = np.mean(np.sqrt(di_df['l1_ratio']))**2
    model = ElasticNet(
        max_iter = 100, 
        random_state = 2, 
        alpha = alpha, 
        l1_ratio = l1_ratio
    )
    di_target_df = key_df.loc[key_df['CTRP_TREAT_ID'] == drugi]
    X = x_pc_df.loc[di_target_df['CCLE_ID']].to_numpy()
    y = di_target_df['AAC'].to_numpy()
    with threadpool_limits(limits = 10, user_api = 'blas'):
        model.fit(X, y)
    
    y_pred_all = model.predict(x_pc)
    preds.append(y_pred_all)

np.any([i.shape[0] < 1293 for i in preds])
preds_arr = np.vstack(preds).T

preds_all_df = pd.DataFrame(
    preds_arr, 
    index = data_dict['cl_exp_rows'], 
    columns = drug_names
)

dr_baseline_path = f"{base_path}../baseline_results/drug_response/ctrp/"
os.makedirs(dr_baseline_path, exist_ok = True)
dr_baseline_preds_fn = f"{dr_baseline_path}pc10_en_preds.csv.gz"
preds_all_df.to_csv(dr_baseline_preds_fn)
