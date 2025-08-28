# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

import importlib.util
import sys

nthreads = 40
pca_model = False

dsc_spec = importlib.util.spec_from_file_location('dataset_collections', 'datasets/cancer_dataset_collections.py')
dsc = importlib.util.module_from_spec(dsc_spec)
sys.modules['dataset_collections'] = dsc
dsc_spec.loader.exec_module(dsc)

bsm_spec = importlib.util.spec_from_file_location('baseline_survival_models', 'baseline_models_surv/baseline_survival_models.py')
bsm = importlib.util.module_from_spec(bsm_spec)
sys.modules['baseline_survival_models'] = bsm
bsm_spec.loader.exec_module(bsm)

#%% 
from sys import platform
from socket import gethostname
if gethostname() == 'teemu-pc':
    base_path = '/home/teemu/research_work/'
    home = True
elif platform == 'linux':
    base_path = '/research/work/rintala/'
    home = False
else:
    base_path = '//research.uefad.uef.fi/workdir/'
    home = False

res_path = f"{base_path}baseline_results/survival/"

#%% Load internal data
#source_dataset = 'scanb'
#target_dataset = 'tcga_brca'
#source_dataset = 'tcga_pancancer'
source_dataset = 'tcga'
target_dataset = None
gene_preselection = False

if source_dataset == 'tcga_brca':
    from dataset_collections import get_tcga_brca_ctrp_ccle_full
    internal_data_dict = get_tcga_brca_ctrp_ccle_full(home = home, gene_preselection = gene_preselection)
    res_path += 'brca/tcga/'
elif source_dataset == 'scanb':
    from dataset_collections import get_scanb_ctrp_ccle_full
    internal_data_dict = get_scanb_ctrp_ccle_full(home = home, gene_preselection = gene_preselection)
    res_path += 'brca/scanb/'
elif source_dataset == 'tcga':
    from dataset_collections import get_tcga_pancan_ctrp_ccle_solid
    internal_data_dict = get_tcga_pancan_ctrp_ccle_solid(
        home = home, 
        gene_preselection = gene_preselection, 
        tissue_classifier = False
    )
    res_path += 'pancan/'
if gene_preselection:
    res_path += 'gene_preselection/'

import os
os.makedirs(res_path, exist_ok = True)

#%% Load external data

if target_dataset == 'scanb':
    from dataset_collections import get_scanb
    external_data_dict = get_scanb(home = home)
elif target_dataset == 'tcga_brca':
    from dataset_collections import get_tcga_brca
    external_data_dict = get_tcga_brca(home = home)
else:
    external_data_dict = None

#%% Baseline model PCA + Cox Elastic net
from sksurv.linear_model import CoxPHSurvivalAnalysis, CoxnetSurvivalAnalysis
from sksurv.metrics import concordance_index_censored, integrated_brier_score
from sklearn.model_selection import GridSearchCV, KFold, RepeatedKFold, RepeatedStratifiedKFold, ShuffleSplit
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.feature_selection import SelectKBest
from sklearn.pipeline import Pipeline
from baseline_survival_models import nested_survival_eval
from copy import copy

eln_model = CoxnetSurvivalAnalysis(max_iter = 100, fit_baseline_model = False)
model = Pipeline(steps = [
    #('scale', StandardScaler(with_mean = True, with_std = True)),
    #('pca', PCA(whiten = False, svd_solver = 'arpack')),
    ('eln', copy(eln_model))])
#model = copy(eln_model)
param_grid = {
    #'pca__n_components' : [10], 
    'eln__alphas' : np.expand_dims(10**np.linspace(-3, -1, 11), (-1)).tolist(), 
    'eln__l1_ratio' : np.power(np.linspace(0.1**4, 0.90**4, 5), 1./4.)}

X = internal_data_dict['patient_exp']
y_pairlist = [(e,t) for e,t in zip(internal_data_dict['survival_event'], internal_data_dict['survival_time'])]
y = np.array(y_pairlist, dtype = [('event', '?'),('time', '<i4')])

cv_outer = RepeatedStratifiedKFold(n_splits = 5, n_repeats = 2, random_state = 0)
#cv_inner = RepeatedStratifiedKFold(n_splits = 5, n_repeats = 1, random_state = 1)
cv_inner = ShuffleSplit(n_splits = 1, test_size = 0.25, random_state = 1)

#pca_instance = PCA(n_components = pca_npc, whiten = True)

if np.logical_not(internal_data_dict['survival_mask']).sum() > 10:
    cvg_pairlist = [(e,t) for e,t in zip(internal_data_dict['survival_event'], internal_data_dict['survival_mask'])]
    cv_outer_args = {'groups' : cvg_pairlist}
else:
    cv_outer_args = {}

#cvg_pairlist = [(e,t) for e,t in zip(internal_data_dict['survival_event'], internal_data_dict['survival_mask'])]
#cv_groups = np.char.add(y[])

N_patient = internal_data_dict.get('patient_exp', np.array([])).shape[0]
N_cl = internal_data_dict.get('cl_exp', np.array([])).shape[0]

strat_var = np.full((N_patient,),'')
sevent = internal_data_dict['survival_event']
stimes = internal_data_dict['survival_time']
smask = internal_data_dict['survival_mask']
stq = np.quantile(stimes[smask], np.arange(1, 10) / 10)
stbins = np.full(smask.shape, -1, dtype = 'int64')
stbins[smask] = np.digitize(stimes[smask], stq)
#stbins[stimes < 0] = -1 # Set all ambiguous/NaNs to -1 (own bin)
strat_var = np.char.add(strat_var, sevent.astype('str'))
strat_var = np.char.add(strat_var, stbins.astype('str'))

X_covar = copy(internal_data_dict['survival_covariates'])
X_covar_vcounts = [np.unique(X_covar[:,i], return_counts = True)[1] for i in np.arange(X_covar.shape[1])]

# Remove last onehot

#%% Inner test
if pca_model:
    import warnings
    warnings.filterwarnings("ignore", category=UserWarning)
    pca10_res = nested_survival_eval(
        model = model, 
        param_grid = param_grid, 
        X = X, 
        y = y, 
        y_mask = internal_data_dict['survival_mask'],
        X_covar = X_covar, 
        strat_var = strat_var, 
        scaler_instance = StandardScaler(),
        covar_scaler_instance = StandardScaler(),
        pca_instance = PCA(n_components = 10, whiten = True, svd_solver = 'arpack'), 
        #pca_instance = None, 
        cv_outer = copy(cv_outer), 
        cv_inner = copy(cv_inner), 
        thread_limit = nthreads)
    
    pca10_train_c = pca10_res['train_c'].mean()
    pca10_test_c = pca10_res['test_c'].mean()
    
    print(f"train C: {pca10_train_c}, test C: {pca10_test_c}")
    
    #param_grid['pca__n_components'] = [100]
    
    pca100_res = nested_survival_eval(
        model = model, 
        param_grid = param_grid, 
        X = X, 
        y = y, 
        y_mask = internal_data_dict['survival_mask'],
        X_covar = X_covar, 
        strat_var = strat_var, 
        scaler_instance = StandardScaler(),
        covar_scaler_instance = StandardScaler(),
        pca_instance = PCA(n_components = 100, whiten = True, svd_solver = 'arpack'), 
        #pca_instance = None, 
        cv_outer = copy(cv_outer), 
        cv_inner = copy(cv_inner), 
        thread_limit = nthreads)
    
    pca100_train_c = pca100_res['train_c'].mean()
    pca100_test_c = pca100_res['test_c'].mean()
    
    print(f"train C: {pca100_train_c}, test C: {pca100_test_c}")
    
    res_df_10pc = pd.DataFrame({
        'npc' : 10, 
        'fold' : pca10_res['fold'], 
        f"{source_dataset}_train_c" : pca10_res['train_c'], 
        f"{source_dataset}_test_c" : pca10_res['test_c']})
    res_df_100pc = pd.DataFrame({
        'npc' : 100, 
        'fold' : pca100_res['fold'], 
        f"{source_dataset}_train_c" : pca100_res['train_c'], 
        f"{source_dataset}_test_c" : pca100_res['test_c']})
    
    res_df = pd.concat([res_df_10pc, res_df_100pc], axis = 0)
    res_df.to_csv(res_path + 'dr_model_internal_cross_validation.csv')

#%% External data
if external_data_dict is not None:
    from dataset_collections import match_genes
    
    X_ext = external_data_dict['patient_exp']
    X_ext = match_genes(X_ext, external_data_dict['gene_ids'], internal_data_dict['gene_ids'])
    
    y_ext_pairlist = [(e,t) for e,t in zip(external_data_dict['survival_event'], external_data_dict['survival_time'])]
    y_ext = np.array(y_ext_pairlist, dtype = [('event', '?'),('time', '<i4')])
    y_ext_mask = external_data_dict['survival_mask']
    
    X_ext_masked = X_ext[y_ext_mask, :]
    y_ext_masked = y_ext[y_ext_mask]
    
    y_masked = y[internal_data_dict['survival_mask']]
    
    scaler = StandardScaler()
    ext_scaler = StandardScaler()
    covar_scaler = StandardScaler()
    
    # Separate scaling for omics
    X_scaled = scaler.fit_transform(X)
    X_ext_scaled = ext_scaler.fit_transform(X_ext)
    X_ext_scaled_masked = X_ext_scaled[y_ext_mask, :]
    
    # Shared scaling for covariates
    X_covar_scaled = covar_scaler.fit_transform(X_covar[internal_data_dict['survival_mask']])
    X_ext_covar_scaled = covar_scaler.transform(external_data_dict['survival_covariates'][y_ext_mask])


#%% PCA
if pca_model and external_data_dict is not None:
    from baseline_survival_models import grid_search
    ext_res = []
    for npc in [10,100]:
        pca = pca_instance = PCA(n_components = npc, whiten = True, svd_solver = 'arpack')
        X_pcs = pca.fit_transform(X_scaled)
        X_pcs_masked = X_pcs[internal_data_dict['survival_mask'], :]
        X_concat_masked = np.concatenate((X_pcs_masked, X_covar_scaled), axis = 1)
        
        X_ext_pcs_masked = pca.transform(X_ext_scaled_masked)
        X_ext_concat_masked = np.concatenate((X_ext_pcs_masked, X_ext_covar_scaled), axis = 1)
        
        # Final grid
        grid_res = grid_search(
            model = model, 
            X = X_concat_masked, 
            y = y_masked, 
            param_grid = param_grid, 
            cv = cv_outer, 
            cv_args = {'groups' : [i for i, masked in zip(cvg_pairlist, internal_data_dict['survival_mask']) if masked]}, 
            thread_limit = nthreads)
        
        # Best params and final model
        best_i = np.argmax(grid_res['test_score'])
        best_par_dict = dict([(k, grid_res[k][best_i]) for k in grid_res.keys()])
        best_par_train_score = best_par_dict.pop('train_score')
        best_par_test_score = best_par_dict.pop('test_score')
        
        final_model = copy(model)
        final_model = final_model.set_params(**best_par_dict)
        final_model = final_model.fit(X_concat_masked, y_masked)
        
        # eval
        internal_score = final_model.score(X_concat_masked, y_masked)
        external_score = final_model.score(X_ext_concat_masked, y_ext_masked)
        print(f"{source_dataset} C: {internal_score}, {target_dataset} C: {external_score}")
        ext_res.append(pd.DataFrame({
            'npc' : npc, 
            f"{source_dataset}_c" : internal_score, 
            f"{target_dataset}_c" : external_score}, 
            index = pd.RangeIndex(0,1)))
    
    ext_res_df = pd.concat(ext_res, axis = 0)
    ext_res_df.to_csv(res_path + 'dr_model_external_validation.csv')

#%% Baseline FS + Cox Elastic Net
from multiprocessing import Pool

def univariate_coxph(X, y, j, m):
    Xj = X[:, j : j + 1]
    if np.var(Xj) > 0. and (Xj > 0.).mean() > 0.1:
        m.fit(Xj, y)
        score = m.score(Xj, y)
    else:
        score = 0.
    return (j, score)

def fit_and_score_features(X, y):
    n_features = X.shape[1]
    
    m = CoxPHSurvivalAnalysis()
    scores = []
    if True:
        with Pool(processes = nthreads) as mp:
            for res in mp.starmap(
                univariate_coxph, 
                [(X, y, j, copy(m)) for j in range(n_features)]
            ):
                scores.append(res)
    else:
        for j in range(n_features):
            res = univariate_coxph(X, y, j, copy(m)) 
            scores.append(res)
    scores_ordered = np.empty(n_features)
    for i,v in scores:
        scores_ordered[i] = v
    return scores_ordered

if False:
    fs_scores = fit_and_score_features(X, y)
    fs_df = pd.DataFrame({'gene' : internal_data_dict['gene_ids'], 'score' : fs_scores})
    fs_df.to_csv(f"{res_path}gene_concordance_index.csv")

eln_model = CoxnetSurvivalAnalysis(max_iter = 100)
model = Pipeline(steps = [
    ('fs', SelectKBest(score_func = fit_and_score_features)), 
    #('scale', StandardScaler()),
    ('eln', copy(eln_model))])
param_grid = {
    'fs__k' : np.logspace(2, 3, 4).astype('int64'), 
    'eln__alphas' : np.expand_dims(10**np.linspace(-3, -1, 4), (-1)).tolist(), 
    'eln__l1_ratio' : np.power(np.linspace(0.1**4, 0.99**4, 3), 1./4.)}

#%%

fs_res = nested_survival_eval(
    model = model, 
    param_grid = param_grid, 
    X = X, 
    y = y, 
    y_mask = internal_data_dict['survival_mask'],
    X_covar = X_covar, 
    strat_var = strat_var, 
    scaler_instance = StandardScaler(),
    covar_scaler_instance = StandardScaler(),
    pca_instance = None, 
    cv_outer = copy(cv_outer), 
    cv_inner = copy(cv_inner), 
    thread_limit = nthreads)

fs_train_c = fs_res['train_c'].mean()
fs_test_c = fs_res['test_c'].mean()

print(f"train C: {fs_train_c}, test C: {fs_test_c}")

fs_res_df = pd.DataFrame({
    'fold' : fs_res['fold'], 
    f"{source_dataset}_train_c" : fs_res['train_c'], 
    f"{source_dataset}_test_c" : fs_res['test_c']})
fs_res_df.to_csv(res_path + 'fs_model_internal_cross_validation.csv')

#%% Final grid
if external_data_dict is not None:
    from baseline_survival_models import grid_search
    
    X_scaled_concat =  np.concatenate(
        (X_scaled[internal_data_dict['survival_mask'], :], 
         X_covar_scaled), 
        axis = 1)
    
    grid_res = grid_search(
        model = model, 
        X = X_scaled_concat, 
        y = y_masked, 
        param_grid = param_grid, 
        cv = cv_outer, 
        cv_args = {'groups' : [i for i, masked in zip(cvg_pairlist, internal_data_dict['survival_mask']) if masked]}, 
        thread_limit = nthreads)

#%% Best params and final model
if external_data_dict is not None:
    best_i = np.argmax(grid_res['test_score'])
    best_par_dict = dict([(k, grid_res[k][best_i]) for k in grid_res.keys()])
    best_par_train_score = best_par_dict.pop('train_score')
    best_par_test_score = best_par_dict.pop('test_score')
    
    final_model = copy(model)
    final_model = final_model.set_params(**best_par_dict)
    final_model = final_model.fit(X_scaled_concat, y_masked)

#%% eval
if external_data_dict is not None:
    X_ext_scaled_concat =  np.concatenate(
        (X_ext_scaled_masked, 
         X_ext_covar_scaled), 
        axis = 1)
    
    internal_score = final_model.score(X_scaled_concat, y_masked)
    external_score = final_model.score(X_ext_scaled_concat, y_ext_masked)
    print(f"{source_dataset} C: {internal_score}, {target_dataset} C: {external_score}")
    
    ext_fs_res_df = pd.DataFrame({
        f"{source_dataset}_c" : internal_score, 
        f"{target_dataset}_c" : external_score}, 
        index = pd.RangeIndex(0,1))
    ext_fs_res_df.to_csv(res_path + 'fs_model_external_validation.csv')


#%% Only covriates

model = CoxPHSurvivalAnalysis()

age_res = nested_survival_eval(
    model = model, 
    param_grid = {}, 
    X = X_covar, 
    y = y, 
    y_mask = internal_data_dict['survival_mask'],
    strat_var = strat_var, 
    scaler_instance = StandardScaler(),
    pca_instance = None, 
    cv_outer = copy(cv_outer), 
    cv_inner = copy(cv_inner), 
    thread_limit = nthreads)

final_model = copy(model)
final_model = final_model.fit(X_covar_scaled, y_masked)

internal_score = final_model.score(X_covar_scaled, y_masked)
age_res_dict = {
    f"{source_dataset}_train_c" : age_res['train_c'].mean(), 
    f"{source_dataset}_c" : internal_score
}
if external_data_dict is not None:
    external_score = final_model.score(X_ext_covar_scaled, y_ext_masked)
    age_res_dict[f"{source_dataset}_test_c"] = age_res['test_c'].mean()
    age_res_dict[f"{target_dataset}_c"] = external_score

age_res_df = pd.DataFrame(age_res_dict, index = pd.RangeIndex(0,1))
age_res_df.to_csv(res_path + 'covar_model_external_validation.csv')