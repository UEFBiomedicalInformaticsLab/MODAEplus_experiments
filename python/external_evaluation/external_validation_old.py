# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 11:35:06 2022

@author: rintala
"""

import tensorflow as tf
import numpy as np
import pandas as pd
import time
import os
import warnings

import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC
import sklearn.metrics as metrics
if False:
    from variational_auto_encoder2 import VAE
    from training_utilities2 import train_aecl_with_pretraining, get_model_losses, process_data_depending_on_model, get_survival_metrics, get_classifier_metrics, train
    from data_utilities import rdata_loader
else:
    from sae.autoencoder_model import superAE
    from sae.training_utilities import train_aecl_with_pretraining, get_model_losses, process_data_depending_on_model, get_survival_metrics, get_classifier_metrics, train
    from sae.data_utilities import rdata_loader
from sys import platform

max_cv_repeats = 20
data_standardize = True
omics_layer = 'mrna'
#task_id = 1
external_validation = True
external_validation_classification = False
external_validation_survival = True
external_validation_combination = False
fail_settings = False
ga_top_settings = True

#%% Load data
if platform == 'linux':
    data_path = '/research/work/rintala/tcga/brca/intermediary_files/mrna_pw_mkkm_mr/'
    res_path = '/research/work/rintala/tcga/brca/intermediary_files/mrna/autoencoder/sr_only_no_es/'
else:
    data_path = '//research/workdir/tcga/brca/intermediary_files/mrna_pw_mkkm_mr/'
    res_path = '//research/workdir/tcga/brca/intermediary_files/mrna/autoencoder/sr_only_no_es/'
if False:
    data_path = '/home/teemu/research_work/tcga/brca/intermediary_files/mrna_pw_mkkm_mr/'
    res_path = '/home/teemu/research_work/tcga/brca/intermediary_files/mrna/autoencoder/sr_only_no_es/'

rdata_path = data_path + 'data_store.rda.gz'
#cv_path = data_path + 'cv_index.csv'
cv_path = data_path + 'stratified_cv_index.csv'
data_dict = rdata_loader(rdata_path = rdata_path, omics_layers = [omics_layer],
                             var_cutoff_rank = -1, #1000,
                             subtype_column = 'BRCA_Subtype_PAM50',
                             survival_data = True,
                             survival_event_column = 'OS',
                             survival_time_column = 'OS.time',
                             survival_id_column = 'sample_id',
                             survival_covariate_names = ['age_at_initial_pathologic_diagnosis'],
                             survival_covariates_onehot = [False])

cv = pd.read_csv(cv_path, header = 0, index_col = 0)

runs = np.unique(cv['run'])
folds = np.flip(np.unique(cv['fold'])) # Put last (usually reference) fold first

runs = runs[:min(runs.shape[0], max_cv_repeats)]

#%% unpack data
x = data_dict['x']
ind = data_dict['index']
y = data_dict['y']

survival_event = data_dict['survival_event']
survival_time = data_dict['survival_time']
survival_covariates = data_dict['survival_covariates']
survival_covariates_categorical = data_dict['survival_covariates_categorical']

#%% External validation settings
def layer_parser(input_str):
    out = input_str.replace('[', '')
    out = out.replace(']', '')
    return [int(i) for i in out.split(',')]
def weight_parser(input_str):
    out = input_str.replace('[', '')
    out = out.replace(']', '')
    return [float(i) for i in out.split('\r\n ')]

if external_validation:
    if False:
        import json
        handle_hof = open("//research/workdir/AESR_HPO/20230321-201127_ga_hof.json", "r")
        ga_hof = json.load(handle_hof)
        handle_hof.close()
        [i for i in ga_hof[0]['model_args'].keys() if ga_hof[0]['model_args'][i] != ga_hof[1]['model_args'][i]] # should not be indentical
        [i['model_args']['encoder_layers'] for i in ga_hof]
    if external_validation_classification:
        model_args = {}
        train_args = {}
    elif external_validation_survival:
        if fail_settings:
            fail_settings_table = pd.read_csv(res_path + 'fail_params.csv', header = 0, index_col = 0)
            faili = 0
            el = layer_parser(fail_settings_table['encoder_layers'].to_numpy()[faili])
            dl = layer_parser(fail_settings_table['decoder_layers'].to_numpy()[faili])
            cl = layer_parser(fail_settings_table['classifier_layers'].to_numpy()[faili])
            sl = layer_parser(fail_settings_table['survival_model_layers'].to_numpy()[faili])
            
            ra = fail_settings_table['reg_a'].to_numpy()[faili]
            fv = fail_settings_table['fix_var'].to_numpy()[faili]
            
            model_args = {'input_dim' : x[omics_layer].shape[1], 
                            'encoder_layers' : layer_parser(fail_settings_table['encoder_layers'].to_numpy()[faili]), 
                            'decoder_layers' : layer_parser(fail_settings_table['decoder_layers'].to_numpy()[faili]), 
                            'classifier_layers' : layer_parser(fail_settings_table['classifier_layers'].to_numpy()[faili]),
                            'survival_model_layers' : layer_parser(fail_settings_table['survival_model_layers'].to_numpy()[faili]), 
                            'reg_a' : fail_settings_table['reg_a'].to_numpy()[faili], 
                            'reg_type' : fail_settings_table['reg_type'].to_numpy()[faili], 
                            'fix_var' : fail_settings_table['fix_var'].to_numpy()[faili], 
                            'variational' : fail_settings_table['variational'].to_numpy()[faili], 
                            'input_noise' : fail_settings_table['input_noise'].to_numpy()[faili], 
                            'noise_sd' : fail_settings_table['noise_sd'].to_numpy()[faili], 
                            'dropout_input' : fail_settings_table['dropout_input'].to_numpy()[faili], 
                            'dropout_autoencoder' : fail_settings_table['dropout_autoencoder'].to_numpy()[faili], 
                            'dropout_survival' : fail_settings_table['dropout_survival'].to_numpy()[faili], 
                            'dropout_classifier' : fail_settings_table['dropout_classifier'].to_numpy()[faili], 
                            'supervised' : fail_settings_table['supervised'].to_numpy()[faili],
                            'survival_model' : fail_settings_table['survival_model'].to_numpy()[faili],
                            'survival_covariate_n' : fail_settings_table['survival_covariate_n'].to_numpy()[faili],
                            'objective_weights' : weight_parser(fail_settings_table['objective_weights'].to_numpy()[faili])}
            train_args = {'mini_batches' : fail_settings_table['mini_batches'].to_numpy()[faili],
                          'shuffle' : fail_settings_table['shuffle'].to_numpy()[faili],
                          'batch_size' : fail_settings_table['batch_size'].to_numpy()[faili],
                          'early_stopping' : fail_settings_table['early_stopping'].to_numpy()[faili],
                          'patience' : fail_settings_table['patience'].to_numpy()[faili], 
                          'return_losses' : True, 
                          'supervised_metric' : fail_settings_table['supervised_metric'].to_numpy()[faili]}
            stepwise_training_args = {'max_epochs' : fail_settings_table['max_epochs'].to_numpy()[faili],
                                      'max_epochs_pre_ae' : fail_settings_table['max_epochs_pre_ae'].to_numpy()[faili],
                                      'max_epochs_pre_cl' : fail_settings_table['max_epochs_pre_cl'].to_numpy()[faili],
                                      'max_epochs_pre_sr' : fail_settings_table['max_epochs_pre_sr'].to_numpy()[faili],
                                      'learning_rate' : fail_settings_table['learning_rate'].to_numpy()[faili],
                                      'pre_learning_rate_ae' : fail_settings_table['pre_learning_rate_ae'].to_numpy()[faili],
                                      'pre_learning_rate_cl' : fail_settings_table['pre_learning_rate_cl'].to_numpy()[faili],
                                      'pre_learning_rate_sr' : fail_settings_table['pre_learning_rate_sr'].to_numpy()[faili],
                                      'early_stopping_set_denominator' : fail_settings_table['early_stopping_set_denominator'].to_numpy()[faili],
                                      'early_stopping_set_stratify' : fail_settings_table['early_stopping_set_stratify'].to_numpy()[faili]}
            data_args = {'data_standardize' : fail_settings_table['data_standardize'].to_numpy()[faili]}
        elif ga_top_settings:
            model_args = {'input_dim' : x[omics_layer].shape[1], 
                            'encoder_layers' : [230, 48, 55], 
                            'decoder_layers' : [53, 154], 
                            'classifier_layers' : [], 
                            'reg_a' : 4.676354e-5, 
                            'fix_var' : True, 
                            'variational' : False, 
                            'input_noise' : False, 
                            'dropout_input' : 0., 
                            'dropout_autoencoder' : 0.03674864, 
                            'dropout_survival' : 0.07425686, 
                            'dropout_classifier' : 0.3551114, 
                            'reg_type' : 'L2', 
                            'supervised' : False,
                            'survival_model' : True,
                            'survival_model_layers' : [16], 
                            'survival_covariate_n' : 1,
                            'objective_weights' : [0.9273858818991182,1.024663356512016,0.36180380772786386]}
            train_args = {'mini_batches' : True,
                'shuffle' : True,
                'batch_size' : 256,
                'early_stopping' : False,
                'patience' : np.inf, #100, 
                'return_losses' : True, 
                'supervised_metric' : 'cross_entropy'}
            stepwise_training_args = {'max_epochs' : 1000,
                                      'max_epochs_pre_ae' : 1000,
                                      'max_epochs_pre_sr' : 1000,
                                      'learning_rate' : 7.21857e-05,
                                      'pre_learning_rate_ae' : 1.516784e-05,
                                      'pre_learning_rate_sr' : 2.169469e-04,
                                      'early_stopping_set_denominator' : 10,
                                      'early_stopping_set_stratify' : True}
            data_args = {'data_standardize' : True}
        else:
            model_args = {'input_dim' : x[omics_layer].shape[1], 
                            'encoder_layers' : [120, 63, 16], 
                            'decoder_layers' : [35, 91], 
                            'classifier_layers' : [], 
                            'reg_a' : 0.00004611802, 
                            'fix_var' : True, 
                            'variational' : False, 
                            'input_noise' : False, 
                            'dropout_input' : 0.441455, 
                            'dropout_autoencoder' : 0.1708341, 
                            'dropout_survival' : 0.0001738632, 
                            'dropout_classifier' : 0.07487409, 
                            'reg_type' : 'L2', 
                            'supervised' : False,
                            'survival_model' : True,
                            'survival_model_layers' : [15], 
                            'survival_covariate_n' : 1,
                            'objective_weights' : [1.37950554,0.85343513,0.76705933]}
            train_args = {'mini_batches' : True,
                'shuffle' : True,
                'batch_size' : 256,
                'early_stopping' : False,
                'patience' : np.inf, #100, 
                'return_losses' : True, 
                'supervised_metric' : 'cross_entropy'}
            stepwise_training_args = {'max_epochs' : 1000,
                                      'max_epochs_pre_ae' : 1000,
                                      'max_epochs_pre_sr' : 1000,
                                      'learning_rate' : 1.190456e-05,
                                      'pre_learning_rate_ae' : 1.667102e-05,
                                      'pre_learning_rate_sr' : 9.237329e-05,
                                      'early_stopping_set_denominator' : 5,
                                      'early_stopping_set_stratify' : True}
            data_args = {'data_standardize' : True}
    elif external_validation_combination:
        model_args = {'input_dim' : x[omics_layer].shape[1], 
                        'encoder_layers' : [79,52,7],#[151, 19, 16], 
                        'decoder_layers' : [51,101],#[29, 107], 
                        'reg_a' : 7.133015e-4,#1.130552e-5, 
                        'fix_var' : True, 
                        'variational' : False, 
                        'input_noise' : False, 
                        'dropout_input' : 0.68369,#0.4646343, 
                        'dropout_autoencoder' : 0.2106352,#0.006667702, 
                        'dropout_survival' : 0.06934017,#0.03302206, 
                        'dropout_classifier' : 0.03746057,#0.006328716, 
                        'reg_type' : 'L2', 
                        'supervised' : True,
                        'classifier_layers' : [16],#[14], 
                        'survival_model' : True,
                        'survival_model_layers' : [2],#[14], 
                        'survival_covariate_n' : 1,
                        'objective_weights' : [1.88854389,0.76794585,0.34351027]}#[0.7089602,0.82628461,1.46475519]}
        train_args = {'mini_batches' : True,
            'shuffle' : True,
            'batch_size' : 256,
            'early_stopping' : True,
            'patience' : 100, 
            'return_losses' : True, 
            'supervised_metric' : 'cross_entropy'}
        stepwise_training_args = {'max_epochs' : 1000,#10000,
                                  'max_epochs_pre_ae' : 1000,#10000,
                                  'max_epochs_pre_sr' : 1000,#10000,
                                  'learning_rate' : 9.730402e-5,#9.234168e-05,
                                  'pre_learning_rate_ae' : 6.837278e-5,#4.83561e-05,
                                  'pre_learning_rate_cl' : 2.662722e-5,#1.464779e-05,
                                  'pre_learning_rate_sr' : 2.293557e-5,#7.84269e-05,
                                  'early_stopping_set_denominator' : 5,
                                  'early_stopping_set_stratify' : True}
        data_args = {'data_standardize' : True}
    else:
        pass
else:
    model_args = {}
    train_args = {}
if False:
    import pickle
    handle_trajectory = open("//research/workdir/AESR_HPO/L2L/simulation/trajectories/trajectory_0_0.bin", "rb")
    trajectory = pickle.load(handle_trajectory)
    handle_trajectory.close()
    
    DropOut_scale = .01
    Layer_scale = 1.
    LogRegularization_scale = .1
    LogLearningRate_scale = .1
    ObjectiveWeight_scale = .01
    FailFitness = -1000.
    def translate_params(model_args, stepwise_training_args, model_param_dict):
        # Layer sizes (round)
        el_sizes = np.rint(model_param_dict['encoder_layers'] * Layer_scale)
        model_args['encoder_layers'] = el_sizes[el_sizes > 0].astype('int')
        dl_sizes = np.rint(model_param_dict['decoder_layers'] * Layer_scale)
        model_args['decoder_layers'] = dl_sizes[dl_sizes > 0].astype('int')
        cl_sizes = np.rint(model_param_dict['classifier_layers'] * Layer_scale)
        model_args['classifier_layers'] = cl_sizes[cl_sizes > 0].astype('int')
        sr_sizes = np.rint(model_param_dict['survival_model_layers'] * Layer_scale)
        model_args['survival_model_layers'] = sr_sizes[sr_sizes > 0].astype('int')
        
        # Dropout rates (model_args)
        do_in = model_param_dict['dropout_input'] * DropOut_scale
        model_args['dropout_input'] = np.clip(do_in, 0., 0.75)
        do_ae = model_param_dict['dropout_autoencoder'] * DropOut_scale
        model_args['dropout_autoencoder'] = np.clip(do_ae, 0., 0.75)
        do_sr = model_param_dict['dropout_survival'] * DropOut_scale
        model_args['dropout_survival'] = np.clip(do_sr, 0., 0.75)
        do_cl = model_param_dict['dropout_classifier'] * DropOut_scale
        model_args['dropout_classifier'] = np.clip(do_cl, 0., 0.75)
        
        # Objective weights (vector of length 3) (model_args)
        model_args['objective_weights'] = np.clip(model_param_dict['objective_weights'] * ObjectiveWeight_scale, 0, np.inf)
        
        # Regularization 
        model_args['reg_a'] = 10.**(model_param_dict['reg_a'] * LogRegularization_scale)
        
        # Learning rates (search_args)
        stepwise_training_args['learning_rate'] = 10.**(model_param_dict['learning_rate'] * LogLearningRate_scale)
        stepwise_training_args['pre_learning_rate_ae'] = 10.**(model_param_dict['pre_learning_rate_ae'] * LogLearningRate_scale)
        stepwise_training_args['pre_learning_rate_cl'] = 10.**(model_param_dict['pre_learning_rate_cl'] * LogLearningRate_scale)
        stepwise_training_args['pre_learning_rate_sr'] = 10.**(model_param_dict['pre_learning_rate_sr'] * LogLearningRate_scale)
        
        #return model_args, stepwise_training_args # not necessary since this is a procedure
    
    translate_params(model_args, stepwise_training_args, trajectory.individual)

np.random.seed(1) # first run seed
model_args['encoder_init_seeds'] = np.random.randint(2**31, size = len(model_args['encoder_layers']), dtype=np.int32)
model_args['decoder_init_seeds'] = np.random.randint(2**31, size = len(model_args['decoder_layers']), dtype=np.int32)
model_args['recon_init_seed'] = np.random.randint(2**31, size = 1, dtype=np.int32)[0]
model_args['classifier_init_seeds'] = np.random.randint(2**31, size = len(model_args['classifier_layers']), dtype=np.int32)
model_args['classifier_final_init_seed'] = np.random.randint(2**31, size = 1, dtype=np.int32)[0]
model_args['survival_model_init_seeds'] = np.random.randint(2**31, size = len(model_args['survival_model_layers']), dtype=np.int32)
model_args['survival_model_final_init_seed'] = np.random.randint(2**31, size = 1, dtype=np.int32)[0]
#%% Process (standardize omics + covariates)
dat_processed = process_data_depending_on_model(x, 
                                                cv_train_ind = np.full(x['mrna'].shape[0], True), 
                                                cv_test_ind = np.full(x['mrna'].shape[0], False), 
                                                supervised = model_args['supervised'],
                                                survival_model = model_args['survival_model'],
                                                y = y, 
                                                survival_event = survival_event,
                                                survival_time = survival_time, 
                                                survival_covariates = survival_covariates,
                                                survival_covariates_categorical = [False],
                                                data_standardize = data_args['data_standardize'])
dat_processed['x_train'] = dat_processed['x_train'][omics_layer]
dat_processed['x_test'] = dat_processed['x_test'][omics_layer]
#dat_processed['rownames_train'] = ind.index[cv_train_ind]
#dat_processed['rownames_test'] = ind.index[cv_test_ind]
dat_processed['labels'] = np.arange(dat_processed['y_n'])
#%% External data
if external_validation:
    ext_file_path = '//research/groups/fortino/luca/preprocessed scanb/With 0.7 zeros filter/'
    sweco_x = pd.read_csv(ext_file_path + 'mrna.csv', index_col = 0, header = 0)
    sweco_y = pd.read_csv(ext_file_path + 'outcome.csv', index_col = 0, header = 0)
    
    genes = data_dict['cols'][omics_layer]
    ext_genes = sweco_x.columns.to_numpy()
    
    genes_not_in = np.setdiff1d(genes, ext_genes)
    genes_in = np.intersect1d(genes, ext_genes)
    genes_not_in_tcga = np.setdiff1d(ext_genes, genes)
    
    ext_x = {'mrna' : sweco_x.loc[:,genes_in].to_numpy(dtype = np.float32)}
    
    ext_y_labels, ext_y = np.unique(sweco_y['Pam50'].to_numpy(), return_inverse=True)
    ext_survival_event = sweco_y['OverallSurv'].to_numpy()
    ext_survival_event = ext_survival_event == 1
    ext_survival_time = sweco_y['SurvDays'].to_numpy()
    
    ext_survival_covariates = sweco_y.loc[:,['Age']].to_numpy()
    
    np.all(ext_y_labels == data_dict['y_labels'][:ext_y_labels.shape[0]]) # same order
    
    # Use processing function to do data scaling
    ext_processed = process_data_depending_on_model(ext_x, 
                                                    cv_train_ind = np.full(ext_x['mrna'].shape[0], True), 
                                                    cv_test_ind = np.full(ext_x['mrna'].shape[0], False), 
                                                    supervised = model_args['supervised'],
                                                    survival_model = model_args['survival_model'],
                                                    y = ext_y, 
                                                    survival_event = ext_survival_event,
                                                    survival_time = ext_survival_time, 
                                                    survival_covariates = ext_survival_covariates,
                                                    survival_covariates_categorical = [False],
                                                    data_standardize = True)
    # Manually apply age scaling (so that it matches with TCGA)
    ext_processed['survival_covariates_train'] = dat_processed['survival_std_scalers'][0].transform(ext_survival_covariates).astype('float32')
    ext_processed['x_train'] = ext_processed['x_train'][omics_layer]
    ext_processed['x_test'] = ext_processed['x_test'][omics_layer]
    #ext_processed['rownames_train'] = ind.index[cv_train_ind]
    #ext_processed['rownames_test'] = ind.index[cv_test_ind]
    ext_processed['labels'] = np.arange(ext_processed['y_n'])

#%% Train model using TCGA

stepwise_training_args['x'] = dat_processed['x_train']
stepwise_training_args['y'] = dat_processed['y_train']
stepwise_training_args['events'] = dat_processed['survival_event_train']
stepwise_training_args['times'] = dat_processed['survival_time_train']
stepwise_training_args['covariates'] = dat_processed['survival_covariates_train']


#model_args['objective_weights'][0] = model_args['objective_weights'][0] * dat_processed['x_train'].shape[1] # issue?
#model_args['activation'] = 'elu'
stepwise_training_args['model_args'] = model_args
train_args['print_period'] = 1
stepwise_training_args['train_args'] = train_args

# Limit epochs for testing
if False:
    stepwise_training_args['max_epochs'] = 10
    stepwise_training_args['max_epochs_pre_ae'] = 10
    stepwise_training_args['max_epochs_pre_sr'] = 10
    stepwise_training_args['max_epochs_pre_cl'] = 10

res = train_aecl_with_pretraining(**stepwise_training_args)

model = res['model']
z_mean = model.encode(dat_processed['x_train']).numpy()
x_recon = model(dat_processed['x_train']).numpy()
train_mse = metrics.mean_squared_error(dat_processed['x_train'], x_recon)
if model_args['survival_model']:
    train_survival_metrics = get_survival_metrics(model = model, 
                                                  z_mean = z_mean, 
                                                  covariates = dat_processed['survival_covariates_train'],
                                                  surv_event = dat_processed['survival_event_train'],
                                                  surv_time = dat_processed['survival_time_train'],
                                                  z_mean_train = None,
                                                  covariates_train = None, 
                                                  surv_event_train = None, 
                                                  surv_time_train = None, 
                                                  brier_times = np.arange(365,2400+1), 
                                                  efron = True)
if model_args['supervised']:
    train_classification_metrics = get_classifier_metrics(model = model, 
                                                          z_mean = z_mean, 
                                                          y = dat_processed['y_train'], 
                                                          labels = dat_processed['labels'])
    y_pred = model.classifier_net(z_mean)
    if False:
        pd.DataFrame(y_pred.numpy()).to_csv('C:/Users/rintala/Documents/vae_scanb_pred.csv')
        pd.DataFrame(z_mean).to_csv('C:/Users/rintala/Documents/vae_scanb_zmean.csv')
        pd.DataFrame(ext_y).to_csv('C:/Users/rintala/Documents/vae_scanb_ytrue.csv')

#%% Plot diagnostics
def plotstuff(y, color_label):
    fig = plt.figure(figsize=(10,10))
    for i in np.arange(y.shape[1] -1):
        for j in np.arange(i+1, y.shape[1]):
            ax = fig.add_subplot(y.shape[1]-1, y.shape[1]-1, i * (y.shape[1]-1) + j)
            ax.scatter(y[:, i], y[:, j], c=color_label, marker = '+')
    plt.show()

if False:
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(res['train_losses_ae']['reconstruction_loss'].to_numpy())
    ax.set_title('AE pre training reconstruction loss')
    
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(res['valid_losses_ae']['reconstruction_loss'].to_numpy())
    ax.set_title('AE pre validation reconstruction loss')
    
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(res['valid_losses_ae']['regularization'].to_numpy())
    ax.set_title('AE pre regularization loss')
    
    if model.supervised:
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(res['train_losses_cl']['reconstruction_loss'].to_numpy())
        ax.set_title('CL pre training reconstruction loss')
        
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(res['valid_losses_cl']['reconstruction_loss'].to_numpy())
        ax.set_title('CL pre validation reconstruction loss')
        
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(res['valid_losses_cl']['regularization'].to_numpy())
        ax.set_title('CL pre regularization loss')
        
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(res['train_losses_cl']['cross_entropy'].to_numpy())
        ax.set_title('CL pre training cross-entropy')
        
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(res['valid_losses_cl']['cross_entropy'].to_numpy())
        ax.set_title('CL pre validation cross-entropy')
    
    if model.survival_model:
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(res['train_losses_sr']['survival_log_likelihood'].to_numpy())
        ax.set_title('SR pre training Cox log-likelihood')
        
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(res['valid_losses_sr']['survival_log_likelihood'].to_numpy())
        ax.set_title('SR pre validation Cox log-likelihood')
    
    
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(res['train_losses']['reconstruction_loss'].to_numpy())
    ax.set_title('Training reconstruction loss')
    
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(res['valid_losses']['reconstruction_loss'].to_numpy())
    ax.set_title('Validation reconstruction loss')
    
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(res['valid_losses']['regularization'].to_numpy())
    ax.set_title('Regularization loss')
    
    if model.supervised:
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(res['train_losses']['cross_entropy'].to_numpy())
        ax.set_title('Training cross-entropy')
        
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(res['valid_losses']['cross_entropy'].to_numpy())
        ax.set_title('Validation cross-entropy')
    
    if model.survival_model:
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(res['train_losses']['survival_log_likelihood'].to_numpy())
        ax.set_title('Training Cox log-likelihood')
        
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(res['valid_losses']['survival_log_likelihood'].to_numpy())
        ax.set_title('Validation Cox log-likelihood')
    
    #plotstuff(z.numpy(), y)
#%% Bad weights analysis
if False:
    model_weights = model.get_weights()
    model.inference_net.summary()
    model.generative_net.summary()
    model.classifier_net.summary()
    model.survival_risk_net.summary()


#%% Deconstruct the training loop

if False:
    return_pre_trained_model_weights = False
    #model_args['input_dim'] = dat_processed['x_train'].shape[1]
    if model_args['supervised']:
        y_u = np.unique(y)
        y_n = y_u.shape[0]
        if np.any(y_u < 0):
            y_n -= 1
        model_args['class_number'] = y_n
    
    model = VAE(**model_args)
    
    train_args['x_train'] = dat_processed['x_train']
    train_args['x_valid'] = None
    
    train_args['y_train'] = dat_processed['y_train']
    train_args['y_valid'] = None
    
    if model.survival_model:
        train_args['event_train'] = dat_processed['survival_event_train']
        train_args['event_valid'] = None
        train_args['time_train'] = dat_processed['survival_time_train']
        train_args['time_valid'] = None
        train_args['covariates_train'] = dat_processed['survival_covariates_train']
        train_args['covariates_valid'] = None
    
    train_args['model'] = model
    train_args['verbose'] = False
    
    if True:
        pre_train_ae_optimizer = tf.keras.optimizers.Adam(stepwise_training_args['pre_learning_rate_ae'])
        train_args_ae = train_args.copy()
        train_args_ae['optimizer'] = pre_train_ae_optimizer
        train_args_ae['epochs'] = stepwise_training_args['max_epochs_pre_ae']
        train_args_ae['train_ae_only'] = True
        
        train_losses_ae, valid_losses_ae = train(**train_args_ae)
        if return_pre_trained_model_weights:
            model_ae_pre_weights = model.get_weights()
        else:
            model_ae_pre_weights = None
    else:
        train_losses_ae = None
        valid_losses_ae = None
        model_ae_pre_weights = None
    
    if model.supervised and True:
        pre_train_cl_optimizer = tf.keras.optimizers.Adam(stepwise_training_args['pre_learning_rate_cl'])
        train_args_cl = train_args.copy()
        train_args_cl['optimizer'] = pre_train_cl_optimizer
        train_args_cl['epochs'] = stepwise_training_args['max_epochs_pre_cl']
        train_args_cl['train_classifier_only'] = True
        
        train_losses_cl, valid_losses_cl = train(**train_args_cl)
        if return_pre_trained_model_weights:
            model_cl_pre_weights = model.get_weights()
        else:
            model_cl_pre_weights = None
    else:
        train_losses_cl = None 
        valid_losses_cl = None
        model_cl_pre_weights = None
    
    from copy import copy
    model_ae_pre_weights = copy(model.get_weights())
    model.set_weights(model_ae_pre_weights)
    
    '''
    Deconstruct survival model training loop
    '''
    train_ae_only = False
    train_cl_only = False
    train_sr_only = True
    if model.survival_model and True:
        pre_train_sr_optimizer = tf.keras.optimizers.Adam(stepwise_training_args['pre_learning_rate_sr'])
        train_args_sr = train_args.copy()
        train_args_sr['optimizer'] = pre_train_sr_optimizer
        train_args_sr['epochs'] = stepwise_training_args['max_epochs_pre_sr']
        train_args_sr['train_survival_only'] = True
        '''
        Training loop for survival pre training
        '''
        supervised_metric_sign = {'cross_entropy' : 1., 'auroc' : -1., 'aupr' : -1., 'acc' : -1.}
        # TODO: do something smarter? This just duplicates the work now
        from training_utilities2 import training_feeder
        if train_args_sr['x_valid'] is None or train_args_sr['x_valid'].shape[0] == 0:
            train_args_sr['x_valid'] = train_args_sr['x_train']
        train_dataset_x = training_feeder(train_args_sr['x_train'], 
                                          train_args_sr['mini_batches'], 
                                          train_args_sr['batch_size'], 
                                          train_args_sr['shuffle'], 
                                          seed = 1)
        valid_dataset_x = training_feeder(train_args_sr['x_valid'], 
                                          train_args_sr['mini_batches'], 
                                          train_args_sr['batch_size'], 
                                          train_args_sr['shuffle'], 
                                          seed = 1)
        
        if model.supervised and not (False or True) :
            # TODO: do something smarter? This just duplicates the work now
            if train_args_sr['y_valid'] is None or train_args_sr['y_valid'].shape[0] == 0:
                train_args_sr['y_valid'] = train_args_sr['y_train']
            train_dataset_y = training_feeder(train_args_sr['y_train'], 
                                              train_args_sr['mini_batches'], 
                                              train_args_sr['batch_size'], 
                                              train_args_sr['shuffle'], 
                                              seed = 1)
            valid_dataset_y = training_feeder(train_args_sr['y_valid'], 
                                              train_args_sr['mini_batches'], 
                                              train_args_sr['batch_size'], 
                                              train_args_sr['shuffle'], 
                                              seed = 1)
        else:
            train_dataset_y = [None]
            valid_dataset_y = [None]
        
        if model.survival_model and not (False or False):
            # TODO: make some more checks here
            # TODO: do something smarter? This just duplicates the work now
            if train_args_sr['event_valid'] is None or train_args_sr['event_valid'].shape[0] == 0:
                train_args_sr['event_valid'] = train_args_sr['event_train']
            train_dataset_event = training_feeder(train_args_sr['event_train'], 
                                              train_args_sr['mini_batches'], 
                                              train_args_sr['batch_size'], 
                                              train_args_sr['shuffle'], 
                                              seed = 1)
            valid_dataset_event = training_feeder(train_args_sr['event_valid'], 
                                              train_args_sr['mini_batches'], 
                                              train_args_sr['batch_size'], 
                                              train_args_sr['shuffle'], 
                                              seed = 1)
            if train_args_sr['time_valid'] is None or train_args_sr['time_valid'].shape[0] == 0:
                train_args_sr['time_valid'] = train_args_sr['time_train']
            train_dataset_time = training_feeder(train_args_sr['time_train'], 
                                              train_args_sr['mini_batches'], 
                                              train_args_sr['batch_size'], 
                                              train_args_sr['shuffle'], 
                                              seed = 11)
            valid_dataset_time = training_feeder(train_args_sr['time_valid'], 
                                              train_args_sr['mini_batches'], 
                                              train_args_sr['batch_size'], 
                                              train_args_sr['shuffle'], 
                                              seed = 1)
            if train_args_sr['covariates_valid'] is None or train_args_sr['covariates_valid'].shape[0] == 0:
                train_args_sr['covariates_valid'] = train_args_sr['covariates_train']
            train_dataset_covariates = training_feeder(train_args_sr['covariates_train'], 
                                              train_args_sr['mini_batches'], 
                                              train_args_sr['batch_size'], 
                                              train_args_sr['shuffle'], 
                                              seed = 1)
            valid_dataset_covariates = training_feeder(train_args_sr['covariates_valid'], 
                                              train_args_sr['mini_batches'], 
                                              train_args_sr['batch_size'], 
                                              train_args_sr['shuffle'], 
                                              seed = 1)
        else:
            train_dataset_event = [None]
            valid_dataset_event = [None]
            train_dataset_time = [None]
            valid_dataset_time = [None]
            train_dataset_covariates = [None]
            valid_dataset_covariates = [None]
        
        if train_args_sr['early_stopping']:
            best_loss = np.Inf
            disappointments = 0
        
        if train_args_sr['return_losses']:
            train_losses = []
            valid_losses = []
        
        if train_args_sr['verbose']:
            if model.variational:
                loss_string = 'ELBO'
            else:
                loss_string = 'MSE'
        
        from copy import copy
        if True: #early_stopping: # NaN protection is always useful
            # Initialize (used if error is nan or always increases, i.e. in case of a bug ...)
            best_model_weights = copy(model.get_weights())
        from itertools import zip_longest
        for epoch in range(1, train_args_sr['epochs'] + 1):
            if train_args_sr['verbose']:
                start_time = time.time()
            # Training for one epoch
            nan_counter = 0
            for x, y, events, times, covariates in zip_longest(train_dataset_x, 
                                                               train_dataset_y, 
                                                               train_dataset_event, 
                                                               train_dataset_time, 
                                                               train_dataset_covariates): 
                if train_ae_only:
                    z = model.ae_pre_train_step(x, train_args_sr['optimizer'], return_z = True)
                elif train_cl_only:
                    z = model.classifier_pre_train_step(x, train_args_sr['optimizer'], y, return_z = True)
                elif train_sr_only:
                    #z = model.survival_model_pre_train_step(x, train_args_sr['optimizer'], y, times, events, covariates, efron = True, return_z = True)
                    model.inference_net.trainable = False
                    model.generative_net.trainable = False
                    if model.supervised:
                        model.classifier_net.trainable = False
                    if model.variational:
                        z_mean, z_logvar = model.encode(x)
                        z = model.reparameterize(z_mean, z_logvar)
                    else:
                        z = model.encode(x, training = False)
                    with tf.GradientTape() as tape:
                        surv_loglikelihood, n_events = model.compute_survival_loss_coxph(z = z, 
                                                                                        times = times, 
                                                                                        events = events, 
                                                                                        covariates = covariates,
                                                                                        efron = True, 
                                                                                        training = True)
                        loss = -tf.reduce_sum(surv_loglikelihood) / n_events * model.objective_weights[2]
                        loss += tf.math.add_n(model.losses) # regularization
                    gradients = tape.gradient(loss, model.trainable_variables)
                    if not tf.math.reduce_all(tf.math.is_finite(surv_loglikelihood)):
                        nan_counter += 1
                        break
                        if False:
                            #print(times.shape) 
                            #print(events.shape)
                            # Filter out missing data
                            time_filter = times >= 0 # negative are assumed to be missing (due to pre-processing)
                            #time_filter = tf.cast(time_filter, dtype = 'float32')
                            
                            z = tf.boolean_mask(z, time_filter, axis = 0)
                            times = tf.boolean_mask(times, time_filter, axis = 0)
                            events = tf.boolean_mask(events, time_filter, axis = 0)
                            if not covariates is None:
                                covariates = tf.boolean_mask(covariates, time_filter, axis = 0)
                            
                            # Create input (x) from AE embedding + covariates
                            if model.survival_covariate_n > 0:
                                z_a = tf.concat([z, covariates], axis = 1)
                            else:
                                z_a = z
                            # Compute log hazards
                            log_hazard = model.survival_risk_net(z_a, training = True)
                            theta = tf.math.exp(log_hazard)
                            
                            # Calculate at risk indicator matrix (tj >= ti)
                            # Each column indicates at risk group for sample i
                            at_risk = tf.expand_dims(times, 1) >= tf.expand_dims(times, 0)
                            #at_risk = times >= tf.transpose(times)
                            at_risk = tf.cast(at_risk, dtype = 'float32')
                            
                            # Implement boolean mask as multiplication
                            #events_mask = tf.transpose(events)
                            events_mask = tf.expand_dims(events, 0)
                            events_mask = tf.cast(events_mask, dtype = 'float32')
                            
                            # At risk group total hazard for each sample
                            total_risk = tf.matmul(tf.transpose(theta), at_risk) + 1e-16
                            
                            # Extract events only
                            total_risk_events = total_risk * events_mask
                            
                            # Sum of event log hazards is always needed for first term
                            log_hazard_events = tf.transpose(log_hazard) * events_mask
                            #total_log_hazard = tf.reduce_sum(log_hazard_events)
                            
                            if True:
                                # Similarly for tied events (tj == ti)
                                #tied = times == tf.transpose(times)
                                tied =  tf.expand_dims(times, 1) == tf.expand_dims(times, 0)
                                tied = tf.cast(tied, dtype = 'float32')
                                
                                # Risk of tied events for Efron's adjustment
                                theta_events = tf.transpose(theta) * events_mask #tf.boolean_mask(theta, events)
                                tied_events = tied * events_mask #tf.boolean_mask(tied, events, axis = 0)
                                tied_events = tied_events * tf.transpose(events_mask) #tf.boolean_mask(tied_events, events, axis = 1)
                                # To avoid issue with ratio calculations we set the diagonal to 1 (results in 0 nominator for non-events)
                                #tied_events = tf.linalg.set_diag(tied_events, tf.ones_like(times))
                                tied_events = tied_events + tf.linalg.tensor_diag(1-events_mask[0,:])
                                total_tied_risk = tf.matmul(theta_events, tied_events)[0,:] + 1e-16
                                
                                # Total risk and total tied risk are equal for all tied events.
                                # And in practice want to sum all of them, but tied events need to be
                                # adjusted based on the number of tied events. Since total_tied_risk 
                                # contains the |h| repeats of the sum, we can make a vector for the 
                                # nominator and denominator for all events. Multiplying them element-
                                # wise and summing we end up with the partial log likelihood.
                                
                                # The denominator is |h| which can be computed from tied_events matrix
                                denominator = tf.reduce_sum(tied_events, axis = 0)
                                # Nominator can be computed using the upper triangular part of the matrix - 1
                                nominator = tf.reduce_sum(tf.linalg.band_part(tied_events, 0, -1), axis = 0) - 1
                                
                                #term2 = tf.reduce_sum(tf.math.log(total_risk_events[0,:] + 1 - events_mask[0,:] - total_tied_risk * nominator/denominator)) # 1 - events_mask ensures log(1) = 0 for non-events
                                #log_likelihood = total_log_hazard - term2
                                term2 = tf.math.log(total_risk_events[0,:] + (1 - events_mask[0,:]) - total_tied_risk * nominator/denominator) # 1 - events_mask ensures log(1) = 0 for non-events
                                log_likelihood = log_hazard_events - term2
                            else:
                                #log_likelihood = total_log_hazard - tf.reduce_sum(tf.math.log(total_risk_events))
                                log_likelihood = log_hazard_events - tf.math.log(total_risk_events)
                            n_events = tf.reduce_sum(events_mask)
                    if tf.math.reduce_any([tf.math.reduce_any(tf.math.is_nan(i)) for i in gradients]):
                        nan_counter += 1
                        break
                    if not tf.math.reduce_all([tf.math.reduce_all(tf.math.is_finite(i)) for i in gradients]):
                        nan_counter += 1
                        break
                    train_args_sr['optimizer'].apply_gradients(zip(gradients, model.trainable_variables))
                    model.inference_net.trainable = True
                    model.generative_net.trainable = True
                    if model.supervised:
                        model.classifier_net.trainable = True
                else:
                    z = model.train_step(x, train_args_sr['optimizer'], y, times, events, covariates, efron = True, return_z = True)
                if tf.math.reduce_any(tf.math.is_nan(z)):
                    nan_counter += 1
            if train_args_sr['verbose']:
                end_time = time.time()
            
            if ((train_args_sr['verbose'] or train_args_sr['return_losses']) and epoch % train_args_sr['print_period'] == 0) or train_args_sr['early_stopping']:
                #loss_W = tf.math.add_n(model.losses)
                #valid_loss = tf.keras.metrics.Mean()
                '''
                valid_kwargs = {'x_data' : valid_dataset_x, 'y_data' : valid_dataset_y, 
                                }
                '''
                metrics_valid = get_model_losses(model, x_data = valid_dataset_x, 
                                                 y_data = valid_dataset_y, 
                                                 survival_event_data = valid_dataset_event,
                                                 survival_time_data = valid_dataset_time, 
                                                 survival_covariate_data = valid_dataset_covariates,
                                                 no_supervised = train_ae_only | train_sr_only, 
                                                 no_survival = train_ae_only | train_cl_only,
                                                 get_metrics = train_args_sr['supervised_metric'] != 'cross_entropy')
                # Print validation set metrics
                if train_args_sr['verbose'] and epoch % train_args_sr['print_period'] == 0:
                    monitor_string = ('Epoch: {}, Test set {}: {}, W: {}').format(epoch,
                                                                                  loss_string,
                                                                                  metrics_valid['reconstruction_loss'],
                                                                                  metrics_valid['regularization'])
                    if model.supervised and not (train_ae_only or train_sr_only):
                        monitor_string += (', predictions {}: {}').format(train_args_sr['supervised_metric'], 
                                                                            metrics_valid[train_args_sr['supervised_metric']])
                    if model.survival_model and not (train_ae_only or train_cl_only):
                        monitor_string += (', survival log-likelihood: {}').format(metrics_valid['survival_log_likelihood'])
                    monitor_string += (', time elapsed for current epoch {}').format(end_time - start_time)
                    print(monitor_string)
                if train_args_sr['return_losses']:
                    metrics_train = get_model_losses(model, x_data = train_dataset_x, 
                                                 y_data = train_dataset_y, 
                                                 survival_event_data = train_dataset_event,
                                                 survival_time_data = train_dataset_time, 
                                                 survival_covariate_data = train_dataset_covariates,
                                                 no_supervised = train_ae_only | train_sr_only, 
                                                 no_survival = train_ae_only | train_cl_only,
                                                 get_metrics = train_args_sr['supervised_metric'] != 'cross_entropy')
                    train_losses.append(pd.DataFrame(metrics_train, index = [0]))
                    valid_losses.append(pd.DataFrame(metrics_valid, index = [0]))
                    '''
                    losses_train_i = [metrics_train['reconstruction_loss'], loss_W]
                    losses_valid_i = [metrics_valid['reconstruction_loss'], loss_W]
                    if model.supervised and not (train_ae_only or train_survival_only):
                        losses_train_i += [metrics_train[supervised_metric]]
                        losses_valid_i += [metrics_valid[supervised_metric]]
                    if model.survival_model and not (train_ae_only or train_classifier_only):
                        losses_train_i += [metrics_train['survival_log_likelihood']]
                        losses_valid_i += [metrics_valid['survival_log_likelihood']]
                    train_losses.append(losses_train_i)
                    valid_losses.append(losses_valid_i)
                    #train_losses = pd.concat(train_losses, pd.DataFrame.from_dict(metrics_train), axis = 0)
                    '''
                if train_args_sr['early_stopping']:
                    total_loss = metrics_valid['reconstruction_loss'] * model.objective_weights[0]
                    if model.supervised and not (train_ae_only or train_sr_only):
                        total_loss += metrics_valid[train_args_sr['supervised_metric']] * supervised_metric_sign[train_args_sr['supervised_metric']] * model.objective_weights[1] # Will be different from gradient update if not cross-entropy
                    if model.survival_model and not (train_ae_only or train_cl_only):
                        total_loss -= metrics_valid['survival_log_likelihood'] * model.objective_weights[2]
                    if total_loss < best_loss:
                        best_loss = total_loss
                        best_losses = metrics_valid
                        best_model_weights = copy(model.get_weights())
                        disappointments = 0
                    elif disappointments < train_args_sr['patience']:
                        disappointments += 1
                    else:
                        if train_args_sr['verbose']:
                            monitor_string = ('Epoch: {}, Test set {}: {}, W: {}').format(epoch,
                                                                                          loss_string,
                                                                                          best_losses['reconstruction_loss'],
                                                                                          best_losses['regularization'])
                            if model.supervised and not (train_ae_only or train_sr_only):
                                monitor_string += (', predictions {}: {}').format(train_args_sr['supervised_metric'], 
                                                                                    best_losses[train_args_sr['supervised_metric']])
                            if model.survival_model and not (train_ae_only or train_cl_only):
                                monitor_string += (', survival log-likelihood: {}').format(best_losses['survival_log_likelihood'])
                            monitor_string += (', time elapsed for current epoch {}').format(end_time - start_time)
                            print(monitor_string)
                        break
                elif nan_counter == 0:
                    # Save model weights at every successfull iteration
                    best_model_weights = copy(model.get_weights())
                if nan_counter > 0:
                    break
        if train_args_sr['early_stopping'] or nan_counter > 0:
            #model.set_weights(best_model_weights)
            pass
        if train_args_sr['return_losses']:
            train_losses = pd.concat(train_losses, axis = 0)
            valid_losses = pd.concat(valid_losses, axis = 0)
            train_losses['iteration'] = np.arange(train_losses.shape[0]) + 1
            valid_losses['iteration'] = np.arange(valid_losses.shape[0]) + 1
        else:
            train_losses = None
            valid_losses = None
        #train_losses_sr, valid_losses_sr = train(**train_args_sr)
        #if return_pre_trained_model_weights:
        #    model_sr_pre_weights = model.get_weights()
        #else:
        #    model_sr_pre_weights = None
    else:
        train_losses_sr = None 
        valid_losses_sr = None
        model_sr_pre_weights = None
    
    # Train final model with both objectives
    optimizer = tf.keras.optimizers.Adam(stepwise_training_args['learning_rate'])
    train_args['optimizer'] = optimizer
    train_args['epochs'] = stepwise_training_args['max_epochs']
    
    train_losses, valid_losses = train(**train_args)

#%% Test on SCAN-B
# Pad missing genes with 0
ext_padded_data = np.full((ext_processed['x_train'].shape[0],) + genes.shape, fill_value = 0., dtype = 'float32')
ext_padded_data = pd.DataFrame(ext_padded_data, columns = genes)
ext_padded_data.loc[:,genes_in] = ext_processed['x_train']

z_mean = model.encode(ext_padded_data.to_numpy()).numpy()
x_recon = model(ext_padded_data.to_numpy()).numpy()

recon_unpadded_data = pd.DataFrame(x_recon, columns = genes)
test_mse = metrics.mean_squared_error(ext_processed['x_train'], recon_unpadded_data.loc[:,genes_in].to_numpy())
if model_args['survival_model']:
    test_survival_metrics = get_survival_metrics(model = model, 
                                                 z_mean = z_mean, 
                                                 covariates = ext_processed['survival_covariates_train'],
                                                 surv_event = ext_processed['survival_event_train'],
                                                 surv_time = ext_processed['survival_time_train'],
                                                 z_mean_train = None,
                                                 covariates_train = None, 
                                                 surv_event_train = None, 
                                                 surv_time_train = None, 
                                                 brier_times = np.arange(365,2400+1), 
                                                 efron = True)
if model_args['supervised']:
    test_classification_metrics = get_classifier_metrics(model = model, 
                                                         z_mean = z_mean, 
                                                         y = ext_processed['y_train'], 
                                                         labels = ext_processed['labels'])
#%% Baseline model PCA + Cox Elastic net
from sksurv.linear_model import CoxPHSurvivalAnalysis, CoxnetSurvivalAnalysis
from sksurv.metrics import concordance_index_censored, integrated_brier_score
from sklearn.model_selection import GridSearchCV, KFold, RepeatedKFold
from sklearn.decomposition import PCA

pca = PCA(10)
x_common = pd.DataFrame(dat_processed['x_train'], columns = genes)
x_pca = pca.fit_transform(x_common.loc[:,genes_in].to_numpy())
ext_x_pca = pca.transform(ext_processed['x_train'])

mse_train = metrics.mean_squared_error(x_common.loc[:,genes_in].to_numpy(), 
                                       pca.inverse_transform(x_pca))
mse_test = metrics.mean_squared_error(ext_processed['x_train'], 
                                      pca.inverse_transform(ext_x_pca))

x_combined = np.concatenate([x_pca, dat_processed['survival_covariates_train']], axis = 1)
ext_x_combined = np.concatenate([ext_x_pca, ext_processed['survival_covariates_train']], axis = 1)

y = np.array([(dat_processed['survival_event_train'][i], dat_processed['survival_time_train'][i]) for i in range(dat_processed['survival_event_train'].shape[0])],
             dtype = [('event', '?'),('time', '<i4')])
ext_y = np.array([(ext_processed['survival_event_train'][i], ext_processed['survival_time_train'][i]) for i in range(ext_processed['survival_event_train'].shape[0])],
                 dtype = [('event', '?'),('time', '<i4')])

# Filter out missing times
time_filter = y['time'] >= 0
ext_time_filter = ext_y['time'] >= 0
x_combined = x_combined[time_filter,]
ext_x_combined = ext_x_combined[ext_time_filter,]
y = y[time_filter,]
ext_y = ext_y[ext_time_filter,]

alphas = 10**np.linspace(-2,0, 3*10+1)
cox_model = CoxnetSurvivalAnalysis(l1_ratio = 0.9, max_iter = 1000, fit_baseline_model = True)
inner_repeated_cv = RepeatedKFold(n_splits=5, n_repeats = 1, random_state=0)
gcv_cox = GridSearchCV(
    cox_model,
    param_grid={'alphas' : [[a] for a in alphas]},
    cv=inner_repeated_cv,
    error_score=0.5,
    n_jobs=5,
    return_train_score=True)
gcv_cox = gcv_cox.fit(x_combined, y)
best_a_index = np.argwhere(alphas == gcv_cox.best_params_['alphas'][0])[0,0]
cox_train_c = gcv_cox.cv_results_['mean_train_score'][best_a_index]

risk_pred = gcv_cox.predict(x_combined)
train_c = concordance_index_censored(y['event'], y['time'], risk_pred)

ext_risk_pred = gcv_cox.predict(ext_x_combined)
test_c = concordance_index_censored(ext_y['event'], ext_y['time'], ext_risk_pred)

survs = gcv_cox.best_estimator_.predict_survival_function(x_combined)
times = np.arange(365, 2400)
preds = np.asarray([[fn(t) for t in times] for fn in survs])
score = integrated_brier_score(y, y, preds, times)

survs = gcv_cox.best_estimator_.predict_survival_function(ext_x_combined)
times = np.arange(365, 2400)
preds = np.asarray([[fn(t) for t in times] for fn in survs])


#%% Baseline model PCA + SVM
from sklearn.model_selection import GridSearchCV, KFold, RepeatedKFold
from sklearn.decomposition import PCA
from sklearn.svm import SVC

pca = PCA(10)
x_common = pd.DataFrame(dat_processed['x_train'], columns = genes)
x_pca = pca.fit_transform(x_common.loc[:,genes_in].to_numpy())
ext_x_pca = pca.transform(ext_processed['x_train'])

mse_train = metrics.mean_squared_error(x_common.loc[:,genes_in].to_numpy(), 
                                       pca.inverse_transform(x_pca))
mse_test = metrics.mean_squared_error(ext_processed['x_train'], 
                                      pca.inverse_transform(ext_x_pca))

model = SVC()
gs = GridSearchCV(model, 
                  param_grid = {'C' : np.logspace(-1, 3, num = 21), 
                                'gamma' : np.logspace(-5, -1, num = 21)},
                  scoring = 'balanced_accuracy',
                  n_jobs = 6)
y_notmissing = dat_processed['y_train'] >= 0
gs.fit(x_pca[y_notmissing, :], dat_processed['y_train'][y_notmissing])

ext_y_notmissing = ext_processed['y_train'] >= 0
y_pred = gs.predict(x_pca[y_notmissing,:])
ext_y_pred = gs.predict(ext_x_pca[ext_y_notmissing,:])

train_f1_macro = metrics.f1_score(dat_processed['y_train'][y_notmissing], y_pred, average = 'macro')
test_f1_macro = metrics.f1_score(ext_processed['y_train'][ext_y_notmissing], ext_y_pred, average = 'macro')
train_bacc = metrics.balanced_accuracy_score(dat_processed['y_train'][y_notmissing], y_pred)
test_bacc = metrics.balanced_accuracy_score(ext_processed['y_train'][ext_y_notmissing], ext_y_pred)
