# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import re

from sae.parsing_utilities import search_arg_generator

def get_kwargs(parameter_file):
    parameters = pd.read_csv(parameter_file, header = 0)
    parameters.drop(
        columns = [
            'serialized_data', 
            'ps_test_sets', 
            'ps_validation_sets', 
            'survival_evaluation_brier_times', 
            'nthreads_interop', 
            'nthreads', 
            'gpu_memory', 
            'parallel', 
            'file_name_prefix', 
            'data_dict'
        ], 
        inplace = True
    )
    parameters_dict = parameters.to_dict()
    parameters_dict = dict([(k, v[0]) for k, v in parameters_dict.items()])
    gym_args = parameters_dict.pop('gym_args', None)
    if gym_args is not None:
        gym_args_pairlist = []
        for i in re.sub('^{|}$', '', gym_args).split(', '):
            k, v = i.split(': ')
            if v == 'True':
                v = True
            elif v == 'False':
                v = False
            gym_args_pairlist.append((re.sub('^\'|\'$', '', k), v))
        gym_args_dict = dict(gym_args_pairlist)
        parameters_dict = {**parameters_dict, **gym_args_dict}
    data_batch_args = parameters_dict.pop('data_batch_args', None)
    if data_batch_args is not None:
        data_batch_args_pairlist = []
        for i in re.sub('^{|}$', '', data_batch_args).split(', '):
            k, v = i.split(': ')
            if v == 'True':
                v = True
            elif v == 'False':
                v = False
            data_batch_args_pairlist.append((re.sub('^\'|\'$', '', k), v))
        data_batch_args_dict = dict(data_batch_args_pairlist)
        parameters_dict = {**parameters_dict, **data_batch_args_dict}
    optimizer_args = parameters_dict.pop('optimizer_args', None) # Currently all identical
    if optimizer_args is not None:
        optimizer_args_pairlist = []
        for i in re.sub('^{|}$', '', optimizer_args).split(', '):
            k, v = i.split(': ')
            if v == 'True':
                v = True
            elif v == 'False':
                v = False
            k_ = re.sub('^\'|\'$', '', k)
            k__ = re.sub('_', '', k_)
            k___ = f"adam_{k__}"
            optimizer_args_pairlist.append((k___, v))
        optimizer_args_dict = dict(optimizer_args_pairlist)
        parameters_dict = {**parameters_dict, **optimizer_args_dict}
    adv_opt_class =  parameters_dict.pop('adversarial_optimizer_class', None)
    if adv_opt_class is not None:
        if 'optimizers.adam.Adam' in adv_opt_class:
            parameters_dict['adversarial_pre_training_optimizer'] = 'adam'
    cor_opt_class =  parameters_dict.pop('correction_optimizer_class', None)
    if cor_opt_class is not None:
        if 'optimizers.adam.Adam' in cor_opt_class:
            parameters_dict['batch_correction_pre_training_optimizer'] = 'adam'
    
    drug_response_model_output_activation = parameters_dict.get('drug_response_model_output_activation', '')
    if isinstance(drug_response_model_output_activation, float):
        # it must be nan since None converts to pandas.NA converts to np.nan
        drug_response_model_output_activation = ''
    
    '''
    Copied from main function of the superAE training
    Purpose is to parse args from the file into the approriate format and load the 
    appropriate dataset for training.
    '''
    # Re-format all args to CLI inputs ...
    arg_dict = {
        'encoder_layers' : re.sub('\\[|\\]', '', parameters_dict['encoder_layers']),
        'decoder_layers' : re.sub('\\[|\\]', '', parameters_dict['decoder_layers']),
        'classifier_layers' : re.sub('\\[|\\]', '', parameters_dict['classifier_layers']),
        'survival_model_layers' : re.sub('\\[|\\]', '', parameters_dict['survival_model_layers']),
        'batch_detector_layers' : re.sub('\\[|\\]', '', parameters_dict['batch_adversarial_model_layers']),
        'drug_response_model_layers' : re.sub('\\[|\\]', '', parameters_dict['drug_response_model_layers']),
        'drug_response_model_drugwise_layers' : re.sub('\\[|\\]', '', parameters_dict.get('drug_response_model_drugwise_layers', '')),
        'reg_weights' : format(parameters_dict['reg_a'], '.20f'),
        'mse' : parameters_dict['fix_var'], 
        'variational' : parameters_dict['variational'], 
        'reg_noise_sd' : format(parameters_dict['noise_sd'], '.20f'), 
        'reg_dropout_rate_input' : format(parameters_dict['dropout_input'], '.20f'), 
        'reg_dropout_rate_autoencoder' : format(parameters_dict['dropout_autoencoder'], '.20f'), 
        'reg_dropout_rate_survival' : format(parameters_dict['dropout_survival'], '.20f'), 
        'reg_dropout_rate_classifier' : format(parameters_dict['dropout_classifier'], '.20f'), 
        'reg_dropout_rate_batch' : format(parameters_dict['dropout_batch'], '.20f'), 
        'reg_dropout_rate_drug_response' : format(parameters_dict['dropout_drug_response'], '.20f'), 
        'reg_weights_type' : parameters_dict['reg_type'], 
        'classifier' : parameters_dict.get('classifier', False), 
        'survival_model' : parameters_dict['survival_model'], 
        'batch_correction' : parameters_dict['batch_adversarial_model'], 
        'batch_loss' : parameters_dict['batch_adversarial_loss_function'], 
        'deconfounder_layers_per_batch' : str(parameters_dict['deconfounder_layers_per_batch']), 
        'deconfounder_norm_penalty' : format(parameters_dict['deconfounder_norm_penalty'], '.20f'), 
        'deconfounder_centered_alignment' : parameters_dict['deconfounder_centered_alignment'], 
        'batch_adversarial_gradient_penalty' : format(parameters_dict['batch_adversarial_gradient_penalty'], '.20f'), 
        'drug_response_model' : parameters_dict['drug_response_model'], 
        'drug_response_model_output_activation' : drug_response_model_output_activation, 
        'reconstruction_weight' : format(float(parameters_dict.get('reconstruction_weight', '1.')), '.20f'), 
        'classifier_weight' : format(float(parameters_dict.get('classifier_weight', '1.')), '.20f'), 
        'survival_risk_weight' : format(float(parameters_dict.get('survival_risk_weight', '1.')), '.20f'), 
        'domain_adaptation_weight' : format(float(parameters_dict.get('domain_adaptation_weight', '1.')), '.20f'), 
        'drug_regression_weight' : format(float(parameters_dict.get('drug_regression_weight', '1.')), '.20f'), 
        'hidden_activation' : parameters_dict['hidden_activation'], 
        'fast_adaptive_multitask_optimization' : parameters_dict.get('fast_adaptive_multitask_optimization', False), 
        'decoder_split_loss' : parameters_dict.get('decoder_split_loss', False),
        'cache_data' : parameters_dict['cache'], 
        'prefetch_data' : parameters_dict['prefetch'], 
        'no_minibatches' : not parameters_dict['mini_batches'], # Applies to shuffle as well
        'minibatchsize' : int(parameters_dict['batch_size']), 
        'earlystop' : bool(parameters_dict['early_stopping']), 
        'print_period' : int(parameters_dict['print_period']), 
        'patience' : str(parameters_dict['patience']), 
        'classifier_metric' : parameters_dict.get('supervised_metric', 'cross-entropy'), 
        'debug_training' : False, 
        'debug_gym' : False, 
        'epochs' : int(parameters_dict['max_epochs']), 
        'epochs_pre_ae' : int(parameters_dict['max_epochs_pre_ae']), 
        'epochs_pre_cl' : int(parameters_dict['max_epochs_pre_cl']), 
        'epochs_pre_sr' : int(parameters_dict['max_epochs_pre_sr']), 
        'epochs_pre_bd' : int(parameters_dict['max_epochs_pre_bd']), 
        'epochs_pre_bc' : int(parameters_dict['max_epochs_pre_bc']), 
        'epochs_pre_dr' : int(parameters_dict['max_epochs_pre_dr']), 
        'learning_rate' : format(float(parameters_dict['learning_rate']), '.20f'), 
        'pre_learning_rate_ae' : format(float(parameters_dict['pre_learning_rate_ae']), '.20f'), 
        'pre_learning_rate_cl' : format(float(parameters_dict['pre_learning_rate_cl']), '.20f'), 
        'pre_learning_rate_sr' : format(float(parameters_dict['pre_learning_rate_sr']), '.20f'), 
        'pre_learning_rate_bd' : format(float(parameters_dict['pre_learning_rate_bd']), '.20f'), 
        'pre_learning_rate_bc' : format(float(parameters_dict['pre_learning_rate_bc']), '.20f'), 
        'pre_learning_rate_dr' : format(float(parameters_dict['pre_learning_rate_dr']), '.20f'), 
        'no_pre_train_ae' : not parameters_dict['pre_train_ae'], 
        'no_pre_train_cl' : not parameters_dict['pre_train_cl'], 
        'no_pre_train_sr' : not parameters_dict['pre_train_sr'], 
        'no_pre_train_bd' : not parameters_dict['pre_train_bd'], 
        'no_pre_train_bc' : not parameters_dict['pre_train_bc'], 
        'no_pre_train_dr' : not parameters_dict['pre_train_dr'], 
        'adam_beta1' : format(float(parameters_dict.get('adam_beta1', '0.9')), '.20f'), 
        'adam_beta2' : format(float(parameters_dict.get('adam_beta2', '0.999')), '.20f'), 
        'adversarial_pre_training_optimizer' : parameters_dict.get('adversarial_pre_training_optimizer', 'sgd'), 
        'batch_correction_pre_training_optimizer' : parameters_dict.get('batch_correction_pre_training_optimizer', 'sgd'), 
        'early_stopping_set_denominator' : int(parameters_dict['early_stopping_set_denominator']), 
        'early_stopping_set_stratify' : parameters_dict['early_stopping_set_stratify'], 
        'adversarial_learning_rate_multi' : format(
            float(parameters_dict['adversarial_learning_rate']) / 
            float(parameters_dict['learning_rate']), '.20f'),
        'cv_runs' : int(parameters_dict['nruns']), 
        'cv_folds' : int(parameters_dict['nfolds']), 
        'data_standardize' : parameters_dict['data_standardize'], # NOTE: not necessarily the same as used in setup
        'task_id' : int(parameters_dict['task_id']), 
        'o' : '',
        'p' : '', 
        'parallel' : False, 
        'gpu_memory' : -1, 
        'nthreads' : 8, 
        'nthreads_interop' : 2, 
        'survival_evaluation_brier_times' : '', 
        'no_ps_validation' : True, 
        'no_ps_test' : True
    }
    
    if 'objective_weights' in parameters_dict.keys():
        obj_str = re.sub('\s+', ',', re.sub('\\[|\\]', '', parameters_dict['objective_weights']).strip())
        obj_split = obj_str.split(',')
        arg_dict['reconstruction_weight'] = obj_split[0]
        arg_dict['classifier_weight'] = obj_split[1]
        arg_dict['survival_risk_weight'] = obj_split[2]
        arg_dict['domain_adaptation_weight'] = obj_split[3]
        arg_dict['drug_regression_weight'] = obj_split[4]
    
    if 'supervised' in parameters_dict.keys():
        arg_dict['classifier'] = parameters_dict['supervised']
    if 'supervised_metric' in parameters_dict.keys():
        arg_dict['supervised_metric'] = parameters_dict['supervised_metric']
    
    search_kwargs = search_arg_generator(
        args = arg_dict, 
        task_id = arg_dict['task_id'], 
        namespace = False)
    search_kwargs['data_standardize'] = True
    
    return search_kwargs

def result_file_dict(result_files):
    file_dict = {}
    for i in np.arange(len(result_files)):
        taski = re.sub('.*task', '', result_files[i])
        taski = re.sub('\\.csv(\\.gz)', '', taski)
        task_files = file_dict.get(taski, [])
        task_files.append(result_files[i])
        file_dict[taski] = task_files
    return file_dict

def load_results(result_files):
    result_list = []
    for f in result_files:
        try:
            result_list.append(pd.read_csv(f, header = 0))
        except Exception as e:
            print(f"Failed to read {f}:{str(e)}")
    
    for i in np.arange(len(result_list)):
        result_list[i].rename(columns = {'Unnamed: 0' : 'id'}, inplace = True)
    
    result_dict = {}
    for i in np.arange(len(result_list)):
        task_col = result_list[i].get('task', None)
        if task_col is None:
            taski = re.sub('.*task', '', result_files[i])
            taski = re.sub('\\.csv(\\.gz)', '', taski)
        else:
            taski = task_col[0]
        dfi = result_list[i]
        task_df = result_dict.get(taski, None)
        if task_df is None:
            result_dict[taski] = dfi
        else:
            result_dict[taski] = pd.concat((result_dict[taski], dfi), axis = 0)
    
    return result_dict