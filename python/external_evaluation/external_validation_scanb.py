#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 11:05:33 2023

@author: teemu
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

bp_spec = importlib.util.spec_from_file_location('best_parameters', 'external_evaluation/best_parameters.py')
bp = importlib.util.module_from_spec(bp_spec)
sys.modules['best_parameters'] = bp
bp_spec.loader.exec_module(bp)

import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

#%% Initialize tensorflow
import tensorflow as tf

gpu_memory = None
nthreads_interop = 2
nthreads = 8

gpus = tf.config.list_physical_devices('GPU')
if len(gpus) and gpu_memory is not None:
    print('Found {} GPUs, using 1 logical GPU with {} MB memory.'.format(len(gpus), gpu_memory))
    '''Fixed allocation, less overhead and easy to calculate'''
    tf.config.set_logical_device_configuration(
            gpus[0], # assume only one gpu
            [tf.config.LogicalDeviceConfiguration(memory_limit = gpu_memory)]
    )
    logical_gpus = tf.config.list_logical_devices('GPU')
else:
    tf.config.threading.set_inter_op_parallelism_threads(nthreads_interop)
    tf.config.threading.set_intra_op_parallelism_threads(nthreads)

#%% Parameters
from sys import platform
if platform == 'linux':
    #res_path = '/research/work/rintala/superAE_HPO/20230821_random_search/brca_test_noclfilter/'
    #res_path = '/research/work/rintala/superAE_HPO/20231027_random_search/brca_test_noclfilter/'
    #res_path = '/research/work/rintala/superAE_HPO/20231101_random_search/brca_test_noclfilter/'
    #res_path = '/research/work/rintala/superAE_HPO/20231103_random_search/brca_test_noclfilter/'
    #res_path = '/research/work/rintala/superAE_HPO/20231120_random_search/brca_test_noclfilter/'
    res_path = '/research/work/rintala/superAE_HPO/20231121_random_search/brca_test_noclfilter/'
else:
    #res_path = '//research/workdir/superAE_HPO/20230821_random_search/brca_test_noclfilter/'
    #res_path = '//research/workdir/superAE_HPO/20231027_random_search/brca_test_noclfilter/'
    #res_path = '//research/workdir/superAE_HPO/20231101_random_search/brca_test_noclfilter/'
    #res_path = '//research/workdir/superAE_HPO/20231103_random_search/brca_test_noclfilter/'
    #res_path = '//research/workdir/superAE_HPO/20231120_random_search/brca_test_noclfilter/'
    res_path = '//research/workdir/superAE_HPO/20231121_random_search/brca_test_noclfilter/'
if False:
    #res_path = '/home/teemu/research_work/superAE_HPO/20230821_random_search/brca_test_noclfilter/'
    #res_path = '/home/teemu/research_work/superAE_HPO/20231103_random_search/brca_test_noclfilter/'
    #res_path = '/home/teemu/research_work/superAE_HPO/20231120_random_search/brca_test_noclfilter/'
    res_path = '/home/teemu/research_work/superAE_HPO/20231121_random_search/brca_test_noclfilter/'

with open(res_path + '../plots/brca_test_noclfilter_best_task.txt', 'r') as f:
    best_task = int(f.readlines()[0])

best_parameter_file = glob.glob(f"{res_path}*parameters_task{best_task}.csv")[0]

from dataset_collections import get_tcga_brca_ctrp_ccle_full

data_dict = get_tcga_brca_ctrp_ccle_full(home = False)

from best_parameters import get_kwargs

search_kwargs = get_kwargs(
    best_parameter_file = best_parameter_file, 
    best_task = best_task, 
    data_dict = data_dict)

if False:
    import seaborn as sns
    
    sns.histplot(data_dict['patient_exp'].flatten(), bins = 30)
    sns.histplot(data_dict['cl_exp'].flatten(), bins = 30)

#%% Serialization
from sae.data_utilities import (
    serialize_dataset_to_file, 
    process_patient_data, 
    process_cl_data, 
    FeatureSpecEncoder, 
    JSONFeatureSpecDecoder
)
import json 
from dataset_processing import process_and_serialize

if platform == 'linux':
    tcga_ccle_serialized_data_path = f"{res_path}external_evaluation/tcga_ccl/"
else:
    tcga_ccle_serialized_data_path = f"{res_path}external_evaluation/tcga_ccl/"
if False:
    tcga_ccle_serialized_data_path = '/home/teemu/Documents/superAE_temp/external_evaluation/tcga_ccl/'

spec_file = f"{tcga_ccle_serialized_data_path}serialized_data_spec.json"
if False and os.path.exists(spec_file):
    handle = open(spec_file, 'r')
    serialized_data_spec = json.load(handle)
    handle.close()
else:
    os.makedirs(tcga_ccle_serialized_data_path, exist_ok = True)
    model_args = search_kwargs['model_args']
    serialized_data_spec = process_and_serialize(
        model_args = model_args, 
        data_dict = data_dict, 
        omics_layer = search_kwargs['omics_layer'], 
        data_standardize = search_kwargs['data_standardize'], 
        serialized_data_path = tcga_ccle_serialized_data_path)

serialized_data = JSONFeatureSpecDecoder(serialized_data_spec)

#%% Setup training
model_args = search_kwargs['model_args']

np.random.seed(0)
if model_args.get('encoder_layers', None):
    model_args['encoder_init_seeds'] = np.random.randint(2**31, size = len(model_args.get('encoder_layers')), dtype=np.int32)
if model_args.get('decoder_layers', None):
    model_args['decoder_init_seeds'] = np.random.randint(2**31, size = len(model_args.get('decoder_layers')), dtype=np.int32)
model_args['recon_init_seed'] = np.random.randint(2**31, size = 1, dtype=np.int32)[0]
if model_args.get('classifier_layers', None):
    model_args['classifier_init_seeds'] = np.random.randint(2**31, size = len(model_args.get('classifier_layers')), dtype=np.int32)
model_args['classifier_final_init_seed'] = np.random.randint(2**31, size = 1, dtype=np.int32)[0]
if model_args.get('survival_model_layers', None):
    model_args['survival_model_init_seeds'] = np.random.randint(2**31, size = len(model_args.get('survival_model_layers')), dtype=np.int32)
model_args['survival_model_final_init_seed'] = np.random.randint(2**31, size = 1, dtype=np.int32)[0]
if model_args.get('batch_adversarial_model_layers', None):
    model_args['batch_adversarial_model_init_seeds'] = np.random.randint(2**31, size = len(model_args.get('batch_adversarial_model_layers')), dtype=np.int32)
model_args['batch_adversarial_model_pred_seed'] = np.random.randint(2**31, size = 1, dtype=np.int32)[0]
if model_args.get('drug_response_model_layers', None):
    model_args['drug_response_model_init_seeds'] = np.random.randint(2**31, size = len(model_args.get('drug_response_model_layers')), dtype=np.int32)
model_args['drug_response_model_pred_seed'] = np.random.randint(2**31, size = 1, dtype=np.int32)[0]

#%% Load serialized data
from sae.data_utilities import parse_serialized_dataset_from_file

data_instance = serialized_data

tcga_patient_dataset = parse_serialized_dataset_from_file(
    serialized_file = data_instance['patient']['filename'],
    feature_spec = data_instance['patient']['feature_spec'])
tcga_patient_rows = np.array(data_instance['patient']['sample_info']['rownames'])
tcga_patient_model_spec = data_instance['patient']['model_spec']
tcga_patient_datasize = len(data_instance['patient']['sample_info']['rownames'])

ccle_cl_dataset = parse_serialized_dataset_from_file(
    serialized_file = data_instance['cl']['filename'],
    feature_spec = data_instance['cl']['feature_spec'])
ccle_cl_rows = np.array(data_instance['cl']['sample_info']['rownames'])
ccle_cl_model_spec = data_instance['cl']['model_spec']
ccle_cl_datasize = len(data_instance['cl']['sample_info']['rownames'])

cl_in = ccle_cl_model_spec.pop('input_dim')
if cl_in != tcga_patient_model_spec['input_dim']:
    raise ValueError('Input dimensions for serialized patient and cell-line data do not match.')

data_model_spec = {**tcga_patient_model_spec, **ccle_cl_model_spec}

#%% Training
from sae.training_utilities import train_aecl_with_pretraining

model_args = copy(search_kwargs['model_args'])
data_batch_args = copy(search_kwargs['data_batch_args'])
train_args = copy(search_kwargs['train_args'])
gym_args = copy(search_kwargs['gym_args'])

if False:
    # Shorter training for testing
    gym_args['max_epochs'] = 10
    gym_args['max_epochs_pre_ae'] = 10
    gym_args['max_epochs_pre_cl'] = 10
    gym_args['max_epochs_pre_sr'] = 10
    gym_args['max_epochs_pre_bd'] = 10
    gym_args['max_epochs_pre_bc'] = 10
    gym_args['max_epochs_pre_dr'] = 10
    # Test variational survival
    #model_args['survival_variational'] = False
    model_args['variational'] = False

result = train_aecl_with_pretraining(
    patient_serialized_dataset = tcga_patient_dataset, 
    cl_serialized_dataset = ccle_cl_dataset, 
    model_args = {**model_args, **data_model_spec},
    data_batch_args = data_batch_args, 
    train_args = train_args, 
    patient_datasize = tcga_patient_datasize,
    cl_datasize = ccle_cl_datasize, 
    **gym_args, 
    return_pre_trained_model_weights = False)
data_batch_args['repeat'] = False

#%% Evaluation
from sae.data_utilities import dataset_batch_setup
from sae.evaluation import get_model_losses

tcga_patient_dataset_batched = dataset_batch_setup(tcga_patient_dataset, **data_batch_args)

tcga_res = get_model_losses(
    result['model'], 
    dataset = tcga_patient_dataset_batched, 
    get_metrics = True, 
    tuple_dataset = False, 
    manual_batch_flag = False)

#%% SCAN-B dataset new
from dataset_collections import get_scanb

scanb_data_dict = get_scanb(home = False)

# Filter to common genes
scanb_original_expression = scanb_data_dict['patient_exp']

reordering_map = dict(zip(data_dict['gene_ids'], np.arange(data_dict['gene_ids'].shape[0])))
new_scanb_gex_mat = np.full(
    (scanb_data_dict['patient_exp'].shape[0], 
    data_dict['gene_ids'].shape[0]), 0.)
for i,g in enumerate(scanb_data_dict['gene_ids']):
    tcga_ind = reordering_map.get(g, None)
    if tcga_ind:
        new_scanb_gex_mat[:,tcga_ind] = scanb_data_dict['patient_exp'][:,i]
#new_scanb_gex_mat[]
scanb_data_dict['patient_exp'] = new_scanb_gex_mat# * 0.

#%% Serialize SCAN-B
from sae.data_utilities import (
    serialize_dataset_to_file, 
    process_patient_data, 
    process_cl_data, 
    FeatureSpecEncoder, 
    JSONFeatureSpecDecoder
)
import json 

scanb_serialized_data_path = tcga_ccle_serialized_data_path + '../scanb/'
scanb_spec_file = f"{scanb_serialized_data_path}serialized_data_spec.json"
if False and os.path.exists(scanb_spec_file):
    handle = open(scanb_spec_file, 'r')
    scanb_serialized_data_spec = json.load(handle)
    handle.close()
else:
    os.makedirs(scanb_serialized_data_path, exist_ok = True)
    model_args = search_kwargs['model_args']
    scanb_serialized_data_spec = process_and_serialize(
        model_args = model_args, 
        data_dict = scanb_data_dict, 
        omics_layer = search_kwargs['omics_layer'], 
        data_standardize = search_kwargs['data_standardize'], 
        serialized_data_path = scanb_serialized_data_path)

serialized_external_data = JSONFeatureSpecDecoder(scanb_serialized_data_spec)

#%% Load serialized data
from sae.data_utilities import parse_serialized_dataset_from_file

external_data_instance = serialized_external_data

scanb_patient_dataset = parse_serialized_dataset_from_file(
    serialized_file = external_data_instance['patient']['filename'],
    feature_spec = external_data_instance['patient']['feature_spec'])
scanb_patient_rows = np.array(external_data_instance['patient']['sample_info']['rownames'])
scanb_patient_model_spec = external_data_instance['patient']['model_spec']
scanb_patient_datasize = len(external_data_instance['patient']['sample_info']['rownames'])

#%% Evaluation
from sae.data_utilities import dataset_batch_setup
from sae.evaluation import get_model_losses

scanb_patient_dataset_batched = dataset_batch_setup(scanb_patient_dataset, **data_batch_args)

scanb_res = get_model_losses(
    result['model'], 
    dataset = scanb_patient_dataset_batched, 
    get_metrics = True, 
    tuple_dataset = False, 
    manual_batch_flag = False)

#%% Save results

result_df = pd.DataFrame({
    'dataset' : ('TCGA_BRCA', 'SCANB'),
    'reconstruction_mse' : (tcga_res['reconstruction_loss_dataset1'], 
                            scanb_res['reconstruction_loss_dataset1']), 
    'regularization' : (tcga_res['regularization'], 
                        scanb_res['regularization']), 
    'survival_log_likelihood' : (tcga_res['survival_log_likelihood'], 
                                 scanb_res['survival_log_likelihood']), 
    'survival_concordance' : (tcga_res['surv_c'], 
                              scanb_res['surv_c']), 
    'confounder_alignment_norm' : (tcga_res['confounder_alignment_norm'], 
                                   scanb_res['confounder_alignment_norm']) 
    })
result_df.to_csv(f"{res_path}external_evaluation/external_survival_validation_res.csv")

