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
    #res_path = '/research/work/rintala/superAE_HPO/20231124_random_search/pancan_test/'
else:
    #res_path = '//research/workdir/superAE_HPO/20230821_random_search/brca_test_noclfilter/'
    #res_path = '//research/workdir/superAE_HPO/20231027_random_search/brca_test_noclfilter/'
    #res_path = '//research/workdir/superAE_HPO/20231101_random_search/brca_test_noclfilter/'
    #res_path = '//research/workdir/superAE_HPO/20231103_random_search/brca_test_noclfilter/'
    #res_path = '//research/workdir/superAE_HPO/20231120_random_search/brca_test_noclfilter/'
    res_path = '//research/workdir/superAE_HPO/20231121_random_search/brca_test_noclfilter/'
    #res_path = '//research/workdir/superAE_HPO/20231124_random_search/pancan_test/'
if False:
    #res_path = '/home/teemu/research_work/superAE_HPO/20230821_random_search/brca_test_noclfilter/'
    #res_path = '/home/teemu/research_work/superAE_HPO/20231103_random_search/brca_test_noclfilter/'
    #res_path = '/home/teemu/research_work/superAE_HPO/20231120_random_search/brca_test_noclfilter/'
    res_path = '/home/teemu/research_work/superAE_HPO/20231121_random_search/brca_test_noclfilter/'
    #res_path = '/home/teemu/research_work/superAE_HPO/20231124_random_search/pancan_test/'

# Should select by name
fn = glob.glob(res_path + '../plots/*best_task.txt')
with open(fn[0], 'r') as f:
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
    from sklearn.preprocessing import StandardScaler
    patient = data_dict['patient_exp']
    ss = StandardScaler()
    patient_scaled = ss.fit_transform(patient)
    
    cl = data_dict['cl_exp']
    cl_scaled = ss.transform(cl)
    
    patient_mean = patient.mean(axis = 0)
    cl_mean = cl.mean(axis = 0)
    
    sns.scatterplot(x  = patient_mean, y = cl_mean)
    sns.scatterplot(x  = patient_scaled.max(axis = 0), y = cl_scaled.max(axis = 0))
    sns.scatterplot(x  = patient_scaled.min(axis = 0), y = cl_scaled.min(axis = 0))
    
    ss_cl = StandardScaler()
    cl_scaled = ss_cl.fit_transform(cl)
    patient_scaled = ss_cl.transform(patient)
    
    sns.scatterplot(x  = patient_scaled.max(axis = 0), y = cl_scaled.max(axis = 0))
    sns.scatterplot(x  = patient_scaled.min(axis = 0), y = cl_scaled.min(axis = 0))
    
    from dataset_collections import get_bruna_pdtc, get_bruna_pdtx

    pdtc_data_dict = get_bruna_pdtc(home = False)
    pdtx_data_dict = get_bruna_pdtx(home = False)

    # Filter to common genes
    for dd in [pdtc_data_dict, pdtx_data_dict]:
        og_expression = dd['cl_exp']
        
        reordering_map = dict(zip(data_dict['gene_ids'], np.arange(data_dict['gene_ids'].shape[0])))
        new_gex_mat = np.full(
            (dd['cl_exp'].shape[0], 
            data_dict['gene_ids'].shape[0]), 0.)
        for i,g in enumerate(dd['gene_ids']):
            tcga_ind = reordering_map.get(g, None)
            if tcga_ind:
                new_gex_mat[:,tcga_ind] = dd['cl_exp'][:,i]
        dd['cl_exp'] = new_gex_mat# * 0.
    
    pdtc = pdtc_data_dict['cl_exp']
    
    ss_cl = StandardScaler()
    cl_scaled = ss_cl.fit_transform(cl)
    ss_pdtc = StandardScaler()
    pdtc_scaled = ss_pdtc.fit_transform(pdtc)
    
    from sklearn.decomposition import PCA
    
    sns.scatterplot(x  = cl.mean(axis = 0), y = pdtc.mean(axis = 0))


if False:
    from dataset_collections import get_tcga_pancan_ctrp_ccle_solid
    data_dict = get_tcga_pancan_ctrp_ccle_solid(home = False)
    
    import seaborn as sns
    
    sns.histplot(data_dict['patient_exp'].flatten(), bins = 30)
    sns.histplot(data_dict['cl_exp'].flatten(), bins = 30)
    
    # Play around with scaling
    from sklearn.preprocessing import StandardScaler
    test = data_dict['patient_exp']
    ss = StandardScaler()
    test_ss = ss.fit_transform(test)
    
    sns.histplot(test_ss.max(axis = 0), bins = 30)
    sns.histplot(test_ss.min(axis = 0), bins = 30)
    
    test_mean = test.mean(axis = 0)
    test_std = test.std(axis = 0)
    
    sns.scatterplot(x  = test_mean, y = test_std)
    
    sns.scatterplot(x  = np.sqrt(test_mean), y = test_std)
    
    sns.scatterplot(x = test_std, y = test_ss.max(axis = 0))
    sns.scatterplot(x = test_std, y = np.abs(test_ss).max(axis = 0))
    sns.scatterplot(x = test_mean, y = np.abs(test_ss).max(axis = 0))
    
    sns.histplot(test_mean)
    sns.histplot(test_std)
    
    test_poiss = (test - test_mean) / np.expand_dims(np.sqrt(test_mean) + 1, axis = 0)
    
    sns.histplot(test_poiss.max(axis = 0), bins = 30)
    sns.histplot(test_poiss.min(axis = 0), bins = 30)
    sns.histplot(test_poiss.flatten(), bins = 30)
    sns.scatterplot(x = test_std, y = np.abs(test_poiss).max(axis = 0))
    
    sns.scatterplot(x = test_std[test_mean < 1.], y = test.max(axis = 0)[test_mean < 1.])
    
    test_std_offset = (test - np.expand_dims(test_mean, axis = 0)) / np.expand_dims(test_std + 1, axis = 0)
    sns.histplot(test_std_offset.flatten(), bins = 30)
    sns.scatterplot(x = test_std, y = test_std_offset.max(axis = 0))
    sns.scatterplot(x = test_std, y = np.abs(test_std_offset).max(axis = 0))
    sns.scatterplot(x = test_mean, y = np.abs(test_std_offset).max(axis = 0))
    
    
    test2 = data_dict['cl_exp']
    test2_mean = test2.mean(axis = 0)
    test2_std = test2.std(axis = 0)
    
    ss2 = StandardScaler()
    test2_ss = ss2.fit_transform(test)
    
    sns.scatterplot(x = test_mean, y = test2_mean)
    sns.scatterplot(x = test_std, y = test2_std)
    
    sns.scatterplot(x = test2_std, y = np.abs(test2_ss).max(axis = 0))
    sns.scatterplot(x = test2_mean, y = np.abs(test2_ss).max(axis = 0))
    sns.histplot(test2_mean)
    sns.histplot(test2_std)
    
    test_filter = test_mean > 1.
    test2_filter = test2_mean > 1.
    
    comb_filter = np.logical_and(test_filter, test2_filter)
    
    sns.scatterplot(x = test_mean[comb_filter], y = test2_mean[comb_filter])
    sns.scatterplot(x = test_std[comb_filter], y = test2_std[comb_filter])
    
    
    test_poiss = (test - test_mean) / np.expand_dims(np.log2(np.sqrt(np.power(2, test_mean)) + 1), axis = 0)
    sns.scatterplot(x = test_std, y = np.abs(test_poiss).max(axis = 0))
    sns.scatterplot(x = test_mean, y = np.abs(test_poiss).max(axis = 0))
    sns.scatterplot(x = test_poiss.mean(axis = 0), y = test_poiss.std(axis = 0))
    
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

if False:
    search_kwargs['data_standardize'] = False

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

if True:
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

if False:
    data_batch_args['repeat'] = True
    
    from sae.data_batch import batch_to_step_args
    from sae.data_utilities import dataset_batch_setup
    import sys
    
    train_args['patient_data'] = train_args['train_dataset']
    train_args['cell_line_data'] = train_args['train_dataset2']
    
    patient_datasize = data_dict['patient_exp'].shape[0]
    cl_datasize = data_dict['cl_exp'].shape[0]
    
    # Only one dataset can be longest, shorter ones are repeated
    patient_repeat = patient_datasize < cl_datasize
    cl_repeat = patient_datasize > cl_datasize
    
    data_batch_args['repeat'] = patient_repeat
    patient_sd = dataset_batch_setup(tcga_patient_dataset, **data_batch_args)
    
    data_batch_args['repeat'] = cl_repeat
    cl_sd = dataset_batch_setup(ccle_cl_dataset, **data_batch_args)
    
    
    train_args['train_dataset'] = patient_sd
    train_args['train_dataset2'] = cl_sd
    
    dataset = train_args['train_dataset']
    dataset2 = train_args['train_dataset2']
    
    tuple_dataset = False
    
    if dataset2 is not None:
        iterator = tf.data.Dataset.zip((dataset, dataset2))
    else:
        iterator = tf.data.Dataset.zip((dataset,))
    
    shape_str = ''
    a = 0
    for batch in iterator:
        step_args2 = batch_to_step_args(
            batch, 
            survival = True, 
            tuple_dataset = tuple_dataset)
        step_args1 = batch_to_step_args(
            batch, 
            batches = True, 
            tuple_dataset = tuple_dataset)
        step_args = batch_to_step_args(
            batch, 
            tuple_dataset = tuple_dataset,
            drugs = True)
        step_arg_tensor_shapes = dict([(k, i.shape) for k,i in step_args2.items() if isinstance(i, tf.Tensor)])
        for k,i in step_arg_tensor_shapes.items():
            print('Survival evaluation step args \'' + str(k) + '\' shape: ' + str(i), file = sys.stderr)
        step_arg_tensor_shapes = dict([(k, i.shape) for k,i in step_args1.items() if isinstance(i, tf.Tensor)])
        for k,i in step_arg_tensor_shapes.items():
            print('AE and batch evaluation step args \'' + str(k) + '\' shape: ' + str(i), file = sys.stderr)
        step_arg_tensor_shapes = dict([(k, i.shape) for k,i in step_args.items() if isinstance(i, tf.Tensor)])
        for k,i in step_arg_tensor_shapes.items():
            print('Drug evaluation step args \'' + str(k) + '\' shape: ' + str(i), file = sys.stderr)
        #break
        shape_str += 'Shapes: '
        for i in np.arange(len(batch)):
            shape_str += 'd' + str(i) + '.shape: ' + str(batch[i]['exp'].shape) + '; '
        shape_str += '\n'
        a += 1
        '''
        if a == 2:
            break
        print(shape_str)
        sys.stdout.write(shape_str)
        sys.stdout.flush()
        '''
    print(shape_str)
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

#%% Predictions and embeddings

from sae.data_utilities import (
    parse_serialized_dataset_from_file, 
    dataset_batch_setup
)
from sae.evaluation import get_model_losses
from sae.model_utilities import get_embeddings, get_predictions


dataset = parse_serialized_dataset_from_file(
    serialized_file = serialized_data['cl']['filename'],
    feature_spec = serialized_data['cl']['feature_spec'])
rows = np.array(serialized_data['cl']['sample_info']['rownames'])
model_spec = serialized_data['cl']['model_spec']
datasize = len(serialized_data['cl']['sample_info']['rownames'])

dataset_batched = dataset_batch_setup(dataset, **data_batch_args)

z = get_embeddings(result['model'], dataset_batched, tuple_dataset = False)
z = pd.DataFrame(
    z, 
    index = rows,
    columns = ['z{}'.format(i+1) for i in range(z.shape[1])])
z['dataset'] = 'ccle'

p = get_predictions(result['model'], dataset_batched, tuple_dataset = False)
for pred, key in zip(p.values(), p.keys()):
    if pred.shape[1]:
        cols = [key + '_' + str(i) for i in np.arange(pred.shape[1])]
    else:
        cols = key
    p[key] = pd.DataFrame(p[key], index = rows, columns = cols)
p = pd.concat(p.values(), axis = 1)
p['dataset'] = 'ccle'

z.to_csv(
    f"{res_path}external_evaluation/internal_drug_response_validation_embeddings.csv.gz", 
    na_rep = 'NA', 
    header = True, 
    index = True)
p.to_csv(
    f"{res_path}external_evaluation/internal_drug_response_validation_predictions.csv.gz", 
    na_rep = 'NA', 
    header = True, 
    index = True)

#%% Bruna PDTC and PDTX gene-expression datasets
from dataset_collections import get_bruna_pdtc, get_bruna_pdtx

pdtc_data_dict = get_bruna_pdtc(home = False)
pdtx_data_dict = get_bruna_pdtx(home = False)

# Filter to common genes
for dd in [pdtc_data_dict, pdtx_data_dict]:
    og_expression = dd['cl_exp']
    
    reordering_map = dict(zip(data_dict['gene_ids'], np.arange(data_dict['gene_ids'].shape[0])))
    new_gex_mat = np.full(
        (dd['cl_exp'].shape[0], 
        data_dict['gene_ids'].shape[0]), 0.)
    for i,g in enumerate(dd['gene_ids']):
        tcga_ind = reordering_map.get(g, None)
        if tcga_ind:
            new_gex_mat[:,tcga_ind] = dd['cl_exp'][:,i]
    dd['cl_exp'] = new_gex_mat# * 0.

#%% Serialize external datasets

serialized_ext_data_list = []
for dd, sub_path in zip([pdtc_data_dict, pdtx_data_dict], ['bruna_pdtc/', 'bruna_pdtx/']):
    serialized_ext_data_path = tcga_ccle_serialized_data_path + '../' + sub_path
    ext_spec_file = f"{serialized_ext_data_path}serialized_data_spec.json"
    if False and os.path.exists(ext_spec_file):
        handle = open(ext_spec_file, 'r')
        serialized_ext_data_spec = json.load(handle)
        handle.close()
    else:
        os.makedirs(serialized_ext_data_path, exist_ok = True)
        model_args = search_kwargs['model_args']
        serialized_ext_data_spec = process_and_serialize(
            model_args = model_args, 
            data_dict = dd, 
            omics_layer = search_kwargs['omics_layer'], 
            data_standardize = search_kwargs['data_standardize'], 
            serialized_data_path = serialized_ext_data_path)
    
    serialized_ext_data_list.append(JSONFeatureSpecDecoder(serialized_ext_data_spec))

#%% Load serialized data
from sae.data_utilities import (
    parse_serialized_dataset_from_file, 
    dataset_batch_setup
)
from sae.evaluation import (
    get_model_losses, 
    get_embeddings, 
    get_predictions
)

embeddings = []
predictions = []
for edi, name in zip(serialized_ext_data_list, ['bruna_pdtc', 'bruna_pdtx']):
    ext_dataset = parse_serialized_dataset_from_file(
        serialized_file = edi['cl']['filename'],
        feature_spec = edi['cl']['feature_spec'])
    ext_rows = np.array(edi['cl']['sample_info']['rownames'])
    ext_model_spec = edi['cl']['model_spec']
    ext_datasize = len(edi['cl']['sample_info']['rownames'])
    
    ext_dataset_batched = dataset_batch_setup(ext_dataset, **data_batch_args)
    
    z = get_embeddings(result['model'], ext_dataset_batched, tuple_dataset = False)
    z = pd.DataFrame(
        z, 
        index = ext_rows,
        columns = ['z{}'.format(i+1) for i in range(z.shape[1])])
    z['dataset'] = name
    embeddings.append(z)
    
    p = get_predictions(result['model'], ext_dataset_batched, tuple_dataset = False)
    for pred, key in zip(p.values(), p.keys()):
        if pred.shape[1]:
            cols = [key + '_' + str(i) for i in np.arange(pred.shape[1])]
        else:
            cols = key
        p[key] = pd.DataFrame(p[key], index = ext_rows, columns = cols)
    p = pd.concat(p.values(), axis = 1)
    p['dataset'] = name
    predictions.append(p)

embeddings = pd.concat(embeddings, axis = 0)
predictions = pd.concat(predictions, axis = 0)
embeddings.to_csv(
    f"{res_path}external_evaluation/external_drug_response_validation_embeddings.csv.gz", 
    na_rep = 'NA', 
    header = True, 
    index = True)
predictions.to_csv(
    f"{res_path}external_evaluation/external_drug_response_validation_predictions.csv.gz", 
    na_rep = 'NA', 
    header = True, 
    index = True)
