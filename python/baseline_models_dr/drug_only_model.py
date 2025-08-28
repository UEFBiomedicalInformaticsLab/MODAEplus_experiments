#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 15:11:58 2023

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

et_spec = importlib.util.spec_from_file_location('evaluation_tools', 'utilities/evaluation_tools.py')
et = importlib.util.module_from_spec(et_spec)
sys.modules['evaluation_tools'] = et
et_spec.loader.exec_module(et)

import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

pd.set_option('display.max_columns', None)

#%% Initialize tensorflow
import tensorflow as tf

gpu_memory = None
nthreads_interop = 2
nthreads = 4

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
#res_path = base_path + '20240301_random_search/scanb_test_noclfilter/'
#res_path = base_path + '20240302_random_search/scanb_test_noclfilter/'
#res_path = base_path + '20240303_random_search/scanb_test_noclfilter/'
res_path = base_path + '20240305_2_random_search/scanb_test_noclfilter/'

os.makedirs(f"{res_path}external_evaluation/", exist_ok = True)

# Should select by name
fn = glob.glob(res_path + '*best_task.txt')
with open(fn[0], 'r') as f:
    best_task = int(f.readlines()[0])

best_parameter_file = glob.glob(f"{res_path}*parameters_task{best_task}.csv")[0]

gene_preselection = False
from dataset_collections import get_ccle_solid
data_dict = get_ccle_solid(home = home, gene_preselection = gene_preselection)
internal_dataset = None
external_dataset = None

from evaluation_tools import get_kwargs

search_kwargs = get_kwargs(
    parameter_file = best_parameter_file, 
    data_dict = data_dict)

#search_kwargs['data_standardize'] = False

#%% Serialization
from sae.data_utilities import JSONFeatureSpecDecoder
import json, re
from dataset_processing import process_and_serialize

if gethostname() == 'teemu-pc':
    internal_serialized_data_path = '/home/teemu/Documents/superAE_temp/external_evaluation/internal/'
elif platform == 'linux':
    internal_serialized_data_path = f"{res_path}external_evaluation/internal/"
else:
    internal_serialized_data_path = f"{res_path}external_evaluation/internal/"

spec_file = f"{internal_serialized_data_path}serialized_data_spec.json"
if False and os.path.exists(spec_file):
    # Loading does not work with scalers
    handle = open(spec_file, 'r')
    serialized_data_spec = json.load(handle)
    handle.close()
else:
    os.makedirs(internal_serialized_data_path, exist_ok = True)
    model_args = search_kwargs['model_args']
    serialized_data_spec, std_scalers = process_and_serialize(
        model_args = model_args, 
        data_dict = data_dict, 
        omics_layer = search_kwargs['omics_layer'], 
        data_standardize = search_kwargs['data_standardize'], 
        serialized_data_path = internal_serialized_data_path, 
        return_scalers = True)

serialized_data = JSONFeatureSpecDecoder(serialized_data_spec)

if True:
    # Default back to using different scalers
    std_scalers = {}

if False:
    #search_kwargs['file_name_prefix'] = res_path
    spec_fn = glob.glob(res_path + '*_serialized_data_spec.json')[0]
    handle = open(spec_fn, 'r')
    serialized_data_cv = json.load(handle)
    handle.close()
    
    serialized_data_cv = JSONFeatureSpecDecoder(serialized_data_cv)
    
    serialized_data = {
        'patient' : serialized_data_cv['0']['1']['test_cv_']['patient']['train'], 
        'cl' : serialized_data_cv['0']['1']['test_cv_']['cell_line']['train']
    }
    
    serialized_data['patient']['filename'] = re.sub(
        '.*/', res_path, serialized_data['patient']['filename'])
    serialized_data['cl']['filename'] = re.sub(
        '.*/', res_path, serialized_data['cl']['filename'])

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
if model_args.get('drug_response_model_drugwise_layers', None):
    model_args['drug_response_model_drugwise_init_seeds'] = np.random.randint(2**31, size = len(model_args.get('drug_response_model_drugwise_layers')), dtype=np.int32)
model_args['drug_response_model_pred_seed'] = np.random.randint(2**31, size = 1, dtype=np.int32)[0]

#%% Load serialized data
from sae.data_utilities import parse_serialized_dataset_from_file

data_instance = copy(serialized_data)

if data_instance.get('patient', None) is not None:
    internal_patient_dataset = parse_serialized_dataset_from_file(
        serialized_file = data_instance['patient']['filename'],
        feature_spec = data_instance['patient']['feature_spec'])
    internal_patient_rows = np.array(data_instance['patient']['sample_info']['rownames'])
    internal_patient_model_spec = data_instance['patient']['model_spec']
    internal_patient_datasize = len(data_instance['patient']['sample_info']['rownames'])
else:
    internal_patient_dataset = None
    internal_patient_model_spec = {}
    internal_patient_datasize = None

if data_instance.get('cl', None) is not None:
    ccle_cl_dataset = parse_serialized_dataset_from_file(
        serialized_file = data_instance['cl']['filename'],
        feature_spec = data_instance['cl']['feature_spec'])
    ccle_cl_rows = np.array(data_instance['cl']['sample_info']['rownames'])
    ccle_cl_model_spec = data_instance['cl']['model_spec']
    ccle_cl_datasize = len(data_instance['cl']['sample_info']['rownames'])
    
    if internal_patient_model_spec.get('input_dim', None) is not None:
        cl_in = ccle_cl_model_spec.pop('input_dim')
        if cl_in != internal_patient_model_spec['input_dim']:
            raise ValueError('Input dimensions for serialized patient and cell-line data do not match.')
else:
    ccle_cl_dataset = None
    ccle_cl_model_spec = {}
    ccle_cl_datasize = None

data_model_spec = {**internal_patient_model_spec, **ccle_cl_model_spec}

#%% Training
from sae.training_utilities import train_aecl_with_pretraining

model_args = copy(search_kwargs['model_args'])
data_batch_args = copy(search_kwargs['data_batch_args'])
train_args = copy(search_kwargs['train_args'])
gym_args = copy(search_kwargs['gym_args'])

if False:
    # Shorter training for testing
    gym_args['max_epochs'] = 3
    gym_args['max_epochs_pre_ae'] = 3
    gym_args['max_epochs_pre_cl'] = 3
    gym_args['max_epochs_pre_sr'] = 3
    gym_args['max_epochs_pre_bd'] = 3
    gym_args['max_epochs_pre_bc'] = 3
    gym_args['max_epochs_pre_dr'] = 3
    
    gym_args['pre_learning_rate_bc'] = 1e-5
    #gym_args['pre_train_ae'] = False
    #gym_args['pre_train_bd'] = False
    
    model_args['deconfounder_norm_penalty'] = 1.
    #model_args['deconfounder_centered_alignment'] = True
    # Test variational survival
    #model_args['survival_variational'] = False
    #model_args['variational'] = True
    
    data_batch_args['batch_size'] = 128
    #data_batch_args['prefetch'] = True
    train_args['return_losses'] = False
    
    # New drugwise models
    model_args['drug_response_model_layers'] = [20]#[model_args['drug_response_model_layers'][0]]
    model_args['drug_response_model_drugwise_layers'] = [4]
    
    # drug only
    model_args['decoder_model'] = False
    model_args['batch_adversarial_model'] = False
    model_args['survival_model'] = False
    gym_args['pre_train_ae'] = False
    gym_args['pre_train_bd'] = False
    gym_args['pre_train_bc'] = False
    gym_args['pre_train_cl'] = False
    gym_args['pre_train_sr'] = False
    gym_args['pre_train_dr'] = False
    
    

result = train_aecl_with_pretraining(
    patient_serialized_dataset = internal_patient_dataset, 
    cl_serialized_dataset = ccle_cl_dataset, 
    model_args = {**model_args, **data_model_spec},
    data_batch_args = data_batch_args, 
    train_args = train_args, 
    patient_datasize = internal_patient_datasize,
    cl_datasize = ccle_cl_datasize, 
    **gym_args, 
    return_pre_trained_model_weights = True)

if False:
    model_final_weights = result['model'].get_weights()
    [(i.shape, np.min(i), np.max(i)) for i in result['model_bc_pre_weights']]
    result['train_losses_ae']
    result['valid_losses_ae']
    result['train_losses_bd']
    result['valid_losses_bd']
    result['train_losses_bc']
    result['valid_losses_bc']
    result['train_losses_sr']
    result['valid_losses_sr']
    result['train_losses_dr']
    result['valid_losses_dr']
    result['train_losses']
    result['valid_losses']

data_batch_args['repeat'] = False

#%% Diagnostics
diagnostics_keys = [
    'train_losses_ae', 
    'valid_losses_ae',
    'train_losses_cl',
    'valid_losses_cl', 
    'train_losses_sr',
    'valid_losses_sr', 
    'train_losses_dr',
    'valid_losses_dr', 
    'train_losses_bd',
    'valid_losses_bd', 
    'train_losses_bc',
    'valid_losses_bc', 
    'train_losses', 
    'valid_losses']
diagnostics = [result[i] for i in diagnostics_keys if result.get(i, None) is not None]
diagnostics_stage = [i for i in diagnostics_keys if result.get(i, None) is not None]
for d, s in zip(diagnostics, diagnostics_stage):
    d['stage'] = s
if len(diagnostics) > 0:
    diagnostics = pd.concat(diagnostics, axis = 0)
    diagnostics.to_csv(
        f"{res_path}external_evaluation/diagnostics.csv.gz", 
        na_rep = 'NA', 
        header = True, 
        index = False)

#%% Evaluation
from sae.data_utilities import dataset_batch_setup
from sae.evaluation import get_model_losses

internal_patient_dataset_batched = dataset_batch_setup(internal_patient_dataset, **data_batch_args)

internal_res = get_model_losses(
    result['model'], 
    dataset = internal_patient_dataset_batched, 
    get_metrics = True, 
    tuple_dataset = False, 
    manual_batch_flag = False)

#%% Predictions and embeddings 

from sae.model_utilities import get_embeddings, get_predictions

ccle_cl_dataset_batched = dataset_batch_setup(ccle_cl_dataset, **data_batch_args)

model = result['model']

datasets = [
    {
        'dataset' : internal_patient_dataset_batched, 
        'rows' : internal_patient_rows, 
        'name' : internal_dataset, 
        'file_string' : f"{res_path}external_evaluation/internal_survival_validation_"
    },
    {
        'dataset' : ccle_cl_dataset_batched, 
        'rows' : ccle_cl_rows, 
        'name' : 'ccle', 
        'file_string' : f"{res_path}external_evaluation/internal_drug_response_validation_"
    }
]

dataset_rows = [i['rows'] for i in datasets if i is not None]
dataset_names = [i['name'] for i in datasets if i is not None]
dataset_files = [i['file_string'] for i in datasets if i is not None]
datasets = [i['dataset'] for i in datasets if i is not None]

phases_of_interest = ['final', 'ae_pt', 'bc_pt']
weight_key_list = ['model_ae_pre_weights', 'model_cl_pre_weights', 
                   'model_sr_pre_weights', 'model_bd_pre_weights', 
                   'model_bc_pre_weights', 'model_dr_pre_weights']
trained_list = [True]
weight_list = [copy(model.get_weights())]
weight_names = ['final']

trained_list += [result.get(i, None) is not None for i in weight_key_list]
weight_list += [result.get(i, None) for i in weight_key_list]
weight_names += ['ae_pt', 'cl_pt', 'sr_pt', 'bd_pt', 'bc_pt', 'dr_pt']

prediction_model = (model.supervised or model.survival_model or model.drug_response_model)
for trained, weights, weight_name in zip(trained_list, weight_list, weight_names):
    if trained and weight_name in phases_of_interest:
        model.set_weights(weights)
        for dataset, rows, dataset_name, file_prefix in zip(datasets, dataset_rows, dataset_names, dataset_files):
            zi = get_embeddings(model, dataset, tuple_dataset = False)
            zi = pd.DataFrame(
                zi, 
                index = rows,
                columns = ['z{}'.format(i+1) for i in range(zi.shape[1])])
            zi['dataset'] = dataset_name
            zi.to_csv(
                f"{file_prefix}{weight_name}_embeddings.csv.gz", 
                na_rep = 'NA', 
                header = True, 
                index = True)
            if weight_name == 'final' and prediction_model:
                pi = get_predictions(model, dataset, tuple_dataset = False)
                for pred, key in zip(pi.values(), pi.keys()):
                    if pred.shape[1]:
                        cols = [key + '_' + str(i) for i in np.arange(pred.shape[1])]
                    else:
                        cols = key
                    pi[key] = pd.DataFrame(pi[key], index = rows, columns = cols)
                pi = pd.concat(pi.values(), axis = 1)
                pi['dataset'] = dataset_name
                pi.to_csv(
                    f"{file_prefix}{weight_name}_predictions.csv.gz", 
                    na_rep = 'NA', 
                    header = True, 
                    index = True)

model.set_weights(weight_list[0]) # Set back to final weights

#%% external dataset
if external_dataset == 'tcga':
    raise ValueError('Pan-cancer TCGA as external dataset is not implemented.')
elif external_dataset == 'tcga_brca':
    from dataset_collections import get_tcga_brca
    external_data_dict = get_tcga_brca(home = home)
elif external_dataset == 'scanb':
    from dataset_collections import get_scanb
    external_data_dict = get_scanb(home = home)
else:
    raise ValueError('external_dataset does not match a known dataset')

# Filter to common genes
external_original_expression = external_data_dict['patient_exp']

reordering_map = dict(zip(data_dict['gene_ids'], np.arange(data_dict['gene_ids'].shape[0])))
new_external_gex_mat = np.full(
    (external_data_dict['patient_exp'].shape[0], 
    data_dict['gene_ids'].shape[0]), 0.)
for i,g in enumerate(external_data_dict['gene_ids']):
    internal_ind = reordering_map.get(g, None)
    if internal_ind:
        new_external_gex_mat[:,internal_ind] = external_data_dict['patient_exp'][:,i]
external_data_dict['patient_exp'] = new_external_gex_mat# * 0.

#%% PCA plots for applicability domain analysis

if False:
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    pca_model = PCA(n_components = 2)
    internal_std_scaler = StandardScaler()
    internal_std = internal_std_scaler.fit_transform(data_dict['patient_exp'])
    internal_pcs = pca_model.fit_transform(internal_std)
    external_std_scaler = StandardScaler()
    use_same_scaler = True
    if use_same_scaler:
        external_std = internal_std_scaler.transform(external_data_dict['patient_exp'])
    else:
        external_std = external_std_scaler.fit_transform(external_data_dict['patient_exp'])
    external_pcs = pca_model.transform(external_std)
    
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    combined_data = pd.DataFrame(
        np.concatenate((internal_pcs, external_pcs), axis = 0))
    n_internal = internal_pcs.shape[0]
    n_external = external_pcs.shape[0]
    combined_data['dataset'] = np.array(['internal'] * n_internal + ['external'] * n_external)
    
    import re
    model_name = re.sub('/$', '', res_path)
    model_name = re.sub('.*?/', '', model_name)
    scale_string = 'same_scaler' if use_same_scaler else 'different_scaler'
    plot_fn = f"{res_path}../plots/{model_name}_pca_applicability_{scale_string}.png"
    
    fig, ax = plt.subplots(figsize=(6, 6))
    sns.set_style("whitegrid")
    sns.scatterplot(
        combined_data.iloc[np.arange(combined_data.shape[0])[::-1]], 
        x = 0, 
        y = 1, 
        hue = 'dataset', 
        markers = ['+'], 
        style = [0] * (n_internal + n_external))
    plt.grid()    
    plt.tight_layout()
    plt.savefig(plot_fn, dpi = 300)
    
    # Inverse PCA order
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    pca_model = PCA(n_components = 2)
    external_std_scaler = StandardScaler()
    external_std = external_std_scaler.fit_transform(external_data_dict['patient_exp'])
    internal_std_scaler = StandardScaler()
    use_same_scaler = False
    if use_same_scaler:
        internal_std = external_std_scaler.transform(data_dict['patient_exp'])
    else:
        internal_std = internal_std_scaler.fit_transform(data_dict['patient_exp'])
    
    external_pcs = pca_model.fit_transform(external_std)
    internal_pcs = pca_model.transform(internal_std)
    
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    combined_data = pd.DataFrame(
        np.concatenate((internal_pcs, external_pcs), axis = 0))
    n_internal = internal_pcs.shape[0]
    n_external = external_pcs.shape[0]
    combined_data['dataset'] = np.array(['internal'] * n_internal + ['external'] * n_external)
    
    import re
    model_name = re.sub('/$', '', res_path)
    model_name = re.sub('.*?/', '', model_name)
    scale_string = 'same_scaler' if use_same_scaler else 'different_scaler'
    plot_fn = f"{res_path}../plots/{model_name}_pca_applicability_{scale_string}_reverse.png"
    
    fig, ax = plt.subplots(figsize=(6, 6))
    sns.set_style("whitegrid")
    sns.scatterplot(
        combined_data.iloc[np.arange(combined_data.shape[0])[::-1]], 
        x = 0, 
        y = 1, 
        hue = 'dataset', 
        markers = ['+'], 
        style = [0] * (n_internal + n_external))
    plt.grid()    
    plt.tight_layout()
    plt.savefig(plot_fn, dpi = 300)


#%% Serialize external data
external_serialized_data_path = internal_serialized_data_path + '../external/'
external_spec_file = f"{external_serialized_data_path}serialized_data_spec.json"
if False and os.path.exists(external_spec_file):
    handle = open(external_spec_file, 'r')
    external_serialized_data_spec = json.load(handle)
    handle.close()
else:
    os.makedirs(external_serialized_data_path, exist_ok = True)
    model_args = search_kwargs['model_args']
    external_serialized_data_spec = process_and_serialize(
        model_args = model_args, 
        data_dict = external_data_dict, 
        omics_layer = search_kwargs['omics_layer'], 
        data_standardize = search_kwargs['data_standardize'], 
        serialized_data_path = external_serialized_data_path, 
        patient_std_scalers = std_scalers.get('patient', {}))

serialized_external_data = JSONFeatureSpecDecoder(external_serialized_data_spec)

#%% Load serialized data
external_data_instance = serialized_external_data

external_patient_dataset = parse_serialized_dataset_from_file(
    serialized_file = external_data_instance['patient']['filename'],
    feature_spec = external_data_instance['patient']['feature_spec'])
external_patient_rows = np.array(external_data_instance['patient']['sample_info']['rownames'])
external_patient_model_spec = external_data_instance['patient']['model_spec']
external_patient_datasize = len(external_data_instance['patient']['sample_info']['rownames'])

#%% external patient evaluation
external_patient_dataset_batched = dataset_batch_setup(external_patient_dataset, **data_batch_args)

external_res = get_model_losses(
    result['model'], 
    dataset = external_patient_dataset_batched, 
    get_metrics = True, 
    tuple_dataset = False, 
    manual_batch_flag = False)

#%% external patient predictions and embeddings

from sae.evaluation import get_embeddings, get_predictions

z = get_embeddings(
    result['model'], 
    external_patient_dataset_batched, 
    tuple_dataset = False)
z = pd.DataFrame(
    z, 
    index = external_patient_rows,
    columns = ['z{}'.format(i+1) for i in range(z.shape[1])])
z['dataset'] = external_dataset

p = get_predictions(result['model'], external_patient_dataset_batched, tuple_dataset = False)
for pred, key in zip(p.values(), p.keys()):
    if pred.shape[1]:
        cols = [key + '_' + str(i) for i in np.arange(pred.shape[1])]
    else:
        cols = key
    p[key] = pd.DataFrame(p[key], index = external_patient_rows, columns = cols)
p = pd.concat(p.values(), axis = 1)
p['dataset'] = external_dataset

z.to_csv(
    f"{res_path}external_evaluation/external_survival_validation_embeddings.csv.gz", 
    na_rep = 'NA', 
    header = True, 
    index = True)
p.to_csv(
    f"{res_path}external_evaluation/external_survival_validation_predictions.csv.gz", 
    na_rep = 'NA', 
    header = True, 
    index = True)

#%% Save results

dataset_label_map = {
    'tcga' : 'TCGA', 
    'tcga_brca' : 'TCGA_BRCA', 
    'scanb' : 'SCANB' 
}

result_df = pd.DataFrame({
    'dataset' : (
        dataset_label_map[internal_dataset], 
        dataset_label_map[external_dataset]
    ),
    'reconstruction_mse' : (
        internal_res['reconstruction_loss_dataset1'], 
        external_res['reconstruction_loss_dataset1']
    ), 
    'regularization' : (
        internal_res['regularization'], 
        external_res['regularization']
    ), 
    'survival_log_likelihood' : (
        internal_res['survival_log_likelihood'], 
        external_res['survival_log_likelihood']
    ), 
    'survival_concordance' : (
        internal_res['surv_c'], 
        external_res['surv_c']
    ), 
    'confounder_alignment_norm' : (
        internal_res['confounder_alignment_norm'], 
        external_res['confounder_alignment_norm']
    ) 
})
result_df.to_csv(f"{res_path}external_evaluation/external_survival_validation_res.csv")

#%% Bruna PDTC and PDTX gene-expression datasets
from dataset_collections import get_bruna_pdtc, get_bruna_pdtx

pdtc_data_dict = get_bruna_pdtc(home = home)
pdtx_data_dict = get_bruna_pdtx(home = home)

# Filter to common genes
for dd in [pdtc_data_dict, pdtx_data_dict]:
    og_expression = dd['cl_exp']
    
    reordering_map = dict(zip(data_dict['gene_ids'], np.arange(data_dict['gene_ids'].shape[0])))
    new_gex_mat = np.full(
        (dd['cl_exp'].shape[0], 
        data_dict['gene_ids'].shape[0]), 0.)
    for i,g in enumerate(dd['gene_ids']):
        internal_ind = reordering_map.get(g, None)
        if internal_ind:
            new_gex_mat[:,internal_ind] = dd['cl_exp'][:,i]
    dd['cl_exp'] = new_gex_mat# * 0.

#%% Cell-line applicability domain via PCA
if False:
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    pca_model = PCA(n_components = 2)
    ccle_std_scaler = StandardScaler()
    ccle_std = ccle_std_scaler.fit_transform(data_dict['cl_exp'])
    ccle_pcs = pca_model.fit_transform(ccle_std)
    bruna_std_scaler = StandardScaler()
    use_same_scaler = False
    if use_same_scaler:
        bruna_std = ccle_std_scaler.transform(pdtc_data_dict['cl_exp'])
    else:
        bruna_std = bruna_std_scaler.fit_transform(pdtc_data_dict['cl_exp'])
    bruna_pcs = pca_model.transform(bruna_std)
    
    import seaborn as sns
    
    combined_data = pd.DataFrame(
        np.concatenate((ccle_pcs, bruna_pcs), axis = 0))
    n_ccle = ccle_pcs.shape[0]
    n_bruna = bruna_pcs.shape[0]
    combined_data['dataset'] = np.array(['ccle'] * n_ccle + ['bruna'] * n_bruna)
    
    import re
    model_name = re.sub('/$', '', res_path)
    model_name = re.sub('.*?/', '', model_name)
    scale_string = 'same_scaler' if use_same_scaler else 'different_scaler'
    plot_fn = f"{res_path}../plots/{model_name}_pca_applicability_{scale_string}_cl.png"
    
    fig, ax = plt.subplots(figsize=(6, 6))
    sns.set_style("whitegrid")
    sns.scatterplot(
        combined_data, 
        x = 0, 
        y = 1, 
        hue = 'dataset', 
        markers = ['+'], 
        style = [0] * (n_ccle + n_bruna))
    plt.grid()    
    plt.tight_layout()
    plt.savefig(plot_fn, dpi = 300)
    
    # Gao
    from dataset_collections import get_gao_pdx

    gao_data_dict = get_gao_pdx(home = home)

    # Filter to common genes
    og_expression = gao_data_dict['cl_exp']
    
    reordering_map = dict(zip(data_dict['gene_ids'], np.arange(data_dict['gene_ids'].shape[0])))
    new_gex_mat = np.full(
        (gao_data_dict['cl_exp'].shape[0], 
        data_dict['gene_ids'].shape[0]), 0.)
    for i,g in enumerate(gao_data_dict['gene_ids']):
        internal_ind = reordering_map.get(g, None)
        if internal_ind:
            new_gex_mat[:,internal_ind] = gao_data_dict['cl_exp'][:,i]
    gao_data_dict['cl_exp'] = new_gex_mat# * 0.
    
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    pca_model = PCA(n_components = 2)
    ccle_std_scaler = StandardScaler()
    ccle_std = ccle_std_scaler.fit_transform(data_dict['cl_exp'])
    ccle_pcs = pca_model.fit_transform(ccle_std)
    gao_std_scaler = StandardScaler()
    use_same_scaler = False
    if use_same_scaler:
        gao_std = ccle_std_scaler.transform(gao_data_dict['cl_exp'])
    else:
        gao_std = gao_std_scaler.fit_transform(gao_data_dict['cl_exp'])
    gao_pcs = pca_model.transform(gao_std)
    
    import seaborn as sns
    
    combined_data = pd.DataFrame(
        np.concatenate((ccle_pcs, gao_pcs), axis = 0))
    n_ccle = ccle_pcs.shape[0]
    n_gao = gao_pcs.shape[0]
    combined_data['dataset'] = np.array(['ccle'] * n_ccle + ['gao'] * n_gao)
    
    import re
    model_name = re.sub('/$', '', res_path)
    model_name = re.sub('.*?/', '', model_name)
    scale_string = 'same_scaler' if use_same_scaler else 'different_scaler'
    plot_fn = f"{res_path}../plots/{model_name}_pca_applicability_{scale_string}_cl_gao.png"
    
    fig, ax = plt.subplots(figsize=(6, 6))
    sns.set_style("whitegrid")
    sns.scatterplot(
        combined_data, 
        x = 0, 
        y = 1, 
        hue = 'dataset', 
        markers = ['+'], 
        style = [0] * (n_ccle + n_gao))
    plt.grid()    
    plt.tight_layout()
    plt.savefig(plot_fn, dpi = 300)

#%% Serialize external datasets

serialized_ext_data_list = []
for dd, sub_path in zip([pdtc_data_dict, pdtx_data_dict], ['bruna_pdtc/', 'bruna_pdtx/']):
    serialized_ext_data_path = internal_serialized_data_path + '../' + sub_path
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
            serialized_data_path = serialized_ext_data_path, 
            cl_std_scalers = std_scalers.get('cl', {}))
    
    serialized_ext_data_list.append(JSONFeatureSpecDecoder(serialized_ext_data_spec))

#%% Load serialized data
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