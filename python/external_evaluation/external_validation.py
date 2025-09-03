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

dsc_spec = importlib.util.spec_from_file_location(
    'dataset_collections', 
    'datasets/cancer_dataset_collections.py'
)
dsc = importlib.util.module_from_spec(dsc_spec)
sys.modules['dataset_collections'] = dsc
dsc_spec.loader.exec_module(dsc)

dsp_spec = importlib.util.spec_from_file_location(
    'dataset_processing', 
    'datasets/cancer_dataset_processing.py'
)
dsp = importlib.util.module_from_spec(dsp_spec)
sys.modules['dataset_processing'] = dsp
dsp_spec.loader.exec_module(dsp)

et_spec = importlib.util.spec_from_file_location(
    'evaluation_tools', 
    'utilities/evaluation_tools.py'
)
et = importlib.util.module_from_spec(et_spec)
sys.modules['evaluation_tools'] = et
et_spec.loader.exec_module(et)

evp_spec = importlib.util.spec_from_file_location(
    'external_evaluation', 
    'external_evaluation/external_validation_patients.py'
)
evp = importlib.util.module_from_spec(evp_spec)
sys.modules['external_validation_patients'] = evp
evp_spec.loader.exec_module(evp)

import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

pd.set_option('display.max_columns', None)

#%% Initialize tensorflow
import tensorflow as tf

gpu_memory = None
nthreads_interop = 1
nthreads = 10

gpus = tf.config.list_physical_devices('GPU')
if len(gpus) and gpu_memory is not None:
    print(f"Found {len(gpus)} GPUs, using 1 logical GPU with {gpu_memory} MB memory.")
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
output_path = os.environ.get('MODAE_OUTPUT_PATH', default = None)
if output_path is None:
    raise ValueError('Please define MODAE_OUTPUT_PATH')
data_root = os.environ.get('MODAE_DATA_PATH', default = None)
if data_root is None:
    raise ValueError('Please define MODAE_DATA_PATH')

res_path = f"{output_path}20250410_random_search/pancan_test/"
#res_path = f"{output_path}20250410_random_search/pancan_ablation_test/"

save_weights = True
weight_analysis = True
save_original_expression = True

if 'ablation' in res_path:
    #ablation_string = 'no_classifier'
    #ablation_string = 'no_famo'
    #ablation_string = 'no_deconfounding'
    #ablation_string = 'no_deconfounding'
    ablation_string = 'no_private_without_deconfounding' # 'no' is a typo
    best_parameter_file = glob.glob(f"{res_path}*{ablation_string}_parameters_taskNone.json")[0]
    parameter_json = True
    res_path = f"{res_path}{ablation_string}/"
    save_original_expression = False
else:
    # Should select by name
    fn = glob.glob(res_path + '*best_task.txt')
    with open(fn[0], 'r') as f:
        best_task = int(f.readlines()[0])
    try:
        best_parameter_file = glob.glob(f"{res_path}*parameters_task{best_task}.csv")[0]
        parameter_json = False
    except IndexError:
        best_parameter_file = glob.glob(f"{res_path}*parameters_task{best_task}.json")[0]
        parameter_json = True

os.makedirs(f"{res_path}external_evaluation/", exist_ok = True)

'''
Note that these must match the training settings. 
TODO: implement automation for settings
'''
dss_sensitivity = True
include_metas_labels = True
exclude_metas_data = False

if parameter_json:
    from sae.parsing_utilities import modae_args_json_decoder
    import json
    param_handle = open(best_parameter_file)
    search_kwargs = json.load(param_handle)
    param_handle.close()
    search_kwargs = modae_args_json_decoder(search_kwargs)
    _ = search_kwargs.pop('nruns', None)
    _ = search_kwargs.pop('nfolds', None)
    _ = search_kwargs.pop('save', None)
    _ = search_kwargs.pop('file_name_prefix', None)
    _ = search_kwargs.pop('parallel', None)
    _ = search_kwargs.pop('gpu_memory', None)
    _ = search_kwargs.pop('nthreads', None)
    _ = search_kwargs.pop('nthreads_interop', None)
    _ = search_kwargs.pop('survival_evaluation_brier_times', None)
    _ = search_kwargs.pop('ps_validation_sets', None)
    _ = search_kwargs.pop('ps_test_sets', None)
    _ = search_kwargs.pop('serialized_data', None)
    search_kwargs['data_standardize'] = True
else:
    from evaluation_tools import get_kwargs
    search_kwargs = get_kwargs(parameter_file = best_parameter_file)

tissue_classifier = search_kwargs['model_args'].get('classifier', False)
gene_preselection = False
if 'pancan' in res_path:
    if dss_sensitivity:
        from dataset_collections import get_xia_ctrp_data
        data_dict = get_xia_ctrp_data(
            data_root = data_root, 
            tissue_classifier = tissue_classifier, 
            include_metas_labels = include_metas_labels, 
            exclude_metas_data = exclude_metas_data
        )
    else:
        from dataset_collections import get_tcga_pancan_ctrp_ccle_solid
        data_dict = get_tcga_pancan_ctrp_ccle_solid(
            data_root = data_root, 
            gene_preselection = gene_preselection, 
            tissue_classifier = tissue_classifier, 
            include_metas_labels = include_metas_labels
        )
    internal_dataset = 'tcga'
    external_dataset = None
elif 'brca' in res_path:
    from dataset_collections import get_tcga_brca_ctrp_ccle_full
    data_dict = get_tcga_brca_ctrp_ccle_full(
        data_root = data_root, 
        gene_preselection = gene_preselection, 
        include_metas_labels = include_metas_labels
    )
    internal_dataset = 'tcga_brca'
    external_dataset = 'scanb'
elif 'scanb' in res_path:
    from dataset_collections import get_scanb_ctrp_ccle_full
    data_dict = get_scanb_ctrp_ccle_full(
        data_root = data_root, 
        gene_preselection = gene_preselection, 
        include_metas_labels = include_metas_labels
    )
    internal_dataset = 'scanb'
    external_dataset = 'tcga_brca'
elif 'pancan_ccle' in res_path:
    from dataset_collections import get_tcga_pancan_ccle_solid
    data_dict = get_tcga_pancan_ccle_solid(
        data_root = data_root, 
        gene_preselection = gene_preselection, 
        include_metas_labels = include_metas_labels
    )
    internal_dataset = 'tcga'
    external_dataset = None
elif 'brca_ccle' in res_path:
    from dataset_collections import get_tcga_brca_ccle_full
    data_dict = get_tcga_brca_ccle_full(
        data_root = data_root, 
        gene_preselection = gene_preselection, 
        include_metas_labels = include_metas_labels
    )
    internal_dataset = 'tcga_brca'
    external_dataset = 'scanb'
elif 'scanb_ccle' in res_path:
    from dataset_collections import get_scanb_ccle_full
    data_dict = get_scanb_ccle_full(
        data_root = data_root, 
        gene_preselection = gene_preselection, 
        include_metas_labels = include_metas_labels
    )
    internal_dataset = 'scanb'
    external_dataset = 'tcga_brca'
else:
    raise ValueError('res_path format does not match known dataset')

if save_original_expression:
    os.makedirs(f"{res_path}external_evaluation/internal/", exist_ok = True)
    # Save filtered datasets
    patient_exp_df = pd.DataFrame(
        data_dict['patient_exp'], 
        columns = data_dict['gene_ids'], 
        index = data_dict['patient_rows']
    )
    patient_exp_df.to_csv(f"{res_path}external_evaluation/internal/patient_mrna.csv.gz")
    cl_exp_df = pd.DataFrame(
        data_dict['cl_exp'], 
        columns = data_dict['gene_ids'], 
        index = data_dict['cl_exp_rows']
    )
    cl_exp_df.to_csv(f"{res_path}external_evaluation/internal/cl_mrna.csv.gz")
    # Save patient cancer types
    fn = f"{res_path}external_evaluation/internal/patient_types.txt"
    with open(fn, 'w') as f:
        f.write('\n'.join(data_dict['patient_cancer_type'].tolist()))

if tissue_classifier:
    os.makedirs(f"{res_path}external_evaluation/", exist_ok = True)
    class_map_df = pd.DataFrame(data_dict['class_map'].items(), columns = ['name', 'key'])
    class_map_df.to_csv(f"{res_path}external_evaluation/class_map.csv")
    

search_kwargs['data_dict'] = data_dict


survival_time_standardize = False

#%% Serialization
from sae.data_utilities import JSONFeatureSpecDecoder
import json
from dataset_processing import process_and_serialize

internal_serialized_data_path = f"{res_path}external_evaluation/internal/"

spec_file = f"{internal_serialized_data_path}serialized_data_spec.json"

os.makedirs(internal_serialized_data_path, exist_ok = True)
model_args = search_kwargs['model_args']
serialized_data_spec, std_scalers = process_and_serialize(
    model_args = model_args, 
    data_dict = data_dict, 
    omics_layer = search_kwargs['omics_layer'], 
    data_standardize = search_kwargs['data_standardize'], 
    serialized_data_path = internal_serialized_data_path, 
    return_scalers = True, 
    shared_scaler = search_kwargs.get('shared_scaler', False), 
    survival_time_standardize = survival_time_standardize)

serialized_data = JSONFeatureSpecDecoder(serialized_data_spec)

different_std_scalers = False
if different_std_scalers:
    std_scalers = {}

#%% Setup training
model_args = search_kwargs['model_args']

maxint = 2**31
np.random.seed(0)
if model_args.get('encoder_layers', None):
    model_args['encoder_init_seeds'] = np.random.randint(
        maxint, 
        size = len(model_args.get('encoder_layers')), 
        dtype=np.int32
    )
if model_args.get('decoder_layers', None):
    model_args['decoder_init_seeds'] = np.random.randint(
        maxint, 
        size = len(model_args.get('decoder_layers')), 
        dtype=np.int32
    )
model_args['recon_init_seed'] = np.random.randint(
    maxint, 
    size = 1, 
    dtype=np.int32
)[0]
if model_args.get('classifier_layers', None):
    model_args['classifier_init_seeds'] = np.random.randint(
        maxint, 
        size = len(model_args.get('classifier_layers')), 
        dtype=np.int32
    )
model_args['classifier_final_init_seed'] = np.random.randint(
    maxint, 
    size = 1, 
    dtype=np.int32
)[0]
if model_args.get('survival_model_layers', None):
    model_args['survival_model_init_seeds'] = np.random.randint(
        maxint, 
        size = len(model_args.get('survival_model_layers')), 
        dtype=np.int32
    )
model_args['survival_model_final_init_seed'] = np.random.randint(
    maxint, 
    size = 1, 
    dtype=np.int32
)[0]
if model_args.get('batch_adversarial_model_layers', None):
    model_args['batch_adversarial_model_init_seeds'] = np.random.randint(
        maxint, 
        size = len(model_args.get('batch_adversarial_model_layers')), 
        dtype=np.int32
    )
model_args['batch_adversarial_model_pred_seed'] = np.random.randint(
    maxint, 
    size = 1, 
    dtype=np.int32
)[0]
if model_args.get('drug_response_model_layers', None):
    model_args['drug_response_model_init_seeds'] = np.random.randint(
        maxint, 
        size = len(model_args.get('drug_response_model_layers')), 
        dtype=np.int32
    )
if model_args.get('drug_response_model_drugwise_layers', None):
    model_args['drug_response_model_drugwise_init_seeds'] = np.random.randint(
        maxint, 
        size = len(model_args.get('drug_response_model_drugwise_layers')), 
        dtype=np.int32
    )
model_args['drug_response_model_pred_seed'] = np.random.randint(
    maxint, 
    size = 1, 
    dtype=np.int32
)[0]

#%% Load serialized data
from sae.data_utilities import parse_serialized_dataset_from_file

data_instance = serialized_data

if data_instance.get('patient', None) is not None:
    internal_patient_dataset = parse_serialized_dataset_from_file(
        serialized_file = data_instance['patient']['filename'],
        feature_spec = data_instance['patient']['feature_spec'])
    internal_patient_rows = np.array(data_instance['patient']['sample_info']['rownames'])
    internal_patient_model_spec = data_instance['patient']['model_spec']
    internal_patient_datasize = len(data_instance['patient']['sample_info']['rownames'])
else:
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
    ccle_cl_model_spec = {}
    ccle_cl_datasize = None

data_model_spec = {**internal_patient_model_spec, **ccle_cl_model_spec}

#%% Training
from sae.training_utilities import train_aecl_with_pretraining

model_args = copy(search_kwargs['model_args'])
data_batch_args = copy(search_kwargs['data_batch_args'])
train_args = copy(search_kwargs['train_args'])
gym_args = copy(search_kwargs['gym_args'])
# Return losses even if disabled for PS
train_args['return_losses'] = True

if False:
    '''
    For rapid testing of parameter alterations. 
    Or set epochs to zero and load trained weights from file later. 
    '''
    gym_args['max_epochs'] = 0

result = train_aecl_with_pretraining(
    patient_serialized_dataset = internal_patient_dataset, 
    cl_serialized_dataset = ccle_cl_dataset, 
    model_args = {**model_args, **data_model_spec},
    data_batch_args = data_batch_args, 
    train_args = train_args, 
    patient_datasize = internal_patient_datasize,
    cl_datasize = ccle_cl_datasize, 
    **gym_args, 
    return_pre_trained_model_weights = True
)

if False:
    '''
    Manually examine results
    '''
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
    result['model'].loss_index_map
    tf.nn.softmax(result['model'].objective_weights)

data_batch_args['repeat'] = False

#%% Load weights from previous run
'''
Run the following block to test trained models
'''
load_weights = False
if load_weights:
    import gzip
    weight_fn = f"{res_path}external_evaluation/final_weights.json.gz"
    os.path.exists(weight_fn)
    with gzip.open(weight_fn, 'r') as handle:
        final_weights = json.load(handle)
    init_weights = result['model'].get_weights()
    if len(init_weights) == len(final_weights) - 1:
        obj_weights = final_weights.pop()
    #[np.array(i).shape for i in final_weights]
    #[np.array(i).shape for i in init_weights]
    final_weights = [np.array(i) for i in final_weights]
    result['model'].set_weights(final_weights)
    save_weights = False


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
        index = False
    )

#%% Evaluation
from sae.data_utilities import dataset_batch_setup
from sae.evaluation import get_model_losses

internal_patient_dataset_batched = dataset_batch_setup(
    internal_patient_dataset, 
    **data_batch_args
)

internal_res = get_model_losses(
    result['model'], 
    dataset = internal_patient_dataset_batched, 
    get_metrics = True, 
    tuple_dataset = False, 
    manual_batch_flag = False
)

#%% Predictions and embeddings 
import json
import gzip

from sae.model_utilities import get_embeddings, get_predictions, JSONWeightEncoder

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

phases_of_interest = ['final', 'ae_pt', 'bc_pt', 'sr_pt', 'dr_pt']
weight_key_list = ['model_ae_pre_weights', 'model_cl_pre_weights', 
                   'model_sr_pre_weights', 'model_bd_pre_weights', 
                   'model_bc_pre_weights', 'model_dr_pre_weights']
trained_list = [True]
weight_list = [copy(model.get_weights())]
weight_names = ['final']

trained_list += [result.get(i, None) is not None for i in weight_key_list]
weight_list += [result.get(i, None) for i in weight_key_list]
weight_names += ['ae_pt', 'cl_pt', 'sr_pt', 'bd_pt', 'bc_pt', 'dr_pt']

prediction_model = (
    model.classifier or 
    model.survival_model or 
    model.drug_response_model or 
    model.batch_adversarial_model
)
for trained, weights, weight_name in zip(trained_list, weight_list, weight_names):
    if not trained: 
        continue
    if weight_name != 'final': 
        model.set_weights(weights)
    if save_weights:
        weights = model.get_weights()
        weights_encoded = JSONWeightEncoder(weights)
        weight_fn = f"{res_path}external_evaluation/{weight_name}_weights.json.gz"
        handle = gzip.open(weight_fn, 'wt')
        handle.write(json.dumps(weights_encoded, cls = json.JSONEncoder))
        handle.close()
    if weight_name in phases_of_interest:
        dataset_iterator = zip(datasets, dataset_rows, dataset_names, dataset_files)
        for dataset, rows, dataset_name, file_prefix in dataset_iterator:
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
                index = True
            )
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
                    index = True
                )

model.set_weights(weight_list[0]) # Set back to final weights

#%% Weight analysis

if weight_analysis:
    import matplotlib.pyplot as plt
    import seaborn as sns
    # Load weights from previous run
    for trained, weights, weight_name in zip(trained_list, weight_list, weight_names):
        if not trained: 
            continue
        weight_fn = f"{res_path}external_evaluation/{weight_name}_weights.json.gz"
        handle = gzip.open(weight_fn, 'r')
        weight_json = json.load(handle)
        handle.close()
        weights = [np.array(i) for i in weight_json]
        print([i.shape for i in weights])
        
        # Hard coded indices, adjust manually
        if tissue_classifier:
            offset = 28
        else:
            offset = 22
        plt.figure()
        p1 = sns.heatmap(weights[offset][:,:], center = 0.)
        plt.figure()
        p2 = sns.heatmap(weights[offset+1].reshape((-1,1)), center = 0.)
        plt.figure()
        p3 = sns.heatmap(weights[offset+2][:,:], center = 0.)
        plt.figure()
        p4 = sns.heatmap(weights[offset+3].reshape((-1,1)), center = 0.)
        
        weight_iterator = zip(
            (p1,p2,p3,p4), 
            ('drug_layer1_w', 'drug_layer1_b', 'drug_layer2_w', 'drug_layer2_b')
        )
        for pi, pn in weight_iterator:
            p_fn = f"{res_path}external_evaluation/{weight_name}_{pn}.png"
            fig = pi.get_figure()
            fig.savefig(p_fn)

#%% Survival validation with external patient dataset

if external_dataset is not None:
    from external_validation_patients import do_external_patient_validation
    epatient_dict = do_external_patient_validation(
        external_dataset = external_dataset, 
        data_dict = data_dict, 
        external_serialized_data_path = internal_serialized_data_path + '../external/', 
        search_kwargs = search_kwargs, 
        std_scalers = std_scalers,
        internal_result = result, 
        data_root = data_root
    )
    
    epatient_dict['embedding'].to_csv(
        f"{res_path}external_evaluation/external_survival_validation_embeddings.csv.gz", 
        na_rep = 'NA', 
        header = True, 
        index = True)
    epatient_dict['pred'].to_csv(
        f"{res_path}external_evaluation/external_survival_validation_predictions.csv.gz", 
        na_rep = 'NA', 
        header = True, 
        index = True)
    
    external_res = epatient_dict['res']
    
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
from dataset_collections import get_bruna_pdtc, get_bruna_pdtx, get_gao_pdx

pdtc_data_dict = get_bruna_pdtc(data_root = data_root)
pdtx_data_dict = get_bruna_pdtx(data_root = data_root)
gao_data_dict = get_gao_pdx(data_root = data_root)

# Filter to common genes
for dd in [pdtc_data_dict, pdtx_data_dict, gao_data_dict]:
    og_expression = dd['cl_exp']
    
    reordering_map = dict(zip(
        data_dict['gene_ids'], np.arange(data_dict['gene_ids'].shape[0])
    ))
    new_gex_mat = np.full(
        (dd['cl_exp'].shape[0], 
        data_dict['gene_ids'].shape[0]), 0.)
    for i,g in enumerate(dd['gene_ids']):
        internal_ind = reordering_map.get(g, None)
        if internal_ind:
            new_gex_mat[:,internal_ind] = dd['cl_exp'][:,i]
    dd['cl_exp'] = new_gex_mat# * 0.



#%% Serialize external datasets

serialized_ext_data_dict = {}
for dd, name in zip(
        [pdtc_data_dict, pdtx_data_dict, gao_data_dict], 
        ['bruna_pdtc', 'bruna_pdtx', 'gao_pdtx']
):
    serialized_ext_data_path = f"{internal_serialized_data_path}../{name}/"
    ext_spec_file = f"{serialized_ext_data_path}serialized_data_spec.json"
    os.makedirs(serialized_ext_data_path, exist_ok = True)
    model_args = search_kwargs['model_args']
    serialized_ext_data_spec = process_and_serialize(
        model_args = model_args, 
        data_dict = dd, 
        omics_layer = search_kwargs['omics_layer'], 
        data_standardize = search_kwargs['data_standardize'], 
        serialized_data_path = serialized_ext_data_path, 
        cl_std_scalers = std_scalers.get('cl', {})
    )
    
    serialized_ext_data_dict[name] = JSONFeatureSpecDecoder(serialized_ext_data_spec)

#%% Save results
embeddings = []
predictions = []
for name, edi in serialized_ext_data_dict.items():
    ext_dataset = parse_serialized_dataset_from_file(
        serialized_file = edi['cl']['filename'],
        feature_spec = edi['cl']['feature_spec']
    )
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

#%% CTRDB clinical trial data

from sae.data_utilities import dataset_batch_setup
from sae.model_utilities import get_embeddings, get_predictions
from dataset_collections import get_ctrdb_datasets, get_ctrdb

ctrdb_datasets = list(get_ctrdb_datasets(data_root = data_root).keys())

for trial_data_id in ctrdb_datasets:
    trial_data_dict = get_ctrdb(data_root = data_root, dataset = trial_data_id)
    og_expression = trial_data_dict['patient_exp']
    
    reordering_map = dict(zip(
        data_dict['gene_ids'], np.arange(data_dict['gene_ids'].shape[0])
    ))
    new_gex_mat = np.full(
        (
            trial_data_dict['patient_exp'].shape[0], 
            data_dict['gene_ids'].shape[0]
        ),
        0.
    )
    for i,g in enumerate(trial_data_dict['gene_ids']):
        internal_ind = reordering_map.get(g, None)
        if internal_ind:
            new_gex_mat[:,internal_ind] = trial_data_dict['patient_exp'][:,i]
    trial_data_dict['patient_exp'] = new_gex_mat# * 0.
    
    serialized_ext_data_path = f"{res_path}external_evaluation/{trial_data_id}/"
    ext_spec_file = f"{serialized_ext_data_path}serialized_data_spec.json"
    os.makedirs(serialized_ext_data_path, exist_ok = True)
    model_args = search_kwargs['model_args']
    serialized_ext_data_spec_json = process_and_serialize(
        model_args = model_args, 
        data_dict = trial_data_dict, 
        omics_layer = search_kwargs['omics_layer'], 
        data_standardize = search_kwargs['data_standardize'], 
        serialized_data_path = serialized_ext_data_path, 
        cl_std_scalers = std_scalers.get('patient', {})
    )
    
    serialized_ext_data_spec = JSONFeatureSpecDecoder(serialized_ext_data_spec_json)
    
    ext_dataset = parse_serialized_dataset_from_file(
        serialized_file = serialized_ext_data_spec['patient']['filename'],
        feature_spec = serialized_ext_data_spec['patient']['feature_spec'])
    ext_rows = np.array(serialized_ext_data_spec['patient']['sample_info']['rownames'])
    ext_model_spec = serialized_ext_data_spec['patient']['model_spec']
    ext_datasize = len(serialized_ext_data_spec['patient']['sample_info']['rownames'])
    
    ext_dataset_batched = dataset_batch_setup(ext_dataset, **data_batch_args)
    
    z = get_embeddings(result['model'], ext_dataset_batched, tuple_dataset = False)
    z = pd.DataFrame(
        z, 
        index = ext_rows,
        columns = ['z{}'.format(i+1) for i in range(z.shape[1])])
    z['dataset'] = trial_data_id
    
    p = get_predictions(result['model'], ext_dataset_batched, tuple_dataset = False)
    for pred, key in zip(p.values(), p.keys()):
        if pred.shape[1]:
            cols = [key + '_' + str(i) for i in np.arange(pred.shape[1])]
        else:
            cols = key
        p[key] = pd.DataFrame(p[key], index = ext_rows, columns = cols)
    p = pd.concat(p.values(), axis = 1)
    p['dataset'] = trial_data_id
    
    z.to_csv(
        f"{serialized_ext_data_path}embeddings.csv.gz", 
        na_rep = 'NA', 
        header = True, 
        index = True)
    p.to_csv(
        f"{serialized_ext_data_path}predictions.csv.gz", 
        na_rep = 'NA', 
        header = True, 
        index = True)
