# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from copy import copy
import os
import json
import importlib.util
import sys

from modae.data_utilities import JSONFeatureSpecDecoder
from modae.data_utilities import parse_serialized_dataset_from_file
from modae.data_utilities import dataset_batch_setup
from modae.evaluation import get_model_losses

script_path = os.environ.get('MODAE_SCRIPT_PATH', default = None)
if script_path is None:
    raise ValueError('Please define MODAE_SCRIPT_PATH')

dsc_spec = importlib.util.spec_from_file_location(
    'dataset_collections', 
    f"{script_path}/python/utilities/cancer_dataset_collections.py"
)
dsc = importlib.util.module_from_spec(dsc_spec)
sys.modules['dataset_collections'] = dsc
dsc_spec.loader.exec_module(dsc)

dsp_spec = importlib.util.spec_from_file_location(
    'dataset_processing', 
    f"{script_path}/python/utilities/cancer_dataset_processing.py"
)
dsp = importlib.util.module_from_spec(dsp_spec)
sys.modules['dataset_processing'] = dsp
dsp_spec.loader.exec_module(dsp)

from dataset_processing import process_and_serialize

def do_external_patient_validation(
        external_dataset, 
        data_dict, 
        external_serialized_data_path, 
        search_kwargs, 
        std_scalers,
        internal_result, 
        data_root
):
    if external_dataset == 'tcga':
        raise ValueError('Pan-cancer TCGA as external dataset is not implemented.')
    elif external_dataset == 'tcga_brca':
        from dataset_collections import get_tcga_brca
        external_data_dict = get_tcga_brca(data_root = data_root)
    elif external_dataset == 'scanb':
        from dataset_collections import get_scanb
        external_data_dict = get_scanb(data_root = data_root)
    else:
        raise ValueError('external_dataset does not match a known dataset')
    
    # Filter to common genes
    #external_original_expression = external_data_dict['patient_exp']
    
    reordering_map = dict(zip(data_dict['gene_ids'], np.arange(data_dict['gene_ids'].shape[0])))
    new_external_gex_mat = np.full(
        (external_data_dict['patient_exp'].shape[0], 
        data_dict['gene_ids'].shape[0]), 0.)
    for i,g in enumerate(external_data_dict['gene_ids']):
        internal_ind = reordering_map.get(g, None)
        if internal_ind:
            new_external_gex_mat[:,internal_ind] = external_data_dict['patient_exp'][:,i]
    external_data_dict['patient_exp'] = new_external_gex_mat# * 0.
    
    # Serialize external data
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
    #external_patient_model_spec = external_data_instance['patient']['model_spec']
    #external_patient_datasize = len(external_data_instance['patient']['sample_info']['rownames'])
    
    #%% external patient evaluation
    data_batch_args = copy(search_kwargs['data_batch_args'])
    data_batch_args['repeat'] = False
    
    external_patient_dataset_batched = dataset_batch_setup(external_patient_dataset, **data_batch_args)
    
    external_res = get_model_losses(
        internal_result['model'], 
        dataset = external_patient_dataset_batched, 
        get_metrics = True, 
        tuple_dataset = False, 
        manual_batch_flag = False)
    
    #%% external patient predictions and embeddings
    
    from modae.model_utilities import get_embeddings, get_predictions
    
    z = get_embeddings(
        internal_result['model'], 
        external_patient_dataset_batched, 
        tuple_dataset = False)
    z = pd.DataFrame(
        z, 
        index = external_patient_rows,
        columns = ['z{}'.format(i+1) for i in range(z.shape[1])])
    z['dataset'] = external_dataset
    
    p = get_predictions(internal_result['model'], external_patient_dataset_batched, tuple_dataset = False)
    for pred, key in zip(p.values(), p.keys()):
        if pred.shape[1]:
            cols = [key + '_' + str(i) for i in np.arange(pred.shape[1])]
        else:
            cols = key
        p[key] = pd.DataFrame(p[key], index = external_patient_rows, columns = cols)
    p = pd.concat(p.values(), axis = 1)
    p['dataset'] = external_dataset
    
    return {'res' : external_res, 'embedding' : z, 'pred' : p}
    
