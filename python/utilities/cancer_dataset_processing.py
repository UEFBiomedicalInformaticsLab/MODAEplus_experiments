# -*- coding: utf-8 -*-

import json
import os
import numpy as np
from sklearn.preprocessing import StandardScaler
from copy import copy

from modae.data_utilities import (
    process_patient_data, 
    process_cl_data, 
    serialize_dataset_to_file, 
    JSONFeatureSpecEncoder
)

def process_and_serialize(
    model_args, 
    data_dict, 
    omics_layer, 
    data_standardize, 
    serialized_data_path, 
    patient_std_scalers = {}, 
    cl_std_scalers = {}, 
    return_scalers = False, 
    survival_time_standardize = False, 
    shared_scaler = False
):
    spec_dict = {}
    pdata_flag = data_dict.get('patient_exp', None) is not None
    cdata_flag = data_dict.get('cl_exp', None) is not None
    
    if pdata_flag:
        patient_train_ind = np.arange(data_dict['patient_exp'].shape[0])
        patient_test_ind = np.array([], dtype = 'int64')
    if cdata_flag:
        cl_train_ind = np.arange(data_dict['cl_exp'].shape[0])
        cl_test_ind = np.array([], dtype = 'int64')
    
    if pdata_flag and cdata_flag and data_standardize and shared_scaler:
        pscaler = StandardScaler().fit(data_dict['patient_exp'][patient_train_ind])
        cscaler = StandardScaler().fit(data_dict['cl_exp'][cl_train_ind])
        sscaler = copy(pscaler)
        sscaler.mean_ = 0.5 * (pscaler.mean_ + cscaler.mean_)
        sscaler.scale_ = np.sqrt(pscaler.scale_ * cscaler.scale_)
        sscaler.var_ = np.sqrt(pscaler.var_ * cscaler.var_)
        patient_std_scalers = {'mrna' : sscaler}
        cl_std_scalers = {'mrna' : sscaler}
    
    if pdata_flag:
        if survival_time_standardize:
            # TODO: implement this inside process_patient_data
            if data_dict.get('survival_time', None) is not None:
                time_scaler = StandardScaler()
                data_dict['survival_time'] = time_scaler.fit_transform(np.expand_dims(data_dict['survival_time'], axis = 1))[:,0]
        
        patient_data, patient_model_spec = process_patient_data(
            x = {'mrna' : data_dict['patient_exp']}, 
            cv_train_ind = patient_train_ind, 
            cv_test_ind = patient_test_ind, 
            classifier = model_args.get('classifier', False),
            #survival_model = model_args.get('survival_model', False),
            survival_model = model_args.get('survival_model', False) and (data_dict.get('survival_event', None) is not None),
            class_label = data_dict.get('patient_class', None), 
            class_mask = data_dict.get('patient_class_mask', None), 
            survival_event = data_dict.get('survival_event', None),
            survival_time = data_dict.get('survival_time', None), 
            survival_covariates = data_dict.get('survival_covariates', None),
            survival_covariates_categorical = data_dict.get('survival_covariates_categorical', None),
            survival_mask = data_dict.get('survival_mask', None),
            data_standardize = data_standardize, 
            std_scalers = patient_std_scalers)
        
        patient_data['x_train'] = patient_data['x_train'][omics_layer]
        patient_data['x_test'] = patient_data['x_test'][omics_layer]
        patient_data['b_train'] = np.full((patient_data['x_train'].shape[0],), 0, dtype = 'int64')
        patient_data['b_test'] = np.full((patient_data['x_test'].shape[0],), 0, dtype = 'int64')
        patient_model_spec['input_dim'] = patient_model_spec['input_dim'][omics_layer]
        patient_info = {}
        patient_info['rownames_train'] = data_dict['patient_rows'][patient_train_ind]
        patient_info['rownames_test'] = data_dict['patient_rows'][patient_test_ind]
        
        fn = os.path.join(serialized_data_path, 'patient_data.tfrecords')
        patient_feature_spec = serialize_dataset_to_file(
            filename = fn,
            expression = patient_data['x_train'], 
            batch = patient_data['b_train'], 
            survival_time = patient_data.get('survival_time_train', None),
            survival_event = patient_data.get('survival_event_train', None),
            survival_covariates = patient_data.get('survival_covariates_train', None),
            survival_mask = patient_data.get('survival_mask_train', None), 
            class_label = patient_data.get('class_label_train', None), 
            class_mask = patient_data.get('class_mask_train', None))
        patient_spec = {
            'filename' : fn,
            'feature_spec' : patient_feature_spec,
            'model_spec' : patient_model_spec,
            'sample_info' : {'rownames' : patient_info['rownames_train']}}
        spec_dict['patient'] = patient_spec
    if cdata_flag:
        cl_data, cl_model_spec = process_cl_data(
            x = {'mrna' : data_dict['cl_exp']}, 
            cv_train_ind = cl_train_ind, 
            cv_test_ind = cl_test_ind, 
            classifier = model_args.get('classifier', False),
            class_label = data_dict.get('cl_class', None), 
            class_mask = data_dict.get('cl_class_mask', None), 
            drug_response_model = model_args.get('drug_response_model', False) and (data_dict.get('dr_table', None) is not None),
            drug_response = data_dict.get('dr_table', None), 
            drug_response_mask = data_dict.get('dr_table_mask', None), 
            data_standardize = data_standardize, 
            std_scalers = cl_std_scalers)
        
        cl_data['x_train'] = cl_data['x_train'][omics_layer]
        cl_data['x_test'] = cl_data['x_test'][omics_layer]
        cl_data['b_train'] = np.full((cl_data['x_train'].shape[0],), 1, dtype = 'int64')
        cl_data['b_test'] = np.full((cl_data['x_test'].shape[0],), 1, dtype = 'int64')
        cl_model_spec['input_dim'] = cl_model_spec['input_dim'][omics_layer]
        cl_info = {}
        cl_info['rownames_train'] = data_dict['cl_exp_rows'][cl_train_ind]
        cl_info['rownames_test'] = data_dict['cl_exp_rows'][cl_test_ind]
        
        fn = f"{serialized_data_path}cl_model_data.tfrecords"
        cl_model_feature_spec = serialize_dataset_to_file(
            filename = fn,
            expression = cl_data['x_train'], 
            batch = cl_data['b_train'], 
            class_label = cl_data.get('class_label_train', None), 
            class_mask = cl_data.get('class_mask_train', None),
            drug_response = cl_data.get('drug_response_train', None), 
            drug_response_mask = cl_data.get('drug_response_mask_train', None))
        cl_spec = {
            'filename' : fn,
            'feature_spec' : cl_model_feature_spec,
            'model_spec' : cl_model_spec, 
            'sample_info' : {'rownames' : cl_info['rownames_train']}}
        spec_dict['cl'] = cl_spec
    if model_args.get('classifier', False): 
        if pdata_flag:
            p_u = patient_model_spec.pop('unique_classes', np.array([]))
        if cdata_flag:
            cl_u = cl_model_spec.pop('unique_classes', np.array([]))
        if pdata_flag and cdata_flag:
            unique_classes = np.union1d(p_u, cl_u)
            patient_spec['model_spec']['class_number'] = unique_classes.shape[0]
            cl_spec['model_spec']['class_number'] = unique_classes.shape[0]
        elif pdata_flag:
            patient_spec['model_spec']['class_number']  = p_u.shape[0]
        elif cdata_flag:
            cl_spec['model_spec']['class_number'] = cl_u.shape[0]
    serialized_data_spec = JSONFeatureSpecEncoder(spec_dict)
    
    fnw = os.path.join(serialized_data_path, 'serialized_data_spec.json')
    handle = open(fnw, 'w')
    handle.write(json.dumps(serialized_data_spec, cls = json.JSONEncoder))
    handle.close()
    
    if return_scalers and data_standardize:
        scalers = {}
        if pdata_flag:
            scalers['patient'] = patient_data['x_scalers']
        if cdata_flag:
            scalers['cl'] = cl_data['x_scalers']
        return serialized_data_spec, scalers
    
    if return_scalers:
        return serialized_data_spec, {}
    
    return serialized_data_spec