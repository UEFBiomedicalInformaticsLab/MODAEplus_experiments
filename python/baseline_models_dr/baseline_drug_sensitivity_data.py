#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 15:39:58 2023

@author: teemu
"""

import numpy as np
import pandas as pd

def load_sensitivity_data(
        root_dir, 
        drug_features = [], 
        sensitivity_sources = ['CTRP']
):
    cl_metadata = pd.read_csv(
        f"{root_dir}combined_cl_metadata", 
        sep = '\t', 
        index_col = 0
    )
    cl_rnaseq = pd.read_csv(
        f"{root_dir}combined_rnaseq_data", 
        sep = '\t', 
        index_col = 0
    )
    cl_rnaseq = cl_rnaseq.loc[cl_metadata.index]
    cl_mapping = pd.read_csv(
        f"{root_dir}cl_mapping", 
        sep = '\t', 
        header = None, 
        names = ['cl1', 'cl2']
    )
    cl_response = pd.read_csv(
        f"{root_dir}combined_single_response_agg",
        sep = '\t'
    )
    
    drug_feature_flag = False
    if np.isin('dragon_desc', drug_features):
        drug_dragon_desc = pd.read_csv(
            f"{root_dir}Combined_PubChem_dragon7_descriptors.tsv", 
            sep = '\t', 
            na_values = ['na']
        )
        drug_feature_flag = True
    if np.isin('dragon_ecfp', drug_features):
        drug_dragon_ecfp = pd.read_csv(
            f"{root_dir}Combined_PubChem_dragon7_ECFP.tsv", 
            sep = '\s+', 
            skiprows = 1,
            header = None, 
            index_col = 0
        )
        drug_feature_flag = True
    if np.isin('dragon_pfp', drug_features):
        drug_dragon_pfp = pd.read_csv(
            f"{root_dir}Combined_PubChem_dragon7_PFP.tsv", 
            sep = '\s+', 
            skiprows = 1, 
            header = None, 
            index_col = 0
        )
        drug_feature_flag = True
    if np.isin('morderd', drug_features):
        drug_mordred = pd.read_csv(
            f"{root_dir}extended_combined_mordred_descriptors", 
            sep = '\t', 
            index_col = 0
        )
        test = drug_mordred.astype('str')
        def check_isnumeric(x):
            return pd.to_numeric(x, errors = 'coerce').isna()
        test_num = test.apply(check_isnumeric, axis = 0)
        test_sums = test_num.apply(pd.Series.value_counts)
        test_mixed = test_sums.apply(pd.Series.isna).any()
        mixed_string_ind = test_mixed.index[np.logical_not(test_mixed)]
        for i in mixed_string_ind:
            drug_mordred.loc[:,i] = pd.to_numeric(
                drug_mordred.loc[:,i], 
                errors = 'coerce'
            )
    
    if drug_feature_flag:
        drug_info = pd.read_csv(
            root_dir + 'drug_info', 
            sep = '\t', index_col = 0
        )
    
    # Response data
    X_key = cl_response.loc[
        cl_response['SOURCE'].isin(sensitivity_sources), 
        :
    ]
    
    # Filter based on availability
    X_key_nonmissing_exp = X_key['CELL'].isin(cl_rnaseq.index)
    X_key = X_key.loc[X_key_nonmissing_exp]
    
    if drug_feature_flag:
        X_key['drug_pubchem'] = drug_info.loc[
            X_key['DRUG'], 
            'PUBCHEM'
        ].reset_index(drop = True)
        X_key['drug_pubchem'].isna().value_counts() # > 30% missing
        X_key.dropna(how = 'any', inplace = True)
        X_key['drug_pubchem'] = [f"PubChem.CID.{i}" for i in X_key['drug_pubchem'].astype('int64')]
        if np.isin('dragon_ecfp', drug_features):
            X_key_nonmissing_drug_ecfp = X_key['drug_pubchem'].isin(drug_dragon_ecfp.index)
            X_key_nonmissing_drug_ecfp.value_counts() # 0 missing
            X_key = X_key.loc[X_key_nonmissing_drug_ecfp]
        if np.isin('dragon_pfp', drug_features):
            X_key_nonmissing_drug_pfp = X_key['drug_pubchem'].isin(drug_dragon_pfp.index)
            X_key_nonmissing_drug_pfp.value_counts() # 0 missing
            X_key = X_key.loc[X_key_nonmissing_drug_pfp]
    if np.isin('morderd', drug_features):
        X_key_nonmissing_drug_mordred = X_key['DRUG'].isin(drug_mordred.index)
        X_key_nonmissing_drug_mordred.value_counts() # > 10% missing
        X_key = X_key.loc[X_key_nonmissing_drug_mordred]
    
    # Gather features
    #X_cl_exp = cl_rnaseq.loc[X_key['CELL']]
    data_dict = {
        'X_key' : X_key, 
        'cl_metadata' : cl_metadata, 
        'cl_mapping' : cl_mapping, 
        'cl_rnaseq' : cl_rnaseq, 
        }
    if np.isin('dragon_desc', drug_features):
        data_dict['X_drug_desc'] = drug_dragon_desc.loc[
            X_key['drug_pubchem']
        ].dropna(axis = 1).to_numpy()
    if np.isin('dragon_ecfp', drug_features):
        data_dict['X_drug_ecfp'] = drug_dragon_ecfp.loc[
            X_key['drug_pubchem']
        ].dropna(axis = 1).to_numpy()
    if np.isin('dragon_pfp', drug_features):
        data_dict['X_drug_pfp'] = drug_dragon_pfp.loc[
            X_key['drug_pubchem']
        ].dropna(axis = 1).to_numpy()
    if np.isin('morderd', drug_features):
        data_dict['X_drug_mordred'] = drug_mordred.loc[
            X_key['DRUG']
        ].dropna(axis = 1).to_numpy()
    
    return data_dict

def load_old_sensitivity_data(
        cell_line_expression_root_dir, 
        cell_line_drug_response_root_dir, 
        cell_line_expression_file = 'CCLE_expression.csv', 
        cell_line_expression_transpose = False, 
        cell_line_expression_mean_cut = 1., 
        cell_line_drug_response_file = 'screen_aac.csv.gz', 
        cell_line_drug_response_row_info_file = 'screen_rowinfo.csv.gz', 
        cell_line_drug_response_row_map_file = 'ctrp_ccle_clid_map.csv.gz',
        cell_line_drug_response_target_column = 'aac_recomputed',
        cell_line_filter_mapping_file = 'Model_augmented.csv', 
        cell_line_filter_mapping_column = 'solid', 
        cell_line_filter_inclusion_list = ['solid'], 
        cell_line_filter_exclusion_list = []):
    cl_data = pd.read_csv(
        f"{cell_line_expression_root_dir}{cell_line_expression_file}", 
        header = 0, 
        index_col = 0
    )
    if cell_line_expression_transpose:
        cl_exp_id = cl_data.columns
        cl_gene_id = cl_data.index
    else:
        cl_exp_id = cl_data.index
        cl_gene_id = cl_data.columns
    cl_data_ind = cl_exp_id.to_numpy()
    if len(cell_line_filter_inclusion_list) > 0:
        cl_filter_table = pd.read_csv(
            f"{cell_line_expression_root_dir}{cell_line_filter_mapping_file}", 
            header = 0, 
            index_col = 0
        )
        cl_filter_code = cl_filter_table.loc[
            cl_data_ind, 
            cell_line_filter_mapping_column
        ]
        cl_filter = np.logical_not(cl_filter_code.isna().to_numpy())
        cl_filter[cl_filter] = np.isin(
            cl_filter_code[cl_filter].to_numpy(), 
            cell_line_filter_inclusion_list
        )
        if np.sum(cl_filter) == 0:
            raise ValueError('Cell-line filter returned 0 cell-lines.')
        cl_data_ind = cl_data_ind[cl_filter]
        if cell_line_expression_transpose:
            cl_data = cl_data.loc[:,cl_filter]
        else:
            cl_data = cl_data.loc[cl_filter,:]
    if len(cell_line_filter_exclusion_list) > 0:
        cl_filter_table = pd.read_csv(
            f"{cell_line_expression_root_dir}{cell_line_filter_mapping_file}", 
            header = 0, 
            index_col = 0
        )
        cl_filter_code = cl_filter_table.loc[
            cl_data_ind, 
            cell_line_filter_mapping_column
        ]
        cl_filter = np.logical_not(cl_filter_code.isna().to_numpy())
        cl_filter[cl_filter] = np.logical_not(np.isin(
                cl_filter_code[cl_filter].to_numpy(), 
                cell_line_filter_exclusion_list
        ))
        if np.sum(cl_filter) == 0:
            raise ValueError('Cell-line filter returned 0 cell-lines.')
        cl_data_ind = cl_data_ind[cl_filter]
        if cell_line_expression_transpose:
            cl_data = cl_data.loc[:,cl_filter]
        else:
            cl_data = cl_data.loc[cl_filter,:]
    n_cl = cl_data_ind.shape[0]
    if cell_line_expression_mean_cut is not None:
        cl_gene_mean_filter = (
            cl_data.to_numpy().mean(axis = 0) > 
            cell_line_expression_mean_cut
        )
    else:
        cl_gene_mean_filter = np.full((cl_data.shape[1],), True)
    
    cl_data = cl_data.loc[:,cl_gene_mean_filter]
    
    dr_data = pd.read_csv(
        f"{cell_line_drug_response_root_dir}{cell_line_drug_response_file}", 
        header = 0, 
        index_col = 0
    )
    dr_row_info = pd.read_csv(
        f"{cell_line_drug_response_root_dir}{cell_line_drug_response_row_info_file}", 
        header = 0, 
        index_col = 0
    )
    dr_screen_cl_map = pd.read_csv(
        f"{cell_line_drug_response_root_dir}{cell_line_drug_response_row_map_file}", 
        header = 0, 
        index_col = 0
    )
    
    nan_cl_ind = np.any(dr_row_info.isna().to_numpy(), axis = 1)
    dr_cl_id_mapped = dr_screen_cl_map.loc[
        dr_row_info['sampleid'], 
        'CCLE_model_id'
    ]
    dr_cl_id_mapped_not_nan_ind = np.logical_not(dr_cl_id_mapped.isna().to_numpy())
    no_expr_cl = np.setdiff1d(
        dr_cl_id_mapped.to_numpy()[dr_cl_id_mapped_not_nan_ind], 
        cl_data_ind
    )
    no_expr_cl_ind = np.isin(dr_cl_id_mapped.to_numpy(), no_expr_cl)
    
    # exclude dr with mismatched cl ids
    dr_filter = np.logical_and(
        np.logical_not(nan_cl_ind), 
        dr_cl_id_mapped_not_nan_ind
    )
    # exclude dr without cl exp
    dr_filter = np.logical_and(
        dr_filter, 
        np.logical_not(no_expr_cl_ind)
    ) 
    
    # exclude dr with missing AAC
    dr_missing_aac = dr_data[cell_line_drug_response_target_column].isna().to_numpy()
    dr_filter = np.logical_and(
        dr_filter, 
        np.logical_not(dr_missing_aac)
    ) 
    
    # Make one df
    dr_cl_vec = dr_cl_id_mapped[dr_filter].to_numpy(dtype = 'str')
    dr_treat_vec = dr_row_info.loc[dr_filter, 'treatmentid_fixed'].to_numpy(dtype = 'str')
    dr_vec = dr_data[cell_line_drug_response_target_column].to_numpy()[dr_filter]
    dr_long = pd.DataFrame({
        'CCLE_ID' : dr_cl_vec, 
        'CTRP_TREAT_ID' : dr_treat_vec, 
        'AAC' : dr_vec})
    
    data_dict = {
        'X_key' : dr_long, 
        'cl_rnaseq' : cl_data, 
        }
    
    return data_dict