# -*- coding: utf-8 -*-

from sys import platform
import pandas as pd
import numpy as np

from modae.data_utilities import complete_data_loader

#%% TCGA and CCLE data
'''
Customized data loading parameters. 
Please make sure they match with settings used in the hyper-parameter search
'''

def get_tcga_brca_ctrp_ccle_full(
        data_root = './', 
        standardize_early = False, 
        gene_preselection = False, 
        drug_response_maxscale = False, 
        tissue_classifier = False, 
        include_metas_labels = False):
    patient_expression_root_dir = f"{data_root}tcga/"
    cell_line_expression_root_dir = f"{data_root}ccle/"
    cell_line_drug_response_root_dir = f"{data_root}ctrp/"
    
    if standardize_early:
        # Load all for standardization context
        patient_expression_cancer_list = 'ACC,BLCA,BRCA,CESC,CHOL,COAD,DLBC,ESCA,GBM,HNSC,KICH,KIRC,KIRP,LGG,LIHC,LUAD,LUSC,MESO,OV,PAAD,PCPG,PRAD,READ,SARC,SKCM,STAD,TGCT,THCA,THYM,UCEC,UCS,UVM'
        patient_expression_cancer_list = patient_expression_cancer_list.split(sep = ',')
        patient_expression_cancer_subset = ['BRCA']
    else:
        patient_expression_cancer_list = ['BRCA']
        patient_expression_cancer_subset = []
    
    if gene_preselection:
        patient_expression_filter_file = 'survival_genes_selected.txt'
        cell_line_expression_filter_file = 'top10_univariate_genes_combined.txt'
    else:
        patient_expression_filter_file = None
        cell_line_expression_filter_file = None
    
    if tissue_classifier:
        patient_class_file = 'oncotree_level1.csv'
        cell_line_class_file = 'oncotree_level1.csv'
        patient_class_col = 'level_1'
        if include_metas_labels:
            cell_line_class_col = 'level_1'
        else:
            cell_line_class_col = 'primary_level_1'
    else:
        patient_class_file = None
        cell_line_class_file = None
        patient_class_col = None
        cell_line_class_col = None
    
    data_dict = complete_data_loader(
        patient_expression_root_dir = patient_expression_root_dir,
        patient_expression_cancer_list = patient_expression_cancer_list, 
        patient_expression_cancer_subset = patient_expression_cancer_subset, 
        patient_expression_file = 'mrna.csv.gz',
        patient_expression_sample_file = 'sample_info.csv.gz', 
        patient_expression_filter_file = patient_expression_filter_file, 
        patient_gene_mapping_file = 'TCGA_gene_mapping.csv', 
        patient_expression_log2 = True, 
        patient_expression_mean_cut = 1., 
        patient_expression_standardize_early = standardize_early, 
        patient_class_file = patient_class_file, 
        patient_class_col = patient_class_col, 
        patient_survival_file = 'survival.csv.gz', 
        patient_survival_event_col = 'OS', 
        patient_survival_time_col = 'OS.time', 
        patient_survival_covar_cols = ['age_at_initial_pathologic_diagnosis'],
        patient_survival_covar_onehot = [False], 
        cell_line_expression_root_dir = cell_line_expression_root_dir,
        cell_line_expression_file = 'CCLE_expression.csv',
        cell_line_expression_filter_file = cell_line_expression_filter_file, 
        cell_line_gene_mapping_file = 'CCLE_gene_mapping.csv', 
        cell_line_expression_log2 = False, 
        cell_line_expression_mean_cut = 1., 
        cell_line_expression_standardize_early = standardize_early, 
        gene_harmonization_union = False, 
        cell_line_class_file = cell_line_class_file, 
        cell_line_class_col = cell_line_class_col, 
        cell_line_drug_response_root_dir = cell_line_drug_response_root_dir,
        cell_line_drug_response_file = 'screen_aac.csv.gz',
        cell_line_drug_response_row_info_file = 'screen_rowinfo.csv.gz', 
        cell_line_drug_response_row_map_file = 'ctrp_ccle_clid_map.csv.gz',
        cell_line_drug_response_target_column = 'aac_recomputed',
        cell_line_drug_response_maxscale = drug_response_maxscale, 
        cell_line_filter_mapping_file = 'Model_augmented.csv', 
        cell_line_filter_mapping_column = 'solid', 
        cell_line_filter_inclusion_list = ['solid'])
    
    return data_dict

def get_tcga_pancan_ctrp_ccle_solid(
        data_root = './', 
        gene_preselection = False, # TODO: create survival list for pan-cancer?
        drug_response_maxscale = False, 
        tissue_classifier = False, 
        include_metas_labels = False):
    patient_expression_root_dir = f"{data_root}tcga/"
    cell_line_expression_root_dir = f"{data_root}ccle/"
    cell_line_drug_response_root_dir = f"{data_root}ctrp/"
    
    patient_expression_cancer_list = 'ACC,BLCA,BRCA,CESC,CHOL,COAD,DLBC,ESCA,GBM,HNSC,KICH,KIRC,KIRP,LGG,LIHC,LUAD,LUSC,MESO,OV,PAAD,PCPG,PRAD,READ,SARC,SKCM,STAD,TGCT,THCA,THYM,UCEC,UCS,UVM'
    patient_expression_cancer_list = patient_expression_cancer_list.split(sep = ',')
    
    if tissue_classifier:
        patient_class_file = 'oncotree_level1.csv'
        cell_line_class_file = 'oncotree_level1.csv'
        patient_class_col = 'level_1'
        if include_metas_labels:
            cell_line_class_col = 'level_1'
        else:
            cell_line_class_col = 'primary_level_1'
    else:
        patient_class_file = None
        cell_line_class_file = None
        patient_class_col = None
        cell_line_class_col = None
    data_dict = complete_data_loader(
        patient_expression_root_dir = patient_expression_root_dir,
        patient_expression_cancer_list = patient_expression_cancer_list,
        patient_expression_file = 'mrna.csv.gz',
        patient_expression_sample_file = 'sample_info.csv.gz', 
        patient_gene_mapping_file = 'TCGA_gene_mapping.csv', 
        patient_expression_log2 = True, 
        patient_expression_mean_cut = 1., 
        patient_class_file = patient_class_file, 
        patient_class_col = patient_class_col, 
        patient_survival_file = 'survival.csv.gz', 
        patient_survival_event_col = 'OS', 
        patient_survival_time_col = 'OS.time', 
        patient_survival_covar_cols = ['gender','age_at_initial_pathologic_diagnosis','type'],
        patient_survival_covar_onehot = [True, False, True], 
        cell_line_expression_root_dir = cell_line_expression_root_dir,
        cell_line_expression_file = 'CCLE_expression.csv',
        cell_line_gene_mapping_file = 'CCLE_gene_mapping.csv', 
        cell_line_expression_log2 = False, 
        cell_line_expression_mean_cut = 1., 
        gene_harmonization_union = False, 
        cell_line_class_file = cell_line_class_file, 
        cell_line_class_col = cell_line_class_col, 
        cell_line_drug_response_root_dir = cell_line_drug_response_root_dir,
        cell_line_drug_response_file = 'screen_aac.csv.gz',
        cell_line_drug_response_row_info_file = 'screen_rowinfo.csv.gz', 
        cell_line_drug_response_row_map_file = 'ctrp_ccle_clid_map.csv.gz',
        cell_line_drug_response_target_column = 'aac_recomputed',
        cell_line_drug_response_maxscale = drug_response_maxscale, 
        cell_line_filter_mapping_file = 'Model_augmented.csv', 
        cell_line_filter_mapping_column = 'solid', 
        cell_line_filter_inclusion_list = ['solid'])
    
    return data_dict

def get_tcga_brca(
        data_root = './', 
        standardize_early = False, 
        gene_preselection = False, 
        drug_response_maxscale = False):
    patient_expression_root_dir = f"{data_root}tcga/"
    
    if standardize_early:
        # Load all for standardization context
        patient_expression_cancer_list = 'ACC,BLCA,BRCA,CESC,CHOL,COAD,DLBC,ESCA,GBM,HNSC,KICH,KIRC,KIRP,LGG,LIHC,LUAD,LUSC,MESO,OV,PAAD,PCPG,PRAD,READ,SARC,SKCM,STAD,TGCT,THCA,THYM,UCEC,UCS,UVM'
        patient_expression_cancer_list = patient_expression_cancer_list.split(sep = ',')
        patient_expression_cancer_subset = ['BRCA']
    else:
        patient_expression_cancer_list = ['BRCA']
        patient_expression_cancer_subset = []
    
    if gene_preselection:
        patient_expression_filter_file = 'survival_genes_selected.txt'
    else:
        patient_expression_filter_file = None
    data_dict = complete_data_loader(
        patient_expression_root_dir = patient_expression_root_dir,
        patient_expression_cancer_list = patient_expression_cancer_list, 
        patient_expression_cancer_subset = patient_expression_cancer_subset, 
        patient_expression_file = 'mrna.csv.gz',
        patient_expression_sample_file = 'sample_info.csv.gz', 
        patient_expression_filter_file = patient_expression_filter_file, 
        patient_gene_mapping_file = 'TCGA_gene_mapping.csv', 
        patient_expression_log2 = True, 
        patient_expression_mean_cut = 1., 
        patient_expression_standardize_early = standardize_early, 
        patient_survival_file = 'survival.csv.gz', 
        patient_survival_event_col = 'OS', 
        patient_survival_time_col = 'OS.time', 
        patient_survival_covar_cols = ['age_at_initial_pathologic_diagnosis'],
        patient_survival_covar_onehot = [False])
    
    return data_dict

#%% CCLE instead of CTRP
def get_tcga_brca_ccle_full(
        data_root = './', 
        standardize_early = False, 
        gene_preselection = False, 
        drug_response_maxscale = False, 
        tissue_classifier = False, 
        include_metas_labels = False):
    patient_expression_root_dir = f"{data_root}tcga/"
    cell_line_expression_root_dir = f"{data_root}ccle/"
    cell_line_drug_response_root_dir = f"{data_root}ccle/drug_sensitivity/"
    
    if standardize_early:
        # Load all for standardization context
        patient_expression_cancer_list = 'ACC,BLCA,BRCA,CESC,CHOL,COAD,DLBC,ESCA,GBM,HNSC,KICH,KIRC,KIRP,LGG,LIHC,LUAD,LUSC,MESO,OV,PAAD,PCPG,PRAD,READ,SARC,SKCM,STAD,TGCT,THCA,THYM,UCEC,UCS,UVM'
        patient_expression_cancer_list = patient_expression_cancer_list.split(sep = ',')
        patient_expression_cancer_subset = ['BRCA']
    else:
        patient_expression_cancer_list = ['BRCA']
        patient_expression_cancer_subset = []
    
    if gene_preselection:
        patient_expression_filter_file = 'survival_genes_selected.txt'
        cell_line_expression_filter_file = 'top10_univariate_genes_combined.txt'
    else:
        patient_expression_filter_file = None
        cell_line_expression_filter_file = None
    
    if tissue_classifier:
        patient_class_file = 'oncotree_level1.csv'
        cell_line_class_file = 'oncotree_level1.csv'
        patient_class_col = 'level_1'
        if include_metas_labels:
            cell_line_class_col = 'level_1'
        else:
            cell_line_class_col = 'primary_level_1'
    else:
        patient_class_file = None
        cell_line_class_file = None
        patient_class_col = None
        cell_line_class_col = None
    
    data_dict = complete_data_loader(
        patient_expression_root_dir = patient_expression_root_dir,
        patient_expression_cancer_list = patient_expression_cancer_list, 
        patient_expression_cancer_subset = patient_expression_cancer_subset, 
        patient_expression_file = 'mrna.csv.gz',
        patient_expression_sample_file = 'sample_info.csv.gz', 
        patient_expression_filter_file = patient_expression_filter_file, 
        patient_gene_mapping_file = 'TCGA_gene_mapping.csv', 
        patient_expression_log2 = True, 
        patient_expression_mean_cut = 1., 
        patient_expression_standardize_early = standardize_early, 
        patient_class_file = patient_class_file, 
        patient_class_col = patient_class_col, 
        patient_survival_file = 'survival.csv.gz', 
        patient_survival_event_col = 'OS', 
        patient_survival_time_col = 'OS.time', 
        patient_survival_covar_cols = ['age_at_initial_pathologic_diagnosis'],
        patient_survival_covar_onehot = [False], 
        cell_line_expression_root_dir = cell_line_expression_root_dir,
        cell_line_expression_file = 'CCLE_expression.csv',
        cell_line_expression_filter_file = cell_line_expression_filter_file, 
        cell_line_gene_mapping_file = 'CCLE_gene_mapping.csv', 
        cell_line_expression_log2 = False, 
        cell_line_expression_mean_cut = 1., 
        cell_line_expression_standardize_early = standardize_early, 
        gene_harmonization_union = False, 
        cell_line_class_file = cell_line_class_file, 
        cell_line_class_col = cell_line_class_col, 
        cell_line_drug_response_root_dir = cell_line_drug_response_root_dir,
        cell_line_drug_response_file = 'screen_aac.csv.gz',
        cell_line_drug_response_row_info_file = 'screen_rowinfo.csv.gz', 
        cell_line_drug_response_row_info_treat_col = 'treatmentid', 
        cell_line_drug_response_row_map_file = 'ccle_clid_map.csv.gz',
        cell_line_drug_response_target_column = 'aac_recomputed',
        cell_line_drug_response_maxscale = drug_response_maxscale, 
        cell_line_filter_mapping_file = 'Model_augmented.csv', 
        cell_line_filter_mapping_column = 'solid', 
        cell_line_filter_inclusion_list = ['solid'])
    
    return data_dict

def get_tcga_pancan_ccle_solid(
        data_root = './', 
        gene_preselection = False, # TODO: create survival list for pan-cancer?
        drug_response_maxscale = False, 
        tissue_classifier = False, 
        include_metas_labels = False):
    patient_expression_root_dir = f"{data_root}tcga/"
    cell_line_expression_root_dir = f"{data_root}ccle/"
    cell_line_drug_response_root_dir = f"{data_root}ccle/drug_sensitivity/"
    
    patient_expression_cancer_list = 'ACC,BLCA,BRCA,CESC,CHOL,COAD,DLBC,ESCA,GBM,HNSC,KICH,KIRC,KIRP,LGG,LIHC,LUAD,LUSC,MESO,OV,PAAD,PCPG,PRAD,READ,SARC,SKCM,STAD,TGCT,THCA,THYM,UCEC,UCS,UVM'
    patient_expression_cancer_list = patient_expression_cancer_list.split(sep = ',')
    
    if tissue_classifier:
        patient_class_file = 'oncotree_level1.csv'
        cell_line_class_file = 'oncotree_level1.csv'
        patient_class_col = 'level_1'
        if include_metas_labels:
            cell_line_class_col = 'level_1'
        else:
            cell_line_class_col = 'primary_level_1'
    else:
        patient_class_file = None
        cell_line_class_file = None
        patient_class_col = None
        cell_line_class_col = None
    
    data_dict = complete_data_loader(
        patient_expression_root_dir = patient_expression_root_dir,
        patient_expression_cancer_list = patient_expression_cancer_list,
        patient_expression_file = 'mrna.csv.gz',
        patient_expression_sample_file = 'sample_info.csv.gz', 
        patient_gene_mapping_file = 'TCGA_gene_mapping.csv', 
        patient_expression_log2 = True, 
        patient_expression_mean_cut = 1., 
        patient_class_file = patient_class_file, 
        patient_class_col = patient_class_col, 
        patient_survival_file = 'survival.csv.gz', 
        patient_survival_event_col = 'OS', 
        patient_survival_time_col = 'OS.time', 
        patient_survival_covar_cols = ['gender','age_at_initial_pathologic_diagnosis','type'],
        patient_survival_covar_onehot = [True, False, True], 
        cell_line_expression_root_dir = cell_line_expression_root_dir,
        cell_line_expression_file = 'CCLE_expression.csv',
        cell_line_gene_mapping_file = 'CCLE_gene_mapping.csv', 
        cell_line_expression_log2 = False, 
        cell_line_expression_mean_cut = 1., 
        gene_harmonization_union = False, 
        cell_line_class_file = cell_line_class_file, 
        cell_line_class_col = cell_line_class_col, 
        cell_line_drug_response_root_dir = cell_line_drug_response_root_dir,
        cell_line_drug_response_file = 'screen_aac.csv.gz',
        cell_line_drug_response_row_info_file = 'screen_rowinfo.csv.gz', 
        cell_line_drug_response_row_info_treat_col = 'treatmentid', 
        cell_line_drug_response_row_map_file = 'ccle_clid_map.csv.gz',
        cell_line_drug_response_target_column = 'aac_recomputed',
        cell_line_drug_response_maxscale = drug_response_maxscale, 
        cell_line_filter_mapping_file = 'Model_augmented.csv', 
        cell_line_filter_mapping_column = 'solid', 
        cell_line_filter_inclusion_list = ['solid'])
    
    return data_dict

def get_scanb_ccle_full(
        data_root = './', 
        standardize_early = False, 
        gene_preselection = False, 
        drug_response_maxscale = False, 
        tissue_classifier = False, 
        include_metas_labels = False):
    patient_expression_root_dir = f"{data_root}scanb_preprocessed/"
    cell_line_expression_root_dir = f"{data_root}ccle/"
    cell_line_drug_response_root_dir = f"{data_root}ccle/drug_sensitivity/"
    
    if gene_preselection:
        patient_expression_filter_file = 'survival_genes_selected.txt'
        cell_line_expression_filter_file = 'top10_univariate_genes_combined.txt'
    else:
        patient_expression_filter_file = None
        cell_line_expression_filter_file = None
    
    if tissue_classifier:
        patient_class_file = 'oncotree_level1.csv'
        cell_line_class_file = 'oncotree_level1.csv'
        patient_class_col = 'level_1'
        if include_metas_labels:
            cell_line_class_col = 'level_1'
        else:
            cell_line_class_col = 'primary_level_1'
    else:
        patient_class_file = None
        cell_line_class_file = None
        patient_class_col = None
        cell_line_class_col = None
    
    data_dict = complete_data_loader(
        patient_expression_root_dir = patient_expression_root_dir,
        patient_expression_cancer_list = ['.'], 
        patient_expression_cancer_subset = [], 
        patient_expression_file = 'scanb_tpm_uq.csv.gz',
        patient_expression_filter_file = patient_expression_filter_file, 
        patient_gene_mapping_file = 'scanb_gene_mapping.csv', 
        patient_expression_log2 = True, 
        patient_expression_mean_cut = 1., 
        patient_class_file = patient_class_file, 
        patient_class_col = patient_class_col, 
        patient_expression_standardize_early = standardize_early, 
        patient_survival_file = 'scanb_pheno.csv.gz', 
        patient_survival_event_col = 'OS_event', 
        patient_survival_time_col = 'OS_days', 
        patient_id_col = 'GEX.assay', 
        patient_id_tcga_barcode = False, 
        patient_redaction_col = '', 
        patient_survival_covar_cols = ['Age (5-year range, e.g., 35(31-35), 40(36-40), 45(41-45) etc.)'],
        patient_survival_covar_onehot = [False], 
        cell_line_expression_root_dir = cell_line_expression_root_dir,
        cell_line_expression_file = 'CCLE_expression.csv',
        cell_line_expression_filter_file = cell_line_expression_filter_file, 
        cell_line_gene_mapping_file = 'CCLE_gene_mapping.csv', 
        cell_line_expression_log2 = False, 
        cell_line_expression_mean_cut = 1., 
        cell_line_expression_standardize_early = standardize_early, 
        gene_harmonization_union = False, 
        cell_line_class_file = cell_line_class_file, 
        cell_line_class_col = cell_line_class_col, 
        cell_line_drug_response_root_dir = cell_line_drug_response_root_dir,
        cell_line_drug_response_file = 'screen_aac.csv.gz',
        cell_line_drug_response_row_info_file = 'screen_rowinfo.csv.gz', 
        cell_line_drug_response_row_info_treat_col = 'treatmentid', 
        cell_line_drug_response_row_map_file = 'ccle_clid_map.csv.gz',
        cell_line_drug_response_target_column = 'aac_recomputed',
        cell_line_drug_response_maxscale = drug_response_maxscale, 
        cell_line_filter_mapping_file = 'Model_augmented.csv', 
        cell_line_filter_mapping_column = 'solid', 
        cell_line_filter_inclusion_list = ['solid'])
    
    return data_dict

#%% CCLE
def get_ccle_solid(
        data_root = './', 
        gene_preselection = False, # TODO: create survival list for pan-cancer?
        drug_response_maxscale = False):
    cell_line_expression_root_dir = f"{data_root}ccle/"
    cell_line_drug_response_root_dir = f"{data_root}ccle/drug_sensitivity/"
    
    data_dict = complete_data_loader(
        cell_line_expression_root_dir = cell_line_expression_root_dir,
        cell_line_expression_file = 'CCLE_expression.csv',
        cell_line_gene_mapping_file = 'CCLE_gene_mapping.csv', 
        cell_line_expression_log2 = False, 
        cell_line_expression_mean_cut = 1., 
        gene_harmonization_union = False, 
        cell_line_drug_response_root_dir = cell_line_drug_response_root_dir,
        cell_line_drug_response_file = 'screen_aac.csv.gz',
        cell_line_drug_response_row_info_file = 'screen_rowinfo.csv.gz', 
        cell_line_drug_response_row_info_treat_col = 'treatmentid', 
        cell_line_drug_response_row_map_file = 'ccle_clid_map.csv.gz',
        cell_line_drug_response_target_column = 'aac_recomputed',
        cell_line_drug_response_maxscale = drug_response_maxscale, 
        cell_line_filter_mapping_file = 'Model_augmented.csv', 
        cell_line_filter_mapping_column = 'solid', 
        cell_line_filter_inclusion_list = ['solid'])
    return data_dict

#%% scanb
def get_scanb(data_root = './'):
    scanb_path = f"{data_root}scanb_preprocessed/"
    
    scanb_data_dict = complete_data_loader(
        patient_expression_root_dir = scanb_path,
        patient_expression_cancer_list = ['.'],
        patient_expression_file = 'scanb_tpm_uq.csv.gz',
        patient_gene_mapping_file = 'scanb_gene_mapping.csv', 
        patient_expression_log2 = True, 
        patient_expression_mean_cut = None, 
        patient_survival_file = 'scanb_pheno.csv.gz', 
        patient_survival_event_col = 'OS_event', 
        patient_survival_time_col = 'OS_days', 
        patient_id_col = 'GEX.assay', 
        patient_id_tcga_barcode = False, 
        patient_redaction_col = '', 
        patient_survival_covar_cols = ['Age (5-year range, e.g., 35(31-35), 40(36-40), 45(41-45) etc.)'],
        patient_survival_covar_onehot = [False])
    
    return scanb_data_dict

#%% scanb with ccle ctrp

def get_scanb_ctrp_ccle_full(
        data_root = './', 
        standardize_early = False, 
        gene_preselection = False, 
        drug_response_maxscale = False, 
        tissue_classifier = False, 
        include_metas_labels = False):
    patient_expression_root_dir = f"{data_root}scanb_preprocessed/"
    cell_line_expression_root_dir = f"{data_root}ccle/"
    cell_line_drug_response_root_dir = f"{data_root}ctrp/"
    
    if gene_preselection:
        patient_expression_filter_file = 'survival_genes_selected.txt'
        cell_line_expression_filter_file = 'top10_univariate_genes_combined.txt'
    else:
        patient_expression_filter_file = None
        cell_line_expression_filter_file = None
    
    if tissue_classifier:
        patient_class_file = 'oncotree_level1.csv'
        cell_line_class_file = 'oncotree_level1.csv'
        patient_class_col = 'level_1'
        if include_metas_labels:
            cell_line_class_col = 'level_1'
        else:
            cell_line_class_col = 'primary_level_1'
    else:
        patient_class_file = None
        cell_line_class_file = None
        patient_class_col = None
        cell_line_class_col = None
    data_dict = complete_data_loader(
        patient_expression_root_dir = patient_expression_root_dir,
        patient_expression_cancer_list = ['.'], 
        patient_expression_cancer_subset = [], 
        patient_expression_file = 'scanb_tpm_uq.csv.gz',
        patient_expression_filter_file = patient_expression_filter_file, 
        patient_gene_mapping_file = 'scanb_gene_mapping.csv', 
        patient_expression_log2 = True, 
        patient_expression_mean_cut = 1., 
        patient_expression_standardize_early = standardize_early, 
        patient_class_file = patient_class_file, 
        patient_class_col = patient_class_col, 
        patient_survival_file = 'scanb_pheno.csv.gz', 
        patient_survival_event_col = 'OS_event', 
        patient_survival_time_col = 'OS_days', 
        patient_id_col = 'GEX.assay', 
        patient_id_tcga_barcode = False, 
        patient_redaction_col = '', 
        patient_survival_covar_cols = ['Age_group'],
        patient_survival_covar_onehot = [False], 
        cell_line_expression_root_dir = cell_line_expression_root_dir,
        cell_line_expression_file = 'CCLE_expression.csv',
        cell_line_expression_filter_file = cell_line_expression_filter_file, 
        cell_line_gene_mapping_file = 'CCLE_gene_mapping.csv', 
        cell_line_expression_log2 = False, 
        cell_line_expression_mean_cut = 1., 
        cell_line_expression_standardize_early = standardize_early, 
        gene_harmonization_union = False, 
        cell_line_class_file = cell_line_class_file, 
        cell_line_class_col = cell_line_class_col, 
        cell_line_drug_response_root_dir = cell_line_drug_response_root_dir,
        cell_line_drug_response_file = 'screen_aac.csv.gz',
        cell_line_drug_response_row_info_file = 'screen_rowinfo.csv.gz', 
        cell_line_drug_response_row_map_file = 'ctrp_ccle_clid_map.csv.gz',
        cell_line_drug_response_target_column = 'aac_recomputed',
        cell_line_drug_response_maxscale = drug_response_maxscale, 
        cell_line_filter_mapping_file = 'Model_augmented.csv', 
        cell_line_filter_mapping_column = 'solid', 
        cell_line_filter_inclusion_list = ['solid'])
    
    return data_dict

#%% Bruna pdtc and pdtx
def get_bruna_pdtc(data_root = './'):
    bruna_path = f"{data_root}breast_pdtc_bruna/"
    
    bruna_data_dict = complete_data_loader(
        cell_line_expression_root_dir = bruna_path,
        cell_line_expression_file = 'pdtc_exp_imputed.csv.gz',
        cell_line_gene_mapping_file = 'bruna_gene_mapping.csv', 
        cell_line_expression_transpose = True, 
        cell_line_expression_log2 = False, 
        cell_line_expression_mean_cut = None)
    
    return bruna_data_dict

def get_bruna_pdtx(data_root = './'):
    bruna_path = f"{data_root}breast_pdtc_bruna/"
    
    bruna_data_dict = complete_data_loader(
        cell_line_expression_root_dir = bruna_path,
        cell_line_expression_file = 'pdtx_exp_imputed.csv.gz',
        cell_line_gene_mapping_file = 'bruna_gene_mapping.csv', 
        cell_line_expression_transpose = True, 
        cell_line_expression_log2 = False, 
        cell_line_expression_mean_cut = None)
    
    return bruna_data_dict


#%% GAO PDX
def get_gao_pdx(data_root = './'):
    gao_path = f"{data_root}pan_cancer_pdx_gao/"
    
    gao_data_dict = complete_data_loader(
        cell_line_expression_root_dir = gao_path,
        cell_line_expression_file = 'exp_tpm_uq.csv.gz',
        cell_line_gene_mapping_file = 'gao_gene_mapping.csv', 
        cell_line_expression_transpose = True, 
        cell_line_expression_log2 = True, 
        cell_line_expression_mean_cut = None)
    
    return gao_data_dict

#%% Alternative drug-sensitivity
def get_xia_ctrp_data(
        data_root = './', 
        tissue_classifier = False, 
        include_metas_labels = False, 
        exclude_metas_data = False):
    patient_expression_root_dir = f"{data_root}tcga/"
    cell_line_expression_root_dir = f"{data_root}ccle/"
    cell_line_drug_response_root_dir = f"{data_root}drug_sensitivity_xia/"
    
    patient_expression_cancer_list = 'ACC,BLCA,BRCA,CESC,CHOL,COAD,DLBC,ESCA,GBM,HNSC,KICH,KIRC,KIRP,LGG,LIHC,LUAD,LUSC,MESO,OV,PAAD,PCPG,PRAD,READ,SARC,SKCM,STAD,TGCT,THCA,THYM,UCEC,UCS,UVM'
    patient_expression_cancer_list = patient_expression_cancer_list.split(sep = ',')
    
    if tissue_classifier:
        patient_class_file = 'oncotree_level1.csv'
        cell_line_class_file = 'oncotree_level1.csv'
        patient_class_col = 'level_1'
        if include_metas_labels:
            cell_line_class_col = 'level_1'
        else:
            cell_line_class_col = 'primary_level_1'
    else:
        patient_class_file = None
        cell_line_class_file = None
        patient_class_col = None
        cell_line_class_col = None
    
    if exclude_metas_data:
        cl_filter_file = 'Model_augmented.csv'
        cl_filter_column = 'solid_primary'
        cl_filter_include = ['solid_primary']
    else:
        cl_filter_file = 'Model_augmented.csv'
        cl_filter_column = 'solid'
        cl_filter_include = ['solid']
    
    xia_data_dict = complete_data_loader(
        patient_expression_root_dir = patient_expression_root_dir,
        patient_expression_cancer_list = patient_expression_cancer_list,
        patient_expression_file = 'mrna.csv.gz',
        patient_expression_sample_file = 'sample_info.csv.gz', 
        patient_gene_mapping_file = 'TCGA_gene_mapping.csv', 
        patient_expression_log2 = True, 
        patient_expression_mean_cut = 1., 
        patient_class_file = patient_class_file, 
        patient_class_col = patient_class_col, 
        patient_survival_file = 'survival.csv.gz', 
        patient_survival_event_col = 'OS', 
        patient_survival_time_col = 'OS.time', 
        patient_survival_covar_cols = ['gender','age_at_initial_pathologic_diagnosis','type'],
        patient_survival_covar_onehot = [True, False, True], 
        cell_line_expression_root_dir = cell_line_expression_root_dir,
        cell_line_expression_file = 'CCLE_expression.csv',
        cell_line_gene_mapping_file = 'CCLE_gene_mapping.csv', 
        cell_line_expression_log2 = False, 
        cell_line_expression_mean_cut = 1., 
        gene_harmonization_union = False, 
        cell_line_class_file = cell_line_class_file, 
        cell_line_class_col = cell_line_class_col, 
        cell_line_drug_response_root_dir = cell_line_drug_response_root_dir,
        cell_line_drug_response_file = 'xia_ctrp_screen_dss.csv.gz',
        cell_line_drug_response_row_info_file = 'xia_ctrp_screen_info.csv.gz', 
        cell_line_drug_response_row_info_sample_col = 'CELL', 
        cell_line_drug_response_row_info_treat_col = 'DRUG',
        cell_line_drug_response_row_map_file = 'xia_id_ccle_clid_map.csv.gz',
        cell_line_drug_response_row_map_column = 'depmap_id', 
        cell_line_drug_response_target_column = 'DSS1',
        cell_line_drug_response_maxscale = False, 
        cell_line_filter_mapping_file = cl_filter_file, 
        cell_line_filter_mapping_column = cl_filter_column, 
        cell_line_filter_inclusion_list = cl_filter_include)
    
    return xia_data_dict

#%% CTRDB 
import os, json
def get_ctrdb(
        data_root = './', 
        dataset = None
):    
    ctrdb_path = "f{data_root}ctrdb/"
    
    dataset_log_map_file = f"{ctrdb_path}log2_transform_map.json"
    with open(dataset_log_map_file, 'r') as f:
        dataset_log_map = json.load(f)
    ctrdb_data_dict = complete_data_loader(
        patient_expression_root_dir = ctrdb_path,
        patient_expression_cancer_list = ['.'],
        patient_expression_file = f"{dataset}_gex.csv.gz",
        patient_gene_mapping_file = f"{dataset}_gene_mapping.csv",
        patient_expression_log2 = dataset_log_map.get(dataset), 
        patient_expression_mean_cut = None)
    
    return ctrdb_data_dict

def get_ctrdb_datasets(data_root = './'):
    ctrdb_path = "f{data_root}ctrdb/"
    
    ctrb_datasets_file = f"{ctrdb_path}non_tcga_datasets.json"
    with open(ctrb_datasets_file, 'r') as f:
        ctrb_datasets = json.load(f)
    return ctrb_datasets

#%%

def match_genes(x_in, rows_in, rows_target):
    reordering_map = dict(zip(rows_target, np.arange(rows_target.shape[0])))
    x_out = np.full((x_in.shape[0], rows_target.shape[0]), 0.)
    for i,g in enumerate(rows_in):
        target_ind = reordering_map.get(g, None)
        if target_ind:
            x_out[:,target_ind] = x_in[:,i]
    return x_out

    
def get_cell_line_labels(data_root = './'):
    ccle_path = f"{data_root}ccle/"
    
    ccle_labels = pd.read_csv(f"{ccle_path}oncotree_level1.csv", index_col = 0)
    return ccle_labels







