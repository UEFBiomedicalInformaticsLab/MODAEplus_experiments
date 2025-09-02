#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 12 10:24:42 2025

@author: teemu
"""

import gzip
import os
import random

import numpy as np
import pandas as pd
import torch

#import data_config

def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.device_count() > 0:
        torch.cuda.manual_seed_all(seed)

data_root = os.environ.get('MODAE_DATA_PATH', default = None)
if data_root is None:
    raise ValueError('Please define MODAE_DATA_PATH')

codeae_data_folder = f"{data_root}CODE-AE-v1.0/data/"
raw_data_folder = os.path.join(codeae_data_folder, 'raw_dat')
ccle_folder = os.path.join(raw_data_folder, 'CCLE')
ccle_sample_file = os.path.join(ccle_folder, 'sample_info.csv')
xena_folder = os.path.join(raw_data_folder, 'Xena')
xena_sample_file = os.path.join(xena_folder, 'TCGA_phenotype_denseDataOnlyDownload.tsv.gz')

preprocessed_data_folder = os.path.join(codeae_data_folder, 'preprocessed_dat')
gex_feature_file = os.path.join(preprocessed_data_folder, 'uq1000_feature.csv')
gex_features_df = pd.read_csv(gex_feature_file, index_col = 0)
seed = 2020
set_seed(seed)

set_seed(seed)
ccle_sample_info_df = pd.read_csv(ccle_sample_file, index_col=0)
with gzip.open(xena_sample_file) as f:
    xena_sample_info_df = pd.read_csv(f, sep='\t', index_col=0)
xena_samples = xena_sample_info_df.index.intersection(gex_features_df.index)
ccle_samples = gex_features_df.index.difference(xena_samples)
#xena_sample_info_df = xena_sample_info_df.loc[xena_samples]
#ccle_sample_info_df = ccle_sample_info_df.loc[ccle_samples.intersection(ccle_sample_info_df.index)]

#xena_df = gex_features_df.loc[xena_samples]
#ccle_df = gex_features_df.loc[ccle_samples]

np.unique(np.in1d(xena_samples, gex_features_df.index), return_counts = True)
np.unique(np.in1d(ccle_samples, gex_features_df.index), return_counts = True)


print(xena_samples.shape[0] + ccle_samples.shape[0])
print(gex_features_df.shape[0])

np.unique(gex_features_df.index == xena_samples.append(ccle_samples), return_counts = True)

fnw = os.path.join(preprocessed_data_folder, 'xena_samples.txt')
with open(fnw, mode = 'w') as f:
    f.writelines([f"{i}\n" for i in xena_samples])

fnw = os.path.join(preprocessed_data_folder, 'ccle_samples.txt')
with open(fnw, mode = 'w') as f:
    f.writelines([f"{i}\n" for i in ccle_samples])