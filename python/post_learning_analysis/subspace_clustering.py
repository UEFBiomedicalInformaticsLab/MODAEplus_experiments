#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 09:57:47 2023

@author: rintala
"""

import numpy as np
import pandas as pd
import os
import glob
import re

import seaborn as sns

#from sklearn import neighbors
from copy import copy
#from sklearn.metrics import silhouette_score
#from sklearn.metrics.pairwise import pairwise_distances

import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

nthreads = 1

#%% Data
from sys import platform
from socket import gethostname
if gethostname() == 'teemu-pc':
    base_path = '/home/teemu/research_work/'
elif platform == 'linux':
    base_path = '/research/work/rintala/'
else:
    base_path = '//research/workdir/'

#res_path = base_path + 'superAE_HPO/20230821_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231027_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231101_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231103_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/random_search_231113/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/231117_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231120_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231121_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231122_random_search/scanb_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231124_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20231128_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231129_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231130_random_search_nostandard/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231201_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231201_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20231218_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231218_random_search_relu/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231219_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231219_random_search_relu/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231219_random_search_variational/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231219_random_search_variational_relu/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231220_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231220_random_search_no_joint/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231221_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20231221_random_search_no_joint/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240103_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240108_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240115_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240118_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240119_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240123_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240124_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240124_random_search/brca_test_noclfilter_alternative/'
#res_path = base_path + 'superAE_HPO/20240125_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240125_random_search/brca_test_noclfilter_nopre/'
#res_path = base_path + 'superAE_HPO/20240206_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240208_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240211_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240213_random_search/brca_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240213_random_search/brca_test_noclfilter_pretrain/'
#res_path = base_path + 'superAE_HPO/20240215_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20240216_random_search/brca_test_noclfilter_pretrain/'
#res_path = base_path + 'superAE_HPO/20240221_random_search/scanb_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240301_random_search/scanb_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240302_random_search/scanb_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240303_random_search/scanb_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240305_random_search/scanb_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240305_2_random_search/scanb_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240306_random_search/scanb_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240310_random_search/scanb_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240315_random_search/scanb_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240319_random_search/scanb_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240328_random_search/scanb_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240430_random_search/scanb_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240605_random_search/scanb_test_noclfilter/'
res_path = base_path + 'superAE_HPO/20240610_random_search/scanb_test_noclfilter/'
#res_path = base_path + 'superAE_HPO/20240614_random_search/pancan_test/'
#res_path = base_path + 'superAE_HPO/20240729_random_search/pancan_test/'

fn = 'annotated_patient_embeddings_with_predicted_aac_brca_only.csv.gz'
fp = f"{res_path}../plots/scanb_test_noclfilter_drug_responses/{fn}"

data = pd.read_csv(fp, header = 0)

del_cols = [
    'SCANB_id', 
    'subtype', 
    'surv_risk_z_bin', 
    'surv_risk'
]
z_cols = [i for i in data.columns if bool(re.search('z[0-9]+', i))]
drug_cols = np.setdiff1d(data.columns, z_cols)
drug_cols = np.setdiff1d(drug_cols, del_cols)
drug_cols = np.setdiff1d(drug_cols, ['age_group', 'grade', 'surv_risk_z'])

#%%


#%% CLIQUE
from pyclustering.cluster.clique import clique, clique_visualizer

X = data.drop(del_cols + z_cols, axis = 1).to_numpy()

clique_instance = clique(X, amount_intervals = 5, density_threshold = 1)

# start clustering process and obtain results
clique_instance.process()
clusters = clique_instance.get_clusters()  # allocated clusters
noise = clique_instance.get_noise()     # points that are considered as outliers (in this example should be empty)
cells = clique_instance.get_cells()     # CLIQUE blocks that forms grid
print("Amount of clusters:", len(clusters))

u, un = np.unique([len(i) for i in clusters], return_counts = True)


# visualize clustering results
clique_visualizer.show_grid(cells, X)    # show grid that has been formed by the algorithm
clique_visualizer.show_clusters(X, clusters, noise)  # show clustering results

#%% 