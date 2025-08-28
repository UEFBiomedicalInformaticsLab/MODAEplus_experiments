# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 13:17:43 2024

@author: rintala
"""
import numpy as np
from sklearn import neighbors

def score_dbscan(
        data, 
        scores, 
        eps=0.5, 
        min_pts=5, 
        score_threshold=0.1, 
        kdt_leaf_size = 15
):
    kdtree = neighbors.KDTree(data, leaf_size = kdt_leaf_size) 
    neighbor_list = kdtree.query_radius(data, eps)
    neighborhood_size = np.array([i.shape[0] for i in neighbor_list])
    neighborhood_std = np.array([np.std(scores[i]) for i in neighbor_list])
    
    core_point = neighborhood_size >= min_pts
    core_point = np.logical_and(core_point, neighborhood_std <= score_threshold)
    
    n_points = data.shape[0]
    clusters = []
    
    for i in np.arange(n_points)[core_point]:
        clusteri = set([i])
        seeds = set(neighbor_list[i]) - clusteri
        while seeds:
            j = seeds.pop()
            clusteri.update([j])
            if core_point[j]:
                seeds.update(neighbor_list[j])
                seeds -= clusteri
        clusters.append(clusteri)
    return clusters

# Example usage
data = np.random.rand(1000, 2)  # 100 points in 2D space
scores = np.random.rand(1000)  # Random scores between 0 and 1
clusters = score_dbscan(data, scores, eps=0.05, min_pts=5, score_threshold=0.25)
print(f"Found {len(clusters)} clusters.")
for idx, cluster in enumerate(clusters):
    print(f"Cluster {idx+1}: {len(cluster)}")

def score_dbscan_optimal_radius(
        data, 
        scores, 
        max_eps=0.5, 
        min_pts=5, 
        score_threshold=0.1, 
        kdt_leaf_size = 15
):
    kdtree = neighbors.KDTree(data, leaf_size = kdt_leaf_size) 
    neighbor_list, neighbor_dist = kdtree.query_radius(
        data, 
        max_eps, 
        return_distance = True, 
        sort_results = True)
    
    neighborhood_size = np.array([i.shape[0] for i in neighbor_list])
    for i in range(len(neighbor_list)):
        n_list = neighbor_list[i]
        n_cum_std = [np.std(scores[n_list[:i]]) for i in np.arange(n_list.shape[0]) + 1]
    neighborhood_std = np.array([np.std(scores[i]) for i in neighbor_list])
    
    core_point = neighborhood_size >= min_pts
    core_point = np.logical_and(core_point, neighborhood_std <= score_threshold)
    
    n_points = data.shape[0]
    clusters = []
    
    for i in np.arange(n_points)[core_point]:
        clusteri = set([i])
        seeds = set(neighbor_list[i]) - clusteri
        while seeds:
            j = seeds.pop()
            clusteri.update([j])
            if core_point[j]:
                seeds.update(neighbor_list[j])
                seeds -= clusteri
        clusters.append(clusteri)
    return clusters

# %%
import numpy as np
import pandas as pd
import glob

import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

#%% Data
from sys import platform
from socket import gethostname
if gethostname() == 'teemu-pc':
    base_path = '/home/teemu/research_work/superAE_HPO/'
elif platform == 'linux':
    base_path = '/research/work/rintala/superAE_HPO/'
else:
    base_path = '//research/workdir/superAE_HPO/'

#res_path = base_path + '20230821_random_search/brca_test_noclfilter/'
#res_path = base_path + '20231027_random_search/brca_test_noclfilter/'
#res_path = base_path + '20231101_random_search/brca_test_noclfilter/'
#res_path = base_path + '20231103_random_search/brca_test_noclfilter/'
#res_path = base_path + 'random_search_231113/brca_test_noclfilter/'
#res_path = base_path + '231117_random_search/brca_test_noclfilter/'
#res_path = base_path + '20231120_random_search/brca_test_noclfilter/'
#res_path = base_path + '20231121_random_search/brca_test_noclfilter/'
#res_path = base_path + '20231122_random_search/scanb_test_noclfilter/'
#res_path = base_path + '20231124_random_search/pancan_test/'
#res_path = base_path + '20240103_random_search/brca_test_noclfilter/'
#res_path = base_path + '20240108_random_search/brca_test_noclfilter/'
#res_path = base_path + '20240115_random_search/brca_test_noclfilter/'
#res_path = base_path + '20240118_random_search/brca_test_noclfilter/'
#res_path = base_path + '20240125_random_search/brca_test_noclfilter/'
#res_path = base_path + '20240125_random_search/brca_test_noclfilter_nopre/'
#res_path = base_path + '20240206_random_search/brca_test_noclfilter/'
res_path = base_path + '20240211_random_search/brca_test_noclfilter/'

res_path = res_path + 'external_evaluation/'

embedding_files = glob.glob(res_path + '*validation_embeddings.csv.gz')

embeddings_list = []
for f in embedding_files:
    embeddings_list.append(pd.read_csv(f, header = 0))

for i in np.arange(len(embeddings_list)):
    embeddings_list[i].rename(columns = {'Unnamed: 0' : 'id'}, inplace = True)

# Fix dataset column
import re
for i in np.arange(len(embeddings_list)):
    if re.search('internal_survival', embedding_files[i]) is not None:
        embeddings_list[i]['dataset'] = 'tcga'
    if re.search('external_survival', embedding_files[i]) is not None:
        embeddings_list[i]['dataset'] = 'scanb'

embeddings_df = pd.concat(embeddings_list, axis = 0)
embeddings_df.value_counts('dataset')

prediction_files = glob.glob(res_path + '*validation_predictions.csv.gz')

prediction_list = []
for f in prediction_files:
    prediction_list.append(pd.read_csv(f, header = 0))

for i in np.arange(len(prediction_list)):
    prediction_list[i].rename(columns = {'Unnamed: 0' : 'id'}, inplace = True)

# Fix dataset column
import re
for i in np.arange(len(prediction_list)):
    if re.search('internal_survival', prediction_files[i]) is not None:
        prediction_list[i]['dataset'] = 'tcga'
    if re.search('external_survival', prediction_files[i]) is not None:
        prediction_list[i]['dataset'] = 'scanb'

predictions_df = pd.concat(prediction_list, axis = 0)
predictions_df.value_counts('dataset')
