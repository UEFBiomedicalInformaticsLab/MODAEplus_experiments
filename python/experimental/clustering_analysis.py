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
