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