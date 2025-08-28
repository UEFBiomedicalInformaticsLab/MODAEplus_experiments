# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import re

def load_results(result_files):
    result_list = []
    for f in result_files:
        result_list.append(pd.read_csv(f, header = 0))
    
    for i in np.arange(len(result_list)):
        result_list[i].rename(columns = {'Unnamed: 0' : 'id'}, inplace = True)
    
    result_dict = {}
    for i in np.arange(len(result_list)):
        task_col = result_list[i].get('task', None)
        if task_col is None:
            taski = re.sub('.*task', '', result_files[i])
            taski = re.sub('\\.csv(\\.gz)', '', taski)
        else:
            taski = task_col[0]
        dfi = result_list[i]
        task_df = result_dict.get(taski, None)
        if task_df is None:
            result_dict[taski] = dfi
        else:
            result_dict[taski] = pd.concat((result_dict[taski], dfi), axis = 0)
    
    return result_dict