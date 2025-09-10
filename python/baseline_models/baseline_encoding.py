# -*- coding: utf-8 -*-

#%% Baseline reconstruction error
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV, KFold, RepeatedKFold, RepeatedStratifiedKFold
from sklearn.metrics import mean_squared_error
from modae.evaluation import DSC

X = np.concatenate((data_dict['patient_exp'], data_dict['cl_exp']), axis = 0)
group = np.concatenate((
    np.full((data_dict['patient_exp'].shape[0],), 0), 
    np.full((data_dict['cl_exp'].shape[0],), 1)))

pca10 = PCA(n_components = 10, whiten = False, svd_solver = 'arpack')
pca100 = PCA(n_components = 100, whiten = False, svd_solver = 'arpack')

res_mse = []
cv_outer = RepeatedStratifiedKFold(n_splits = 5, n_repeats = 2, random_state = 0)
for train_ind, test_ind in cv_outer.split(X, y = group, groups = group):
    X_train = X[train_ind,:]
    X_test = X[test_ind,:]
    group_train = group[train_ind]
    group_test = group[test_ind]
    
    patient_scaleri = StandardScaler(with_mean = True, with_std = True)
    cl_scaleri = StandardScaler(with_mean = True, with_std = True)
    
    X_patient_train_scaled = patient_scaleri.fit_transform(X_train[group_train == 0,:])
    X_patient_test_scaled = patient_scaleri.transform(X_test[group_test == 0,:])
    X_cl_train_scaled = cl_scaleri.fit_transform(X_train[group_train == 1,:])
    X_cl_test_scaled = cl_scaleri.transform(X_test[group_test == 1,:])
    
    X_train_scaled = np.concatenate((X_patient_train_scaled, X_cl_train_scaled), axis = 0)
    X_test_scaled = np.concatenate((X_patient_test_scaled, X_cl_test_scaled), axis = 0)
    
    group_train = np.concatenate((
        np.full((X_patient_train_scaled.shape[0],), 0), 
        np.full((X_cl_train_scaled.shape[0],), 1)))
    group_test = np.concatenate((
        np.full((X_patient_test_scaled.shape[0],), 0), 
        np.full((X_cl_test_scaled.shape[0],), 1)))
    
    pca10_i = copy(pca10)
    pca100_i = copy(pca100)
    X_train_pca10 = pca10_i.fit_transform(X_train_scaled)
    X_train_pca100 = pca100_i.fit_transform(X_train_scaled)
    X_test_pca10 = pca10_i.transform(X_test_scaled)
    X_test_pca100 = pca100_i.transform(X_test_scaled)
    
    X_train_pca10_recon = pca10_i.inverse_transform(X_train_pca10)
    X_test_pca10_recon = pca10_i.inverse_transform(X_test_pca10)
    X_train_pca100_recon = pca100_i.inverse_transform(X_train_pca100)
    X_test_pca100_recon = pca100_i.inverse_transform(X_test_pca100)
    
    pca10_train_patient_recon_mse = mean_squared_error(X_train_pca10_recon[group_train == 0,:], X_train_scaled[group_train == 0,:], multioutput='uniform_average', squared = True)
    pca10_test_patient_recon_mse = mean_squared_error(X_test_pca10_recon[group_test == 0,:], X_test_scaled[group_test == 0,:], multioutput='uniform_average', squared = True)
    pca100_train_patient_recon_mse = mean_squared_error(X_train_pca100_recon[group_train == 0,:], X_train_scaled[group_train == 0,:], multioutput='uniform_average', squared = True)
    pca100_test_patient_recon_mse = mean_squared_error(X_test_pca100_recon[group_test == 0,:], X_test_scaled[group_test == 0,:], multioutput='uniform_average', squared = True)
    
    pca10_train_cl_recon_mse = mean_squared_error(X_train_pca10_recon[group_train == 1,:], X_train_scaled[group_train == 1,:], multioutput='uniform_average', squared = True)
    pca10_test_cl_recon_mse = mean_squared_error(X_test_pca10_recon[group_test == 1,:], X_test_scaled[group_test == 1,:], multioutput='uniform_average', squared = True)
    pca100_train_cl_recon_mse = mean_squared_error(X_train_pca100_recon[group_train == 1,:], X_train_scaled[group_train == 1,:], multioutput='uniform_average', squared = True)
    pca100_test_cl_recon_mse = mean_squared_error(X_test_pca100_recon[group_test == 1,:], X_test_scaled[group_test == 1,:], multioutput='uniform_average', squared = True)
    
    res_mse.append(pd.DataFrame({
        'fold' : len(res_mse), 
        'pca10_train_patient_recon_mse' : pca10_train_patient_recon_mse, 
        'pca10_test_patient_recon_mse' : pca10_test_patient_recon_mse, 
        'pca100_train_patient_recon_mse' : pca100_train_patient_recon_mse, 
        'pca100_test_patient_recon_mse' : pca100_test_patient_recon_mse, 
        'pca10_train_cl_recon_mse' : pca10_train_cl_recon_mse, 
        'pca10_test_cl_recon_mse' : pca10_test_cl_recon_mse, 
        'pca100_train_cl_recon_mse' : pca100_train_cl_recon_mse, 
        'pca100_test_cl_recon_mse' : pca100_test_cl_recon_mse, 
        'pca10_train_DSC' : DSC(X_train_pca10, group_train).numpy(), 
        'pca10_test_DSC' : DSC(X_test_pca10, group_test).numpy(), 
        'pca100_train_DSC' : DSC(X_train_pca100, group_train).numpy(), 
        'pca100_test_DSC' : DSC(X_test_pca100, group_test).numpy()},
        index = pd.RangeIndex(0,1)))
res_mse_df = pd.concat(res_mse, axis = 0)

res_mse_df.mean()