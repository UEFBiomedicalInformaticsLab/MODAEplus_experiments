import re
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


#%% PCA plots for applicability domain analysis in patient data


pca_model = PCA(n_components = 2)
internal_std_scaler = StandardScaler()
internal_std = internal_std_scaler.fit_transform(data_dict['patient_exp'])
internal_pcs = pca_model.fit_transform(internal_std)
external_std_scaler = StandardScaler()
use_same_scaler = True
if use_same_scaler:
    external_std = internal_std_scaler.transform(external_data_dict['patient_exp'])
else:
    external_std = external_std_scaler.fit_transform(external_data_dict['patient_exp'])
external_pcs = pca_model.transform(external_std)

combined_data = pd.DataFrame(
    np.concatenate((internal_pcs, external_pcs), axis = 0))
n_internal = internal_pcs.shape[0]
n_external = external_pcs.shape[0]
combined_data['dataset'] = np.array(['internal'] * n_internal + ['external'] * n_external)

model_name = re.sub('/$', '', res_path)
model_name = re.sub('.*?/', '', model_name)
scale_string = 'same_scaler' if use_same_scaler else 'different_scaler'
plot_fn = f"{res_path}../plots/{model_name}_pca_applicability_{scale_string}.png"

fig, ax = plt.subplots(figsize=(6, 6))
sns.set_style("whitegrid")
sns.scatterplot(
    combined_data.iloc[np.arange(combined_data.shape[0])[::-1]], 
    x = 0, 
    y = 1, 
    hue = 'dataset', 
    markers = ['+'], 
    style = [0] * (n_internal + n_external))
plt.grid()    
plt.tight_layout()
plt.savefig(plot_fn, dpi = 300)

#%% Inverse PCA order in patient data
pca_model = PCA(n_components = 2)
external_std_scaler = StandardScaler()
external_std = external_std_scaler.fit_transform(external_data_dict['patient_exp'])
internal_std_scaler = StandardScaler()
use_same_scaler = False
if use_same_scaler:
    internal_std = external_std_scaler.transform(data_dict['patient_exp'])
else:
    internal_std = internal_std_scaler.fit_transform(data_dict['patient_exp'])

external_pcs = pca_model.fit_transform(external_std)
internal_pcs = pca_model.transform(internal_std)

combined_data = pd.DataFrame(
    np.concatenate((internal_pcs, external_pcs), axis = 0))
n_internal = internal_pcs.shape[0]
n_external = external_pcs.shape[0]
combined_data['dataset'] = np.array(['internal'] * n_internal + ['external'] * n_external)

model_name = re.sub('/$', '', res_path)
model_name = re.sub('.*?/', '', model_name)
scale_string = 'same_scaler' if use_same_scaler else 'different_scaler'
plot_fn = f"{res_path}../plots/{model_name}_pca_applicability_{scale_string}_reverse.png"

fig, ax = plt.subplots(figsize=(6, 6))
sns.set_style("whitegrid")
sns.scatterplot(
    combined_data.iloc[np.arange(combined_data.shape[0])[::-1]], 
    x = 0, 
    y = 1, 
    hue = 'dataset', 
    markers = ['+'], 
    style = [0] * (n_internal + n_external))
plt.grid()    
plt.tight_layout()
plt.savefig(plot_fn, dpi = 300)

#%% Cell-line applicability domain via PCA

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
pca_model = PCA(n_components = 2)
ccle_std_scaler = StandardScaler()
ccle_std = ccle_std_scaler.fit_transform(data_dict['cl_exp'])
ccle_pcs = pca_model.fit_transform(ccle_std)
bruna_std_scaler = StandardScaler()
use_same_scaler = False
if use_same_scaler:
    bruna_std = ccle_std_scaler.transform(pdtc_data_dict['cl_exp'])
else:
    bruna_std = bruna_std_scaler.fit_transform(pdtc_data_dict['cl_exp'])
bruna_pcs = pca_model.transform(bruna_std)

import seaborn as sns

combined_data = pd.DataFrame(
    np.concatenate((ccle_pcs, bruna_pcs), axis = 0))
n_ccle = ccle_pcs.shape[0]
n_bruna = bruna_pcs.shape[0]
combined_data['dataset'] = np.array(['ccle'] * n_ccle + ['bruna'] * n_bruna)

import re
model_name = re.sub('/$', '', res_path)
model_name = re.sub('.*?/', '', model_name)
scale_string = 'same_scaler' if use_same_scaler else 'different_scaler'
plot_fn = f"{res_path}../plots/{model_name}_pca_applicability_{scale_string}_cl.png"

fig, ax = plt.subplots(figsize=(6, 6))
sns.set_style("whitegrid")
sns.scatterplot(
    combined_data, 
    x = 0, 
    y = 1, 
    hue = 'dataset', 
    markers = ['+'], 
    style = [0] * (n_ccle + n_bruna))
plt.grid()    
plt.tight_layout()
plt.savefig(plot_fn, dpi = 300)

# Gao
from dataset_collections import get_gao_pdx

gao_data_dict = get_gao_pdx(home = home)

# Filter to common genes
og_expression = gao_data_dict['cl_exp']

reordering_map = dict(zip(data_dict['gene_ids'], np.arange(data_dict['gene_ids'].shape[0])))
new_gex_mat = np.full(
    (gao_data_dict['cl_exp'].shape[0], 
    data_dict['gene_ids'].shape[0]), 0.)
for i,g in enumerate(gao_data_dict['gene_ids']):
    internal_ind = reordering_map.get(g, None)
    if internal_ind:
        new_gex_mat[:,internal_ind] = gao_data_dict['cl_exp'][:,i]
gao_data_dict['cl_exp'] = new_gex_mat# * 0.

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
pca_model = PCA(n_components = 2)
ccle_std_scaler = StandardScaler()
ccle_std = ccle_std_scaler.fit_transform(data_dict['cl_exp'])
ccle_pcs = pca_model.fit_transform(ccle_std)
gao_std_scaler = StandardScaler()
use_same_scaler = True
if use_same_scaler:
    gao_std = ccle_std_scaler.transform(gao_data_dict['cl_exp'])
else:
    gao_std = gao_std_scaler.fit_transform(gao_data_dict['cl_exp'])
gao_pcs = pca_model.transform(gao_std)

import seaborn as sns

combined_data = pd.DataFrame(
    np.concatenate((ccle_pcs, gao_pcs), axis = 0))
n_ccle = ccle_pcs.shape[0]
n_gao = gao_pcs.shape[0]
combined_data['dataset'] = np.array(['ccle'] * n_ccle + ['gao'] * n_gao)

import re
model_name = re.sub('/$', '', res_path)
model_name = re.sub('.*?/', '', model_name)
scale_string = 'same_scaler' if use_same_scaler else 'different_scaler'
plot_fn = f"{res_path}../plots/{model_name}_pca_applicability_{scale_string}_cl_gao.png"

fig, ax = plt.subplots(figsize=(6, 6))
sns.set_style("whitegrid")
sns.scatterplot(
    combined_data, 
    x = 0, 
    y = 1, 
    hue = 'dataset', 
    markers = ['+'], 
    style = [0] * (n_ccle + n_gao))
plt.grid()    
plt.tight_layout()
plt.savefig(plot_fn, dpi = 300)