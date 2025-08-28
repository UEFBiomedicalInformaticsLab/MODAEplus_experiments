library(ggplot2)
source("../setup.R")

fn <- paste0(ccle_path, "CCLE_RRBS_TSS_1kb_20180614.txt")
tss_methylation <- readr::read_tsv(fn)
fn <- paste0(ccle_path, "CCLE_RRBS_TSS_CpG_clusters_20180614.txt")
cpg_methylation <- readr::read_tsv(fn)

cl_cols <- colnames(tss_methylation)[-1:-7]
cl_ids <- gsub("_.*$", "", cl_cols)
cl_tissue <- gsub("^.*?_", "", cl_cols)

cpg_cl_cols <- colnames(cpg_methylation)[-1:-5]

fn <- paste0(ccle_path, "Model_augmented.csv")
model_df <- readr::read_csv(fn)

# Stripped id is not 100% reliable
table(cl_ids %in% model_df[["StrippedCellLineName"]])
cl_ids[!cl_ids %in% model_df[["StrippedCellLineName"]]]

table(
  name = cl_cols %in% model_df[["CCLEName"]], 
  id = cl_ids %in% model_df[["StrippedCellLineName"]]
)

meth_id_df <- tibble::tibble(
  CCLEName = cl_cols, 
  StrippedCellLineName = cl_ids
)
cpg_meth_id_df <- tibble::tibble(
  CCLEName = cpg_cl_cols
)

meth_id_df1 <- plyr::join(
  meth_id_df, 
  model_df[,c("StrippedCellLineName", "ModelID")], 
  type = "left", 
  by = "StrippedCellLineName"
)
meth_id_df2 <- plyr::join(
  meth_id_df, 
  model_df[,c("CCLEName", "ModelID")], 
  type = "left", 
  by = "CCLEName"
)
table(meth_id_df1[["ModelID"]] == meth_id_df2[["ModelID"]], useNA = "always")
meth_id_df1[meth_id_df1[["ModelID"]] != meth_id_df2[["ModelID"]], ]

meth_id_df <- meth_id_df2

cpg_meth_id_df <- plyr::join(
  cpg_meth_id_df, 
  model_df[,c("CCLEName", "ModelID")], 
  type = "left", 
  by = "CCLEName"
)

meth_id_df[["idx"]] <- match(meth_id_df[["ModelID"]], model_df[["ModelID"]])
cpg_meth_id_df[["idx"]] <- match(cpg_meth_id_df[["ModelID"]], model_df[["ModelID"]])

# Check distributions of features
met_dat <- as.matrix(tss_methylation[, cl_cols])
table(duplicated(tss_methylation[["gene"]]))
rownames(met_dat) <- tss_methylation[["TSS_id"]]
table(is.na(met_dat)) #  686017
table(met_dat == 0)   # 4312143
table(met_dat == 1)   #  435258

hist(apply(is.na(met_dat), 1, mean, na.rm = TRUE))
hist(apply(met_dat == 0, 1, mean, na.rm = TRUE))
hist(apply(met_dat == 1, 1, mean, na.rm = TRUE))

cpg_met_dat <- as.matrix(cpg_methylation[, cl_cols])
rownames(cpg_met_dat) <- cpg_methylation[["TSS_id"]]
table(is.na(cpg_met_dat)) # 4214652
table(cpg_met_dat == 0)   # 6409687
table(cpg_met_dat == 1)   # 3918806

hist(apply(is.na(cpg_met_dat), 1, mean))
hist(apply(cpg_met_dat == 0, 1, mean, na.rm = TRUE))
hist(apply(cpg_met_dat == 1, 1, mean, na.rm = TRUE))

# Impute missing values
met_dat_imputed <- impute::impute.knn(met_dat, k = 10, rng.seed = 0)
cpg_met_dat_imputed <- impute::impute.knn(cpg_met_dat, k = 10, rng.seed = 0)

met_dat_imputed <- met_dat_imputed[["data"]]
cpg_met_dat_imputed <- cpg_met_dat_imputed[["data"]]

# Visualize
fn <- paste0(ccle_path, "oncotree_level1.csv")
ccle_oncotree_level1 <- read.csv(fn, row.names = 1, header = TRUE)
meth_sample_ids <- union(meth_id_df[["ModelID"]], cpg_meth_id_df[["ModelID"]])
n_labels <- length(unique(ccle_oncotree_level1[meth_sample_ids, "level_1"]))
set.seed(0)
COPS::umap_viz(
  data = t(met_dat_imputed), 
  category = ccle_oncotree_level1[meth_id_df[["ModelID"]], "level_1"], 
  category_label = "tissue", 
  color_scale = scale_color_manual(values = pals::kovesi.rainbow(n_labels)), 
  umap_args = list(init = "normlaplacian", min_dist = 0.5)
)
set.seed(0)
COPS::umap_viz(
  data = t(met_dat_imputed), 
  category = ccle_oncotree_level1[meth_id_df[["ModelID"]], "primary_level_1"], 
  category_label = "tissue", 
  color_scale = scale_color_manual(values = pals::kovesi.rainbow(n_labels)), 
  umap_args = list(init = "normlaplacian", min_dist = 0.5)
)
set.seed(0)
COPS::umap_viz(
  data = t(cpg_met_dat_imputed), 
  category = ccle_oncotree_level1[cpg_meth_id_df[["ModelID"]], "level_1"], 
  category_label = "tissue", 
  color_scale = scale_color_manual(values = pals::kovesi.rainbow(n_labels)), 
  umap_args = list(init = "normlaplacian", min_dist = 0.5)
)
set.seed(0)
COPS::umap_viz(
  data = t(cpg_met_dat_imputed), 
  category = ccle_oncotree_level1[cpg_meth_id_df[["ModelID"]], "primary_level_1"], 
  category_label = "tissue", 
  color_scale = scale_color_manual(values = pals::kovesi.rainbow(n_labels)), 
  umap_args = list(init = "normlaplacian", min_dist = 0.5)
)





