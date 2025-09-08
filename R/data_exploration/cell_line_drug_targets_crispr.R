script_root <- Sys.getenv("MODAE_SCRIPT_PATH")
if (script_root == "") {
  stop("Please define the MODAE_SCRIPT_PATH environment variable.")
}
source(paste0(script_root, "R/setup.R"))

ctrp_dr_df <- get_ctrp_dr_df()

fn <- paste0(ccle_path, "CRISPRGeneDependency.csv")
crispr_dependency <- readr::read_csv(fn)
crispr_dep_rows <- crispr_dependency[[1]]
crispr_dep_cols <- gsub(" ([0-9]+)$", "", colnames(crispr_dependency)[-1])
crispr_dependency <- as.matrix(crispr_dependency[,-1])
rownames(crispr_dependency) <- crispr_dep_rows
colnames(crispr_dependency) <- crispr_dep_cols

if (FALSE) {
  crispr_dependency_imputed <- t(impute::impute.knn(t(crispr_dependency), k = 10)$data)
  set.seed(0)
  COPS::umap_viz(
    crispr_dependency_imputed, 
    factor(1), 
    "cat", 
    umap_neighbors = 20, 
    umap_args = list(init = "normlaplacian", min_dist = 0.5)
  )
}

fn <- paste0(ccle_path, "CRISPRGeneEffect.csv")
crispr_effect <- readr::read_csv(fn)
crispr_eff_rows <- crispr_effect[[1]]
crispr_eff_cols <- gsub(" ([0-9]+)$", "", colnames(crispr_effect)[-1])
crispr_effect <- as.matrix(crispr_effect[,-1])
rownames(crispr_effect) <- crispr_eff_rows
colnames(crispr_effect) <- crispr_eff_cols

if (FALSE) {
  crispr_effect_imputed <- t(impute::impute.knn(t(crispr_effect), k = 10)$data)
  set.seed(0)
  COPS::umap_viz(
    crispr_effect_imputed, 
    factor(1), 
    "cat", 
    umap_neighbors = 20, 
    umap_args = list(init = "normlaplacian", min_dist = 0.5)
  )
}


