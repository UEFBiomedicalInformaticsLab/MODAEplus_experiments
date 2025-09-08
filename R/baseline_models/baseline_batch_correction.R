script_root <- Sys.getenv("MODAE_SCRIPT_PATH")
if (script_root == "") {
  stop("Please define the MODAE_SCRIPT_PATH environment variable.")
}
source(paste0(script_root, "R/setup.R"))
library(ggplot2)

# Get pan-cancer mrna dataset from one of the external validation runs
fn <- paste0(base_dir, "20250202_random_search/pancan_test/external_evaluation/internal/cl_mrna.csv.gz")
cl_exp <- read.csv(fn, row.names = 1, header = TRUE)
fn <- paste0(base_dir, "20250202_random_search/pancan_test/external_evaluation/internal/patient_mrna.csv.gz")
patient_exp <- read.csv(fn, row.names = 1, header = TRUE)
fn <- paste0(base_dir, "20250202_random_search/pancan_test/external_evaluation/internal/patient_types.txt")
patient_types <- readLines(fn)

output_path <- paste0(output_dir, "baseline_results/batch_correction/pan_cancer/")
dir.create(output_path, recursive = TRUE)

#sva::ComBat
# ?batchelor::fastMNN
#batchelor_res <- batchelor::batchCorrect(t(cl_exp), t(patient_exp), PARAM = batchelor::FastMnnParam())
set.seed(0)
batchelor_res <- batchelor::fastMNN(CCLE = t(cl_exp), TCGA = t(patient_exp), k = 20)
plot1 <- COPS::umap_viz(
  SingleCellExperiment::reducedDim(batchelor_res), 
  category = SingleCellExperiment::colData(batchelor_res)[["batch"]], 
  category_label = "dataset", 
  #color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=n_labels)), 
  umap_neighbors = 20, 
  pre_manifold_pca = FALSE, 
  umap_args = list(init = "normlaplacian", min_dist = 0.5)
)
save_figure_safe(
  plot1, 
  png, 
  paste0(output_path, "batchelor_fastMNN.png"), 
  width = plot_width, 
  height = plot_width, 
  res = plot_res, 
  units = plot_units
)

bvar <- rep(c("CCLE", "TCGA"), c(nrow(cl_exp), nrow(patient_exp)))
set.seed(0)
ComBat_seq_res <- sva::ComBat_seq(
  2**cbind(t(cl_exp), t(patient_exp)) - 1, 
  batch = bvar
)
plot2 <- COPS::umap_viz(
  log2(t(ComBat_seq_res)+1), 
  bvar, 
  category_label = "dataset", 
  #color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=n_labels)), 
  umap_neighbors = 20, 
  pre_manifold_pca = 50, 
  umap_args = list(init = "normlaplacian", min_dist = 0.5)
)
save_figure_safe(
  plot2, 
  png, 
  paste0(output_path, "ComBat_seq_no_covar.png"), 
  width = plot_width, 
  height = plot_width, 
  res = plot_res, 
  units = plot_units
)
write.csv(log2(t(ComBat_seq_res)+1), gzfile(paste0(output_path, "ComBat_seq_no_covar.csv.gz")))

set.seed(0)
ComBat_res <- sva::ComBat(
  cbind(t(cl_exp), t(patient_exp)), 
  batch = bvar
)
plot3 <- COPS::umap_viz(
  t(ComBat_res), 
  bvar, 
  category_label = "dataset", 
  #color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=n_labels)), 
  umap_neighbors = 20, 
  pre_manifold_pca = 50, 
  umap_args = list(init = "normlaplacian", min_dist = 0.5)
)
save_figure_safe(
  plot3, 
  png, 
  paste0(output_path, "ComBat_no_covar.png"), 
  width = plot_width, 
  height = plot_width, 
  res = plot_res, 
  units = plot_units
)
write.csv(t(ComBat_res), gzfile(paste0(output_path, "ComBat_no_covar.csv.gz")))

# Use oncotree to map cell lines to tissue of origin
cell_line_oncotree_mapping_file <- 'Model_oncotree.csv'
fn <- paste0(cell_line_expression_root_dir, cell_line_oncotree_mapping_file)
cell_line_oncotree_mappings <- read.csv(fn, header = TRUE, row.names = 1)

# Other cell-line info
cell_line_info_file <- 'Model_augmented.csv'
fn <- paste0(cell_line_expression_root_dir, cell_line_info_file)
cell_line_info <- read.csv(fn, header = TRUE, row.names = 1)

with(cell_line_info[rownames(cl_exp),], table(Sex, useNA = "always"))
with(cell_line_info[rownames(cl_exp),], table(is.na(Age)))
with(cell_line_info[rownames(cl_exp),], table(PrimaryOrMetastasis))

ccle_primary_ptr <- rownames(cell_line_info)[cell_line_info[["PrimaryOrMetastasis"]] == "Primary"]
ccle_primary_ptr <- intersect(ccle_primary_ptr, rownames(cl_exp))

if (FALSE) {
  # Primary CL only
  set.seed(0)
  tissue_plot <- COPS::umap_viz(
    cl_exp[ccle_primary_ptr,],
    category = cell_line_oncotree_mappings[ccle_primary_ptr, "level_1"], 
    category_label = "tissue", 
    color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=n_labels)), 
    umap_neighbors = 20, 
    pre_manifold_pca = TRUE, 
    max_pcs = 50, 
    umap_args = list(init = "normlaplacian", min_dist = 0.5)
  )
  save_figure_safe(
    tissue_plot, 
    png, 
    paste0(output_path, "CCLE_gex_primary_only_OT_level1.png"), 
    width = plot_width * 1.2, 
    height = plot_width * 0.9, 
    res = plot_res, 
    units = plot_units
  )
  # Primary and metastatic CL
  tissue_label <- cell_line_oncotree_mappings[rownames(cl_exp), "level_1"]
  n_labels <- length(unique(tissue_label))
  primary_label <- ifelse(
    rownames(cl_exp) %in% ccle_primary_ptr, 
    "primary", 
    "metastasis/NA"
  )
  set.seed(0)
  tissue_plot <- tissue_visualizer(
    cl_exp, 
    reference_shape_label = "primary",
    labeled_ind = !is.na(tissue_label), 
    color_var = tissue_label, 
    color_name = "tissue",
    shape_var = primary_label, 
    shape_name = "primary",
    umap_neighbors = 20, 
    pre_manifold_pca = TRUE, 
    max_pcs = 50, 
    umap_args = list(init = "normlaplacian", min_dist = 0.5),
    primary_color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=n_labels)), 
    primary_fill_scale = scale_fill_manual(values = pals::kovesi.rainbow(n=n_labels)), 
    secondary_color_scale = scale_color_manual(values = c(NA, "#000000")), 
    primary_shape_scale = scale_shape_manual(values = c(3,NA)),
    secondary_shape_scale = scale_shape_manual(values = c(NA,21)), 
    point_size = 3, 
    point_alpha = 0.5, 
    annotation_size = 6, 
    annotation_force = 10
  )
  save_figure_safe(
    tissue_plot, 
    png, 
    paste0(output_path, "CCLE_gex_primary_v_metastasis_OT_level1.png"), 
    width = plot_width * 1.2*2, 
    height = plot_width * 0.9*2, 
    res = plot_res, 
    units = plot_units
  )
}

# TCGA sample info
cancer_folders <- dir(patient_expression_root_dir)
cancer_folders <- cancer_folders[-grep("\\.", cancer_folders)]

fns <- paste0(patient_expression_root_dir, cancer_folders, "/sample_info.csv.gz")
patient_info <- lapply(fns, read.csv, header = TRUE, row.names = 1)
patient_info <- Reduce(rbind, patient_info)
patient_info_ind <- match(rownames(patient_exp), patient_info[["colname"]])

#table(is.na(patient_info_ind))
#table(rownames(patient_exp) %in% patient_info[["colname"]])
#table(patient_info[patient_info_ind, "type"], useNA = "always")

patient_cancer_type <- patient_info[patient_info_ind, "disease"]
patient_id <- patient_info[patient_info_ind, "primary"]
skin_sample_type_ind <- patient_info[patient_info_ind, "disease"] == "SKCM"
skin_sample_type <- patient_info[
  patient_info_ind[skin_sample_type_ind], 
  "type"
]

fn <- paste0(ccle_path, "oncotree.json")
if (!file.exists(fn)) {
  response <- httr::GET("http://oncotree.mskcc.org/api/tumorTypes")
  json_str <- httr::content(response, as = "text")
  writeLines(json_str, fn)
} else {
  json_str <- readLines(fn)
}

oncotree <- jsonlite::fromJSON(json_str, flatten = TRUE)
oncotree <- jsonlite::fromJSON("http://oncotree.mskcc.org/api/tumorTypes", 
                               flatten = TRUE)
oncotree_parent_map <- oncotree$parent
names(oncotree_parent_map) <- oncotree$code
oncotree_level_map <- oncotree$level
names(oncotree_level_map) <- oncotree$code

sarc_subtypes <- TCGAbiolinks::TCGAquery_subtype("SARC")
table(sarc_subtypes[["histology"]], useNA = "always")
table(sarc_subtypes[sarc_subtypes[["histology"]] == "leiomyosarcoma - gynecologic", "Tumor Site"], useNA = "always")

blca_subtypes <- TCGAbiolinks::TCGAquery_subtype("BLCA")
table(blca_subtypes[["Histologic subtype"]], useNA = "always")

brca_subtypes <- TCGAbiolinks::TCGAquery_subtype("BRCA")
table(brca_subtypes[["BRCA_Pathology"]], useNA = "always")

cesc_subtypes <- TCGAbiolinks::TCGAquery_subtype("CESC")
table(cesc_subtypes[["CLIN:Dx_merged"]], useNA = "always")
cesc_subtypes_map <- c(
  Adenocarcinoma = "CEAD", 
  Adenosquamous = "CEAS", 
  Squamous = "CESC", 
  default = "CERVIX"
)
coad_subtypes <- TCGAbiolinks::TCGAquery_subtype("COAD")
#read_subtypes <- TCGAbiolinks::TCGAquery_subtype("READ")
#identical(coad_subtypes, read_subtypes) # TRUE
table(coad_subtypes[["histological_type"]], useNA = "always")
#table(coad_subtypes[["cancer"]], useNA = "always")
#table(coad_subtypes[["anatomic_organ_subdivision"]], useNA = "always")

esca_subtypes <- TCGAbiolinks::TCGAquery_subtype("ESCA")
table(esca_subtypes[["Histological Type - Oesophagus"]], useNA = "always")
table(esca_subtypes[["Histological Type"]], useNA = "always")
esca_subtypes_map <- c(
  AC = "ESCA", 
  ESCC = "ESCC", 
  default = "ESCA"
)

gbm_subtypes <- TCGAbiolinks::TCGAquery_subtype("GBM")
table(gbm_subtypes[["Histology"]], useNA = "always")

lgg_subtypes <- TCGAbiolinks::TCGAquery_subtype("LGG")
table(lgg_subtypes[["Histology"]], useNA = "always")
lgg_subtypes_map <- c(
  astrocytoma = "ASTR",
  oligoastrocytoma = "OAST", 
  oligodendroglioma = "ODG", 
  default = "DIFG"
)
hnsc_subtypes <- TCGAbiolinks::TCGAquery_subtype("HNSC")
paad_subtypes <- TCGAbiolinks::TCGAquery_subtype("PAAD")
table(paad_subtypes[["Histological type by RHH"]], useNA = "always")

pcpg_subtypes <- TCGAbiolinks::TCGAquery_subtype("PCPG")
table(pcpg_subtypes[["Pheochromocytoma or Paraganglioma"]], useNA = "always")
pcpg_subtypes_map <- c(
  paraganglioma = "PGNG", 
  pheochromocytoma = "PHC", 
  default = "ADRENAL_GLAND"
)

prad_subtypes <- TCGAbiolinks::TCGAquery_subtype("PRAD")
stad_subtypes <- TCGAbiolinks::TCGAquery_subtype("STAD")
table(stad_subtypes[["Lauren.Class"]], useNA = "always")
stad_subtypes_map <- c(
  Diffuse = "DSTAD", 
  Intestinal = "ISTAD", 
  default = "STAD"
)

ucec_subtypes <- TCGAbiolinks::TCGAquery_subtype("UCEC")
table(ucec_subtypes[["histology"]], useNA = "always")
ucec_subtypes_map <- c(
  Endometrioid = "UEC", 
  Serous = "USC", 
  default = "UCEC"
)

ucs_subtypes <- TCGAbiolinks::TCGAquery_subtype("UCS")
#table(ucs_subtypes[["histologic subtype"]], useNA = "always")
#table(ucs_subtypes[["Histologic classification"]], useNA = "always")

uvm_subtypes <- TCGAbiolinks::TCGAquery_subtype("UVM")

tcga_cancer_type_manual_map <- c(
  ACC = "ACC",
  BLCA = "BLCA", # no more detail available
  BRCA = "BRCA", # ductal / lobular, mixed, other, missing
  CESC = "CERVIX", # manual detail from subtypes
  CHOL = "CHOL", 
  COAD = "COAD", # might be able to split into AD and mucinous (MACR)
  DLBC = "DLBCLNOS", 
  ESCA = "ESCA", # squamous cell suubtype not separated
  GBM = "GB", 
  HNSC = "HNSC", 
  KICH = "CHRCC", 
  KIRC = "CCRCC", 
  KIRP = "PRCC", 
  LGG = "DIFG", # Must manually update based on subtypes
  LIHC = "HCC", # hepatocellular carcinoma
  LUAD = "LUAD", 
  LUSC = "LUSC", 
  MESO = "PLMESO", # More detail in https://doi.org/10.1158/2159-8290.CD-18-0804 Supp. Table 1
  OV = "HGSOC", # high-grade serous ovarian cancer
  PAAD = "PAAD", # Too few other types to matter
  PCPG = "ADRENAL_GLAND", # technically PHC and PGNG match different tissues
  PRAD = "PRAD", # by far the most common prostate cancer type
  READ = "READ", # might be able to split into AD and mucinous (MACR)
  SARC = "SOFT_TISSUE", # manual detail from subtypes
  SKCM = "SKCM", 
  STAD = "STAD", # diffuse and intestinal from subtypes
  TGCT = "TESTIS", # https://doi.org/10.1016/j.celrep.2018.05.039 Supp. Table 2
  THCA = "THPA", # papillary thyroid carcinoma
  THYM = "TET", # More detail in https://doi.org/10.1016/j.ccell.2018.01.003 Supp. Table 1
  UCEC = "UCEC", # endometrioid, serous from subtypes, NOTE: UCEC contains UCS
  UCS = "UCS", 
  UVM = "UM"
)

patient_cancer_type_fixed <- tcga_cancer_type_manual_map[patient_cancer_type]

tcga_cancer_type_oncotree <- data.frame(cancer_type = unique(patient_cancer_type))
# Initial type code
codesi <- unique(patient_cancer_type_fixed)
for (i in rev(sort(unique(oncotree$level)))) {
  levelsi <- oncotree_level_map[codesi]
  if (any(is.na(levelsi))) {
    missingsi <- codesi[!codesi %in% names(oncotree_level_map)]
    stop(paste("Could not find code level for", paste(missingsi, collapse = ", ")))
  }
  tcga_cancer_type_oncotree[[paste0("level_", i)]] <- ifelse(levelsi == i, codesi, NA)
  parentsi <- oncotree_parent_map[codesi]
  codesi <- ifelse(levelsi == i, parentsi, codesi)
}

rownames(tcga_cancer_type_oncotree) <- tcga_cancer_type_oncotree[["cancer_type"]]
tcga_cancer_type_oncotree[["cancer_type"]] <- NULL

ccle_level1_labels <- cell_line_oncotree_mappings[rownames(cl_exp), "level_1"]
table(ccle_level1_labels, useNA = "always")
ccle_level1_labels_primary_only <- rep_len(NA, nrow(cl_exp))
names(ccle_level1_labels_primary_only) <- rownames(cl_exp)
ccle_level1_labels_primary_only[ccle_primary_ptr] <- cell_line_oncotree_mappings[ccle_primary_ptr, "level_1"]
table(ccle_level1_labels_primary_only, useNA = "always")
table(ccle_level1_labels == ccle_level1_labels_primary_only, useNA = "always")
tcga_level1_labels <- tcga_cancer_type_oncotree[patient_cancer_type, "level_1"]
table(tcga_level1_labels, useNA = "always")

# Write level1 to file
fnw <- paste0(ccle_path, "oncotree_level1.csv")
ccle_level1_labels_df <- data.frame(
  level_1 = ccle_level1_labels, 
  primary_level_1 = ccle_level1_labels_primary_only)
rownames(ccle_level1_labels_df) <- rownames(cl_exp)
write.csv(ccle_level1_labels_df, fnw)

fnw <- paste0(patient_expression_root_dir, "oncotree_level1.csv")
tcga_level1_labels_df <- data.frame(level_1 = tcga_level1_labels, type = patient_cancer_type)
rownames(tcga_level1_labels_df) <- rownames(patient_exp)
write.csv(tcga_level1_labels_df, fnw)

ccle_labeled_ind <- !is.na(ccle_level1_labels)
tcga_labeled_ind <- !is.na(tcga_level1_labels)

lvar <- factor(c(
  ccle_level1_labels[ccle_labeled_ind],
  tcga_level1_labels[tcga_labeled_ind]
))
label_mod <- model.matrix(~lvar)

# Rerun with covariates
bvar_labeled <- rep(c("CCLE", "TCGA"), c(sum(ccle_labeled_ind), sum(tcga_labeled_ind)))
set.seed(0)
ComBat_seq_res <- sva::ComBat_seq(
  2**cbind(
    t(cl_exp[ccle_labeled_ind,]), 
    t(patient_exp[tcga_labeled_ind,])
  ) - 1, 
  batch = bvar_labeled, 
  covar_mod = label_mod, 
  full_mod = TRUE
)
plot2 <- COPS::umap_viz(
  log2(t(ComBat_seq_res)+1), 
  bvar_labeled, 
  category_label = "dataset", 
  #color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=n_labels)), 
  umap_neighbors = 20, 
  pre_manifold_pca = 50, 
  umap_args = list(init = "normlaplacian", min_dist = 0.5)
)
save_figure_safe(
  plot2, 
  png, 
  paste0(output_path, "ComBat_seq.png"), 
  width = plot_width, 
  height = plot_width, 
  res = plot_res, 
  units = plot_units
)
write.csv(log2(t(ComBat_seq_res)+1), gzfile(paste0(output_path, "ComBat_seq.csv.gz")))

set.seed(0)
ComBat_res <- sva::ComBat(
  cbind(
    t(cl_exp[ccle_labeled_ind,]), 
    t(patient_exp[tcga_labeled_ind,])
  ), 
  batch = bvar_labeled, 
  mod = label_mod
)
plot3 <- COPS::umap_viz(
  t(ComBat_res), 
  bvar_labeled, 
  category_label = "dataset", 
  #color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=n_labels)), 
  umap_neighbors = 20, 
  pre_manifold_pca = 50, 
  umap_args = list(init = "normlaplacian", min_dist = 0.5)
)
save_figure_safe(
  plot3, 
  png, 
  paste0(output_path, "ComBat.png"), 
  width = plot_width, 
  height = plot_width, 
  res = plot_res, 
  units = plot_units
)
write.csv(t(ComBat_res), gzfile(paste0(output_path, "ComBat.csv.gz")))

# Tissue labeled UMAPs
output_path <- paste0(output_dir, "baseline_results/batch_correction/pan_cancer/")
dir.create(output_path, recursive = TRUE)

fns <- list(
  ComBat_seq_no_covar = "ComBat_seq_no_covar.csv.gz", 
  ComBat_no_covar = "ComBat_no_covar.csv.gz", 
  ComBat_seq = "ComBat_seq.csv.gz", 
  ComBat = "ComBat.csv.gz"
)
fns <- lapply(fns, function(x) paste0(output_path, x))

dat_mat_list <- lapply(fns, readr::read_csv, show_col_types = FALSE)
for (i in names(dat_mat_list)) {
  sample_names <- dat_mat_list[[i]][[1]]
  dat_mat_list[[i]] <- data.frame(dat_mat_list[[i]][,-1])
  rownames(dat_mat_list[[i]]) <- sample_names
}

fn <- paste0(ccle_path, "oncotree_level1.csv")
ccle_level1_labels_df <- read.csv(fn, row.names = 1, header = TRUE)
ccle_level1_labels_df[["dataset"]] <- "CCLE"

fn <- paste0(patient_expression_root_dir, "oncotree_level1.csv")
tcga_level1_labels_df <- read.csv(fn, row.names = 1, header = TRUE)
tcga_level1_labels_df[["dataset"]] <- "TCGA"

n_labels <- length(union(
  unique(ccle_level1_labels_df[["level_1"]]), 
  unique(tcga_level1_labels_df[["level_1"]])
))

# Raw data
plot_raw <- TRUE
if (plot_raw) {
  fn <- paste0(base_dir, "20250410_random_search/pancan_test/external_evaluation/internal/cl_mrna.csv.gz")
  cl_exp <- readr::read_csv(fn, show_col_types = FALSE)
  cl_ids <- cl_exp[[1]]
  cl_exp <- data.frame(cl_exp[,-1])
  rownames(cl_exp) <- cl_ids
  fn <- paste0(base_dir, "20250410_random_search/pancan_test/external_evaluation/internal/patient_mrna.csv.gz")
  patient_exp <- readr::read_csv(fn, show_col_types = FALSE)
  patient_ids <- patient_exp[[1]]
  patient_exp <- data.frame(patient_exp[,-1])
  rownames(patient_exp) <- patient_ids
  fn <- paste0(base_dir, "20250410_random_search/pancan_test/external_evaluation/internal/patient_types.txt")
  patient_types <- readLines(fn)
  
  dat_mat_list[["Unadjusted"]] <- rbind(cl_exp, patient_exp)
}

plot_tissue_umap <- TRUE
run_tissue_classifier <- FALSE
metrics_correlation_list <- list()
metrics_euclidean_list <- list()
for (i in names(dat_mat_list)[4:5]) {
  dat_mat <- dat_mat_list[[i]]
  tissue_label <- rep_len(NA, nrow(dat_mat))
  names(tissue_label) <- rownames(dat_mat)
  ccle_ptr <- intersect(rownames(dat_mat), rownames(ccle_level1_labels_df))
  tcga_ptr <- intersect(rownames(dat_mat), rownames(tcga_level1_labels_df))
  tissue_label[ccle_ptr] <- ccle_level1_labels_df[ccle_ptr, "level_1"]
  tissue_label[tcga_ptr] <- tcga_level1_labels_df[tcga_ptr, "level_1"]
  tissue_label <- tolower(gsub("_", " ", tissue_label))
  
  dataset <- rep_len(NA, nrow(dat_mat))
  names(dataset) <- rownames(dat_mat)
  dataset[ccle_ptr] <- ccle_level1_labels_df[ccle_ptr, "dataset"]
  dataset[tcga_ptr] <- tcga_level1_labels_df[tcga_ptr, "dataset"]
  dataset <- factor(dataset, levels = c("TCGA", "CCLE"))
  
  nna_ind <- !is.na(tissue_label)
  if (run_tissue_classifier) {
    metrics_correlation_list[[i]] <- tissue_classifier_evaluation(
      data = dat_mat, 
      dataset_label = dataset, 
      reference_dataset_class = "TCGA" , 
      predicted_dataset_class = "CCLE", 
      tissue_label = tissue_label, 
      labeled_ind = nna_ind, 
      excluded_labels = c("adrenal gland", "testis"),
      scale_datasets_separately = TRUE, 
      knn_method = correlation_knn
    )
  }
  
  if(plot_tissue_umap) {
    set.seed(0)
    tissue_plot <- tissue_visualizer(
      dat_mat, 
      reference_shape_label = "TCGA",
      labeled_ind = nna_ind, 
      color_var = tissue_label, 
      color_name = "tissue",
      shape_var = dataset, 
      shape_name = "dataset",
      umap_neighbors = 20, 
      pre_manifold_pca = FALSE, 
      scale_datasets_separately = TRUE, 
      max_pcs = 50, 
      umap_args = list(init = "normlaplacian", min_dist = 0.5, nn_method = "annoy", metric = "correlation"),
      primary_color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=n_labels)), 
      primary_fill_scale = scale_fill_manual(values = pals::kovesi.rainbow(n=n_labels)), 
      secondary_color_scale = scale_color_manual(values = c(NA, "#000000")), 
      primary_shape_scale = scale_shape_manual(values = c(3,NA)),
      secondary_shape_scale = scale_shape_manual(values = c(NA,21)), 
      point_size = 3, 
      point_alpha = 0.5, 
      annotation_size = 6, 
      annotation_force = 10, 
      knn_accuracy = FALSE,
      knn_predicted_dataset = "CCLE", 
      knn_title_method_name = i
    )
    saveRDS(tissue_plot, file = paste0(output_path, i, "_OT_level1.rds"))
    save_figure_safe(
      with(tissue_plot, tissue_plot), 
      png, 
      paste0(output_path, i, "_OT_level1.png"), 
      width = plot_width * 1.2, 
      height = plot_width * 0.9, 
      res = plot_res, 
      units = plot_units
    )
  }
}

# include only primary cell-lines
cell_line_info_file <- 'Model_augmented.csv'
fn <- paste0(cell_line_expression_root_dir, cell_line_info_file)
cell_line_info <- read.csv(fn, header = TRUE, row.names = 1)
ccle_primary_ptr <- rownames(cell_line_info)[cell_line_info[["PrimaryOrMetastasis"]] == "Primary"]
for (i in names(dat_mat_list)) {
  dat_mat <- dat_mat_list[[i]]
  tissue_label <- rep_len(NA, nrow(dat_mat))
  names(tissue_label) <- rownames(dat_mat)
  ccle_ptr <- intersect(rownames(dat_mat), rownames(ccle_level1_labels_df))
  ccle_ptr <- intersect(ccle_ptr, ccle_primary_ptr)
  tcga_ptr <- intersect(rownames(dat_mat), rownames(tcga_level1_labels_df))
  tissue_label[ccle_ptr] <- ccle_level1_labels_df[ccle_ptr, "level_1"]
  tissue_label[tcga_ptr] <- tcga_level1_labels_df[tcga_ptr, "level_1"]
  dataset <- rep_len(NA, nrow(dat_mat))
  names(dataset) <- rownames(dat_mat)
  dataset[ccle_ptr] <- ccle_level1_labels_df[ccle_ptr, "dataset"]
  dataset[tcga_ptr] <- tcga_level1_labels_df[tcga_ptr, "dataset"]
  dataset <- factor(dataset, levels = c("TCGA", "CCLE"))
  
  set.seed(0)
  tissue_plot <- tissue_visualizer(
    dat_mat, 
    reference_shape_label = "TCGA",
    labeled_ind = nna_ind, 
    color_var = tissue_label, 
    color_name = "tissue",
    shape_var = dataset, 
    shape_name = "dataset",
    umap_neighbors = 20, 
    pre_manifold_pca = TRUE, 
    max_pcs = 50, 
    umap_args = list(init = "normlaplacian", min_dist = 0.5, nn_method = "annoy", metric = "correlation"),
    primary_color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=n_labels)), 
    primary_fill_scale = scale_fill_manual(values = pals::kovesi.rainbow(n=n_labels)), 
    secondary_color_scale = scale_color_manual(values = c(NA, "#000000")), 
    primary_shape_scale = scale_shape_manual(values = c(3,NA)),
    secondary_shape_scale = scale_shape_manual(values = c(NA,21)), 
    point_size = 3, 
    point_alpha = 0.5, 
    annotation_size = 6, 
    annotation_force = 10
  )
  save_figure_safe(
    with(tissue_plot, tissue_plot), 
    png, 
    paste0(output_path, i, "_OT_level1_primary_cl_only.png"), 
    width = plot_width * 2, 
    height = plot_width * 1.5, 
    res = plot_res, 
    units = plot_units
  )
}

# Celligner
fn <- paste0(celligner_path, "Celligner_aligned_data.csv")
celligner_mat <- readr::read_csv(fn, show_col_types = FALSE)
celligner_sample_ids <- celligner_mat[[1]]
celligner_mat <- as.data.frame(celligner_mat[,-1])
rownames(celligner_mat) <- celligner_sample_ids

tcga_level1_labels_df_truncated <- tcga_level1_labels_df
rownames(tcga_level1_labels_df_truncated) <- substr(rownames(tcga_level1_labels_df), 1, 15)

dat_mat <- celligner_mat
tissue_label <- rep_len(NA, nrow(dat_mat))
names(tissue_label) <- rownames(dat_mat)
ccle_ptr <- intersect(rownames(dat_mat), rownames(ccle_level1_labels_df))
tcga_ptr <- intersect(rownames(dat_mat), rownames(tcga_level1_labels_df_truncated))
tissue_label[ccle_ptr] <- ccle_level1_labels_df[ccle_ptr, "level_1"]
tissue_label[tcga_ptr] <- tcga_level1_labels_df_truncated[tcga_ptr, "level_1"]
tissue_label <- tolower(gsub("_", " ", tissue_label))

dataset <- rep_len(NA, nrow(dat_mat))
names(dataset) <- rownames(dat_mat)
dataset[ccle_ptr] <- ccle_level1_labels_df[ccle_ptr, "dataset"]
dataset[tcga_ptr] <- tcga_level1_labels_df_truncated[tcga_ptr, "dataset"]
dataset <- factor(dataset, levels = c("TCGA", "CCLE"))
dataset_ind <- !is.na(dataset)

nna_ind <- !is.na(tissue_label)

metrics_correlation_list[["Celligner"]] <- tissue_classifier_evaluation(
  data = dat_mat[dataset_ind,], 
  dataset_label = dataset[dataset_ind], 
  reference_dataset_class = "TCGA", 
  predicted_dataset_class = "CCLE", 
  tissue_label = tissue_label[dataset_ind], 
  labeled_ind = nna_ind[dataset_ind], 
  excluded_labels = c("adrenal gland", "testis"),
  scale_datasets_separately = FALSE, 
  knn_method = correlation_knn
)

set.seed(0)
tissue_plot_celligner <- tissue_visualizer(
  dat_mat, 
  reference_shape_label = "TCGA",
  labeled_ind = nna_ind, 
  color_var = tissue_label, 
  color_name = "tissue",
  shape_var = dataset, 
  shape_name = "dataset",
  umap_neighbors = 20, 
  pre_manifold_pca = TRUE, 
  max_pcs = 50, 
  umap_args = list(init = "normlaplacian", min_dist = 0.5, nn_method = "annoy", metric = "correlation"),
  primary_color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=n_labels)), 
  primary_fill_scale = scale_fill_manual(values = pals::kovesi.rainbow(n=n_labels)), 
  secondary_color_scale = scale_color_manual(values = c(NA, "#000000")), 
  primary_shape_scale = scale_shape_manual(values = c(3,NA)),
  secondary_shape_scale = scale_shape_manual(values = c(NA,21)), 
  point_size = 3, 
  point_alpha = 0.5, 
  annotation_size = 6, 
  annotation_force = 10, 
  knn_accuracy = FALSE,
  knn_predicted_dataset = "CCLE", 
  knn_title_method_name = "Celligner"
)
saveRDS(tissue_plot_celligner, file = paste0(output_path, "Celligner_OT_level1.rds"))
save_figure_safe(
  with(tissue_plot_celligner, tissue_plot), 
  png, 
  paste0(output_path, "Celligner_OT_level1.png"), 
  width = plot_width * 1.2, 
  height = plot_width * 0.9, 
  res = plot_res, 
  units = plot_units
)

# CODE-AE
codeae_data_path <- paste0(codeae_path, "data/preprocessed_dat/")
fn <- paste0(codeae_data_path, "xena_samples.txt")
tcga_samples <- readLines(fn)
fn <- paste0(codeae_data_path, "ccle_samples.txt")
ccle_samples <- readLines(fn)

codeae_result_path <- paste0(
  codeae_path, 
  "XieResearchGroup-CODE-AE-6dc17a5/", 
  "intermediate_results/encoded_features/"
)
fn <- paste0(codeae_result_path, "codeae_tcga_features.csv")
codae_tcga_mat <- readr::read_csv(fn)
fn <- paste0(codeae_result_path, "codeae_ccle_features.csv")
codae_ccle_mat <- readr::read_csv(fn)

tcga_level1_labels_df_truncated <- tcga_level1_labels_df
rownames(tcga_level1_labels_df_truncated) <- substr(rownames(tcga_level1_labels_df), 1, 15)

library(dplyr)
dat_mat <- bind_rows(codae_tcga_mat[,-1], codae_ccle_mat[,-1])
dat_mat <- as.data.frame(dat_mat)
rownames(dat_mat) <- c(tcga_samples, ccle_samples)

tissue_label <- rep_len(NA, nrow(dat_mat))
names(tissue_label) <- rownames(dat_mat)
ccle_ptr <- intersect(rownames(dat_mat), rownames(ccle_level1_labels_df))
tcga_ptr <- intersect(rownames(dat_mat), rownames(tcga_level1_labels_df_truncated))
tissue_label[ccle_ptr] <- ccle_level1_labels_df[ccle_ptr, "level_1"]
tissue_label[tcga_ptr] <- tcga_level1_labels_df_truncated[tcga_ptr, "level_1"]
tissue_label <- tolower(gsub("_", " ", tissue_label))

dataset <- rep_len(NA, nrow(dat_mat))
names(dataset) <- rownames(dat_mat)
dataset[ccle_ptr] <- ccle_level1_labels_df[ccle_ptr, "dataset"]
dataset[tcga_ptr] <- tcga_level1_labels_df_truncated[tcga_ptr, "dataset"]
dataset <- factor(dataset, levels = c("TCGA", "CCLE"))
dataset_ind <- !is.na(dataset)

nna_ind <- !is.na(tissue_label)

metrics_correlation_list[["CODE-AE"]] <- tissue_classifier_evaluation(
  data = dat_mat[dataset_ind,], 
  dataset_label = dataset[dataset_ind], 
  reference_dataset_class = "TCGA", 
  predicted_dataset_class = "CCLE", 
  tissue_label = tissue_label[dataset_ind], 
  labeled_ind = nna_ind[dataset_ind], 
  excluded_labels = c("adrenal gland", "testis"),
  scale_datasets_separately = FALSE, 
  knn_method = correlation_knn
)

saveRDS(metrics_correlation_list, file = paste0(output_path, "metrics_list.rds"))

set.seed(0)
tissue_plot_codeae <- tissue_visualizer(
  dat_mat, 
  reference_shape_label = "TCGA",
  labeled_ind = nna_ind, 
  color_var = tissue_label, 
  color_name = "tissue",
  shape_var = dataset, 
  shape_name = "dataset",
  umap_neighbors = 20, 
  pre_manifold_pca = FALSE, 
  max_pcs = 50, 
  umap_args = list(init = "normlaplacian", min_dist = 0.5, nn_method = "annoy", metric = "correlation"),
  primary_color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=n_labels)), 
  primary_fill_scale = scale_fill_manual(values = pals::kovesi.rainbow(n=n_labels)), 
  secondary_color_scale = scale_color_manual(values = c(NA, "#000000")), 
  primary_shape_scale = scale_shape_manual(values = c(3,NA)),
  secondary_shape_scale = scale_shape_manual(values = c(NA,21)), 
  point_size = 3, 
  point_alpha = 0.5, 
  annotation_size = 6, 
  annotation_force = 10, 
  knn_accuracy = FALSE,
  knn_predicted_dataset = "CCLE", 
  knn_title_method_name = "CODE-AE"
)
saveRDS(tissue_plot_codeae, file = paste0(output_path, "CODE_AE_OT_level1.rds"))
save_figure_safe(
  with(tissue_plot_codeae, tissue_plot), 
  png, 
  paste0(output_path, "CODE_AE_OT_level1.png"), 
  width = plot_width * 1.2, 
  height = plot_width * 0.9, 
  res = plot_res, 
  units = plot_units
)

# Celligner labels
fn <- paste0(celligner_path, "Celligner_info.csv")
celligner_info <- readr::read_csv(fn, show_col_types = FALSE)

celligner_info

library(dplyr)
ccle_ids <- ccle_level1_labels_df |> filter(dataset == "CCLE") |> rownames()

celligner_ptr <- match(ccle_ids, celligner_info[["sampleID"]])

table(
  ccle_level1_labels_df[["level_1"]], 
  celligner_info[["lineage"]][celligner_ptr], 
  useNA = "always"
)

table(
  ccle_level1_labels_df[["level_1"]], 
  celligner_info[["Yu_et_al_annotation"]][celligner_ptr], 
  useNA = "always"
)

table(
  ccle_level1_labels_df[["level_1"]], 
  celligner_info[["CancerCellNet_annotation"]][celligner_ptr], 
  useNA = "always"
)


