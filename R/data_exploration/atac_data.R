source("../setup.R")

fn <- paste0(data_root, "tcga_atac_seq/aav1898_data_s6.xlsx")

footprint_depth <- readxl::read_xlsx(
  fn, 
  sheet = 1, 
  range = "A18:OR1782", 
  col_names = TRUE
)
flanking_access <- readxl::read_xlsx(
  fn, 
  sheet = 2, 
  range = "A18:OR1782", 
  col_names = TRUE
)

fd_mat <- as.matrix(footprint_depth[,-(1:4)])
fa_mat <- as.matrix(flanking_access[,-(1:4)])

# Tissues
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

fn <- paste0(ccle_path, "oncotree_level1.csv")
ccle_level1_labels_df <- read.csv(fn, row.names = 1, header = TRUE)
ccle_level1_labels_df[["dataset"]] <- "CCLE"

fn <- paste0(patient_expression_root_dir, "oncotree_level1.csv")
tcga_level1_labels_df <- read.csv(fn, row.names = 1, header = TRUE)
tcga_level1_labels_df[["dataset"]] <- "TCGA"

#n_labels <- length(union(
#  unique(ccle_level1_labels_df[["level_1"]]), 
#  unique(tcga_level1_labels_df[["level_1"]])
#))

oncotree_idx <- match(
  colnames(fd_mat), 
  substr(rownames(tcga_level1_labels_df), 1, 12)
)
n_labels <- length(unique(tcga_level1_labels_df[oncotree_idx,"level_1"]))

COPS::umap_viz(
  t(fd_mat), 
  tcga_level1_labels_df[oncotree_idx,"level_1"], 
  category_label = "tissue",
  color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=n_labels)), 
  umap_args = list(
    init = "normlaplacian", 
    n_neighbors = 20, 
    min_dist = 0.5
  )
)

COPS::umap_viz(
  t(fa_mat), 
  tcga_level1_labels_df[oncotree_idx,"level_1"], 
  category_label = "tissue",
  color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=n_labels)), 
  umap_args = list(
    init = "normlaplacian", 
    n_neighbors = 20, 
    min_dist = 0.5
  )
)