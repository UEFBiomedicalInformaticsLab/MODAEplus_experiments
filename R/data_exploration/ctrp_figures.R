script_root <- Sys.getenv("MODAE_SCRIPT_PATH")
if (script_root == "") {
  stop("Please define the MODAE_SCRIPT_PATH environment variable.")
}
source(paste0(script_root, "R/setup.R"))

ctrp_dr_df <- get_ctrp_dr_df()
db_brca_drugs <- readLines(paste0(ctrp_path, "drug_bank_brca_drugs.txt"))
if (FALSE) {
  # Old code from validation_analysis.R
  high_confidence_drugs <- dr_r2_best[dr_r2_best[["cl_test"]] > 0.1, "drug"]
  ext_drugs <- unique(ext_dr_res[["drug"]])
  
  if (scanb || tcga_brca) {
    brca_cl_ids <- unique(ctrp_dr_df[
      ctrp_dr_df[["cl_type"]] == "BRCA" & 
        !is.na(ctrp_dr_df[["cl_type"]]), 
      "cl_id"])
    cl_of_interest <- brca_cl_ids
    drugs_of_interest <- union(db_brca_drugs, ext_drugs)
    table(known_or_ext = ctrp_drugs %in% drugs_of_interest, 
          confidence = ctrp_drugs %in% high_confidence_drugs)
  } else {
    cl_of_interest <- TRUE
    drugs_of_interest <- ctrp_drugs
  }
  
  final_drug_list <- intersect(drugs_of_interest, high_confidence_drugs)
} 
brca_cl_ids <- unique(ctrp_dr_df[
  ctrp_dr_df[["cl_type"]] == "BRCA" & 
    !is.na(ctrp_dr_df[["cl_type"]]), 
  "cl_id"])

cor_dist_order_impute <- function(x, na_mean_cut = 0.5, linkage = "complete") {
  cor_dist <- 1 - cor(x, use = "pairwise.complete.obs")
  if (any(is.na(cor_dist) | is.nan(cor_dist))) {
    exclude_ind <- apply(
      is.na(cor_dist) | is.nan(cor_dist) | cor_dist == 0, 
      1, 
      mean
    ) > na_mean_cut
    cor_dist <- impute::impute.knn(
      cor_dist[!exclude_ind, !exclude_ind], k = 5)$data
  } else {
    exclude_ind <- rep(FALSE, nrow(x))
  }
  hc_order <- hclust(
    dist(cor_dist, method = "euclidean"), 
    method = linkage
  )$order
  
  return(colnames(x)[!exclude_ind][hc_order])
}

ctrp_dr_df_nna <- ctrp_dr_df[
  !is.na(ctrp_dr_df[["aac"]]) & 
    !is.nan(ctrp_dr_df[["aac"]]),
]
ctrp_dr_df_nna_mat <- reshape2::acast(
  ctrp_dr_df_nna, 
  cl_id ~ drug, 
  fun.aggregate = mean, 
  value.var = "aac")

ctrp_plot_path <- paste0(ctrp_path, "figures/")
dir.create(ctrp_plot_path)

hc_cl_order <- cor_dist_order_impute(t(ctrp_dr_df_nna_mat))
hc_drug_order <- cor_dist_order_impute(ctrp_dr_df_nna_mat)

# Annotate breast cancer drugs and cell-lines
col_vec <- pals::kovesi.rainbow(n=2)
names(col_vec) <- c("TRUE", "FALSE")
drug_ha <- ComplexHeatmap::rowAnnotation(
  BRCA_drug = hc_drug_order %in% db_brca_drugs, 
  col = list(BRCA_drug = col_vec))
cl_ha <- ComplexHeatmap::HeatmapAnnotation(
  BRCA_CL = hc_cl_order %in% brca_cl_ids, 
  col = list(BRCA_CL = col_vec))

save_figure_safe(
  ComplexHeatmap::Heatmap(
    t(ctrp_dr_df_nna_mat)[hc_drug_order, hc_cl_order], 
    name = "CTRP AAC", 
    cluster_rows = FALSE, 
    cluster_columns = FALSE, 
    show_row_names = FALSE, 
    show_column_names = FALSE, 
    top_annotation = cl_ha, 
    left_annotation = drug_ha
  ), 
  png, 
  paste0(ctrp_plot_path, "ctrp_aac_heatmap.png"), 
  width = plot_width * 2, 
  height = plot_height * 2, 
  res = plot_res,
  units = plot_units
)

ctrp_dr_df_brca_nna <- ctrp_dr_df[
  ctrp_dr_df[["cl_type"]] == "BRCA" & 
    !is.na(ctrp_dr_df[["aac"]]) &
    !is.nan(ctrp_dr_df[["aac"]]),
]
ctrp_dr_mat_brca <- reshape2::acast(
  ctrp_dr_df_brca_nna, 
  cl_id ~ drug, 
  fun.aggregate = mean, 
  value.var = "aac")

hc_cl_order <- cor_dist_order_impute(t(ctrp_dr_mat_brca[,db_brca_drugs]))
hc_drug_order <- cor_dist_order_impute(ctrp_dr_mat_brca[,db_brca_drugs])

save_figure_safe(
  ComplexHeatmap::Heatmap(
    t(ctrp_dr_mat_brca)[hc_drug_order, hc_cl_order], 
    name = "CTRP AAC", 
    cluster_rows = FALSE, 
    cluster_columns = FALSE
  ),
  png, 
  paste0(ctrp_plot_path, "brca_known_drug_ctrp_aac_heatmap.png"), 
  width = plot_width, 
  height = plot_height, 
  res = plot_res, 
  units = plot_units
)

# Correlations
hm_color_scale <- circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
ctrp_brca_true_cor <- cor(ctrp_dr_mat_brca, use = "pairwise.complete.obs")
ctrp_brca_drug_order <- cor_dist_order_impute(ctrp_dr_mat_brca)
ctrp_brca_true_hm <- ComplexHeatmap::Heatmap(
  ctrp_brca_true_cor[ctrp_brca_drug_order,ctrp_brca_drug_order], 
  col = hm_color_scale, 
  name = "CTRP BRCA", 
  show_row_names = FALSE, 
  show_column_names = FALSE, 
  cluster_rows = FALSE, 
  cluster_columns = FALSE, 
)
save_figure_safe(
  ctrp_brca_true_hm, 
  png, 
  paste0(ctrp_path, "ctrp_brca_true_cor_heatmap.png"), 
  width = plot_width * 2, 
  height = plot_width * 2, 
  res = plot_res, 
  units = plot_units
)
fnw <- paste0(ctrp_path, "brca_drug_correlation_order.txt")
writeLines(ctrp_brca_drug_order, fnw)

ctrp_dr_df_nna <- ctrp_dr_df[!is.na(ctrp_dr_df[["aac"]]),]
ctrp_dr_mat <- reshape2::acast(
  ctrp_dr_df_nna, 
  cl_id ~ drug, 
  fun.aggregate = mean, 
  value.var = "aac")

ctrp_true_cor <- cor(ctrp_dr_mat, use = "pairwise.complete.obs")
ctrp_drug_order <- cor_dist_order_impute(ctrp_dr_mat)
ctrp_true_hm <- ComplexHeatmap::Heatmap(
  ctrp_true_cor[ctrp_drug_order,ctrp_drug_order], 
  col = hm_color_scale, 
  name = "CCLE CTRP", 
  show_row_names = FALSE, 
  show_column_names = FALSE, 
  cluster_rows = FALSE, 
  cluster_columns = FALSE)
save_figure_safe(
  ctrp_true_hm, 
  png, 
  paste0(ctrp_path, "ctrp_true_cor_heatmap.png"), 
  width = plot_width * 2, 
  height = plot_width * 2, 
  res = plot_res, 
  units = plot_units
)
writeLines(ctrp_drug_order, paste0(ctrp_path, "drug_correlation_order.txt"))

# Check baseline predictor correlations for reference
base_dr_dir <- paste0(output_dir, "baseline_results/drug_response/ctrp/")
baseline_dr_pred_file <- paste0(base_dr_dir, "pc10_en_preds.csv.gz")
base_dr_pred <- read.csv(baseline_dr_pred_file, header = TRUE, row.names = 1)
table(is.na(base_dr_pred))
table(apply(base_dr_pred, 2, sd) == 0)

base_dr_pred_cor <- cor(base_dr_pred)

# Fix R index garbage
ctrp_drug_order_strip <- gsub("[^a-zA-Z0-9]", ".", ctrp_drug_order)
ctrp_drug_order_strip <- ifelse(
  grepl("^[0-9]", ctrp_drug_order_strip), 
  paste0("X", ctrp_drug_order_strip), 
  ctrp_drug_order_strip
)

base_cl_pred_hm <- ComplexHeatmap::Heatmap(
  base_dr_pred_cor[ctrp_drug_order_strip,ctrp_drug_order_strip], 
  col = hm_color_scale, 
  name = "CCLE pred", 
  show_row_names = FALSE, 
  show_column_names = FALSE, 
  cluster_rows = FALSE, 
  cluster_columns = FALSE)

save_figure_safe(
  base_cl_pred_hm, 
  png, 
  paste0(ctrp_path, "cl_aac_correlation_heatmap_pc10_en_baseline.png"), 
  width = plot_width * 2, 
  height = plot_width * 2, 
  res = plot_res, 
  units = plot_units
)
