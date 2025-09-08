script_root <- Sys.getenv("MODAE_SCRIPT_PATH")
if (script_root == "") {
  stop("Please define the MODAE_SCRIPT_PATH environment variable.")
}
source(paste0(script_root, "R/setup.R"))

for (save_dir in save_dirs) {
  if(exists("var_list")) try(detach(var_list))
  var_list <- get_paths_and_parameters(save_dir, remove_incomplete = FALSE)
  attach(var_list)
  best_task <- readLines(paste0(save_path, "best_task.txt"))
  param_best_task_ind <- match(best_task, parameters[["task"]])
  shared_embedding_names <- paste0("z", 1:parameters[param_best_task_ind, "bottle_neck"])
  
  # Final predictions
  int_pred <- get_final_internal_predictions()
  internal_predictions_patient <- int_pred[["patient"]]
  internal_predictions_cl <- int_pred[["cl"]]
  
  ctrp_drugs <- get_ctrp_drugs()
  
  hm_color_scale <- circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  if (tcga_brca || scanb) {
    d_order_fn <- paste0(ctrp_path, "brca_drug_correlation_order.txt")
  } else {
    d_order_fn <- paste0(ctrp_path, "drug_correlation_order.txt")
  }
  drug_correlation_order <- readLines(d_order_fn)
  drug_correlation_order <- drug_correlation_order[drug_correlation_order %in% ctrp_drugs]
  patient_dr_cols <- grep("dr_pred_[0-9]+", colnames(internal_predictions_patient))
  patient_pred_cor <- cor(internal_predictions_patient[,patient_dr_cols])
  colnames(patient_pred_cor) <- rownames(patient_pred_cor) <- ctrp_drugs
  patient_pred_hm <- ComplexHeatmap::Heatmap(
    patient_pred_cor[drug_correlation_order,drug_correlation_order], 
    col = hm_color_scale, 
    name = "patient", 
    show_row_names = FALSE, 
    show_column_names = FALSE, 
    cluster_rows = FALSE, 
    cluster_columns = FALSE)
  
  if (FALSE) {
    sens_hm_color_scale <- circlize::colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
    #internal_predictions_patient[,patient_dr_cols]
    #drug_interest_ind <- ctrp_drugs %in% c("docetaxel", "paclitaxel", "doxorubicin", "gemcitabine")
    
    
    drug_name_map <- ctrp_drugs
    names(drug_name_map) <- colnames(internal_predictions_patient)[patient_dr_cols]
    temp <- internal_predictions_patient[,c(patient_dr_cols)]
    colnames(temp) <- drug_name_map[colnames(temp)]
    
    drug_sens_pair_plot <- GGally::ggpairs(
      temp[,c("docetaxel", "paclitaxel", "doxorubicin", "gemcitabine")]
    )
    save_figure_safe(
      drug_sens_pair_plot, 
      png, 
      paste0(best_dr_res_path, "patient_predicted_aac_known_drugs_correlation.png"), 
      width = plot_width * 0.5, 
      height = plot_width * 0.5, 
      res = plot_res, 
      units = plot_units
    )
    
    hist(apply(temp, 1, sd))
    temp_z <- scale(temp, center = TRUE, scale = TRUE)
    
    patient_pred_sens_hm <- ComplexHeatmap::Heatmap(
      temp_z, 
      #col = sens_hm_color_scale, 
      name = "patient AAC z-score", 
      show_row_names = FALSE, 
      show_column_names = FALSE, 
      cluster_rows = TRUE, 
      cluster_columns = TRUE, 
      show_row_dend = FALSE, 
      show_column_dend = FALSE)
    save_figure_safe(
      patient_pred_sens_hm, 
      png, 
      paste0(best_dr_res_path, "patient_aac_z_score_all_drugs_heatmap.png"), 
      width = plot_width * 4, 
      height = plot_width * 2, 
      res = plot_res, 
      units = plot_units
    )
    rownames(temp) <- internal_predictions_patient[["X"]]
    write.csv(temp, gzfile(paste0(best_dr_res_path, "patient_aac_predictions.csv.gz")))
    
    if (FALSE) {
      long_df <- reshape2::melt(
        internal_predictions_patient[,c(1, patient_dr_cols)], 
        id.vars = "X", 
        measure.vars = colnames(internal_predictions_patient)[patient_dr_cols], 
        variable.name = "drug", 
        value.name = "aac"
      )
      
      drug_name_map <- ctrp_drugs
      names(drug_name_map) <- colnames(internal_predictions_patient)[patient_dr_cols]
      long_df[["drug_name"]] <- drug_name_map[long_df[["drug"]]]
      long_df[[""]]
      
      scatter_plots <- ggplot(long_df, aes)
    }
    
  }
  
  cl_dr_cols <- grep("dr_pred_[0-9]+", colnames(internal_predictions_cl))
  cl_pred_cor <- cor(internal_predictions_cl[,cl_dr_cols])
  
  colnames(cl_pred_cor) <- rownames(cl_pred_cor) <- ctrp_drugs
  
  cl_pred_hm <- ComplexHeatmap::Heatmap(
    cl_pred_cor[drug_correlation_order,drug_correlation_order], 
    col = hm_color_scale, 
    name = "cl", 
    show_row_names = FALSE, 
    show_column_names = FALSE, 
    cluster_rows = FALSE, 
    cluster_columns = FALSE)
  
  if (FALSE) {
    # Check difference compared to parameter search
    cor_diff <- full_cor_impute - cl_pred_cor
    table(ctrp_drugs == colnames(cl_pred_cor))
    
    cor_diff_order <- hclust(
      dist(
        cor_diff, 
        method = "euclidean"
      ), 
      method = "complete"
    )$order
    
    cor_diff_hm <- ComplexHeatmap::Heatmap(
      cor_diff[cor_diff_order,cor_diff_order], 
      name = "AAC cor diff", 
      cluster_rows = FALSE, 
      cluster_columns = FALSE, 
      show_column_names = FALSE, 
      show_row_names = FALSE
    )
  }
  
  library(ComplexHeatmap)
  save_figure_safe(
    patient_pred_hm, #%v% ctrp_brca_true_hm, 
    png, 
    paste0(best_dr_res_path, "patient_aac_correlation_heatmap.png"), 
    width = plot_width * 2, 
    height = plot_width * 2, 
    res = plot_res, 
    units = plot_units
  )
  
  save_figure_safe(
    cl_pred_hm, #%v% ctrp_true_hm, 
    png, 
    paste0(best_dr_res_path, "cl_aac_correlation_heatmap.png"), 
    width = plot_width * 2, 
    height = plot_width * 2, 
    res = plot_res, 
    units = plot_units
  )
}