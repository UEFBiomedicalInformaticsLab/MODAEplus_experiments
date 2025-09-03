source("../setup.R")

for (save_dir in save_dirs) {
  if(exists("var_list")) try(detach(var_list))
  var_list <- get_paths_and_parameters(save_dir, remove_incomplete = FALSE)
  attach(var_list)
  best_task <- readLines(paste0(save_path, "best_task.txt"))
  param_best_task_ind <- match(best_task, parameters[["task"]])
  shared_embedding_names <- paste0("z", 1:parameters[param_best_task_ind, "bottle_neck"])
  
  best_dr_res_path <- paste0(save_path, "../plots/", model_name, "_drug_responses/")
  dir.create(best_dr_res_path, recursive = TRUE)
  
  # CV drug performance for best setting
  drug_response_r2 <- read_result_files("drug_response_r2.csv.gz", path = save_path, model_name = model_name)
  drug_response_r2_best <- drug_response_r2[drug_response_r2[["task"]] == best_task,]
  n_drugs <- length(unique(drug_response_r2_best[["drug"]]))
  drug_response_r2_best[["drug"]] <- process_ctrp_drug_names(
    tolower(
      fix_ctrp_drug_names(
        drug_response_r2_best[["drug"]], 
        n_drugs = n_drugs
      )
    )
  )
  dr_r2_best_long <- plyr::ddply(drug_response_r2_best, c("drug", "dataset"), function(x) data.frame(r2 = mean(x$dr_r2)))
  dr_r2_best <- reshape2::dcast(dr_r2_best_long, drug ~ dataset, value.var = "r2")
  #plot(dr_r2_best[["cl_train"]], dr_r2_best[["cl_test"]])
  dr_r2_best[["rank"]] <- rank(-dr_r2_best[["cl_test"]])
  #dr_r2_best[dr_r2_best[["rank"]] < 20,]
  
  plot1 <- ggplot(dr_r2_best_long, aes(pmax(0, r2), fill = dataset)) + 
    geom_density(alpha = 0.4, adjust = 0.5) + 
    theme_bw() + scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
  save_figure_safe(
    plot1, 
    png, 
    paste0(best_dr_res_path, "drug_r2_density_plot.png"), 
    width = plot_width, 
    height = plot_height, 
    res = plot_res, 
    units = plot_units
  )
  
  # External drug performance
  ext_dr_results <- read.csv(
    paste0(save_path, "external_evaluation/external_drug_response_validation_performance.csv"), 
    row.names = NULL, header = TRUE)
  ext_dr_res <- plyr::join(ext_dr_results, dr_r2_best, by = c("drug"), type = "left")
  #ext_dr_res[ext_dr_res[["cl_test"]] > 0.1,]
  
  # Final training embeddings
  fn_ext <- paste0(save_path, "external_evaluation/external_drug_response_validation_embeddings.csv.gz")
  fn_cl <- paste0(save_path, "external_evaluation/internal_drug_response_validation_embeddings.csv.gz")
  fn_p <- paste0(save_path, "external_evaluation/internal_survival_validation_embeddings.csv.gz")
  if (!file.exists(fn_cl)) {
    fn_cl <- paste0(save_path, "external_evaluation/internal_drug_response_validation_final_embeddings.csv.gz")
  }
  if (!file.exists(fn_p)) {
    fn_p <- paste0(save_path, "external_evaluation/internal_survival_validation_final_embeddings.csv.gz")
  }
  final_embeddings_ext <- read.csv(fn_ext, row.names = NULL, header = TRUE)
  final_embeddings_cl <- read.csv(fn_cl, row.names = NULL, header = TRUE)
  final_embeddings_p <- read.csv(fn_p, row.names = NULL, header = TRUE)
  final_dr_embeddings <- rbind(
    final_embeddings_cl, 
    #final_embeddings_ext[final_embeddings_ext[["dataset"]] == "bruna_pdtc",])
    final_embeddings_ext)
  
  # Use oncotree to map cell lines to tissue of origin
  cell_line_oncotree_mapping_file <- 'Model_oncotree.csv'
  cell_line_oncotree_mappings <- read.csv(paste0(cell_line_expression_root_dir, cell_line_oncotree_mapping_file), header = TRUE, row.names = 1)
  cell_line_brca_map <- cell_line_oncotree_mappings[["level_2"]] == "BRCA"
  names(cell_line_brca_map) <- rownames(cell_line_oncotree_mappings)
  
  final_dr_embeddings[["cell_line_level1"]] <- cell_line_oncotree_mappings[final_dr_embeddings[["X"]], "level_1"]
  final_dr_embeddings[["brca_cell_line"]] <- final_dr_embeddings[["X"]] %in% names(which(cell_line_brca_map))
  final_dr_embeddings[final_dr_embeddings[["brca_cell_line"]], "dataset"] <- "cl_brca"
  
  # UMAP of embeddings labelled by cell line type
  if (FALSE) {
    cl_ind <- grepl("ACH", final_dr_embeddings[["X"]])
    set.seed(0)
    COPS::umap_viz(
      final_dr_embeddings[cl_ind, shared_embedding_names], 
      category = final_dr_embeddings[cl_ind, "cell_line_level1"], 
      category_label = "source tissue", pre_manifold_pca = FALSE, 
      umap_args = list(init = "normlaplacian", min_dist = 0.5)) + 
      scale_color_manual(values = pals::kovesi.rainbow(length(table(final_dr_embeddings[cl_ind, "cell_line_level1"]))))
    
    set.seed(0)
    COPS::umap_viz(
      final_dr_embeddings[cl_ind, shared_embedding_names], 
      category = final_dr_embeddings[cl_ind, "cell_line_level1"] == "LYMPH", 
      category_label = "lymph", pre_manifold_pca = FALSE, 
      umap_args = list(init = "normlaplacian", min_dist = 0.5))
  }
  # UMAP of CL expression data
  if (FALSE) {
    ccle_expression <- read.csv(
      paste0(ccle_path, "CCLE_expression.csv"), 
      row.names = 1, header = TRUE)
    ccle_info <- read.csv(
      paste0(ccle_path, "Model_augmented.csv"), 
      row.names = 1, header = TRUE)
    
    cl_filter <- ccle_info[rownames(ccle_expression), "solid"] == "solid"
    
    # Lymph is a prominent cluster
    set.seed(0)
    COPS::umap_viz(
      as.matrix(ccle_expression), 
      category = cell_line_oncotree_mappings[rownames(ccle_expression), "level_1"] == "LYMPH", 
      category_label = "lymph", pre_manifold_pca = TRUE, 
      umap_args = list(init = "normlaplacian", min_dist = 0.5))
    set.seed(0)
    COPS::umap_viz(
      as.matrix(ccle_expression), 
      category = cell_line_oncotree_mappings[rownames(ccle_expression), "level_1"], 
      category_label = "tissue", pre_manifold_pca = TRUE, 
      umap_args = list(init = "normlaplacian", min_dist = 0.5)) + 
      scale_color_manual(values = pals::kovesi.rainbow(length(table(cell_line_oncotree_mappings[rownames(ccle_expression), "level_1"]))))
  }
  # UMAP of embeddings labelled by datasets
  if (FALSE) {
    set.seed(0)
    COPS::umap_viz(
      final_dr_embeddings[,shared_embedding_names], 
      category = final_dr_embeddings[["dataset"]], category_label = "dataset", pre_manifold_pca = FALSE, 
      umap_args = list(init = "normlaplacian", min_dist = 0.5))
    
    brca_only_filter <- final_dr_embeddings[["dataset"]] %in% c("bruna_pdtc", "cl_brca")
    set.seed(0)
    COPS::umap_viz(
      final_dr_embeddings[brca_only_filter,shared_embedding_names], 
      category = final_dr_embeddings[brca_only_filter, "dataset"], 
      category_label = "dataset", pre_manifold_pca = FALSE, 
      umap_args = list(init = "normlaplacian", min_dist = 0.5))
  }
  
  # Final predictions
  int_pred <- get_final_internal_predictions()
  internal_predictions_patient <- int_pred[["patient"]]
  internal_predictions_cl <- int_pred[["cl"]]
  
  # Prediction correlation heatmap
  if (FALSE) {
    patient_dr_cols <- grep("dr_pred_[0-9]+", colnames(internal_predictions_patient))
    #hist(as.matrix(internal_predictions_patient[,patient_dr_cols]))
    patient_pred_cor <- cor(internal_predictions_patient[,patient_dr_cols])
    
    patient_pred_hm <- ComplexHeatmap::Heatmap(
      patient_pred_cor, name = "patient", 
      show_row_names = FALSE, 
      show_column_names = FALSE)
    
    cl_dr_cols <- grep("dr_pred_[0-9]+", colnames(internal_predictions_cl))
    cl_pred_cor <- cor(internal_predictions_cl[,cl_dr_cols])
    
    cl_pred_hm <- ComplexHeatmap::Heatmap(
      cl_pred_cor, name = "cl", 
      show_row_names = FALSE, 
      show_column_names = FALSE)
    
    library(ComplexHeatmap)
    patient_pred_hm %v% cl_pred_hm
  }
  
  # Check average prediction for each drug to find highest sensitivity drugs
  patient_dr_pred_mean <- apply(
    internal_predictions_patient[
      , grep("dr_pred_[0-9]+", colnames(internal_predictions_patient))], 
    2, 
    mean
  )
  
  # Drug names
  ctrp_drugs <- get_ctrp_drugs(n_drugs = length(patient_dr_pred_mean))
  
  patient_dr_pred_mean <- data.frame(
    drug = ctrp_drugs, 
    patient_aac_mean = patient_dr_pred_mean
  )
  
  # Combine with CV performance
  dr_r2_best <- plyr::join(dr_r2_best, patient_dr_pred_mean)
  dr_r2_best[["r2_rank"]] <- dr_r2_best[["rank"]]
  dr_r2_best[["patient_aac_rank"]] <- rank(-dr_r2_best[["patient_aac_mean"]])
  #plot(dr_r2_best$cl_test, dr_r2_best$patient_aac_mean)
  #plot(dr_r2_best$r2_rank, dr_r2_best$patient_aac_rank)
  #best_drugs <- dr_r2_best[dr_r2_best$r2_rank <= 50 & dr_r2_best$patient_aac_rank <= 50, "drug"]
  
  # Get list of BRCA drugs from drugbank
  db_brca_drugs <- readLines(paste0(ctrp_path, "drug_bank_brca_drugs.txt"))
  
  brca_dr_r2_best_db_subset <- dr_r2_best[dr_r2_best[["drug"]] %in% db_brca_drugs, ]
  brca_dr_r2_best_db_subset[["r2_rank"]] <- rank(-brca_dr_r2_best_db_subset[["cl_test"]])
  brca_dr_r2_best_db_subset[["patient_aac_rank"]] <- rank(-brca_dr_r2_best_db_subset[["patient_aac_mean"]])
  
  known_brca_drugs_perf_cleaned <- brca_dr_r2_best_db_subset[
    , c("drug", "cl_test", "r2_rank", "patient_aac_mean", "patient_aac_rank")]
  colnames(known_brca_drugs_perf_cleaned)[c(2,4)] <- c("cl_test_r2", "patient_aac_pred_mean")
  
  # CTRP drug response data
  ctrp_dr_df <- get_ctrp_dr_df()
  
  ctrp_dr_df_brca <- ctrp_dr_df[
    ctrp_dr_df[["cl_type"]] == "BRCA" & !is.na(ctrp_dr_df[["cl_type"]]), ]
  
  ctrp_dr_df_brca_mean <- plyr::ddply(
    ctrp_dr_df_brca, c("drug"), 
    function(x) data.frame(
      brca_cl_aac_mean = mean(x$aac, na.rm = TRUE), 
      brca_cl_aac_drugwise_z_mean = mean(x$aac_drugwise_z_score, na.rm = TRUE)))
  ctrp_dr_df_brca_mean[["drug"]] <- tolower(ctrp_dr_df_brca_mean[["drug"]])
  
  known_brca_drugs_perf_cleaned <- plyr::join(
    known_brca_drugs_perf_cleaned, 
    ctrp_dr_df_brca_mean, type = "left")
  
  # Baseline drug-response prediction models
  {
    #old1_elnet_dr_r2 <- read.csv(paste0(xia_path, "elasticnet_aac1.csv"))
    #old2_elnet_dr_r2 <- read.csv(paste0(xia_path, "elasticnet_aac_old_fixed.csv"))
    #old3_elnet_dr_r2 <- read.csv(paste0(xia_path, "elasticnet_aac_old.csv"))
    
    fs1_elnet_dr_r2 <- read.csv(paste0(xia_path, "fs_elasticnet_aac_old.csv"))
    #fs2_elnet_dr_r2 <- read.csv(paste0(xia_path, "fs_elasticnet_aac_old_fixed.csv"))
    
    pca10_elnet_dr_r2 <- read.csv(paste0(xia_path, "pca10_elasticnet_aac_old.csv"))
    pca100_elnet_dr_r2 <- read.csv(paste0(xia_path, "pca100_elasticnet_aac_old.csv"))
    
    fs1_elnet_dr_r2[["model"]] <- "fs_elnet"
    pca10_elnet_dr_r2[["model"]] <- "pca10_elnet"
    pca100_elnet_dr_r2[["model"]] <- "pca100_elnet"
    
    baseline_dr <- Reduce(COPS::rbind_fill, list(fs1_elnet_dr_r2, pca10_elnet_dr_r2, pca100_elnet_dr_r2))
    baseline_dr_summary <- plyr::ddply(
      baseline_dr, c("model", "drug"), 
      function(x) as.data.frame(lapply(x[grep("_mse$|_r2$", colnames(x))], mean)))
    
    dl_dr_summary <- dr_r2_best[,c("drug", "cl_train", "cl_test")]
    
    colnames(dl_dr_summary) <- c("drug", "train_r2", "test_r2")
    dl_dr_summary[["model"]] <- "MODAE"
    
    combined_dr_table <- COPS::rbind_fill(dl_dr_summary, baseline_dr_summary)
  }
  
  # Known BRCA drug performances
  if (scanb || tcga_brca) {
    # Barplots
    sensitive_brca_drugs <- ctrp_dr_df_brca_mean[rank(-ctrp_dr_df_brca_mean[["brca_cl_aac_mean"]]) <= 10, "drug"]
    resistant_brca_drugs <- ctrp_dr_df_brca_mean[rank(ctrp_dr_df_brca_mean[["brca_cl_aac_mean"]]) <= 10, "drug"]
    plot_selections <- list(
      known_brca_drugs = tolower(combined_dr_table[["drug"]]) %in% db_brca_drugs, 
      brca_sensitive_drug = tolower(combined_dr_table[["drug"]]) %in% tolower(sensitive_brca_drugs),
      brca_resistant_drug = tolower(combined_dr_table[["drug"]]) %in% tolower(resistant_brca_drugs)
    )
    combined_dr_table[["train_rmse"]] <- sqrt(combined_dr_table[["train_mse"]])
    combined_dr_table[["test_rmse"]] <- sqrt(combined_dr_table[["test_mse"]])
    for (selection_name in names(plot_selections)){
      rows <- plot_selections[[selection_name]]
      temp <- combined_dr_table[rows,]
      write.csv(temp, paste0(best_dr_res_path, selection_name, "metrics.csv"))
      for (metric in c("train_r2", "test_r2", "train_rmse", "test_rmse")) {
        plot1 <- ggplot(temp, aes(x = model, y = !!ggplot2::ensym(metric), fill = model)) + 
          geom_bar(stat = "identity") + theme_bw() + theme(axis.text.x = element_blank()) + 
          scale_fill_brewer(palette = "Dark2") + facet_wrap(drug ~., scales = "fixed")
        save_figure_safe(
          plot1, 
          png, 
          paste0(best_dr_res_path, selection_name, "_", metric, ".png"), 
          width = plot_width, 
          height = plot_height, 
          res = plot_res, 
          units = plot_units
        )
      }
      temp <- reshape2::melt(
        temp, id.vars = c("drug", "model"), 
        variable.name = "metric", 
        value.name = "value")
      
      temp <- temp[temp[["metric"]] %in% c("train_rmse", "test_rmse", "train_r2", "test_r2"),]
      temp[["model"]] <- factor(
        as.character(temp[["model"]]), 
        levels = c("MODAE", "fs_elnet", "pca10_elnet", "pca100_elnet"))
      
      plot2_r2 <- ggplot(temp[grep("_r2$", temp[["metric"]]), ], aes(x = metric, y = value, fill = model)) + 
        geom_bar(stat = "identity", position = "dodge") + 
        theme_bw() + #theme(axis.text.x = element_blank()) + 
        scale_fill_brewer(palette = "Dark2") + facet_wrap(drug ~., scales = "fixed")
      save_figure_safe(
        plot2_r2, 
        png, 
        paste0(best_dr_res_path, selection_name, "_r2.png"), 
        width = plot_width, 
        height = plot_height, 
        res = plot_res, 
        units = plot_units
      )
      
      plot2_rmse <- ggplot(temp[grep("_rmse$", temp[["metric"]]), ], aes(x = metric, y = value, fill = model)) + 
        geom_bar(stat = "identity", position = "dodge") + 
        theme_bw() + #theme(axis.text.x = element_blank()) + 
        scale_fill_brewer(palette = "Dark2") + facet_wrap(drug ~., scales = "fixed")
      save_figure_safe(
        plot2_rmse, 
        png, 
        paste0(best_dr_res_path, selection_name, "_rmse.png"), 
        width = plot_width, 
        height = plot_height, 
        res = plot_res, 
        units = plot_units
      )
    }
    # Box + swarm plots from CV results
    dr_r2_best_cv <- reshape2::dcast(drug_response_r2_best, fold + run + drug ~ dataset, value.var = "dr_r2")
    colnames(dr_r2_best_cv) <- c("fold", "run", "drug", "test_r2", "train_r2")
    dr_r2_best_cv[["model"]] <- "MODAE"
    full_combined_dr_table <- COPS::rbind_fill(
      dr_r2_best_cv[,c("model", "drug", "train_r2", "test_r2")], 
      baseline_dr[,c("model", "drug", "train_r2", "test_r2", "train_mse", "test_mse")])
    full_combined_dr_table[["model"]] <- factor(
      as.character(full_combined_dr_table[["model"]]), 
      levels = c("MODAE", "fs_elnet", "pca10_elnet", "pca100_elnet"))
    plot_selections <- list(
      known_brca_drugs = tolower(full_combined_dr_table[["drug"]]) %in% db_brca_drugs, 
      brca_sensitive_drug = tolower(full_combined_dr_table[["drug"]]) %in% tolower(sensitive_brca_drugs),
      brca_resistant_drug = tolower(full_combined_dr_table[["drug"]]) %in% tolower(resistant_brca_drugs)
    )
    full_combined_dr_table[["train_rmse"]] <- sqrt(full_combined_dr_table[["train_mse"]])
    full_combined_dr_table[["test_rmse"]] <- sqrt(full_combined_dr_table[["test_mse"]])
    for (selection_name in names(plot_selections)){
      rows <- plot_selections[[selection_name]]
      temp <- full_combined_dr_table[rows,]
      write.csv(temp, paste0(best_dr_res_path, selection_name, "metrics_cv.csv"))
      for (metric in c("train_r2", "test_r2", "train_rmse", "test_rmse")) {
        plot1 <- ggplot(temp, aes(x = model, y = !!ggplot2::ensym(metric))) + 
          geom_boxplot(outlier.shape = NA) + 
          ggbeeswarm::geom_beeswarm(aes(color = model), cex = 3, size = 5, shape = "+") + 
          theme_bw() + theme(axis.text.x = element_blank()) + 
          scale_fill_brewer(palette = "Dark2") + 
          facet_wrap(drug ~., scales = "fixed")
        save_figure_safe(
          plot1, 
          png, 
          paste0(best_dr_res_path, selection_name, "_", metric, ".png"), 
          width = plot_width, 
          height = plot_height, 
          res = plot_res, 
          units = plot_units
        )
      }
    }
  }
  
  # TODO:Bruna AUC
  
  ## Heatmaps
  # internal_predictions_patient
  # internal_predictions_cl
  # (drug_response_r2_best), (dr_r2_best)
  # ext_dr_results
  # ctrp_dr_df_brca
  
  high_confidence_drugs <- dr_r2_best[dr_r2_best[["cl_test"]] > 0.1, "drug"]
  ext_drugs <- unique(ext_dr_res[["drug"]])
  drugs_of_interest <- ctrp_drugs
  db_brca_drugs <- readLines(paste0(ctrp_path, "drug_bank_brca_drugs.txt"))
  
  # Get the clustered drug order from CTRP AAC correlation
  brca_d_order_fn <- paste0(ctrp_path, "brca_drug_correlation_order.txt")
  brca_drug_correlation_order <- readLines(brca_d_order_fn)
  d_order_fn <- paste0(ctrp_path, "drug_correlation_order.txt")
  drug_correlation_order <- readLines(d_order_fn)
  
  # Omit mismatched / non-unique drugs (tipifarnib P1&2)
  drug_correlation_order <- drug_correlation_order[drug_correlation_order %in% ctrp_drugs]
  brca_drug_correlation_order <- brca_drug_correlation_order[brca_drug_correlation_order %in% db_brca_drugs]
  
  if (tcga_brca || scanb) {
    drug_selection <- brca_drug_correlation_order
  } else {
    drug_selection <- high_confidence_drugs
  }
  
  prediction_mat_patient <- as.matrix(internal_predictions_patient[,grep("dr_pred_[0-9]+", colnames(internal_predictions_patient))])
  colnames(prediction_mat_patient) <- ctrp_drugs
  rownames(prediction_mat_patient) <- internal_predictions_patient[["X"]]
  save_figure_safe(
    ComplexHeatmap::Heatmap(
      t(prediction_mat_patient)[drug_selection,], 
      name = "predicted AAC", 
      cluster_rows = TRUE, 
      clustering_distance_rows = "pearson", 
      clustering_method_rows = "average", 
      show_row_dend = FALSE, 
      show_row_names = length(drug_selection) < 20, 
      cluster_columns = TRUE, 
      clustering_distance_columns = "pearson", 
      clustering_method_columns = "average", 
      show_column_dend = FALSE, 
      show_column_names = FALSE
    ),
    png, 
    paste0(best_dr_res_path, "patient_predicted_aac_heatmap.png"), 
    width = plot_width, 
    height = plot_height, 
    res = plot_res, 
    units = plot_units
  )
  if (tcga_brca || pan_cancer) {
    # TCGA BRCA annotations
    brca_annotations <- read.csv(
      paste0(
        patient_expression_root_dir, 
        "../../../tcga_extra/brca_labels_extended.csv"), 
      header = TRUE, 
      row.names = 1)
    brca_annotations <- as.list(brca_annotations)
    brca_annotations[["ER"]] <- ifelse(brca_annotations[["ER.Status"]] < 0, "-", "+")
    brca_annotations[["PR"]] <- ifelse(brca_annotations[["PR.Status"]] < 0, "-", "+")
    brca_annotations[["HER2"]] <- ifelse(brca_annotations[["HER2.Final.Status"]] < 0, "-", "+")
    brca_annotations[["age group"]] <- cut(
      brca_annotations[["age_at_initial_pathologic_diagnosis"]], 
      breaks = c(20, 40, 50, 60, 70, 90))
    brca_annotations[["pathologic stage"]] <- brca_annotations[["pathologic_stage"]] 
    brca_annotations[["subtype"]] <- brca_annotations[["BRCA_Subtype_PAM50"]] 
    if (FALSE) {
      lapply(
        brca_annotations[
          c("age group", "pathologic stage", "subtype")
        ], 
        table
      )
    }
    brca_annotations_matched_ind <- match(
      substring(rownames(prediction_mat_patient), 1, 12), 
      brca_annotations[["patient"]]
    )
    brca_annotations_matched <- lapply(
      brca_annotations[
        c("age group", "pathologic stage", "subtype", 
          "ER", "PR", "HER2")
      ], 
      function(x) x[brca_annotations_matched_ind]
    )
    brca_annotations_matched <- lapply(brca_annotations_matched, function(x) {x[x == "NA"] <- NA; return(x)})
    #if (all(sapply(brca_annotations_matched, function(x) all(is.na(x))))) brca_annotations_matched <- list()
  }
  if (scanb) {
    # SCANB
    scanb_y_fn <- paste0(scanb_path, "scanb_pheno.csv.gz")
    scanb_y <- read.csv(scanb_y_fn, row.names = 1, header = TRUE)
    scanb_match <- match(rownames(prediction_mat_patient), scanb_y[["GEX.assay"]])
    scanb_col_map <- c(
      Age_group = "age group", 
      NHG = "grade", 
      ER = "ER", 
      PR = "PR", 
      HER2 = "HER2", 
      NCN.PAM50 = "subtype"
    )
    brca_annotations_matched <- scanb_y[scanb_match, names(scanb_col_map)]
    colnames(brca_annotations_matched) <- scanb_col_map
    receptor_class_map <- c(Negative = "-", Positive = "+")
    for (receptor in c("ER", "PR", "HER2")) {
      brca_annotations_matched[[receptor]] <- receptor_class_map[brca_annotations_matched[[receptor]]]
    }
    brca_annotations_matched_ind <- scanb_match
  }
  brca_ind <- !is.na(brca_annotations_matched_ind)
  
  # Continue BRCA specific analysis
  brca_annotations_matched <- lapply(brca_annotations_matched, function(x) x[brca_ind])
  brca_annotations_matched <- lapply(brca_annotations_matched, factor)
  brca_annotations_matched <- lapply(brca_annotations_matched, droplevels)
  survival_risk <- scale(internal_predictions_patient[["survival_risk_0"]])
  if (FALSE) {
    survival_binned <- cut(
      survival_risk, 
      breaks = c(-Inf, quantile(
        survival_risk, 
        probs = c(0.2, 0.4, 0.6, 0.8), 
        na.rm = TRUE), Inf))
  } else {
    survival_binned <- cut(
      survival_risk, 
      breaks = c(-Inf, qnorm(c(0.2, 0.4, 0.6, 0.8)), Inf))
  }
  brca_annotations_matched[["triple negative"]] <- ifelse(
    brca_annotations_matched[["ER"]] == "-" & 
      brca_annotations_matched[["PR"]] == "-" & 
      brca_annotations_matched[["HER2"]] == "-", 
    "yes", 
    "no"
  )
  brca_annotations_matched[["survival risk z bin"]] <- survival_binned[brca_ind]
  annot_nbins <- sapply(brca_annotations_matched, function(x) length(table(x)))
  #if(is.na(annot_nbins[7])) annot_nbins[7] <- annot_nbins[4]
  if (scanb) {
    brca_second_annotation = "grade"
  } else {
    brca_second_annotation = "pathologic stage"
  }
  
  brca_heatmap_annots <- COPS::heatmap_annotations(
    brca_annotations_matched, 
    factor_color_sets = list(
      rev(pals::kovesi.rainbow(annot_nbins["age group"])), 
      rev(RColorBrewer::brewer.pal(annot_nbins[brca_second_annotation], "PiYG")), 
      c("#67A9CF", "#EF8A62"), 
      c("#67A9CF", "#EF8A62"), 
      c("#67A9CF", "#EF8A62"), 
      c("#67A9CF", "#EF8A62"), 
      RColorBrewer::brewer.pal(annot_nbins["subtype"], "Dark2"), 
      rev(RColorBrewer::brewer.pal(annot_nbins["survival risk z bin"], "RdBu"))
    ), 
    annotation_legend_param = list(nrow = 6)
  )
  
  max_drug_aac <- apply(t(prediction_mat_patient)[brca_drug_correlation_order,], 1, max)
  min_drug_aac <- apply(t(prediction_mat_patient)[brca_drug_correlation_order,], 1, min)
  minmax_drug_ind <- max_drug_aac > 0.4
  q90_drug_aac <- apply(t(prediction_mat_patient)[brca_drug_correlation_order,], 1, quantile, probs = 0.9)
  aac_lims <- c(0,0.5,0.8) * max(max_drug_aac)
  
  #mean(ctrp_true_cor[brca_drug_correlation_order[minmax_drug_ind], brca_drug_correlation_order[minmax_drug_ind]], na.rm = TRUE)
  #mean(ctrp_brca_true_cor[brca_drug_correlation_order[minmax_drug_ind], brca_drug_correlation_order[minmax_drug_ind]], na.rm = TRUE)
  #mean(patient_pred_cor[brca_drug_correlation_order[minmax_drug_ind], brca_drug_correlation_order[minmax_drug_ind]], na.rm = TRUE)
  #mean(cl_pred_cor[brca_drug_correlation_order[minmax_drug_ind], brca_drug_correlation_order[minmax_drug_ind]], na.rm = TRUE)
  
  #with(brca_annotations_matched, table(triple_negative, subtype))
  
  brca_heatmap <- ComplexHeatmap::Heatmap(
    t(prediction_mat_patient)[brca_drug_correlation_order, brca_ind][minmax_drug_ind,], 
    name = "predicted AAC", 
    cluster_rows = FALSE, 
    cluster_columns = TRUE, 
    clustering_distance_columns = "pearson", 
    clustering_method_columns = "average", 
    show_column_names = FALSE, 
    show_column_dend = FALSE, 
    col = circlize::colorRamp2(aac_lims, c("blue", "white", "red"))
  )
  
  brca_heatmap_list <- list(brca_heatmap_annots, brca_heatmap)
  
  brca_heatmap_df <- brca_annotations_matched
  brca_heatmap_df[["surv_risk"]] <- internal_predictions_patient[brca_ind, "survival_risk_0"]
  brca_heatmap_df[["surv_risk_z_bin"]] <- brca_heatmap_df[["surv_risk_z"]]
  brca_heatmap_df[["surv_risk_z"]] <- survival_risk[brca_ind,1]
  brca_heatmap_df <- tibble::as_tibble(brca_heatmap_df)
  brca_sens_df <- tibble::as_tibble(
    prediction_mat_patient[brca_ind,], 
    .name_repair = "minimal"
  )
  brca_rownames_df <- tibble::tibble(ID = rownames(prediction_mat_patient)[brca_ind])
  brca_heatmap_df <- tibble::add_column(brca_rownames_df, brca_heatmap_df)
  brca_heatmap_df_all <- tibble::add_column(
    brca_heatmap_df, 
    brca_sens_df, 
    .name_repair = "minimal"
  )
  if (!all(rownames(prediction_mat_patient) == final_embeddings_p[["X"]])) {
    stop("Invalid embedding indexing")
  }
  brca_embeddings_df <- tibble::as_tibble(final_embeddings_p[brca_ind,shared_embedding_names])
  #selected_drugs <- intersect(brca_drugs, high_confidence_drugs)
  #brca_heatmap_df <- tibble::add_column(brca_heatmap_df, brca_embeddings_df)
  #brca_heatmap_df <- tibble::add_column(brca_heatmap_df, brca_sens_df[,selected_drugs])
  brca_heatmap_df_all <- tibble::add_column(
    brca_heatmap_df_all, 
    brca_embeddings_df, 
    .name_repair = "minimal"
  )
  fn <- paste0(
    best_dr_res_path, 
    "annotated_brca_patient_predicted_aac.csv.gz"
  )
  readr::write_csv(brca_heatmap_df, file = gzfile(fn))
  
  if (FALSE) {
    spls_tune <- mixOmics::tune.spls(
      X = brca_embeddings_df, 
      Y = brca_sens_df[,selected_drugs], 
      ncomp = 5,
      test.keepX = 1:10,
      test.keepY = 1:4,
      nrepeat = 2, 
      folds = 10,
      mode = 'regression', 
      measure = 'cor'
    )
    plot(spls_tune)
    spls_final <- mixOmics::spls(
      X = brca_embeddings_df, 
      Y = brca_sens_df[,selected_drugs], 
      ncomp = 5, 
      keepX = spls_tune$choice.keepX,
      keepY = spls_tune$choice.keepY,
      mode = "regression"
    )
    pmo <- mixOmics::plotIndiv(
      spls_final, 
      ind.names = FALSE, 
      rep.space = "XY-variate", 
      group = brca_heatmap_df$subtype, # colour by time group
      pch = as.factor(brca_heatmap_df$surv_risk_z_bin), # select symbol
      col.per.group = mixOmics::color.mixo(1:6),                      # by dose group
      legend = TRUE, 
      legend.title = 'Subtype', 
      legend.title.pch = 'Risk'
    )
  }
  
  save_figure_safe(
    Reduce(ComplexHeatmap::`%v%`, brca_heatmap_list), 
    png, 
    paste0(best_dr_res_path, "brca_known_drug_patient_predicted_aac_heatmap_decorated.png"), 
    width = plot_width * 1.5, 
    height = plot_height * 1.3, 
    res = plot_res, 
    units = plot_units
  )
  
  prediction_mat_cl <- as.matrix(internal_predictions_cl[,grep("dr_pred_[0-9]+", colnames(internal_predictions_cl))])
  colnames(prediction_mat_cl) <- ctrp_drugs
  rownames(prediction_mat_cl) <- internal_predictions_cl[["X"]]
  brca_ind <- cell_line_brca_map[rownames(prediction_mat_cl)]
  brca_ind[is.na(brca_ind)] <- FALSE
  
  save_figure_safe(
    ComplexHeatmap::Heatmap(
      t(prediction_mat_cl)[brca_drug_correlation_order, ], 
      name = "predicted AAC", 
      cluster_rows = TRUE, 
      clustering_distance_rows = "pearson", 
      clustering_method_rows = "average", 
      show_row_dend = FALSE, 
      show_row_names = TRUE, 
      cluster_columns = TRUE, 
      clustering_distance_columns = "pearson", 
      clustering_method_columns = "average", 
      show_column_dend = FALSE, 
      show_column_names = FALSE
    ),
    png, 
    paste0(best_dr_res_path, "brca_known_drug_cl_predicted_aac_heatmap.png"), 
    width = plot_width, 
    height = plot_height, 
    res = plot_res, 
    units = plot_units
  )
  
  temp <- plyr::join(ext_dr_results, dr_r2_best, type = "left")
  plot1 <- ggplot(temp, aes(x = cl_test, y = SpearmanR, color = dataset)) + 
    geom_point(shape = 3, size = 3) + 
    theme_bw() + xlab("test R^2")
  
  save_figure_safe(
    plot1, 
    png, 
    paste0(best_dr_res_path, "external_validation_drug_performance.png"), 
    width = plot_height, 
    height = plot_height, 
    res = plot_res, 
    units = plot_units
  )
  
  # One table
  drug_performance <- plyr::join(dr_r2_best, ext_dr_results, type = "left")
  drug_performance <- drug_performance[,c("drug", "cl_train", "cl_test", "PearsonR", "SpearmanR", "dataset")]
  colnames(drug_performance) <- c("drug", "CTRP_train_R2", "CTRP_test_R2", "ext_PearsonCor", "ext_SpearmanCor", "dataset")
  
  write.csv(drug_performance, paste0(plot_path, "external_drug_response_performance.csv"), row.names = FALSE)
  
  temp <- drug_performance#[drug_performance[["drug"]] %in% drugs_of_interest,]
  plot1 <- ggplot(temp, aes(x = CTRP_test_R2, y = ext_SpearmanCor, label = drug, color = dataset)) + 
    geom_point(shape = "+", size = 3) + 
    ggrepel::geom_text_repel() +
    theme_bw() + 
    xlim(c(-0.5, 0.4))
  save_figure_safe(
    plot1, 
    png, 
    paste0(plot_path, "external_validation_drug_performance_text.png"), 
    width = plot_width, 
    height = 0.7 * plot_width, 
    res = plot_res, 
    units = plot_units
  )
  
  # Final predictions for BRCA cell-lines
  # Use to identify interesting drugs?
  if (FALSE) {
    best_cl_predictions <- internal_predictions_cl
    
    # TODO: remove duplicated
    ctrp_dr <- read.csv(
      paste0(ctrp_path, "screen_aac.csv.gz"), 
      row.names = 1, header = TRUE)
    ctrp_dr_rowinfo <- read.csv(
      paste0(ctrp_path, "screen_rowinfo.csv.gz"), 
      row.names = 1, header = TRUE)
    ctrp_dr_row_cl_id <- read.csv(
      paste0(ctrp_path, "ctrp_ccle_clid_map.csv.gz"), 
      row.names = 1, header = TRUE)
    
    ccle_brca_ind <- cell_line_brca_map[rownames(best_cl_predictions)]
    ccle_brca_ind[is.na(ccle_brca_ind)] <- FALSE
    best_brca_cl_ccle_pred <- best_cl_predictions[ccle_brca_ind,]
    
    dr_pred_drug_names <- readLines(paste0(ctrp_path, "drug_names.txt"))
    dr_pred_drug_ind <- as.numeric(gsub("dr_pred_", "", colnames(best_cl_predictions)))
    dr_pred_drug_names_matched <- dr_pred_drug_names[dr_pred_drug_ind + 1]
    
    ctrp_dr_drugwise_z <- plyr::ddply(
      cbind(ctrp_dr, ctrp_dr_rowinfo), c("treatmentid_fixed"), 
      function(x) data.frame(sample_id = x[["sampleid"]], z = scale(x[["aac_recomputed"]])))
    
    ctrp_dr_drugwise_z[["ccle_id"]] <- ctrp_dr_row_cl_id[ctrp_dr_drugwise_z[["sample_id"]], "CCLE_model_id"]
    ctrp_dr_drugwise_z[["brca"]] <- cell_line_brca_map[ctrp_dr_drugwise_z[["ccle_id"]]]
    ctrp_dr_drugwise_z[["brca"]][is.na(ctrp_dr_drugwise_z[["brca"]])] <- FALSE
    
    ctrp_dr_drugwise_brca_z_mean <- plyr::ddply(
      ctrp_dr_drugwise_z[ctrp_dr_drugwise_z[["brca"]], ], c("treatmentid_fixed"), 
      function(x) data.frame(z_mean = mean(x$z, na.rm = TRUE)))
    brca_specific_ind <- which(ctrp_dr_drugwise_brca_z_mean$z_mean > 0.5)
    ctrp_dr_drugwise_brca_z_mean[brca_specific_ind, "treatmentid_fixed"]
    
    # How to identify best drugs predicted for BRCA?
    # How to quantify agreement with BRCA CL sensitivities in CTRP?
    # * Much easier if we have pan-cancer predictions
  }
}