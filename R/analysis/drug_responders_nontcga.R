source("../setup.R")

for (save_dir in save_dirs) {
  if(exists("var_list")) try(detach(var_list))
  var_list <- get_paths_and_parameters(save_dir, remove_incomplete = FALSE)
  attach(var_list)
  best_task <- readLines(paste0(save_path, "best_task.txt"))
  param_best_task_ind <- match(best_task, parameters[["task"]])
  shared_embedding_names <- paste0("z", 1:parameters[param_best_task_ind, "bottle_neck"])
  
  fn <- paste0(save_path, "external_evaluation/internal_survival_validation_final_embeddings.csv.gz")
  tcga_embeddings <- readr::read_csv(fn, show_col_types = FALSE)
  colnames(tcga_embeddings)[1] <- "id"
  
  fn <- paste0(patient_expression_root_dir, "oncotree_level1.csv")
  patient_ot_level1_labels <- read.csv(fn, row.names = 1, header = TRUE)
  patient_ot_level1_labels[["level_1"]] <- tolower(gsub("_", " ", patient_ot_level1_labels[["level_1"]]))
  tcga_tissue_label <- rep_len(NA, nrow(tcga_embeddings))
  names(tcga_tissue_label) <- tcga_embeddings[["id"]]
  
  tcga_tissue_label[rownames(patient_ot_level1_labels)] <- patient_ot_level1_labels[["level_1"]]
  tcga_nna_ind <- !is.na(tcga_tissue_label)
  
  class_map_df <- readr::read_csv(paste0(save_path, "external_evaluation/class_map.csv"))
  class_map <- class_map_df[["name"]]
  names(class_map) <- paste0("class_pred_", class_map_df[["key"]])
  
  classification_list <- list()
  
  for (study_id in names(ctrdb_datasets)) {
    fn <- paste0(ctrdb_path, study_id, "_response.csv")
    ctrdb_response <- readr::read_csv(fn, show_col_types = FALSE)
    
    prediction_path <- paste0(save_path, "external_evaluation/", study_id, "/")
    figure_str <- ""
    
    fn <- paste0(prediction_path, "predictions.csv.gz")
    predictions <- readr::read_csv(fn, show_col_types = FALSE)
    
    fn <- paste0(prediction_path, "embeddings.csv.gz")
    embeddings <- readr::read_csv(fn, show_col_types = FALSE)
    colnames(embeddings)[1] <- "id"
    
    tissue_label <- rep_len(NA, nrow(embeddings))
    names(tissue_label) <- embeddings[["id"]]
    
    response_patient_z_ptr <- match(embeddings[["id"]], ctrdb_response[["Sample_id"]])
    tissue_label <- ctrdb_response[["Cancer_type_level1"]][response_patient_z_ptr]
    tissue_label <- ctrdb_types_to_oncotree[tissue_label]
    
    n_labels <- length(unique(c(tcga_tissue_label, tissue_label)))
    
    dataset <- factor(
      rep(
        c("TCGA", study_id), 
        c(nrow(tcga_embeddings), nrow(embeddings))
      ), 
      levels = c("TCGA", study_id)
    )
    
    set.seed(0)
    tissue_plot <- tissue_visualizer(
      dplyr::bind_rows(
        tcga_embeddings[,shared_embedding_names], 
        embeddings[,shared_embedding_names]
      ),
      labeled_ind = c(tcga_nna_ind, !is.na(tissue_label)), 
      color_var = c(tcga_tissue_label, tissue_label), 
      color_name = "tissue",
      shape_var = dataset, 
      shape_name = "dataset",
      reference_shape_label = "TCGA", 
      umap_neighbors = 20, 
      pre_manifold_pca = FALSE, 
      scale_datasets_separately = FALSE, 
      umap_args = list(init = "normlaplacian", min_dist = 0.5),
      primary_color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=n_labels)), 
      primary_fill_scale = scale_fill_manual(values = pals::kovesi.rainbow(n=n_labels)), 
      secondary_color_scale = scale_color_manual(values = c(NA, "#000000")), 
      primary_shape_scale = scale_shape_manual(values = c(3,NA)),
      secondary_shape_scale = scale_shape_manual(values = c(NA,21)), 
      point_size = 3, 
      point_alpha = 0.5, 
      annotation_size = 6, 
      annotation_force = 10, 
      knn_accuracy = TRUE,
      knn_predicted_dataset = study_id, 
      knn_title_method_name = "MODAE"
    )
    source_cancer <- paste(
      tolower(unique(ctrdb_response[["Cancer_type_level1"]])), 
      collapse = ","
    )
    canceri_path <- paste0(plot_path, "treatment_response/", source_cancer)
    dir.create(canceri_path, recursive = TRUE)
    save_figure_safe(
      with(tissue_plot, tissue_plot), 
      png, 
      paste0(canceri_path, "/", study_id, "_embeddings.png"), 
      width = plot_width, 
      height = plot_width, 
      res = plot_res, 
      units = plot_units
    )
    
    class_predictions <- predictions[, grep("^class_pred_[0-9]+$", colnames(predictions))]
    class_predictions_order <- t(apply(-class_predictions, 1, order))
    
    ctrp_drugs <- get_ctrp_drugs(n_drugs = length(grep("^dr_pred_[0-9]+$", colnames(predictions))))
    ctrp_drug_processed <- process_ctrp_drug_names(ctrp_drugs)
    
    pred_patient_id <- predictions[[1]]
    response_patient_pred_ptr <- match(ctrdb_response[["Sample_id"]], pred_patient_id)
    
    classification_list[[study_id]] <- data.frame(
      ctrdb = ctrdb_response[["Cancer_type_level1"]], 
      prediction = class_map[colnames(class_predictions)[class_predictions_order[response_patient_pred_ptr,1]]]
    )
    
    process_ctrdb_drug_list <- function(x) {
      x <- tolower(x)
      x <- strsplit(x, split = "/")
      x <- lapply(x, function(x) gsub(".*\\(|\\)$", "", x))
      x <- sapply(x, paste, collapse = "+")
      return(x)
    }
    drug_comb <- process_ctrdb_drug_list(ctrdb_response[["Drug_list"]])
    drug_comb_match <- match(drug_comb, ctrp_drug_processed)
    match_ind <- !is.na(response_patient_pred_ptr) & !is.na(drug_comb_match)
    direct_match_table <- with(
      ctrdb_response[match_ind,], 
      table(
        Source, 
        Drug_ptr = drug_comb_match[match_ind], 
        Response
      )
    )
    direct_match_table_total <- apply(direct_match_table, c(1,2), sum)
    direct_match_table_min <- apply(direct_match_table, c(1,2), min)
    
    direct_match_groups <- which(direct_match_table_min > 2 & direct_match_table_total >= 4, arr.ind = TRUE)
    if (nrow(direct_match_groups)>0) for (i in 1:nrow(direct_match_groups)) {
      match_table_drug_ptr <- as.integer(colnames(direct_match_table_total)[direct_match_groups[i,"Drug_ptr"]])
      drugi_name <- ctrp_drug_processed[match_table_drug_ptr]
      match_table_drug_ptr <- paste0("dr_pred_", match_table_drug_ptr-1)
      match_table_source <- rownames(direct_match_table_total)[direct_match_groups[i,"Source"]]
      response_source_ind <- ctrdb_response[["Source"]] == match_table_source & drug_comb == drugi_name
      response_source_ind <- response_source_ind & !is.na(response_patient_pred_ptr)
      matched_patients_ptr <- response_patient_pred_ptr[response_source_ind]
      predsi <- predictions[matched_patients_ptr, match_table_drug_ptr][[1]]
      canceri <- paste0(match_table_source, " (", source_cancer, ")")
      
      plot_df <- data.frame(
        drug_response = ctrdb_response[["Response"]][response_source_ind], 
        predicted_sensitivity = predsi
      )
      plot_df[["drug_response"]] <- factor(
        plot_df[["drug_response"]], 
        levels = c("Non_response", "Response")
      )
      p1 <- ggpubr::ggboxplot(
        plot_df, 
        x = "drug_response", 
        y = "predicted_sensitivity", 
        fill = "drug_response"
      ) + 
        theme_bw() + 
        ggpubr::stat_compare_means() + 
        ggtitle(paste(drugi_name, "response in\n", canceri))
      
      canceri_path <- paste0(plot_path, "treatment_response/", source_cancer)
      dir.create(canceri_path, recursive = TRUE)
      save_figure_safe(
        p1, 
        png, 
        paste0(canceri_path, "/", match_table_source, "_", drugi_name, "_direct_comparison.png"), 
        width = plot_width / 2, 
        height = plot_width / 2, 
        res = plot_res, 
        units = plot_units
      )
    }
    
    # Combined matches
    drug_lists <- strsplit(drug_comb, split = "\\+")
    drug_lists <- lapply(drug_lists, tolower)
    drug_counts <- sapply(drug_lists, length)
    response_drugs_in_ctrp_ind <- sapply(
      drug_lists, 
      function(x) all(x %in% ctrp_drug_processed)
    )
    comb_match_ind <- !is.na(response_patient_pred_ptr) & response_drugs_in_ctrp_ind & drug_counts > 1
    comb_match_table <- with(
      ctrdb_response[comb_match_ind,], 
      table(
        Source, 
        Drug_list, 
        Response
      )
    )
    comb_match_table_total <- apply(comb_match_table, c(1,2), sum)
    comb_match_table_min <- apply(comb_match_table, c(1,2), min)
    
    comb_match_groups <- which(comb_match_table_min > 2 & comb_match_table_total >= 4, arr.ind = TRUE)
    
    if (nrow(comb_match_groups)>0) for (i in 1:nrow(comb_match_groups)) {
      drug_listi <- colnames(comb_match_table_total)[comb_match_groups[i,"Drug_list"]]
      drug_listi_separated <- strsplit(process_ctrdb_drug_list(drug_listi), split = "\\+")[[1]]
      drug_ptrs <- paste0("dr_pred_", match(drug_listi_separated, ctrp_drug_processed) - 1)
      source_name <- rownames(comb_match_table_total)[comb_match_groups[i,"Source"]]
      
      response_source_ind <- ctrdb_response[["Source"]] == source_name & ctrdb_response[["Drug_list"]] == drug_listi
      response_source_ind <- response_source_ind & !is.na(response_patient_pred_ptr)
      matched_patients_ptr <- response_patient_pred_ptr[response_source_ind]
      predsi <- predictions[matched_patients_ptr, drug_ptrs]
      colnames(predsi) <- drug_listi_separated
      canceri <- paste0(source_name, " (", source_cancer, ")")
      plot_df <- data.frame(
        drug_response = ctrdb_response[["Response"]][response_source_ind], 
        predsi
      )
      plot_df[["drug_response"]] <- factor(
        plot_df[["drug_response"]], 
        levels = c("Non_response", "Response")
      )
      p1 <- GGally::ggpairs(
        plot_df, 
        aes(color = drug_response, shape = drug_response), 
        columns = 2:ncol(plot_df),
        legend = c(ncol(predsi), 1), 
        diag = list(continuous = GGally::wrap("densityDiag", alpha = 0.5)),
        upper = list(continuous = "blank"), 
        lower = list(continuous = GGally::wrap("points", size = 3))) + 
        scale_color_brewer(palette = "Dark2") + 
        scale_fill_brewer(palette = "Dark2") + 
        theme_bw() + 
        ggtitle(paste(drug_listi, "sensitivity in\n", canceri))
      canceri_path <- paste0(plot_path, "treatment_response/", source_cancer)
      dir.create(canceri_path, recursive = TRUE)
      save_figure_safe(
        p1, 
        png, 
        paste0(canceri_path, "/", source_name, "_", drug_listi, "_scatter.png"), 
        width = plot_width, 
        height = plot_width / 2, 
        res = plot_res, 
        units = plot_units
      )
    }
  }
  
  unique_ctrdb_types <- Reduce(union, lapply(classification_list, function(x) unique(x[["ctrdb"]])))
  
  classification_results <- list()
  classification_accuracy <- c()
  # Classification accuracy
  for (i in names(classification_list)) {
    y_pred <- classification_list[[i]][["prediction"]]
    y_true <- classification_list[[i]][["ctrdb"]]
    y_true <- ctrdb_types_to_oncotree[y_true]
    classification_accuracy[i] <- mean(y_pred == y_true, na.rm = TRUE)
    classification_results[[i]] <- list(
      predicted = unique(y_pred), 
      actual = unique(y_true), 
      accuracy = classification_accuracy[i]
    )
  }
  classification_accuracy["GSE109211"]
  classification_accuracy["GSE206501"]
}