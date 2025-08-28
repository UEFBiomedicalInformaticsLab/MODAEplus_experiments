source("../setup.R")

for (save_dir in save_dirs) {
  if(exists("var_list")) try(detach(var_list))
  var_list <- get_paths_and_parameters(save_dir, remove_incomplete = FALSE)
  attach(var_list)
  best_task <- readLines(paste0(save_path, "best_task.txt"))
  param_best_task_ind <- match(best_task, parameters[["task"]])
  shared_embedding_names <- paste0("z", 1:parameters[param_best_task_ind, "bottle_neck"])
  
  fn <- paste0(ctrdb_path, "tcga_responses.csv")
  tcga_response <- readr::read_csv(fn)
  
  # Final predictions
  prediction_paths <- paste0(save_path, "external_evaluation/")
  if (dir.exists(paste0(save_path, "reverse_external_evaluation"))) {
    prediction_paths[2] <- paste0(save_path, "reverse_external_evaluation/")
  }
  for (prediction_path in prediction_paths) {
    if (grepl("reverse", prediction_path)) {
      figure_str <- "_reverse"
    } else {
      figure_str <- ""
    }
    
    internal_predictions <- read.csv(
      paste0(prediction_path, "internal_survival_validation_final_predictions.csv.gz"), 
      header = TRUE, row.names = 1)
    ctrp_drugs <- get_ctrp_drugs(n_drugs = length(grep("^dr_pred_[0-9]+$", colnames(internal_predictions))))
    ctrp_drug_processed <- process_ctrp_drug_names(ctrp_drugs)
    
    pred_patient_id <- substr(rownames(internal_predictions), 1, 12)
    tcga_response_patient_pred_ptr <- match(tcga_response[["Sample_id"]], pred_patient_id)
    tcga_drug_comb <- tolower(tcga_response[["Drug_list"]])
    tcga_drug_comb_match <- match(tcga_drug_comb, ctrp_drug_processed)
    match_ind <- !is.na(tcga_response_patient_pred_ptr) & !is.na(tcga_drug_comb_match)
    direct_match_table <- with(
      tcga_response[match_ind,], 
      table(
        Source, 
        Drug_ptr = tcga_drug_comb_match[match_ind], 
        Response
      )
    )
    direct_match_table_total <- apply(direct_match_table, c(1,2), sum)
    direct_match_table_min <- apply(direct_match_table, c(1,2), min)
    
    direct_match_groups <- which(direct_match_table_min > 2 & direct_match_table_total >= 10, arr.ind = TRUE)
    if (nrow(direct_match_groups)>0) for (i in 1:nrow(direct_match_groups)) {
      match_table_drug_ptr <- as.integer(colnames(direct_match_table_total)[direct_match_groups[i,"Drug_ptr"]])
      drugi_name <- ctrp_drug_processed[match_table_drug_ptr]
      match_table_drug_ptr <- paste0("dr_pred_", match_table_drug_ptr-1)
      match_table_source <- rownames(direct_match_table_total)[direct_match_groups[i,"Source"]]
      tcga_response_source_ind <- tcga_response[["Source"]] == match_table_source & tcga_drug_comb == drugi_name
      tcga_response_source_ind <- tcga_response_source_ind & !is.na(tcga_response_patient_pred_ptr)
      matched_patients_ptr <- tcga_response_patient_pred_ptr[tcga_response_source_ind]
      predsi <- internal_predictions[matched_patients_ptr, match_table_drug_ptr]
      canceri <- gsub("TCGA_", "", match_table_source)
      
      plot_df <- data.frame(
        drug_response = tcga_response[["Response"]][tcga_response_source_ind], 
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
        ggtitle(paste(drugi_name, "response in", canceri))
      
      canceri_path <- paste0(plot_path, "treatment_response/", canceri)
      dir.create(canceri_path, recursive = TRUE)
      save_figure_safe(
        p1, 
        png, 
        paste0(canceri_path, "/", canceri, "_", drugi_name, "_direct_comparison.png"), 
        width = plot_width / 2, 
        height = plot_width / 2, 
        res = plot_res, 
        units = plot_units
      )
    }
    
    # Combined matches
    tcga_drug_lists <- strsplit(tcga_response[["Drug_list"]], split = "\\+")
    tcga_drug_lists <- lapply(tcga_drug_lists, tolower)
    tcga_drug_counts <- sapply(tcga_drug_lists, length)
    response_drugs_in_ctrp_ind <- sapply(
      tcga_drug_lists, 
      function(x) all(x %in% ctrp_drug_processed)
    )
    comb_match_ind <- !is.na(tcga_response_patient_pred_ptr) & response_drugs_in_ctrp_ind & tcga_drug_counts > 1
    comb_match_table <- with(
      tcga_response[comb_match_ind,], 
      table(
        Source, 
        Drug_list, 
        Response
      )
    )
    comb_match_table_total <- apply(comb_match_table, c(1,2), sum)
    comb_match_table_min <- apply(comb_match_table, c(1,2), min)
    
    comb_match_groups <- which(comb_match_table_min > 2 & comb_match_table_total >= 10, arr.ind = TRUE)
    
    if (nrow(comb_match_groups)>0) for (i in 1:nrow(comb_match_groups)) {
      drug_listi <- colnames(comb_match_table_total)[comb_match_groups[i,"Drug_list"]]
      drug_listi_separated <- tolower(strsplit(drug_listi, split = "\\+")[[1]])
      drug_ptrs <- paste0("dr_pred_", match(drug_listi_separated, ctrp_drug_processed) - 1)
      source_name <- rownames(comb_match_table_total)[comb_match_groups[i,"Source"]]
      
      tcga_response_source_ind <- tcga_response[["Source"]] == source_name & tcga_response[["Drug_list"]] == drug_listi
      tcga_response_source_ind <- tcga_response_source_ind & !is.na(tcga_response_patient_pred_ptr)
      matched_patients_ptr <- tcga_response_patient_pred_ptr[tcga_response_source_ind]
      predsi <- internal_predictions[matched_patients_ptr, drug_ptrs]
      colnames(predsi) <- drug_listi_separated
      canceri <- gsub("TCGA_", "", source_name)
      plot_df <- data.frame(
        drug_response = tcga_response[["Response"]][tcga_response_source_ind], 
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
        ggtitle(paste(drug_listi, "sensitivity in", canceri))
      canceri_path <- paste0(plot_path, "treatment_response/", canceri)
      dir.create(canceri_path, recursive = TRUE)
      save_figure_safe(
        p1, 
        png, 
        paste0(canceri_path, "/", canceri, "_", drug_listi, "_scatter.png"), 
        width = plot_width, 
        height = plot_width / 2, 
        res = plot_res, 
        units = plot_units
      )
    }
  }
}