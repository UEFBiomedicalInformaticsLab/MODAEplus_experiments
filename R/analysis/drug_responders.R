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
        Drug_ptr = tcga_drug_comb_match[match_ind], 
        Response
      )
    )
    direct_match_table_total <- apply(direct_match_table, c(1), sum)
    direct_match_table_min <- apply(direct_match_table, c(1), min)
    
    direct_auroc <- c()
    direct_match_groups <- which(direct_match_table_min > 2 & direct_match_table_total >= 10)
    if (length(direct_match_groups)>0) for (i in 1:length(direct_match_groups)) {
      match_table_drug_ptr <- as.integer(names(direct_match_table_total)[direct_match_groups[i]])
      drugi_name <- ctrp_drug_processed[match_table_drug_ptr]
      match_table_drug_ptr <- paste0("dr_pred_", match_table_drug_ptr-1)
      #match_table_source <- rownames(direct_match_table_total)[direct_match_groups[i,"Source"]]
      tcga_response_ind <- tcga_drug_comb == drugi_name
      tcga_response_ind <- tcga_response_ind & !is.na(tcga_response_patient_pred_ptr)
      matched_patients_ptr <- tcga_response_patient_pred_ptr[tcga_response_ind]
      predsi <- internal_predictions[matched_patients_ptr, match_table_drug_ptr]
      #canceri <- gsub("TCGA_", "", match_table_source)
      
      plot_df <- data.frame(
        drug_response = tcga_response[["Response"]][tcga_response_ind], 
        predicted_sensitivity = predsi
      )
      plot_df[["drug_response"]] <- factor(
        plot_df[["drug_response"]], 
        levels = c("Non_response", "Response")
      )
      try(direct_auroc[drugi_name] <- with(
        plot_df, 
        pROC::auc(pROC::roc(
          drug_response, 
          predicted_sensitivity, 
          levels = levels(drug_response), 
          direction = "<"
        ))
      ))
      p1 <- ggpubr::ggboxplot(
        plot_df, 
        x = "drug_response", 
        y = "predicted_sensitivity", 
        fill = "drug_response"
      ) + 
        theme_bw() + 
        ggbeeswarm::geom_beeswarm() +
        ggpubr::stat_compare_means() + 
        ggtitle(paste(drugi_name, "response \nAUROC:", direct_auroc[drugi_name]))
      
      treat_path <- paste0(plot_path, "treatment_response/")
      dir.create(treat_path, recursive = TRUE)
      save_figure_safe(
        p1, 
        png, 
        paste0(treat_path, drugi_name, "_direct_comparison.png"), 
        width = plot_width / 2, 
        height = plot_width / 2, 
        res = plot_res, 
        units = plot_units
      )
    }
    
    # Partial matches
    tcga_drug_lists <- strsplit(tcga_response[["Drug_list"]], split = "\\+")
    tcga_drug_lists <- lapply(tcga_drug_lists, tolower)
    tcga_drug_lists <- lapply(tcga_drug_lists, only_alphanumericals)
    tcga_drug_unique <- Reduce(union, tcga_drug_lists)
    tcga_drug_unique <- intersect(
      tcga_drug_unique, 
      only_alphanumericals(ctrp_drug_processed)
    )
    partial_auroc <- c()
    for (drugi in tcga_drug_unique) {
      drugi_idx <- paste0("dr_pred_", which(only_alphanumericals(ctrp_drug_processed) == drugi)-1)
      tcga_response_ind <- sapply(tcga_drug_lists, function(x) drugi %in% x)
      tcga_response_ind <- tcga_response_ind & !is.na(tcga_response_patient_pred_ptr)
      matched_patients_ptr <- tcga_response_patient_pred_ptr[tcga_response_ind]
      predsi <- internal_predictions[matched_patients_ptr, drugi_idx]
      
      plot_df <- data.frame(
        drug_response = tcga_response[["Response"]][tcga_response_ind], 
        predicted_sensitivity = predsi
      )
      plot_df[["drug_response"]] <- factor(
        plot_df[["drug_response"]], 
        levels = c("Non_response", "Response")
      )
      
      try(partial_auroc[drugi] <- with(
        plot_df, 
        pROC::auc(pROC::roc(
          drug_response, 
          predicted_sensitivity, 
          levels = levels(drug_response), 
          direction = "<"
        ))
      ))
      
      p1 <- ggpubr::ggboxplot(
        plot_df, 
        x = "drug_response", 
        y = "predicted_sensitivity", 
        fill = "drug_response"
      ) + 
        theme_bw() + 
        ggbeeswarm::geom_beeswarm() +
        ggpubr::stat_compare_means() + 
        ggtitle(paste(drugi, "response \nAUROC:", partial_auroc[drugi]))
      
      treat_path <- paste0(plot_path, "treatment_response/")
      dir.create(treat_path, recursive = TRUE)
      save_figure_safe(
        p1, 
        png, 
        paste0(treat_path, drugi, "_partial_comparison.png"), 
        width = plot_width / 2, 
        height = plot_width / 2, 
        res = plot_res, 
        units = plot_units
      )
    }
    
    # Combined matches
    if (FALSE) {
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
          Drug_list, 
          Response
        )
      )
      comb_match_table_total <- apply(comb_match_table, c(1), sum)
      comb_match_table_min <- apply(comb_match_table, c(1), min)
      
      comb_match_groups <- which(comb_match_table_min > 2 & comb_match_table_total >= 10)
      
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
}