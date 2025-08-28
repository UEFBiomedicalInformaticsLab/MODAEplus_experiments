source("../setup.R")

codeae_comparison <- TRUE

for (save_dir in save_dirs) {
  if(exists("var_list")) try(detach(var_list))
  var_list <- get_paths_and_parameters(save_dir, remove_incomplete = FALSE)
  attach(var_list)
  best_task <- readLines(paste0(save_path, "best_task.txt"))
  param_best_task_ind <- match(best_task, parameters[["task"]])
  shared_embedding_names <- paste0("z", 1:parameters[param_best_task_ind, "bottle_neck"])
  
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
    
    tcga_surv_data <- read.csv(
      paste0(patient_expression_root_dir, "survival.csv.gz"), 
      header = TRUE, row.names = 1)
    
    internal_match_tcga <- match(
      substr(rownames(internal_predictions), 1, 12), 
      tcga_surv_data[["bcr_patient_barcode"]])
    
    if (!all(is.na(internal_match_tcga))) {
      tcga_predictions <- internal_predictions
    } else if (all(!is.na(external_match_tcga))) {
      tcga_predictions <- external_predictions
    } else {
      stop("No predictions rows match all TCGA ids.")
    }
    
    tcga_surv_ind <- match(
      substr(rownames(tcga_predictions), 1, 12), 
      tcga_surv_data[["bcr_patient_barcode"]])
    
    tcga_surv_type <- tcga_surv_data[tcga_surv_ind, "type"]
    tcga_surv_type[is.na(tcga_surv_type)] <- "NA"
    #hist(tcga_predictions[["survival_risk_0"]])
    #hist(scanb_predictions[["survival_risk_0"]])
    
    dr_dir <- paste0(
      base_dir, 
      "../CODE-AE-v1.0/data/tcga/"
    )
    first_treatment <- TRUE
    if (first_treatment) {
      tcga_new_tumour_table_raw <- read.csv(
        paste0(dr_dir, "tcga_drug_first_response.csv"), 
        header = TRUE
      )
      tcga_dr_table_raw <- read.csv(
        paste0(dr_dir, "tcga_drug_first_treatment.csv"), 
        header = TRUE
      )
      figure_str <- paste0(figure_str, "_first")
    } else {
      tcga_new_tumour_table_raw <- read.csv(
        paste0(dr_dir, "tcga_drug_response.csv"), 
        header = TRUE
      )
      tcga_dr_table_raw <- read.csv(
        paste0(dr_dir, "tcga_drug_treatment.csv"), 
        header = TRUE
      )
    }
    
    tcga_dr_joined <- dplyr::inner_join(
      tcga_dr_table_raw, 
      tcga_new_tumour_table_raw, 
      by = c("bcr_patient_barcode", "tcga_project")
    )
    tcga_dr_joined <- tcga_dr_joined[
      tcga_dr_joined[["bcr_patient_barcode"]] %in% substr(rownames(internal_predictions), 1, 12), 
    ]
    tcga_dr_labeled <- plyr::ddply(
      tcga_dr_joined, 
      c("pharmaceutical_therapy_drug_name"), 
      function(x) with(x, data.frame(
        bcr_patient_barcode = bcr_patient_barcode, 
        responder = days_to_new_tumor_event_after_initial_treatment > median(days_to_new_tumor_event_after_initial_treatment)
      ))
    )
    
    codeae_pred_path <- paste0(
      codeae_path, 
      "XieResearchGroup-CODE-AE-6dc17a5/", 
      "intermediate_results/tcga_prediction/"
    )
    fn <- paste0(codeae_pred_path, "tcga_scores.csv")
    codeae_predictions <- readr::read_csv(fn, show_col_types = FALSE)
    codeae_drug_cols <- colnames(codeae_predictions)[-1:-2]
    
    codeae_drugs_in_tcga <- c(intersect(
      codeae_drug_cols, 
      unique(tcga_dr_labeled[["pharmaceutical_therapy_drug_name"]])
    ), "Sorafenib")
    codae_auc <- c()
    for (drugi in codeae_drugs_in_tcga) {
      drugi_ind <- tcga_dr_labeled[["pharmaceutical_therapy_drug_name"]] == drugi
      tcga_labeli <- tcga_dr_labeled[drugi_ind, c("bcr_patient_barcode", "responder")]
      tcga_predi_ptr <- match(
        tcga_labeli[["bcr_patient_barcode"]], 
        substr(codeae_predictions[["Sample"]], 1, 12)
      )
      if (drugi == "Sorafenib") {
        for (drugj in c("Sorafenib...20", "Sorafenib...51")) {
          tcga_labeli[["pred"]] <- codeae_predictions[[drugj]][tcga_predi_ptr]
          codae_auc[drugj] <- with(tcga_labeli, pROC::auc(pROC::roc(
            responder, 
            pred, 
            direction = "<"
          )))
        }
        next
      }
      tcga_labeli[["pred"]] <- codeae_predictions[[drugi]][tcga_predi_ptr]
      codae_auc[drugi] <- with(tcga_labeli, pROC::auc(pROC::roc(
        responder, 
        pred, 
        direction = "<"
      )))
    }
  }
}