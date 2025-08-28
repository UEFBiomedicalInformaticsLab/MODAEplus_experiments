source("../setup.R")
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
    internal_embeddings <- read.csv(
      paste0(prediction_path, "internal_survival_validation_final_embeddings.csv.gz"), 
      header = TRUE, row.names = 1)
    
    if (tcga_brca || scanb) {
      p_ext_fn <- paste0(prediction_path, "external_survival_validation_predictions.csv.gz") 
      if (file.exists(p_ext_fn)) {
        external_predictions <- read.csv(p_ext_fn, header = TRUE, row.names = 1)
      } else {
        external_predictions <- data.frame()
      }
      tcga_surv_data <- read.csv(
        paste0(patient_expression_root_dir, "BRCA/survival.csv.gz"), 
        header = TRUE, row.names = 1)
      scanb_surv_data <- read.csv(
        paste0(scanb_path, "scanb_pheno.csv.gz"), 
        header = TRUE, row.names = 1)
      external_match_tcga <- match(
        substr(rownames(external_predictions), 1, 12), 
        tcga_surv_data[["bcr_patient_barcode"]])
    } else {
      tcga_surv_data <- read.csv(
        paste0(patient_expression_root_dir, "survival.csv.gz"), 
        header = TRUE, row.names = 1)
    }
    
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
    if (tcga_brca || scanb) {
      internal_match_scanb <- match(
        rownames(internal_predictions), 
        scanb_surv_data[["GEX.assay"]])
      external_match_scanb <- match(
        rownames(external_predictions), 
        scanb_surv_data[["GEX.assay"]])
      if (!all(is.na(internal_match_scanb))) {
        scanb_predictions <- internal_predictions
      } else if (all(!is.na(external_match_scanb))) {
        scanb_predictions <- external_predictions
      } else {
        stop("No predictions rows match all SCANB ids.")
      }
      
      tcga_subtypes <- TCGAbiolinks::TCGAquery_subtype("BRCA")
      col_intersect <- intersect(colnames(tcga_surv_data), colnames(tcga_subtypes))
      tcga_surv_data[["patient"]] <- tcga_surv_data[["bcr_patient_barcode"]]
      tcga_surv_data <- plyr::join(
        tcga_surv_data, 
        tcga_subtypes[, !colnames(tcga_subtypes) %in% col_intersect], 
        type = "left")
      tcga_surv_data[["simplified_stage"]] <- gsub(
        "[A-C]+$", "", tcga_surv_data[["ajcc_pathologic_tumor_stage"]])
      
      scanb_surv_ind <- match(
        rownames(scanb_predictions), 
        scanb_surv_data[["GEX.assay"]])
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
    first_treatment <- FALSE
    if (first_treatment) {
      tcga_dr_table_raw <- read.csv(
        paste0(dr_dir, "tcga_drug_first_treatment.csv"), 
        header = TRUE
      )
      figure_str <- paste0(figure_str, "_first")
    } else {
      tcga_dr_table_raw <- read.csv(
        paste0(dr_dir, "tcga_drug_treatment.csv"), 
        header = TRUE
      )
    }
    
    # Omit treatments given prior to biopsy
    tcga_dr_table <- tcga_dr_table_raw[tcga_dr_table_raw[["pharmaceutical_tx_started_days_to"]]>=0,]
    # Fix drug name whitespace
    tcga_dr_table[["pharmaceutical_therapy_drug_name"]] <- gsub("\\s", "", tcga_dr_table[["pharmaceutical_therapy_drug_name"]])
    # Sometimes drug regimen continue for a while, but we must simplify for survival analysis
    # Sort the table so that later timepoints are identified as duplicates 
    tcga_dr_table <- tcga_dr_table[order(tcga_dr_table[["pharmaceutical_tx_started_days_to"]]),]
    duped_treats <- duplicated(tcga_dr_table[c("bcr_patient_barcode", "pharmaceutical_therapy_drug_name")])
    #with(tcga_dr_table[duped_treats, ], table(tcga_project))
    #dup_rows <- duplicated(tcga_dr_table)
    #tcga_dr_table <- tcga_dr_table[!dup_rows,]
    tcga_dr_table <- tcga_dr_table[!duped_treats,]
    tcga_dr_table[["pharmaceutical_therapy_drug_name"]] <- tolower(tcga_dr_table[["pharmaceutical_therapy_drug_name"]])
    
    if (tcga_brca || scanb) {
      tcga_dr_table <- tcga_dr_table[tcga_dr_table[["tcga_project"]] == "BRCA",]
    }
    tcga_drugs <- unique(tcga_dr_table[["pharmaceutical_therapy_drug_name"]])
    ctrp_drugs <- get_ctrp_drugs(n_drugs = length(grep("^dr_pred_[0-9]+$", colnames(internal_predictions))))
    ctrp_drug_processed <- process_ctrp_drug_names(ctrp_drugs)
    ctrp_drug_processed[grep("tipifarnib", ctrp_drug_processed)] <- "tipifarnib"
    
    # Create drug prediction table and ensure correct column naming
    tcga_drug_preds <- tibble::tibble(
      tcga_predictions[,grep("dr_pred_[0-9]+", colnames(tcga_predictions))])
    pred_col_idx <- as.integer(gsub("dr_pred_", "", colnames(tcga_drug_preds)))+1
    colnames(tcga_drug_preds)[pred_col_idx] <- paste0(process_ctrp_drug_names(ctrp_drugs), "_sensitivity")
    
    tcga_output_table <- tcga_drug_preds
    tcga_output_table[["survival_risk"]] <- tcga_predictions[["survival_risk_0"]]
    tcga_output_table <- tibble::tibble(
      id = rownames(tcga_predictions), 
      tcga_surv_data[tcga_surv_ind,-1], 
      internal_embeddings[rownames(tcga_predictions), shared_embedding_names], 
      tcga_output_table, 
      .name_repair = "unique"
    )
    tcga_dr_table_wide <- tidyr::pivot_wider(
      tcga_dr_table, 
      id_cols = "bcr_patient_barcode", 
      names_from = "pharmaceutical_therapy_drug_name", 
      values_from = "pharmaceutical_tx_started_days_to"
    )
    colnames(tcga_dr_table_wide)[-1] <- paste0(
      tolower(colnames(tcga_dr_table_wide)[-1]), 
      "_tx_started_days_to")
    
    tcga_final_output_table <- dplyr::left_join(
      tcga_output_table, 
      tcga_dr_table_wide, 
      by = "bcr_patient_barcode"
    )
    
    # NRF2 score and KEAP1, NFE2L2 mutations
    fn <- paste0(harkonen_path, "1-s2.0-S2213231723000459-mmc3.xlsx")
    nrf2_muts <- readxl::read_xlsx(fn, sheet = 1)
    nrf2_muts_ptr <- match(
      tcga_final_output_table[["bcr_patient_barcode"]],
      nrf2_muts[["barcode"]]
    )
    tcga_final_output_table[["NRF2_score"]] <- nrf2_muts[["Nrf2_Score"]][nrf2_muts_ptr]
    tcga_final_output_table[["NRF2_mutations"]] <- nrf2_muts[["Mutation"]][nrf2_muts_ptr]
    
    # Subtypes
    tcga_subtypes <- TCGAbiolinks::PanCancerAtlas_subtypes()
    subtype_ptr <- match(
      tcga_final_output_table[["bcr_patient_barcode"]], 
      substr(tcga_subtypes[["pan.samplesID"]], 1, 12)
    )
    subtype_cols <- grep("subtype", colnames(tcga_subtypes), ignore.case = TRUE)
    tcga_final_output_table <- data.frame(tcga_final_output_table, tcga_subtypes[subtype_ptr, subtype_cols])
    
    fnw <- paste0(prediction_path, "internal/unified_patient_output_table.csv.gz")
    readr::write_csv(tcga_final_output_table, file = gzfile(fnw))
  }
  
  # DR performance
  best_dr_res_path <- paste0(save_path, "../plots/", model_name, "_drug_responses/")
  dir.create(best_dr_res_path, recursive = TRUE)
  
  drug_response_r2 <- read_result_files(
    "drug_response_r2.csv.gz", 
    path = save_path, 
    model_name = model_name)
  drug_response_r2_best <- drug_response_r2[drug_response_r2[["task"]] == best_task,]
  dr_r2_best_long <- plyr::ddply(
    drug_response_r2_best, 
    c("drug", "dataset"), 
    function(x) data.frame(
      r2_mean = mean(x$dr_r2), 
      r2_sd = sd(x$dr_r2)))
  dr_r2_mean_best <- reshape2::dcast(
    dr_r2_best_long, 
    drug ~ dataset, 
    value.var = "r2_mean")
  mean_col_map <- c(drug = "drug", cl_test = "r2_test_mean", cl_train = "r2_train_mean")
  colnames(dr_r2_mean_best) <- mean_col_map[colnames(dr_r2_mean_best)]
  dr_r2_sd_best <- reshape2::dcast(
    dr_r2_best_long, 
    drug ~ dataset, 
    value.var = "r2_sd")
  sd_col_map <- c(drug = "drug", cl_test = "r2_test_sd", cl_train = "r2_train_sd")
  colnames(dr_r2_sd_best) <- sd_col_map[colnames(dr_r2_sd_best)]
  dr_r2_formatted <- plyr::join(dr_r2_mean_best, dr_r2_sd_best, by = "drug")
  
  drug_response_mse <- read_result_files(
    "drug_response_mse.csv.gz", 
    path = save_path, 
    model_name = model_name)
  drug_response_mse_best <- drug_response_mse[drug_response_mse[["task"]] == best_task,]
  dr_mse_best_long <- plyr::ddply(
    drug_response_mse_best, 
    c("drug", "dataset"), 
    function(x) data.frame(
      mse_mean = mean(x$dr_mse), 
      mse_sd = sd(x$dr_mse)))
  dr_mse_mean_best <- reshape2::dcast(
    dr_mse_best_long, 
    drug ~ dataset, 
    value.var = "mse_mean")
  mean_col_map <- c(drug = "drug", cl_test = "mse_test_mean", cl_train = "mse_train_mean")
  colnames(dr_mse_mean_best) <- mean_col_map[colnames(dr_mse_mean_best)]
  dr_mse_sd_best <- reshape2::dcast(
    dr_mse_best_long, 
    drug ~ dataset, 
    value.var = "mse_sd")
  sd_col_map <- c(drug = "drug", cl_test = "mse_test_sd", cl_train = "mse_train_sd")
  colnames(dr_mse_sd_best) <- sd_col_map[colnames(dr_mse_sd_best)]
  dr_mse_formatted <- plyr::join(dr_mse_mean_best, dr_mse_sd_best, by = "drug")
  
  dr_perf_table <- plyr::join(dr_r2_formatted, dr_mse_formatted, by = "drug")
  dr_perf_table[["drug"]] <- process_ctrp_drug_names(
    fix_ctrp_drug_names(
      dr_perf_table[["drug"]], 
      n_drugs = length(grep("^dr_pred_[0-9]+$", colnames(internal_predictions)))
    )
  )
  dr_perf_table[["drug"]] <- tolower(dr_perf_table[["drug"]])
  dr_perf_table[["in TCGA"]] <- dr_perf_table[["drug"]] %in% tolower(tcga_drugs)
  readr::write_csv(dr_perf_table, paste0(plot_path, "best_model_drugwise_performance.csv"))
}