script_root <- Sys.getenv("MODAE_SCRIPT_PATH")
if (script_root == "") {
  stop("Please define the MODAE_SCRIPT_PATH environment variable.")
}
source(paste0(script_root, "R/setup.R"))

codeae_comparison <- TRUE
baseline_comparison <- TRUE
save_plots <- FALSE

for (save_dir in save_dirs) {
  if(exists("var_list")) try(detach(var_list))
  var_list <- get_paths_and_parameters(save_dir, remove_incomplete = FALSE)
  attach(var_list)
  best_task <- readLines(paste0(save_path, "best_task.txt"))
  param_best_task_ptr <- match(best_task, parameters[["task"]])
  shared_embedding_names <- paste0(
    "z", 
    1:parameters[param_best_task_ptr, "bottle_neck"]
  )
  
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
    
    fn <- paste0(
      prediction_path, 
      "internal_survival_validation_final_predictions.csv.gz"
    )
    internal_predictions <- read.csv(fn, header = TRUE, row.names = 1)
    
    if (tcga_brca || scanb) {
      p_ext_fn <- paste0(
        prediction_path, 
        "external_survival_validation_predictions.csv.gz"
      ) 
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
      col_intersect <- intersect(
        colnames(tcga_surv_data), 
        colnames(tcga_subtypes)
      )
      tcga_surv_data[["patient"]] <- tcga_surv_data[["bcr_patient_barcode"]]
      tcga_surv_data <- plyr::join(
        tcga_surv_data, 
        tcga_subtypes[, !colnames(tcga_subtypes) %in% col_intersect], 
        type = "left")
      tcga_surv_data[["simplified_stage"]] <- gsub(
        "[A-C]+$", "", tcga_surv_data[["ajcc_pathologic_tumor_stage"]])
      
      scanb_surv_ptr <- match(
        rownames(scanb_predictions), 
        scanb_surv_data[["GEX.assay"]])
    }
    tcga_surv_ptr <- match(
      substr(rownames(tcga_predictions), 1, 12), 
      tcga_surv_data[["bcr_patient_barcode"]])
    
    tcga_surv_type <- tcga_surv_data[tcga_surv_ptr, "type"]
    tcga_surv_type[is.na(tcga_surv_type)] <- "NA"
    
    dr_dir <- paste0(
      data_root, 
      "CODE-AE-v1.0/data/tcga/"
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
    
    # Omit prior treatments
    tcga_dr_table <- tcga_dr_table_raw[
      tcga_dr_table_raw[["pharmaceutical_tx_started_days_to"]]>=0,
    ]
    # Fix drug name whitespace
    tcga_dr_table[["pharmaceutical_therapy_drug_name"]] <- gsub(
      "\\s", 
      "", 
      tcga_dr_table[["pharmaceutical_therapy_drug_name"]]
    )
    # Sometimes drug regimen continue for a while, but we must simplify 
    # Sort the table so that later timepoints are identified as duplicates 
    tcga_dr_table <- tcga_dr_table[
      order(tcga_dr_table[["pharmaceutical_tx_started_days_to"]]),
    ]
    duped_treats <- duplicated(tcga_dr_table[
      c("bcr_patient_barcode", "pharmaceutical_therapy_drug_name")
    ])
    tcga_dr_table <- tcga_dr_table[!duped_treats,]
    tcga_dr_single_drug_patients <- names(
      which(table(tcga_dr_table[["bcr_patient_barcode"]]) == 1)
    )
    
    if (tcga_brca || scanb) {
      tcga_dr_table <- tcga_dr_table[tcga_dr_table[["tcga_project"]] == "BRCA",]
    }
    tcga_drugs <- unique(tcga_dr_table[["pharmaceutical_therapy_drug_name"]])
    n_drug_cols <- length(grep("^dr_pred_[0-9]+$", colnames(internal_predictions)))
    ctrp_drugs <- get_ctrp_drugs(n_drugs = n_drug_cols)
    ctrp_drug_processed <- process_ctrp_drug_names(ctrp_drugs)
    ctrp_drug_processed[grep("tipifarnib", ctrp_drug_processed)] <- "tipifarnib"
    tcga_drug_in_ctrp <- match(tolower(tcga_drugs), tolower(ctrp_drug_processed))
    
    codeae_pred_path <- paste0(
      codeae_path, 
      "XieResearchGroup-CODE-AE-6dc17a5/", 
      "intermediate_results/tcga_prediction/"
    )
    fn <- paste0(codeae_pred_path, "tcga_scores.csv")
    codeae_predictions <- readr::read_csv(fn, show_col_types = FALSE)
    codeae_drug_cols <- colnames(codeae_predictions)[-1:-2]
    
    baseline_models <- c(
      "pca10_elasticnet",
      "pca100_elasticnet",
      "fs_elasticnet"
    )
    baseline_model_preds <- list()
    for (i in baseline_models) {
      fn <- file.path(
        output_dir,
        "baseline_results",
        "drug_response",
        "ctrp",
        paste0(i, "_predictions.csv.gz")
      )
      baseline_model_preds[[i]] <- readr::read_csv(fn)
      colnames(baseline_model_preds[[i]]) <- c("id", paste0("dr_pred_", 0:(n_drug_cols-1)))
    }
    
    drug_cols <- list()
    for (i in (1:length(tcga_drugs))[!is.na(tcga_drug_in_ctrp)]) {
      drug_cols <- c(
        drug_cols, 
        list(list(
          name = tcga_drugs[i], 
          column = paste0("dr_pred_", tcga_drug_in_ctrp[i] - 1), 
          drug_adjust = TRUE, 
          tcga_quantile = TRUE,
          tcga_all_projects = FALSE, 
          treated_subset = TRUE, 
          covariates = TRUE, 
          surv_event_col = "PFI", 
          surv_time_col = "PFI.time"
        )
      ))
      drug_cols <- c(
        drug_cols, 
        list(list(
          name = tcga_drugs[i], 
          column = paste0("dr_pred_", tcga_drug_in_ctrp[i] - 1), 
          drug_adjust = TRUE, 
          tcga_quantile = TRUE,
          tcga_all_projects = FALSE, 
          treated_subset = TRUE, 
          covariates = TRUE, 
          surv_event_col = "OS", 
          surv_time_col = "OS.time"
        )
        ))
    }
    pred_cols <- c(
      drug_cols, 
      list(list(
        name = "survival_risk", 
        column = "survival_risk_0", 
        drug_adjust = FALSE, 
        tcga_quantile = FALSE, 
        tcga_all_projects = TRUE, 
        treated_subset = FALSE, 
        covariates = FALSE, 
        surv_event_col = "OS", 
        surv_time_col = "OS.time"
      ))
    )
    pred_cols <- c(
      drug_cols, 
      list(list(
        name = "survival_risk", 
        column = "survival_risk_0", 
        drug_adjust = FALSE, 
        tcga_quantile = FALSE, 
        tcga_all_projects = TRUE, 
        treated_subset = FALSE, 
        covariates = FALSE, 
        surv_event_col = "PFI", 
        surv_time_col = "PFI.time"
      ))
    )
    cancer_drug_counts <- with(
      tcga_dr_table, 
      table(
        tcga_project, 
        drug = tolower(pharmaceutical_therapy_drug_name)
      )
    )
    cancer_single_drug_counts <- with(
      tcga_dr_table[
        tcga_dr_table$bcr_patient_barcode %in% tcga_dr_single_drug_patients,
      ], 
      table(tcga_project, drug = tolower(pharmaceutical_therapy_drug_name))
    )
    
    cox_results <- list()
    km_results <- list()
    concordance_results <- list()
    codeae_cox_results <- list()
    codeae_km_results <- list()
    codeae_concordance_results <- list()
    baseline_cox_results <- list()
    baseline_km_results <- list()
    baseline_concordance_results <- list()
    treated_pairs_n <- list()
    
    for (i in 1:length(pred_cols)) {
      pred_col <- pred_cols[[i]][["column"]]
      pred_name <- pred_cols[[i]][["name"]]
      pred_time_adjust <- pred_cols[[i]][["drug_adjust"]]
      pred_tcga_quantile <- pred_cols[[i]][["tcga_quantile"]]
      treated_subset <- pred_cols[[i]][["treated_subset"]]
      covariates <- pred_cols[[i]][["covariates"]]
      surv_event_col <- pred_cols[[i]][["surv_event_col"]]
      surv_time_col <- pred_cols[[i]][["surv_time_col"]]
      
      if (pred_cols[[i]][["tcga_all_projects"]]) {
        cancers <- unique(tcga_surv_type)
      } else {
        cancers <- as.list(names(which(cancer_drug_counts[,tolower(pred_name)] > 5)))
      }
      for (canceri in cancers) {
        canceri_path <- paste0(plot_path, "treatment_survival_km/", canceri)
        if (canceri == "NA") next
        dir.create(canceri_path, recursive = TRUE)
        cancer_ind <- tcga_surv_type %in% canceri
        cancer_sample_ids <- rownames(tcga_predictions)[cancer_ind]
        cancer_surv_ptr <- tcga_surv_ptr[cancer_ind]
        tcga_risk_pred <- tcga_predictions[cancer_sample_ids, pred_col]
        tcga_risk_quantiles <- mean(tcga_risk_pred, na.rm = TRUE)
        if (tcga_brca || scanb) {
          scanb_risk_quantiles <- quantile(
            scanb_predictions[cancer_sample_ids, pred_col], 
            probs = 1/2, 
            na.rm = TRUE)
          use_quants <- ifelse(
            pred_tcga_quantile, 
            tcga_risk_quantiles, 
            scanb_risk_quantiles)
          scanb_risk_bin <- cut(
            scanb_predictions[cancer_sample_ids, pred_col], 
            breaks = c(-Inf, use_quants, Inf), 
            labels = c("low", "high"))
        } else {
          use_quants <- tcga_risk_quantiles
        }
        tcga_risk_bin <- cut(
          tcga_risk_pred, 
          breaks = c(-Inf, use_quants, Inf), 
          labels = c("low", "high"))
        survival_cox_f <- function(predictors, event, time, maxtime = 3000) {
          event[time > maxtime] <- FALSE
          if (!any(event)) return(NULL)
          if (all(event)) return(NULL)
          time[time > maxtime] <- maxtime
          nna_ind <- !is.na(event)
          surv_df <- data.frame(predictors, event, time)[nna_ind,]
          sf <- survival::survfit(
            survival::Surv(time, event, type = "right") ~ group, 
            data = surv_df, 
            model = TRUE)
          attributes(sf)$survival_pvalue <- survminer::surv_pvalue(
            sf, 
            data = surv_df
          )
          return(sf)
        }
        survival_f <- function(group, event, time, maxtime = 3000) {
          event[time > maxtime] <- FALSE
          if (!any(event)) return(NULL)
          if (all(event)) return(NULL)
          time[time > maxtime] <- maxtime
          nna_ind <- !is.na(event)
          surv_df <- data.frame(group, event, time)[nna_ind,]
          sf <- survival::survfit(
            survival::Surv(time, event, type = "right") ~ group, 
            data = surv_df, 
            model = TRUE)
          attributes(sf)$survival_pvalue <- survminer::surv_pvalue(
            sf, 
            data = surv_df
          )
          return(sf)
        }
        tcga_nna_ptr <- which(!(
          is.na(tcga_surv_data[cancer_surv_ptr, surv_event_col]) | 
            is.na(tcga_surv_data[cancer_surv_ptr, surv_time_col]) | 
            is.na(tcga_risk_bin)
        ))
        # Can't skip missing or invariantt predictions for performance comparison
        #if (length(tcga_nna_ptr) == 0) next
        #if (sd(tcga_risk_pred[tcga_nna_ptr]) == 0) next
        tcga_covar_names <- c(
          age_at_initial_pathologic_diagnosis = "age", 
          #gender = "sex", 
          #race = "ancestry", 
          ajcc_pathologic_tumor_stage = "stage"
        )
        cancer_covars <- tcga_surv_data[
          cancer_surv_ptr[tcga_nna_ptr], 
          names(tcga_covar_names), 
          drop = FALSE
        ]
        colnames(cancer_covars) <- tcga_covar_names[colnames(cancer_covars)]
        na_vals <- c(
          "[Discrepancy]", 
          "[Not Applicable]", 
          "[Not Available]", 
          "[Unknown]", 
          "[Not Evaluated]"
        )
        cancer_covars[is.na(cancer_covars)] <- NA
        cancer_covars <- lapply(
          cancer_covars, 
          function(x) ifelse(x %in% na_vals, NA, x)
        )
        cancer_covar_missing <- sapply(
          cancer_covars, 
          function(x) mean(is.na(x))
        ) > 0.5
        cancer_response_predictors <- data.frame(
          cancer_covars[!cancer_covar_missing],
          sensitivity = scale(tcga_risk_pred[tcga_nna_ptr])
        )
        # Simplify stage information
        if (!is.null(cancer_response_predictors[["stage"]])) {
          cancer_response_predictors[["stage"]] <- simplify_cancer_stages(
            cancer_response_predictors[["stage"]]
          )
        }
        rownames(cancer_response_predictors) <- cancer_sample_ids[tcga_nna_ptr]
        if (treated_subset) {
          treat_ind <- tcga_dr_table[["pharmaceutical_therapy_drug_name"]] == pred_name
          treat_id <- tcga_dr_table[treat_ind, "bcr_patient_barcode"]
          if (any(duplicated(treat_id))) stop("Duplicated treatments")
          treat_time <- tcga_dr_table[
            treat_ind, 
            "pharmaceutical_tx_started_days_to"
          ]
          
          #treat_ind <- which(treat_ind)[treat_time >= 0]
          treat_id <- treat_id[treat_time >= 0]
          treat_time <- treat_time[treat_time >= 0]
          names(treat_time) <- treat_id
          
          # Index of binned values with valid survival data and matching treatment
          bin_treat_ptr <- intersect(
            which(substr(cancer_sample_ids, 1, 12) %in% treat_id), 
            tcga_nna_ptr
          )
          bin_treat_id <- substr(cancer_sample_ids, 1, 12)[bin_treat_ptr]
          bin_treat_ptr <- match(bin_treat_id, substr(cancer_sample_ids, 1, 12))
          if (any(bin_treat_ptr != bin_treat_ptr)) {
            stop("Indexing assumption is incorrect.") 
            # meaning equivalent code in single drug subset
          }
          # Match binned index to treatment start time index
          time_ptr <- match(bin_treat_id, treat_id)
          
          time_adjusted <- (
            tcga_surv_data[
              cancer_surv_ptr[bin_treat_ptr], 
              surv_time_col
            ] - 
              treat_time[time_ptr] * ifelse(pred_time_adjust, 1, 0)
          )
          positive_time_ind <- time_adjusted > 0
          treated_pairs_n <- c(treated_pairs_n, list(data.frame(
            drug = pred_name, 
            cancer = canceri, 
            survival = surv_event_col, 
            n = sum(positive_time_ind)
          )))
          # Skip if not at least 2 patients with positive time
          if (sum(positive_time_ind)<2) next
          
          cox_table <- data.frame(
            cancer_response_predictors[bin_treat_id,][positive_time_ind,], 
            event = tcga_surv_data[cancer_surv_ptr[bin_treat_ptr][positive_time_ind], surv_event_col] == 1, 
            time = time_adjusted[positive_time_ind]
          )
          if (nrow(cox_table) > 1) {
            unitary_factors <- sapply(
              lapply(cox_table, table, useNA = "no"), 
              length
            ) == 1
            cox_table <- cox_table[
              , 
              !sapply(cox_table, is.character) | !unitary_factors, 
              drop = FALSE
            ]
            cox_formula <- as.formula(
              paste0(
                "survival::Surv(time, event) ~ ", 
                paste(
                  c(
                    intersect(
                      tcga_covar_names[!cancer_covar_missing], 
                      colnames(cox_table)
                    ), 
                    "sensitivity"
                    ), 
                  collapse = " + "
                )
              )
            )
            try({
              cox_model <- survival::coxph(cox_formula, data = cox_table)
              test <- summary(cox_model)
              temp <- data.frame(
                drug = pred_name, 
                cancer = canceri, 
                survival = surv_event_col, 
                var = rownames(test$coefficients), 
                coef = test$coefficients[,"coef"], 
                se_coef = test$coefficients[,"se(coef)"],
                z = test$coefficients[,"z"],
                p = test$coefficients[,"Pr(>|z|)"]
              )
              cox_results <- c(cox_results, list(temp))
            })
            # Investigate NA coefficients
            if (length(cox_results)>0) {
              coef_only_na <- with(
                cox_results[[length(cox_results)]], 
                is.na(coef[var=="sensitivity"]) & !is.na(p[var=="sensitivity"]) 
              )
              if (coef_only_na) stop("Unexplained Cox NA.")
            }
            try({
              cox_plot <- survminer::ggforest(
                cox_model, data = cox_table, fontsize = 0.6, main = NULL, 
                cpositions = c(0.02, 0.25, 0.5))
              fn <- paste0(
                canceri_path, 
                "/", 
                canceri, 
                "_", 
                pred_name, 
                figure_str, 
                "_", 
                surv_event_col, 
                "_coxph.png"
              )
              if (save_plots) save_figure_safe(
                cox_plot, 
                png, 
                fn, 
                width = plot_width, 
                height = plot_height, 
                res = plot_res, 
                units = plot_units
              )
            })
            if (codeae_comparison) {
              codeae_ptr <- grep(
                tolower(only_alphanumericals(pred_name)),  
                tolower(only_alphanumericals(codeae_drug_cols))
              )
              for (codeae_i in codeae_ptr) {
                codeae_pred_name <- codeae_drug_cols[codeae_i]
                codeae_predi <- dplyr::left_join(
                  data.frame(Sample = substr(rownames(cox_table), 1, 15)), 
                  codeae_predictions[,c("Sample", codeae_pred_name)], 
                  by = "Sample"
                )
                cox_table[["sensitivity"]] <- codeae_predi[[2]]
                try({
                  cox_model <- survival::coxph(cox_formula, data = cox_table)
                  test <- summary(cox_model)
                  temp <- data.frame(
                    drug = pred_name, 
                    drug_codeae = codeae_pred_name, 
                    cancer = canceri, 
                    survival = surv_event_col, 
                    var = rownames(test$coefficients), 
                    coef = test$coefficients[,"coef"], 
                    se_coef = test$coefficients[,"se(coef)"],
                    z = test$coefficients[,"z"],
                    p = test$coefficients[,"Pr(>|z|)"]
                  )
                  codeae_cox_results <- c(codeae_cox_results, list(temp))
                })
                try({
                  cox_plot <- survminer::ggforest(
                    cox_model, data = cox_table, fontsize = 0.6, main = NULL, 
                    cpositions = c(0.02, 0.25, 0.5))
                  fn <- paste0(
                    canceri_path, 
                    "/", 
                    canceri, 
                    "_CODE-AE_", 
                    codeae_pred_name, 
                    figure_str, 
                    "_", 
                    surv_event_col, 
                    "_coxph.png"
                  )
                  if (save_plots) save_figure_safe(
                    cox_plot, 
                    png, 
                    fb, 
                    width = plot_width, 
                    height = plot_height, 
                    res = plot_res, 
                    units = plot_units
                  )
                })
              }
            }
            if (baseline_comparison) {
              for (blmod in names(baseline_model_preds)) {
                baseline_predi <- dplyr::left_join(
                  data.frame(id = rownames(cox_table)), 
                  baseline_model_preds[[blmod]][,c("id", pred_col)], 
                  by = "id"
                )
                cox_table[["sensitivity"]] <- baseline_predi[[2]]
                try({
                  cox_model <- survival::coxph(cox_formula, data = cox_table)
                  test <- summary(cox_model)
                  temp <- data.frame(
                    baseline_model = blmod,
                    drug = pred_name, 
                    cancer = canceri, 
                    survival = surv_event_col, 
                    var = rownames(test$coefficients), 
                    coef = test$coefficients[,"coef"], 
                    se_coef = test$coefficients[,"se(coef)"],
                    z = test$coefficients[,"z"],
                    p = test$coefficients[,"Pr(>|z|)"]
                  )
                  baseline_cox_results <- c(baseline_cox_results, list(temp))
                })
                try({
                  cox_plot <- survminer::ggforest(
                    cox_model, data = cox_table, fontsize = 0.6, main = NULL, 
                    cpositions = c(0.02, 0.25, 0.5))
                  fn <- paste0(
                    canceri_path, 
                    "/", 
                    canceri, 
                    "_",
                    blmod,
                    "_", 
                    codeae_pred_name, 
                    figure_str, 
                    "_", 
                    surv_event_col, 
                    "_coxph.png"
                  )
                  if (save_plots) save_figure_safe(
                    cox_plot, 
                    png, 
                    fb, 
                    width = plot_width, 
                    height = plot_height, 
                    res = plot_res, 
                    units = plot_units
                  )
                })
              }
            }
          }
          # Patients treated with only one drug
          single_drug_treat_id <- intersect(
            bin_treat_id, 
            tcga_dr_single_drug_patients
          )
          if (length(single_drug_treat_id) > 1) {
            single_drug_bin_treat_ptr <- match(
              single_drug_treat_id, 
              substr(cancer_sample_ids, 1, 12)
            )
            single_drug_time_ptr <- match(single_drug_treat_id, treat_id)
            sd_cox_table <- data.frame(
              cancer_response_predictors[single_drug_treat_id,], 
              event = tcga_surv_data[
                cancer_surv_ptr[single_drug_bin_treat_ptr], 
                surv_event_col
              ] == 1, 
              time = (
                tcga_surv_data[
                  cancer_surv_ptr[single_drug_bin_treat_ptr], 
                  surv_time_col
                ] - (
                  treat_time[single_drug_time_ptr] * 
                    ifelse(pred_time_adjust, 1, 0)
                )
              )
            )
            # Remove any observations of with 0 or negative time
            sd_cox_table <- sd_cox_table[
              sd_cox_table[["time"]] > 0,, 
              drop= FALSE
            ]
            unitary_factors <- sapply(
              lapply(sd_cox_table, table, useNA = "no"), 
              length
            ) == 1
            sd_cox_table <- sd_cox_table[
              , 
              !sapply(sd_cox_table, is.character) | !unitary_factors, 
              drop = FALSE
            ]
            cox_formula <- as.formula(
              paste0(
                "survival::Surv(time, event) ~ ", 
                paste(
                  c(
                    intersect(
                      tcga_covar_names[!cancer_covar_missing], 
                      colnames(sd_cox_table)
                    ), 
                    "sensitivity"
                  ), 
                  collapse = " + "
                )
              )
            )
            try({
              cox_model <- survival::coxph(cox_formula, data = sd_cox_table)
              cox_plot <- survminer::ggforest(
                cox_model, data = sd_cox_table, fontsize = 0.6, main = NULL, 
                cpositions = c(0.02, 0.25, 0.5))
              fn <- paste0(
                canceri_path, 
                "/", 
                canceri, 
                "_", 
                pred_name, 
                figure_str, 
                "_single_drug_", 
                surv_event_col, 
                "_coxph.png"
              )
              if (save_plots) save_figure_safe(
                cox_plot, 
                png, 
                fn, 
                width = plot_width, 
                height = plot_height, 
                res = plot_res, 
                units = plot_units
              )
            })
            if (codeae_comparison) {
              codeae_ptr <- grep(
                tolower(only_alphanumericals(pred_name)),  
                tolower(only_alphanumericals(codeae_drug_cols))
              )
              for (codeae_i in codeae_ptr) {
                codeae_pred_name <- codeae_drug_cols[codeae_i]
                codeae_predi <- dplyr::left_join(
                  data.frame(Sample = substr(rownames(sd_cox_table), 1, 15)), 
                  codeae_predictions[,c("Sample", codeae_pred_name)], 
                  by = "Sample"
                )
                sd_cox_table[["sensitivity"]] <- codeae_predi[[2]]
                sd_cox_table <- sd_cox_table[!is.na(sd_cox_table[["sensitivity"]]),]
                try({
                  cox_model <- survival::coxph(cox_formula, data = sd_cox_table)
                  cox_plot <- survminer::ggforest(
                    cox_model, data = sd_cox_table, fontsize = 0.6, main = NULL, 
                    cpositions = c(0.02, 0.25, 0.5))
                  fn <- paste0(
                    canceri_path, 
                    "/", 
                    canceri, 
                    "_CODE-AE_", 
                    codeae_pred_name, 
                    figure_str, 
                    "_single_drug_", 
                    surv_event_col, 
                    "_coxph.png"
                  )
                  if (save_plots) save_figure_safe(
                    cox_plot, 
                    png, 
                    fn, 
                    width = plot_width, 
                    height = plot_height, 
                    res = plot_res, 
                    units = plot_units
                  )
                })
              }
            }
            if (baseline_comparison) {
              for (blmod in names(baseline_model_preds)) {
                baseline_predi <- dplyr::left_join(
                  data.frame(id = rownames(sd_cox_table)), 
                  baseline_model_preds[[blmod]][,c("id", pred_col)], 
                  by = "id"
                )
                sd_cox_table[["sensitivity"]] <- baseline_predi[[2]]
                sd_cox_table <- sd_cox_table[!is.na(sd_cox_table[["sensitivity"]]),]
                try({
                  cox_model <- survival::coxph(cox_formula, data = sd_cox_table)
                  cox_plot <- survminer::ggforest(
                    cox_model, data = sd_cox_table, fontsize = 0.6, main = NULL, 
                    cpositions = c(0.02, 0.25, 0.5))
                  fn <- paste0(
                    canceri_path, 
                    "/", 
                    canceri, 
                    "_",
                    blmod,
                    "_", 
                    pred_name, 
                    figure_str, 
                    "_single_drug_", 
                    surv_event_col, 
                    "_coxph.png"
                  )
                  if (save_plots) save_figure_safe(
                    cox_plot, 
                    png, 
                    fn, 
                    width = plot_width, 
                    height = plot_height, 
                    res = plot_res, 
                    units = plot_units
                  )
                })
              }
            }
          }
          
          if(any(tcga_surv_data[cancer_surv_ptr[bin_treat_ptr], "bcr_patient_barcode"] != bin_treat_id)) {
            stop("Mis-matched IDs in treatment-survival tables.")
          }
          
          surv_df <- data.frame(
            event = tcga_surv_data[
              cancer_surv_ptr[bin_treat_ptr], 
              surv_event_col
            ][positive_time_ind] == 1, 
            time = time_adjusted[positive_time_ind],
            sensitivity = tcga_risk_pred[bin_treat_ptr][positive_time_ind]
          )
          surv_c <- survival::concordance(
            survival::Surv(time, event) ~ sensitivity, 
            data = surv_df
          )
          temp <- data.frame(
            drug = pred_name, 
            cancer = canceri, 
            survival = surv_event_col, 
            concordance = surv_c$concordance, 
            n = surv_c$n
          )
          concordance_results <- c(concordance_results, list(temp))
          
          tcga_sf <- survival_f(
            group = tcga_risk_bin[bin_treat_ptr][positive_time_ind], 
            event = tcga_surv_data[
              cancer_surv_ptr[bin_treat_ptr], 
              surv_event_col
            ][positive_time_ind] == 1, 
            time = time_adjusted[positive_time_ind]
          )
          try({
            temp <- data.frame(
              drug = pred_name, 
              cancer = canceri, 
              survival = surv_event_col, 
              p = attributes(tcga_sf)$survival_pvalue$pval
            )
            km_results <- c(km_results, list(temp))
          })
          if (codeae_comparison) {
            codeae_ptr <- grep(
              tolower(only_alphanumericals(pred_name)),  
              tolower(only_alphanumericals(codeae_drug_cols))
            )
            for (codeae_i in codeae_ptr) {
              codeae_pred_name <- codeae_drug_cols[codeae_i]
              codeae_predi <- dplyr::left_join(
                data.frame(
                  Sample = substr(cancer_sample_ids[bin_treat_ptr], 1, 15)
                ), 
                codeae_predictions[,c("Sample", codeae_pred_name)], 
                by = "Sample"
              )
              surv_df <- data.frame(
                event = tcga_surv_data[
                  cancer_surv_ptr[bin_treat_ptr], 
                  surv_event_col
                ][positive_time_ind] == 1, 
                time = time_adjusted[positive_time_ind],
                sensitivity = codeae_predi[[codeae_pred_name]][positive_time_ind]
              )
              surv_c <- survival::concordance(
                survival::Surv(time, event) ~ sensitivity, 
                data = surv_df
              )
              temp <- data.frame(
                drug = pred_name, 
                drug_codeae = codeae_pred_name, 
                cancer = canceri, 
                survival = surv_event_col, 
                concordance = surv_c$concordance, 
                n = surv_c$n
              )
              codeae_concordance_results <- c(
                codeae_concordance_results, 
                list(temp)
              )
              codeae_tcga_sf <- survival_f(
                group = factor(
                  ifelse(codeae_predi[[codeae_pred_name]]>0.5, "high", "low"), 
                  levels = c("low", "high")
                )[positive_time_ind], 
                event = tcga_surv_data[
                  cancer_surv_ptr[bin_treat_ptr], 
                  surv_event_col
                ][positive_time_ind] == 1, 
                time = time_adjusted[positive_time_ind]
              )
              try({
                temp <- data.frame(
                  drug = pred_name, 
                  drug_codeae = codeae_pred_name, 
                  cancer = canceri, 
                  survival = surv_event_col, 
                  p = attributes(codeae_tcga_sf)$survival_pvalue$pval
                )
                codeae_km_results <- c(codeae_km_results, list(temp))
              })
              if (!is.null(codeae_tcga_sf)) {
                group_colors <- c("#4DAF4A", "#377EB8", "#E41A1C")
                names(group_colors) <- c("low", "medium", "high")
                km_color_scale <- scale_color_manual(values = group_colors)
                plot_title <- paste(
                  canceri, 
                  codeae_pred_name, 
                  "survival", 
                  attributes(codeae_tcga_sf)$survival_pvalue["pval.txt"]
                )
                codeae_tcga_km <- GGally::ggsurv(
                  codeae_tcga_sf, 
                  cens.shape = 3, 
                  order.legend = FALSE
                ) + 
                  theme_bw() + 
                  ylim(c(0,1)) + 
                  km_color_scale + 
                  ggtitle(plot_title)
                fn <- paste0(
                  canceri_path, 
                  "/", 
                  canceri, 
                  "_CODE-AE_", 
                  codeae_pred_name, 
                  figure_str, 
                  "_", 
                  surv_event_col, 
                  ".png"
                )
                if (save_plots) save_figure_safe(
                  codeae_tcga_km, 
                  png, 
                  fn, 
                  width = plot_width * 0.6, 
                  height = plot_height, 
                  res = plot_res, 
                  units = plot_units
                )
              }
            }
          }
          if (baseline_comparison) {
            for (blmod in names(baseline_model_preds)) {
              baseline_predi <- dplyr::left_join(
                data.frame(id = cancer_sample_ids[bin_treat_ptr]), 
                baseline_model_preds[[blmod]][,c("id", pred_col)], 
                by = "id"
              )
              surv_df <- data.frame(
                event = tcga_surv_data[
                  cancer_surv_ptr[bin_treat_ptr], 
                  surv_event_col
                ][positive_time_ind] == 1, 
                time = time_adjusted[positive_time_ind],
                sensitivity = baseline_predi[[pred_col]][positive_time_ind]
              )
              surv_c <- survival::concordance(
                survival::Surv(time, event) ~ sensitivity, 
                data = surv_df
              )
              temp <- data.frame(
                drug = pred_name, 
                baseline_model = blmod, 
                cancer = canceri, 
                survival = surv_event_col, 
                concordance = surv_c$concordance, 
                n = surv_c$n
              )
              baseline_concordance_results <- c(
                baseline_concordance_results, 
                list(temp)
              )
              baseline_tcga_risk_pred <- dplyr::left_join(
                data.frame(id = cancer_sample_ids), 
                baseline_model_preds[[blmod]][,c("id", pred_col)], 
                by = "id"
              )
              baseline_group = cut(
                baseline_tcga_risk_pred[[pred_col]], 
                breaks = c(-Inf, mean(baseline_tcga_risk_pred[[pred_col]], na.rm = TRUE), Inf), 
                labels = c("low", "high"))
              baseline_tcga_sf <- survival_f(
                group = baseline_group[bin_treat_ptr][positive_time_ind], 
                event = tcga_surv_data[
                  cancer_surv_ptr[bin_treat_ptr], 
                  surv_event_col
                ][positive_time_ind] == 1, 
                time = time_adjusted[positive_time_ind]
              )
              try({
                temp <- data.frame(
                  drug = pred_name, 
                  baseline_model = blmod, 
                  cancer = canceri, 
                  survival = surv_event_col, 
                  p = attributes(baseline_tcga_sf)$survival_pvalue$pval
                )
                baseline_km_results <- c(baseline_km_results, list(temp))
              })
              if (!is.null(baseline_tcga_sf)) {
                group_colors <- c("#4DAF4A", "#377EB8", "#E41A1C")
                names(group_colors) <- c("low", "medium", "high")
                km_color_scale <- scale_color_manual(values = group_colors)
                plot_title <- paste(
                  blmod,
                  canceri,
                  pred_name,
                  "survival", 
                  attributes(baseline_tcga_sf)$survival_pvalue["pval.txt"]
                )
                baseline_tcga_km <- GGally::ggsurv(
                  baseline_tcga_sf, 
                  cens.shape = 3, 
                  order.legend = FALSE
                ) + 
                  theme_bw() + 
                  ylim(c(0,1)) + 
                  km_color_scale + 
                  ggtitle(plot_title)
                fn <- paste0(
                  canceri_path, 
                  "/", 
                  canceri, 
                  "_",
                  blmod,
                  "_", 
                  pred_name, 
                  figure_str, 
                  "_", 
                  surv_event_col, 
                  ".png"
                )
                if (save_plots) save_figure_safe(
                  baseline_tcga_km, 
                  png, 
                  fn, 
                  width = plot_width * 0.6, 
                  height = plot_height, 
                  res = plot_res, 
                  units = plot_units
                )
              }
            }
          }
        } else {
          surv_df <- data.frame(
            event = tcga_surv_data[
              cancer_surv_ptr[tcga_nna_ptr], 
              surv_event_col
            ] == 1, 
            time = tcga_surv_data[
              cancer_surv_ptr[tcga_nna_ptr], 
              surv_time_col
            ],
            sensitivity = tcga_risk_pred[tcga_nna_ptr]
          )
          surv_c <- survival::concordance(
            survival::Surv(time, event) ~ sensitivity, 
            data = surv_df
          )
          temp <- data.frame(
            drug = pred_name, 
            cancer = canceri, 
            survival = surv_event_col, 
            concordance = surv_c$concordance, 
            n = surv_c$n
          )
          concordance_results <- c(concordance_results, list(temp))
          
          tcga_sf <- survival_f(
            group = tcga_risk_bin[tcga_nna_ptr], 
            event = tcga_surv_data[
              cancer_surv_ptr[tcga_nna_ptr], 
              surv_event_col
            ] == 1, 
            time = tcga_surv_data[
              cancer_surv_ptr[tcga_nna_ptr], 
              surv_time_col
            ]
          )
        }
        
        if (is.null(tcga_sf)) next
        
        group_colors <- c(
          low = "#4DAF4A", 
          medium = "#377EB8", 
          high = "#E41A1C"
        )
        km_color_scale <- scale_color_manual(values = group_colors)
        plot_title <- paste(
          canceri, 
          pred_name, 
          "survival", 
          attributes(tcga_sf)$survival_pvalue["pval.txt"]
        )
        tcga_km <- GGally::ggsurv(
          tcga_sf, 
          cens.shape = 3, 
          order.legend = FALSE
        ) + 
          theme_bw() + 
          ylim(c(0,1)) + 
          km_color_scale + 
          ggtitle(plot_title)
        
        # TODO: add ggplot2::annotate for p-value and group sizes
        
        if (tcga_brca || scanb) {
          scanb_sf <- survival_f(
            group = scanb_risk_bin, 
            event = scanb_surv_data[scanb_surv_ptr, "OS_event"] == 1, 
            time = scanb_surv_data[scanb_surv_ptr, "OS_days"])
          scanb_km <- GGally::ggsurv(scanb_sf, cens.shape = NA) + 
            theme_bw() + ylim(c(0,1)) + 
            km_color_scale
          model_legend <- cowplot::get_legend(tcga_km)
          add_theme <-  theme(
            legend.position = "none", 
            plot.caption = element_text(hjust = 0., face = "bold", size = 10))
          combined_km_plot <- gridExtra::grid.arrange(
            tcga_km + add_theme + labs(caption = "A"), 
            scanb_km + add_theme + labs(caption = "B"), 
            model_legend, 
            ncol = 3, 
            widths = c(4,4,1))
          fn <- paste0(canceri_path, "/", pred_name, figure_str, ".pdf")
          pdf(fn, width = plot_width / 25.4 * 1, height = plot_height / 25.4 * 0.5)
          grid::grid.draw(combined_km_plot)
          dev.off()
        } else {
          fn <- paste0(
            canceri_path, 
            "/", 
            canceri, 
            "_", 
            pred_name, 
            figure_str, 
            "_", 
            surv_event_col, 
            ".png"
          )
          if (save_plots) save_figure_safe(
            tcga_km, 
            png, 
            fn, 
            width = plot_width * 0.6, 
            height = plot_height, 
            res = plot_res, 
            units = plot_units
          )
        }
        if (length(single_drug_treat_id) > 1) {
          time_adjusted <- (
            tcga_surv_data[
              cancer_surv_ptr[single_drug_bin_treat_ptr], 
              surv_time_col
            ] - (
              treat_time[single_drug_time_ptr] * 
                ifelse(pred_time_adjust, 1, 0)
            )
          )
          positive_time_ind <- time_adjusted > 0
          single_drug_tcga_sf <- survival_f(
            group = tcga_risk_bin[single_drug_bin_treat_ptr][positive_time_ind], 
            event = tcga_surv_data[
              cancer_surv_ptr[single_drug_bin_treat_ptr], 
              surv_event_col
            ][positive_time_ind] == 1, 
            time = time_adjusted[positive_time_ind]
          )
          group_colors <- c(
            low = "#4DAF4A", 
            medium = "#377EB8", 
            high = "#E41A1C"
          )
          km_color_scale <- scale_color_manual(values = group_colors)
          plot_title <- paste(
            canceri, 
            pred_name, 
            "survival", 
            attributes(single_drug_tcga_sf)$survival_pvalue["pval.txt"]
          )
          single_drug_tcga_km <- GGally::ggsurv(
            single_drug_tcga_sf, 
            cens.shape = 3, 
            order.legend = FALSE, 
            plot.cens = FALSE
          ) + 
            theme_bw() + 
            ylim(c(0,1)) + 
            km_color_scale + 
            ggtitle(plot_title)
          fn <- paste0(
            canceri_path, 
            "/", 
            canceri, 
            "_", 
            pred_name, 
            figure_str, 
            "_single_drug_", 
            surv_event_col, 
            ".png"
          )
          if (save_plots) save_figure_safe(
            single_drug_tcga_km, 
            png, 
            fn, 
            width = plot_width * 0.6, 
            height = plot_height, 
            res = plot_res, 
            units = plot_units
          )
        }
        
        # Hazard ratio plot with risk groups and clinicals
        if (tcga_brca || scanb) {
          clinical_surv_list <- list(
            tcga = list(
              group = factor(tcga_risk_bin, levels = c("mid", "low", "high")), 
              data = tcga_surv_data[cancer_surv_ptr,], 
              vars = c(
                "age_at_initial_pathologic_diagnosis", 
                #"ajcc_pathologic_tumor_stage" # Too many subdivisions
                #"pathologic_stage" 
                "simplified_stage" # Has "Stage X"
              ), 
              var_names = c("age_z", "stage"), 
              var_std = c(TRUE, FALSE), 
              #time_var = "OS.time", 
              #event_var = "OS"
              time_var = "PFI.time", 
              event_var = "PFI"
            ), 
            scanb = list(
              group = factor(scanb_risk_bin, c("mid", "low", "high")), 
              data = scanb_surv_data[scanb_surv_ptr,], 
              vars = c(
                "Age_group", 
                "NHG", 
                "Size.mm"
              ), 
              var_names = c("age_z", "NH_grade", "tumor_size_z"), 
              var_std = c(TRUE, FALSE, TRUE), 
              #time_var = "OS_days", 
              #event_var = "OS_event"
              time_var = "RFi_days", 
              event_var = "RFi_event"
            )
          )
          
          sub_na_values <- function(x) ifelse(
            x %in% c("NA", "[Not Available]", "[Discrepancy]", "Stage X"), 
            NA, 
            x
          )
          
          hr_list <- list()
          for (surv_data in clinical_surv_list) {
            temp <- lapply(
              surv_data[["data"]][,surv_data[["vars"]]], 
              sub_na_values
            )
            temp <- data.frame(temp)
            for (j in (1:ncol(temp))[surv_data[["var_std"]]]) {
              if (
                FALSE && # TODO: verify if this still works for BRCA
                surv_data[["var_names"]][j] == "age_z" && 
                pred_time_adjust
              ) {
                # Adjust age by delay in treatment
                temp[[j]] <- temp[[j]] + (treat_time[time_ptr] / 365)
              }
              temp[[j]] <- scale(temp[[j]], center = TRUE, scale = TRUE)
            }
            colnames(temp) <- surv_data[["var_names"]]
            surv_df <- data.frame(
              group = surv_data[["group"]], 
              time = surv_data[["data"]][,surv_data[["time_var"]]],
              event = surv_data[["data"]][,surv_data[["event_var"]]],
              temp)
            if ("stage" %in% colnames(surv_df)) {
              surv_df[["stage"]] <- factor(
                gsub("Stage ", "", surv_df[["stage"]]), 
                levels = c("I", "II", "III", "IV"))
            }
            if ("NH_grade" %in% colnames(surv_df)) {
              surv_df[["NH_grade"]] <- factor(
                surv_df[["NH_grade"]], 
                levels = 1:3)
            }
            cox_formula <- as.formula(
              paste(
                "survival::Surv(time, event == 1) ~", 
                paste(
                  c("group", surv_data[["var_names"]]), 
                  collapse = " + "
                )
              )
            )
            coxph <- survival::coxph(cox_formula, data = surv_df) 
            hr_list[[length(hr_list)+1]] <- survminer::ggforest(
              coxph, data = surv_df, fontsize = 0.6, main = NULL, 
              cpositions = c(0.02, 0.25, 0.5))
          }
          
          hr_add_theme <- theme(
            plot.caption = element_text(hjust = 0.05, face = "bold", size = 10), 
            plot.margin = margin(0, 6, 6, 0, "mm"))
          combined_hr_plot <- gridExtra::grid.arrange(
            hr_list[[1]] + hr_add_theme + labs(caption = "C"), 
            hr_list[[2]] + hr_add_theme + labs(caption = "D"), 
            ncol = 2, 
            widths = c(4.5,4.5))
          combined_plot <- gridExtra::grid.arrange(
            combined_km_plot, 
            combined_hr_plot, 
            ncol = 1, 
            heights = c(0.5,0.5))
          fn <- paste0(
            canceri_path, 
            "/survival_figure_", 
            pred_name, 
            figure_str, 
            ".pdf"
          )
          pdf(fn, width = plot_width / 25.4 * 1, height = plot_height / 25.4 * 1.5)
          grid::grid.draw(combined_plot)
          dev.off()
          fn <- paste0(
            canceri_path, 
            "/tcga_hazards_", 
            pred_name, 
            figure_str, 
            ".pdf"
          )
          pdf(fn, width = plot_width / 25.4 * 0.5, height = plot_height / 25.4 * 0.75)
          print(hr_list[[1]] + hr_add_theme)
          dev.off()
          fn <- paste0(
            canceri_path, 
            "/scanb_hazards_", 
            pred_name, 
            figure_str, 
            ".pdf"
          )
          pdf(fn, width = plot_width / 25.4 * 0.5, height = plot_height / 25.4 * 0.75)
          print(hr_list[[2]] + hr_add_theme)
          dev.off()
        }
      }
    }
    
    # Generate summary tables
    library(dplyr)
    cox_results_df <- bind_rows(cox_results)
    km_results_df <- bind_rows(km_results)
    concordance_results_df <- bind_rows(concordance_results)
    codeae_cox_results_df <- bind_rows(codeae_cox_results)
    codeae_km_results_df <- bind_rows(codeae_km_results)
    codeae_concordance_results_df <- bind_rows(codeae_concordance_results)
    baseline_cox_results_df <- bind_rows(baseline_cox_results)
    baseline_km_results_df <- bind_rows(baseline_km_results)
    baseline_concordance_results_df <- bind_rows(baseline_concordance_results)
    treated_pairs_n <- bind_rows(treated_pairs_n)
    
    # How many drugs?
    length(tcga_drugs[!is.na(tcga_drug_in_ctrp)])
    sum(apply(cancer_drug_counts[,tolower(tcga_drugs[!is.na(tcga_drug_in_ctrp)])] > 5,2,any))
    length(unique(treated_pairs_n$drug))
    
    cox_cleaner <- function(x) {
      out <- x |> filter(var == "sensitivity")
      sel_cols <- c("drug")
      if ("drug_codeae" %in% colnames(x)) {
        sel_cols <- c(sel_cols, "drug_codeae")
      }
      if ("baseline_model" %in% colnames(x)) {
        sel_cols <- c(sel_cols, "baseline_model")
      }
      sel_cols <- c(
        sel_cols, 
        "cancer", 
        "survival", 
        "coef", 
        "se_coef", 
        "p"
      )
      out <- out[, sel_cols]
      colnames(out)[ncol(out) - (2:0)] <- c(
        "cox_coef", 
        "cox_coef_se", 
        "cox_coef_p"
      )
      return(out)
    }
    km_cleaner <- function(x) {
      sel_cols <- c("drug")
      if ("drug_codeae" %in% colnames(x)) {
        sel_cols <- c(sel_cols, "drug_codeae")
      }
      if ("baseline_model" %in% colnames(x)) {
        sel_cols <- c(sel_cols, "baseline_model")
      }
      sel_cols <- c(sel_cols, "cancer", "survival", "p")
      out <- x[, sel_cols]
      colnames(out)[ncol(out)] <- "km_logrank_p"
      return(out)
    }
    
    # MODAE tables
    survival_res_df <- dplyr::full_join(
      cox_cleaner(cox_results_df), 
      km_cleaner(km_results_df), 
      by = c("drug", "cancer", "survival")
    )
    survival_res_df <- dplyr::full_join(
      survival_res_df, 
      concordance_results_df, 
      by = c("drug", "cancer", "survival")
    )
    survival_res_df[["method"]] <- "MODAE"
    
    # CODE-AE tables
    codeae_survival_res_df <- dplyr::full_join(
      cox_cleaner(codeae_cox_results_df), 
      km_cleaner(codeae_km_results_df), 
      by = c("drug", "drug_codeae", "cancer", "survival")
    )
    codeae_survival_res_df <- dplyr::full_join(
      codeae_survival_res_df, 
      codeae_concordance_results_df, 
      by = c("drug", "drug_codeae", "cancer", "survival")
    )
    codeae_survival_res_df[["method"]] <- "CODE-AE"
    
    # Baseline tables
    baseline_survival_res_df <- dplyr::full_join(
      cox_cleaner(baseline_cox_results_df), 
      km_cleaner(baseline_km_results_df), 
      by = c("drug", "baseline_model", "cancer", "survival")
    )
    baseline_survival_res_df <- dplyr::full_join(
      baseline_survival_res_df, 
      baseline_concordance_results_df, 
      by = c("drug", "baseline_model", "cancer", "survival")
    )
    baseline_survival_res_df[["method"]] <- baseline_survival_res_df[["baseline_model"]]
    baseline_survival_res_df[["baseline_model"]] <- NULL
    
    # Check patient overlap between CODE-AE predictions and our predictions
    mean(
      substr(rownames(tcga_predictions), 1, 15) %in% 
        codeae_predictions[["Sample"]]
    )
    
    # Check patient overlaps by cancer
    fn <- paste0(prediction_path, "internal/unified_patient_output_table.csv.gz")
    tcga_final_output_table <- readr::read_csv(fn)
    table(
      substr(tcga_final_output_table[["id"]], 1, 15) %in% 
        codeae_predictions[["Sample"]],
      tcga_final_output_table[["type"]], 
      useNA = "always"
    )
    
    # Which Sorafenib prediction performs better for CODE-AE?
    codeae_survival_res_df |> filter(grepl("Sorafenib", drug_codeae))
    
    # Comprehensive comparison table
    full_surv_comp_table_long <- bind_rows(
      survival_res_df |> filter(drug != "survival_risk"),
      codeae_survival_res_df |>
        filter(drug_codeae != "Sorafenib...51") |>
        select(-drug_codeae),
      baseline_survival_res_df
    )
    fn <- paste0(plot_path, "treatment_survival_comparison_long.csv")
    readr::write_csv(full_surv_comp_table_long, fn)
    full_surv_comp_table_wide <- tidyr::pivot_wider(
      full_surv_comp_table_long,
      id_cols = c("cancer", "survival", "drug"),
      names_from = "method",
      values_from = c("cox_coef", "cox_coef_p", "km_logrank_p", "concordance", "n")
    )
    fn <- paste0(plot_path, "treatment_survival_comparison_wide.csv")
    readr::write_csv(full_surv_comp_table_wide, fn)
    
    # Compact summary
    library(magrittr)
    library(dplyr)
    full_surv_comp_table_long %>%
      group_by(method, survival) %>%
      summarise(
        ci_mean=mean(ifelse(is.na(concordance), 0.5, concordance)),
        km_sig_rate=mean(ifelse(is.na(km_logrank_p), 1, km_logrank_p) < 0.05),
        cox_sig_rate=mean(ifelse(is.na(cox_coef_p), 1, cox_coef_p) < 0.05),
        n=n()
      )
    
    shared_obs_concordance <- plyr::ddply(
      full_surv_comp_table_wide,
      c("drug", "survival"),
      function(x) {
        shared_obs <- with(x, !(is.na(concordance_MODAE) | is.na(`concordance_CODE-AE`)))
        if (!any(shared_obs)) return(data.frame())
        out <- with(
          x,
          data.frame(
            MODAE_mean = mean(concordance_MODAE[shared_obs]),
            MODAE_sd = sd(concordance_MODAE[shared_obs]),
            `CODE-AE_mean` = mean(`concordance_CODE-AE`[shared_obs]),
            `CODE-AE_sd` = sd(`concordance_CODE-AE`[shared_obs])
          )
        )
        return(out)
      }
    )
    fn <- paste0(plot_path, "treatment_survival_comparison_shared_obs_concordance.csv")
    readr::write_csv(shared_obs_concordance, fn)
    
    shared_obs_km <- plyr::ddply(
      full_surv_comp_table_wide,
      c("drug", "survival"),
      function(x) {
        shared_obs <- with(x, !(is.na(km_logrank_p_MODAE) | is.na(`km_logrank_p_CODE-AE`)))
        if (!any(shared_obs)) return(data.frame())
        out <- with(
          x,
          data.frame(
            MODAE_rate = mean(km_logrank_p_MODAE[shared_obs] < 0.05),
            `CODE-AE_rate` = mean(`km_logrank_p_CODE-AE`[shared_obs] < 0.05)
          )
        )
        return(out)
      }
    )
    fn <- paste0(plot_path, "treatment_survival_comparison_shared_obs_km.csv")
    readr::write_csv(shared_obs_km, fn)
    
    shared_obs_cox <- plyr::ddply(
      full_surv_comp_table_wide,
      c("drug", "survival"),
      function(x) {
        shared_obs <- with(x, !(is.na(cox_coef_p_MODAE) | is.na(`cox_coef_p_CODE-AE`)))
        if (!any(shared_obs)) return(data.frame())
        out <- with(
          x,
          data.frame(
            MODAE_rate = mean(cox_coef_p_MODAE[shared_obs] < 0.05),
            `CODE-AE_rate` = mean(`cox_coef_p_CODE-AE`[shared_obs] < 0.05)
          )
        )
        return(out)
      }
    )
    fn <- paste0(plot_path, "treatment_survival_comparison_shared_obs_cox.csv")
    readr::write_csv(shared_obs_cox, fn)
    
    # Survival across cancers
    # TODO: implement
    
    # Test Cox PH assumptions
    pi <- tcga_predictions
    
    tcga_event <- tcga_surv_data[tcga_surv_ptr, "OS"] == 1
    tcga_time <- tcga_surv_data[tcga_surv_ptr, "OS.time"]
    
    tcga_event[tcga_time > 3000] <- FALSE
    tcga_time[tcga_time > 3000] <- 3000
    
    hazard_cols <- colnames(pi)[grep("survival_hazard_[0-9]+", colnames(pi))]
    hazards <- pi[tcga_surv_ptr, hazard_cols]
    
    hazard_formula <- as.formula(
      paste0(
        "survival::Surv(tcga_time, tcga_event) ~ ", 
        paste(hazard_cols, collapse = " + ")))
    cox_model <- survival::coxph(hazard_formula, data = hazards)
    test_ph <- survival::cox.zph(cox_model)
    
    fn <- paste0(plot_path, "tcga_survival_assumptions", figure_str, ".pdf")
    pdf(fn, width = plot_width / 25.4 * 4, height = plot_height / 25.4 * 6)
    print(survminer::ggcoxzph(test_ph))
    dev.off()
    
    if (tcga_brca || scanb) {
      pi <- scanb_predictions
      
      scanb_event <- scanb_surv_data[scanb_surv_ptr, "OS_event"] == 1
      scanb_time <- scanb_surv_data[scanb_surv_ptr, "OS_days"]
      
      scanb_event[scanb_time > 3000] <- FALSE
      scanb_time[scanb_time > 3000] <- 3000
      
      hazard_cols <- colnames(pi)[grep("survival_hazard_[0-9]+", colnames(pi))]
      hazards <- pi[scanb_surv_ptr, hazard_cols]
      
      hazard_formula <- as.formula(
        paste0(
          "survival::Surv(scanb_time, scanb_event) ~ ", 
          paste(hazard_cols, collapse = " + ")))
      cox_model <- survival::coxph(hazard_formula, data = hazards)
      test_ph <- survival::cox.zph(cox_model)
      
      fn <- paste0(plot_path, "scanb_survival_assumptions", figure_str, ".pdf")
      pdf(fn, width = plot_width / 25.4 * 4, height = plot_height / 25.4 * 6)
      print(survminer::ggcoxzph(test_ph))
      dev.off()
    }
  }
}