source("../setup.R")

for (save_dir in save_dirs) {
  if(exists("var_list")) try(detach(var_list))
  var_list <- get_paths_and_parameters(save_dir, remove_incomplete = FALSE)
  attach(var_list)
  #table(gsub("^.*?task|\\.csv\\.gz$|\\.csv$", "", incomplete_files))
  
  ps_files <- all_files[grep(paste0(model_name, "_paramsearch|_ps_cv"), all_files)]
  t_files <- all_files[grep(paste0(model_name, "_testing|_test_cv"), all_files)]
  
  metrics_list <- list()
  for (file_list in list(ps_files, t_files)) {
    metrics_files <- file_list[grep("metrics", file_list)]
    
    metrics <- lapply(metrics_files, read_result_files, path = save_path, model_name = model_name)
    metrics <- Reduce(COPS::rbind_fill, metrics)
    
    # Fix early bugs
    if (!"combined_train_batch_loss" %in% colnames(metrics))
      metrics[["combined_train_batch_loss"]] <- as.numeric(gsub("^tf.Tensor\\(|,.*$", "", metrics[["cl_test_batch_loss"]]))
    if (!"combined_train_batch_dsc" %in% colnames(metrics))
      metrics[["combined_train_batch_dsc"]] <- as.numeric(gsub("^tf.Tensor\\(|,.*$", "", metrics[["cl_test_batch_dsc"]]))
    if (!"combined_test_batch_loss" %in% colnames(metrics))
      metrics[["combined_test_batch_loss"]] <- as.numeric(gsub("^tf.Tensor\\(|,.*$", "", metrics[["cl_test_batch_loss.1"]]))
    if (!"combined_test_batch_dsc" %in% colnames(metrics))
      metrics[["combined_test_batch_dsc"]] <- as.numeric(gsub("^tf.Tensor\\(|,.*$", "", metrics[["cl_test_batch_dsc.1"]]))
    
    metrics_list <- c(metrics_list, list(metrics))
  }
  names(metrics_list) <- c("ps", "test")
  
  additional_metrics_list <- list()
  stnr_cut <- 100
  signal_to_noise_files <- all_files[grep("neighborhood_signal_to_noise", all_files)]
  mean_diff_files <- all_files[grep("neighborhood_mean_diff", all_files)]
  predicted_signal_to_noise_files <- signal_to_noise_files[grep("^predicted_", signal_to_noise_files)]
  predicted_mean_diff_files <- mean_diff_files[grep("^predicted_", mean_diff_files)]
  signal_to_noise_files <- signal_to_noise_files[!signal_to_noise_files %in% predicted_signal_to_noise_files]
  mean_diff_files <- mean_diff_files[!mean_diff_files %in% predicted_mean_diff_files]
  metrics_from_raw_data <- list()
  if (length(signal_to_noise_files) > 0) metrics_from_raw_data[["stnr"]] <- signal_to_noise_files
  if (length(mean_diff_files) > 0) metrics_from_raw_data[["mdiff"]] <- mean_diff_files
  for (file_listi in names(metrics_from_raw_data)) {
    processed_listi <- list()
    raw_listi <- lapply(metrics_from_raw_data[[file_listi]], read_result_files, path = save_path, model_name = model_name)
    # Treat patient and CL metrics differently, because drug STNR should be summarised
    settings <- lapply(
      raw_listi, 
      function(x) list(
        summarise = x[1, "sample_type"] == "cl")) # Should all be same
    for (i in 1:length(raw_listi)) {
      # TODO: update for new format, done?
      metric_df <- reshape2::dcast(
        raw_listi[[i]], task + run + fold + sample_type ~ col, 
        value.var = file_listi)
      metric_cols <- colnames(metric_df)[!colnames(metric_df) %in% c("task", "run", "fold", "sample_type")]
      if (settings[[i]][["summarise"]]) {
        metric_mat <- as.matrix(metric_df[,metric_cols])
        metric_quantiles <- apply(metric_mat, 1, quantile, probs = seq(0,1,0.1), na.rm = TRUE)
        metric_quantiles <- as.data.frame(t(metric_quantiles))
        metric_cut_mean <- apply(metric_mat > stnr_cut, 1, mean, na.rm = TRUE)
        metric_mat[!is.finite(metric_mat)] <- NA
        metric_mean <- apply(metric_mat, 1, mean, na.rm = TRUE)
        colnames(metric_quantiles) <- paste0("dr_", file_listi, "_Q", gsub("\\%", "", colnames(metric_quantiles)))
        processed_listi[[i]] <- metric_df[, c("task", "run", "fold", "sample_type")]
        processed_listi[[i]][[paste0("dr_", file_listi, "_mean")]] <- metric_mean
        processed_listi[[i]][[paste0("dr_", file_listi, "_cut_mean")]] <- metric_cut_mean
        processed_listi[[i]] <- data.frame(processed_listi[[i]], metric_quantiles)
      } else {
        colnames(metric_df)[colnames(metric_df) %in% metric_cols] <- paste0(metric_cols, "_", file_listi) # same order by definition
        processed_listi[[i]] <- metric_df
      }
    }
    additional_metrics_list[[file_listi]] <- Reduce(outer_join, processed_listi)
  }
  
  metrics_from_predictions <- list()
  if (length(signal_to_noise_files) > 0) metrics_from_predictions[["stnr"]] <- predicted_signal_to_noise_files
  if (length(predicted_mean_diff_files) > 0) metrics_from_predictions[["mdiff"]] <- predicted_mean_diff_files
  for (file_listi in names(metrics_from_predictions)) {
    processed_listi <- list()
    raw_listi <- lapply(metrics_from_predictions[[file_listi]], read_result_files, path = save_path, model_name = model_name)
    # Treat patient and CL metrics differently, because drug STNR should be summarised
    settings <- lapply(
      raw_listi, 
      function(x) list(
        summarise = x[1, "sample_type"] == "cl")) # Should all be same
    for (i in 1:length(raw_listi)) {
      metric_cols <- colnames(raw_listi[[i]])[grepl("_[0-9]+$", colnames(raw_listi[[i]]))]
      # Check if metric should be summarised (survival should not be summarised)
      if (settings[[i]][["summarise"]]) {
        metric_mat <- as.matrix(raw_listi[[i]][,metric_cols])
        metric_quantiles <- apply(metric_mat, 1, quantile, probs = seq(0,1,0.1), na.rm = TRUE)
        metric_quantiles <- as.data.frame(t(metric_quantiles))
        metric_cut_mean <- apply(metric_mat > stnr_cut, 1, mean, na.rm = TRUE)
        metric_mat[!is.finite(metric_mat)] <- NA
        metric_mean <- apply(metric_mat, 1, mean, na.rm = TRUE)
        colnames(metric_quantiles) <- paste0("dr_pred_", file_listi, "_Q", gsub("\\%", "", colnames(metric_quantiles)))
        processed_listi[[i]] <- raw_listi[[i]][, c("task", "run", "fold", "sample_type")]
        processed_listi[[i]][[gsub("_[0-9]+$", paste0("_", file_listi, "_mean"), metric_cols[1])]] <- metric_mean
        processed_listi[[i]][[gsub("_[0-9]+$", paste0("_", file_listi, "_cut_mean"), metric_cols[1])]] <- metric_cut_mean
        processed_listi[[i]] <- data.frame(processed_listi[[i]], metric_quantiles)
      } else {
        metric_df <- raw_listi[[i]][, c("task", "run", "fold", "sample_type", metric_cols)]
        colnames(metric_df)[colnames(metric_df) %in% metric_cols] <- paste0(metric_cols, "_", file_listi)
        processed_listi[[i]] <- metric_df
      }
    }
    additional_metrics_list[[paste0("predicted_", file_listi)]] <- Reduce(outer_join, processed_listi)
  }
  
  silhouette_files <- all_files[grep("silhouette", all_files)]
  #predicted_silhouette_files <- silhouette_files[grep("predicted", silhouette_files)]
  #silhouette_files <- silhouette_files[silhouette_files %in% predicted_silhouette_files]
  if (length(silhouette_files) > 0) {
    silhouette_list <- lapply(silhouette_files, read_result_files, path = save_path, model_name = model_name)
    silhouette_list_processed <- list()
    # Treat patient and CL metrics differently, because drug STNR should be summarised
    settings <- lapply(
      silhouette_list, 
      function(x) list(
        summarise = x[1, "sample_type"] == "cl")) # Should all be same
    for (i in 1:length(silhouette_list)) {
      value_name <- ifelse(grepl("predicted", silhouette_files[i]), "predicted_silhouette", "silhouette")
      metric_df <- reshape2::dcast(
        silhouette_list[[i]], task + run + fold + sample_type ~ col, 
        value.var = "silhouette")
      metric_cols <- colnames(metric_df)[!colnames(metric_df) %in% c("task", "run", "fold", "sample_type")]
      if (settings[[i]][["summarise"]]) {
        metric_mat <- as.matrix(metric_df[,metric_cols])
        metric_quantiles <- apply(metric_mat, 1, quantile, probs = seq(0,1,0.1), na.rm = TRUE)
        metric_quantiles <- as.data.frame(t(metric_quantiles))
        metric_mat[!is.finite(metric_mat)] <- NA
        metric_mean <- apply(metric_mat, 1, mean, na.rm = TRUE)
        colnames(metric_quantiles) <- paste0("dr_", value_name, "_Q", gsub("\\%", "", colnames(metric_quantiles)))
        silhouette_list_processed[[i]] <- metric_df[, c("task", "run", "fold", "sample_type")]
        silhouette_list_processed[[i]][[paste0("dr_", value_name, "_mean")]] <- metric_mean
        silhouette_list_processed[[i]] <- data.frame(silhouette_list_processed[[i]], metric_quantiles)
      } else {
        colnames(metric_df)[colnames(metric_df) %in% metric_cols] <- paste0(metric_cols, "_", value_name) # same order by definition
        silhouette_list_processed[[i]] <- metric_df
      }
    }
    additional_metrics_list[["silhouette"]] <- Reduce(outer_join, silhouette_list_processed)
  }
  
  drug_response_r2_processed <- NULL
  drug_response_r2_files <- all_files[grep("drug_response_r2.csv.gz$",  all_files)]
  #if ("drug_response_r2.csv.gz" %in% all_files) {
  if (length(drug_response_r2_files) > 0) {
    settings <- list(
      label = ifelse(grepl("^ps_", drug_response_r2_files), "parameter_search_", ""), 
      label2 = ifelse(grepl("^adjusted_", drug_response_r2_files), "adjusted_", ""), 
      summarise = !grepl("overall_", drug_response_r2_files)
    )
    drug_response_r2_processed <- list()
    for (i in 1:length(drug_response_r2_files)) {
      drug_response_r2 <- read_result_files(drug_response_r2_files[i], path = save_path, model_name = model_name)
      if (settings[["summarise"]][i]) {
        metric_df <- reshape2::dcast(
          drug_response_r2, task + run + fold + dataset ~ drug, 
          value.var = "dr_r2")
        #hist(drug_response_r2[drug_response_r2$dataset == "cl_train", "dr_r2"])
        out <- list()
        r2_list_raw <- split(metric_df, f = metric_df$dataset)
        for(dati in names(r2_list_raw)) {
          metric_cols <- colnames(metric_df)[!colnames(metric_df) %in% c("task", "run", "fold", "dataset")]
          metric_mat_raw <- as.matrix(r2_list_raw[[dati]][,metric_cols])
          metric_mat <- matrix(NA, nrow = nrow(metric_mat_raw), ncol = ncol(metric_mat_raw))
          metric_mat[] <- as.numeric(gsub("\\[|\\]", "", metric_mat_raw))
          metric_quantiles <- apply(metric_mat, 1, quantile, probs = seq(0,1,0.1), na.rm = TRUE)
          metric_quantiles <- as.data.frame(t(metric_quantiles))
          metric_mat[!is.finite(metric_mat)] <- NA
          metric_mat[metric_mat < 0] <- 0 # ignore sub zero R^2
          metric_mean <- apply(metric_mat, 1, mean, na.rm = TRUE)
          colnames(metric_quantiles) <- paste0(
            "dr_", dati, "_", settings[["label"]][i], settings[["label2"]][i], "R2_Q", 
            gsub("\\%", "", colnames(metric_quantiles)))
          
          out[[dati]] <- r2_list_raw[[dati]][, c("task", "run", "fold")]
          out[[dati]][[paste0("dr_", dati, "_", settings[["label"]][i], settings[["label2"]][i], "R2_mean")]] <- metric_mean
          out[[dati]] <- data.frame(out[[dati]], metric_quantiles)
        }
        drug_response_r2_processed[[i]] <- Reduce(outer_join, out)
      } else {
        out <- list()
        r2_list_raw <- split(drug_response_r2, f = drug_response_r2$dataset)
        for(dati in names(r2_list_raw)) {
          out[[dati]] <- r2_list_raw[[dati]][,c("task", "run", "fold")]
          out[[dati]][[paste0("dr_", dati, "_", settings[["label"]][i], settings[["label2"]][i], "overall_R2")]] <- r2_list_raw[[dati]][["dr_r2"]]
        }
        drug_response_r2_processed[[i]] <- Reduce(outer_join, out)
      }
    }
    drug_response_r2_processed <- Reduce(outer_join, drug_response_r2_processed)
    drug_response_r2_processed[["sample_type"]] <- "cl"
    additional_metrics_list[["r2"]] <- drug_response_r2_processed
  }
  batch_prediction <- NULL
  batch_prediction_files <- all_files[grep("batch_predictions\\.csv$", all_files)]
  if (length(batch_prediction_files) > 0) {
    batch_prediction <- list()
    for (i in 1:length(batch_prediction_files)) {
      name_first <- strsplit(batch_prediction_files[i], split = "_")[[1]][1]
      if (name_first == "batch") {
        name_first <- ""
      } else {
        name_first <- paste0(name_first, "_")
      }
      batch_prediction_i <- read_result_files(
        batch_prediction_files[i], path = save_path, model_name = model_name)
      batch_prediction_i[["X"]] <- NULL
      res_cols <- colnames(batch_prediction_i)
      res_ind <- !res_cols %in% c("task", "run", "fold", "name")
      res_cols <- res_cols[res_ind]
      colnames(batch_prediction_i)[res_ind] <- paste0(
        name_first, "batch_prediction_", colnames(batch_prediction_i)[res_ind])
      
      batch_prediction[[i]] <- batch_prediction_i
    }
    batch_prediction <- Reduce(outer_join, batch_prediction)
    #plot(batch_prediction[["batch_prediction_dsc"]], batch_prediction[["batch_prediction_svm_auc"]])
    #plot(batch_prediction[["batch_prediction_svm_auc"]], batch_prediction[["batch_prediction_rfc_auc"]])
  }
  tissue_classifier_processed <- NULL
  tissue_classifier_files <- all_files[grep("tissue_classifier_performance.csv.gz$",  all_files)]
  if (length(tissue_classifier_files) > 0) {
    label <- ifelse(grepl("^ps_", tissue_classifier_files), "parameter_search_", "")
    tissue_classifier_processed <- list()
    for (i in 1:length(tissue_classifier_files)) {
      tissue_classifier <- read_result_files(tissue_classifier_files[i], path = save_path, model_name = model_name)
      out <- list()
      out_list_raw <- split(tissue_classifier, f = tissue_classifier$dataset)
      for(dati in names(out_list_raw)) {
        out[[dati]] <- out_list_raw[[dati]][,c("task", "run", "fold")]
        perf_cols <- intersect(colnames(out_list_raw[[dati]]), c("apr", "f1", "bacc", "acc", "auroc"))
        out[[dati]][paste0("tissue_class_", dati, "_", label[i], perf_cols)] <- out_list_raw[[dati]][perf_cols]
      }
      tissue_classifier_processed[[i]] <- Reduce(outer_join, out)
    }
    tissue_classifier_processed <- Reduce(outer_join, tissue_classifier_processed)
    #additional_metrics_list[["tissue_class_perf"]] <- tissue_classifier_processed
  }
  tissue_classifier_stability_processed <- NULL
  tissue_classifier_stability_files <- all_files[grep("tissue_classifier_stability_within_run.csv.gz$",  all_files)]
  if (length(tissue_classifier_stability_files) > 0) {
    label <- ifelse(grepl("^ps_", tissue_classifier_stability_files), "parameter_search_", "")
    tissue_classifier_stability_processed <- list()
    for (i in 1:length(tissue_classifier_stability_files)) {
      tissue_classifier_stability <- read_result_files(tissue_classifier_stability_files[i], path = save_path, model_name = model_name)
      out <- list()
      out_list_raw <- split(tissue_classifier_stability, f = tissue_classifier_stability$dataset)
      for(dati in names(out_list_raw)) {
        out[[dati]] <- out_list_raw[[dati]][,c("task", "run", "fold1", "fold2")]
        perf_cols <- intersect(colnames(out_list_raw[[dati]]), c("ari", "soft_max_mean_dot"))
        out[[dati]][paste0("tissue_class_", dati, "_", label[i], perf_cols)] <- out_list_raw[[dati]][perf_cols]
      }
      tissue_classifier_stability_processed[[i]] <- Reduce(outer_join, out)
    }
    tissue_classifier_stability_processed <- Reduce(outer_join, tissue_classifier_stability_processed)
    tissue_classifier_stability_processed <- plyr::ddply(
      tissue_classifier_stability_processed[,!colnames(tissue_classifier_stability_processed) %in% c("fold1", "fold2")], 
      c("task", "run"), 
      function(x) data.frame(lapply(x, mean))
    )
    #additional_metrics_list[["tissue_class_perf"]] <- tissue_classifier_stability_processed
  }
  
  if (length(additional_metrics_list) > 0 | 
      !is.null(drug_response_r2_processed) | 
      !is.null(batch_prediction) |
      !is.null(tissue_classifier_processed) |
      !is.null(tissue_classifier_stability_processed) ) {
    if (length(additional_metrics_list) > 0) {
      additional_metrics <- Reduce(outer_join, additional_metrics_list)
    } else {
      additional_metrics <- NULL
    }
    
    final_metrics <- metrics_list[["test"]]
    final_metrics <- final_metrics[final_metrics[["stage"]] == "final",]
    
    final_metrics_ps <- metrics_list[["ps"]]
    final_metrics_ps <- final_metrics_ps[final_metrics_ps[["stage"]] == "final",]
    colnames(final_metrics_ps) <- paste0("ps_", colnames(final_metrics_ps))
    final_metrics_ps[["stage"]] <- final_metrics_ps[["ps_stage"]]
    final_metrics_ps[["task"]] <- final_metrics_ps[["ps_task"]]
    final_metrics_ps[["run"]] <- final_metrics_ps[["ps_run"]]
    final_metrics_ps[["fold"]] <- final_metrics_ps[["ps_fold"]]
    final_metrics_ps[["name"]] <- final_metrics_ps[["ps_name"]]
    final_metrics_ps[["ps_stage"]] <- final_metrics_ps[["ps_task"]] <- final_metrics_ps[["ps_run"]] <- final_metrics_ps[["ps_fold"]] <- final_metrics_ps[["ps_name"]] <- NULL
    
    final_metrics <- plyr::join(final_metrics, final_metrics_ps)
    if (!is.null(drug_response_r2_processed)) {
      final_metrics <- plyr::join(final_metrics, drug_response_r2_processed)
    }
    if (!is.null(batch_prediction)) {
      final_metrics <- plyr::join(final_metrics, batch_prediction)
    }
    if (!is.null(tissue_classifier_processed)) {
      final_metrics <- plyr::join(final_metrics, tissue_classifier_processed)
    }
    if (!is.null(tissue_classifier_stability_processed)) {
      final_metrics <- plyr::join(final_metrics, tissue_classifier_stability_processed)
    }
    
    unnecessary_cols <- c("stage", "name") # already filtered
    unnecessary_cols <- c(unnecessary_cols, "run", "fold") # to be summarised
    mean_fun <- function(x) as.data.frame(lapply(x, function(y) if(all(is.na(y))) return(NA) else return(mean(y, na.rm = TRUE))))
    final_metrics_mean <- plyr::ddply(
      final_metrics[,!colnames(final_metrics) %in% unnecessary_cols], 
      c("task"), mean_fun)
    
    if (!is.null(additional_metrics)) {
      split_join_by_sample_type <- function(x) {
        type_ind <- which(colnames(x) == "sample_type")
        patient_ametrics <- x[x$sample_type == "patient", -type_ind]
        cl_ametrics <- x[x$sample_type == "cl", -type_ind]
        metric_ind <- which(!colnames(patient_ametrics) %in% c("task", "run", "fold"))
        colnames(patient_ametrics)[metric_ind] <- paste0("patient_", colnames(patient_ametrics)[metric_ind])
        colnames(cl_ametrics)[metric_ind] <- paste0("cl_", colnames(cl_ametrics)[metric_ind])
        
        additional_metrics_widened <- outer_join(patient_ametrics, cl_ametrics)
        nnacols <- !apply(is.na(additional_metrics_widened), 2, all)
        additional_metrics_widened <- additional_metrics_widened[, nnacols]
        return(additional_metrics_widened)
      }
      additional_metrics_widened <- split_join_by_sample_type(additional_metrics)
      
      colinds <- !colnames(additional_metrics_widened) %in% unnecessary_cols
      additional_metrics_widened_mean <- plyr::ddply(
        additional_metrics_widened[,colinds], 
        c("task"), mean_fun)
      
      metrics_table <- outer_join(additional_metrics_widened_mean, final_metrics_mean)
      metrics_table_full <- outer_join(additional_metrics_widened, final_metrics)
    } else {
      metrics_table <- final_metrics_mean
      metrics_table_full <- final_metrics
    }
    
    if (!inherits(parameters, "list")) {
      
    }
    parameter_performance_table <- outer_join(metrics_table, parameters)
    
    # Objective balance
    shrink_f <- function(x, a, b) {
      x <- ifelse(x < a, a, x)
      x <- ifelse(x > b, b, x)
      return(x)
    }
    rescale_metrics_f <- function(x) {
      for (i in 1:9) {
        for (j in c("dr_cl_train_R2_Q", 
                    "dr_cl_test_R2_Q", 
                    "dr_cl_train_adjusted_R2_Q", 
                    "dr_cl_test_adjusted_R2_Q", 
                    "dr_cl_train_parameter_search_R2_Q", 
                    "dr_cl_test_parameter_search_R2_Q")) {
          try({
            qi <- paste0(j, i, "0")
            dr_test_r2_clean <- shrink_f(x[[qi]], 0, 1)
            x[[paste0(qi, "_clean")]] <- dr_test_r2_clean
          })
        }
      }
      for (i in list(c("patient_", "_recon"), 
                     c("cl_", "_recon"), 
                     c("ps_patient_", "_recon"), 
                     c("ps_cl_", "_recon"))) {
        for (j in c("train", "test")) {
          unscaled_name <- paste0(i[1], j, i[2], "struction_loss")
          scaled_name <- paste0(i[1], j, i[2], "_nlogloss")
          x[[scaled_name]] <- -log10(x[[unscaled_name]])
        }
      }
      return(x)
    }
    parameter_performance_table <- rescale_metrics_f(parameter_performance_table)
    metrics_table_full <- rescale_metrics_f(metrics_table_full)
    
    #parameter_performance_table[["combined_test_batch_ndsc"]] <- -parameter_performance_table[["combined_test_batch_ndsc"]]
    
    objectives <- c(
      #"dr_cl_test_adjusted_R2_Q90",
      #"dr_cl_test_adjusted_R2_Q90_clean",
      "dr_cl_test_R2_Q90_clean",
      #"dr_cl_test_R2_mean", 
      "patient_test_surv_c",
      #"combined_test_batch_dsc",
      #"batch_prediction_svm_auc",
      "batch_prediction_rfc_bacc",
      "tissue_class_cl_test_bacc", 
      "patient_test_recon_nlogloss", 
      "cl_test_recon_nlogloss")
    objective_names <- c(
      "drug sens. R^2", 
      "survival C-index", 
      #"batch DSC", 
      #"batch SVM AUROC",
      #"batch RFC AUROC",
      "batch RFC BACC",
      "CL tissue BACC", 
      "PT rec. -log10MSE", 
      "CL rec. -log10MSE"
      )
    objective_comparators <- list(
      .Primitive(">"), 
      .Primitive(">"), 
      .Primitive("<"), 
      .Primitive(">"), 
      .Primitive(">"), 
      .Primitive(">"))
    
    objectives_in_table <- objectives %in% colnames(parameter_performance_table)
    objectives <- objectives[objectives_in_table]
    objective_names <- objective_names[objectives_in_table]
    objective_comparators <- objective_comparators[objectives_in_table]
    
    naind_list <- lapply(
      objectives, 
      function(i) is.na(parameter_performance_table[[i]]))
    nnaind <- !Reduce("|", naind_list)
    if (!any(nnaind)) {
      parameter_performance_table[is.na(parameter_performance_table)] <- -1
      metrics_table_full[is.na(metrics_table_full)] <- -1
      nnaind <- !nnaind
    }
    parameter_performance_table <- parameter_performance_table[nnaind,]
    
    parameter_performance_table[["pareto_front"]] <- COPS::pareto_fronts(
      parameter_performance_table, 
      objectives, 
      objective_comparators)
    
    pareto_table <- parameter_performance_table
    obj_col <- match(objectives, colnames(pareto_table))
    colnames(pareto_table)[obj_col] <- objective_names
    
    #png(paste0(plot_path, "objectives_plot.png"), width = plot_width * 2, height = plot_height * 2, res = plot_res, units = plot_units)
    #pdf(paste0(plot_path, "objectives_plot.pdf"), width = plot_width / 25.4 * 1.25, height = plot_height / 25.4 * 2.25)
    save_figure_safe(
      COPS::pareto_plot(
        pareto_table, 
        metrics = objective_names, 
        color_var = NULL, 
        shape_var = NULL, 
        size_var = NULL, 
        plot_pareto_front = TRUE, 
        front_color = "red", 
        metric_comparators = objective_comparators, 
        point_args = list(shape = "+", size = 4)), 
      pdf, 
      paste0(plot_path, "objectives_plot.pdf"), 
      width = plot_width / 25.4 * 1.25, 
      height = plot_height / 25.4 * 2.25
    )
    
    objective_weights <- c(1,1,-1,1,0.25,0.25)
    names(objective_weights) <- objectives
    weighted_objective_list <- lapply(
      objectives, 
      function(i) objective_weights[i] * parameter_performance_table[[i]])
    
    weighted_objective <- Reduce("+", weighted_objective_list)
    parameter_performance_table[["selected"]] <- rank(-weighted_objective) == 1
    parameter_performance_table[["weighted_objective"]] <- weighted_objective
    
    ps_objectives <- c(
      "dr_cl_test_parameter_search_R2_Q90_clean",
      #"dr_cl_test_parameter_search_R2_mean", 
      "ps_patient_test_surv_c",
      #"ps_combined_test_batch_dsc",
      #"ps_batch_prediction_svm_auc",
      "ps_batch_prediction_rfc_auc",
      "tissue_class_cl_test_parameter_search_bacc", 
      "ps_patient_test_recon_nlogloss", 
      "ps_cl_test_recon_nlogloss")
    ps_objectives <- ps_objectives[objectives_in_table]
    ps_objective_weights <- objective_weights
    names(ps_objective_weights) <- ps_objectives
    ps_weighted_objective_list <- lapply(
      ps_objectives, 
      function(i) ps_objective_weights[i] * parameter_performance_table[[i]])
    
    ps_weighted_objective <- Reduce("+", ps_weighted_objective_list)
    parameter_performance_table[["ps_weighted_objective"]] <- ps_weighted_objective
    
    plot(parameter_performance_table[["weighted_objective"]], parameter_performance_table[["ps_weighted_objective"]])
    
    # Predict weighted sum
    pp_table <- parameter_performance_table
    pp_table[["log10_alpha"]] <- log10(pp_table[["reg_a"]])
    pp_table[["log10_learning_rate"]] <- log10(pp_table[["learning_rate"]])
    pp_table[["log10_adversarial_learning_rate"]] <- log10(pp_table[["adversarial_learning_rate"]])
    pp_table[["log10_pre_learning_rate_ae"]] <- log10(pp_table[["pre_learning_rate_ae"]])
    pp_table[["log10_pre_learning_rate_sr"]] <- log10(pp_table[["pre_learning_rate_sr"]])
    pp_table[["log10_pre_learning_rate_cl"]] <- log10(pp_table[["pre_learning_rate_cl"]])
    pp_table[["log10_pre_learning_rate_dr"]] <- log10(pp_table[["pre_learning_rate_dr"]])
    pp_table[["log10_pre_learning_rate_bd"]] <- log10(pp_table[["pre_learning_rate_bd"]])
    pp_table[["log10_pre_learning_rate_bc"]] <- log10(pp_table[["pre_learning_rate_bc"]])
    pp_table[["log10_pre_learning_rate_bc"]] <- log10(pp_table[["pre_learning_rate_bc"]])
    pp_relevant <- c(
      paste0("encoder_layer", 1:2), "bottle_neck", paste0("decoder_layer", 1:2), 
      paste0("survival_model_layer", 1:2), #paste0("classifier_layer", 1:2), 
      paste0("drug_response_model_layer", 1:2), paste0("batch_adversarial_model_layer", 1:2), 
      paste0("dropout_", c("input", "autoencoder", "survival", "drug_response", "batch")), #"classifier", 
      paste0("log10_", c("", rep("pre_", 5)), "learning_rate", c("", "_ae", "_sr", "_dr", "_bd", "_bc")), #"_cl", 
      paste0(c("reconstruction", "survival", "drug", "batch"), "_weight"), #"classifier", 
      "log10_alpha", "deconfounder_norm_penalty", "deconfounder_layers_per_batch", 
      "batch_adversarial_gradient_penalty", "log10_adversarial_learning_rate",
      "weight_optimizer_args.weight_decay", "drug_response_model_drugwise_layers"
    )
    #colnames(pp_table)[grep("learning_rate", colnames(pp_table))]
    
    # Should probably work in python instead
    write.csv(pp_table[, pp_relevant[pp_relevant %in% colnames(pp_table)]], paste0(save_path, "objective_x.csv"))
    write.csv(pp_table[, c(objectives, "weighted_objective")], paste0(save_path, "objective_y.csv"))
    write.csv(pp_table[, c(ps_objectives, "ps_weighted_objective")], paste0(save_path, "ps_objective_y.csv"))
    
    plot(pp_table[["encoder_layer2"]], pp_table[["ps_combined_test_batch_dsc"]])
    
    # Parameter search test set performance
    # (prefix, postfix, sign, set indicator)
    performance_criteria <- list(
      list("patient_", "_recon_nlogloss", 0, TRUE), # 0.25, TRUE), 
      list("cl_", "_recon_nlogloss", 0, TRUE), # 0.25, TRUE),
      list("patient_", "_surv_c", 1, TRUE),
      list("dr_cl_", "_R2_Q90", 1, TRUE), 
      #list("dr_cl_", "_R2_mean", 1, TRUE), 
      list("tissue_class_cl_", "_bacc", 1, TRUE), 
      #list("batch_prediction_svm_auc", "", -1, FALSE), 
      #list("batch_prediction_rfc_auc", "", -1, FALSE))#-1, FALSE))
      list("batch_prediction_rfc_bacc", "", -1, FALSE))#-1, FALSE))
      #list("batch_prediction_svm_auc", "", -0.5, FALSE))#-1, FALSE))
      #list("combined_", "_batch_dsc", -1))
    ps_col_ind <- grep("^ps_|_parameter_search_", colnames(metrics_table_full))
    metrics_list_copy <- list()
    metrics_list_copy[["ps"]] <- metrics_table_full[,ps_col_ind]
    colnames(metrics_list_copy[["ps"]]) <- gsub("^ps_|_parameter_search", "", colnames(metrics_list_copy[["ps"]]))
    metrics_list_copy[["ps"]][["run"]] <- metrics_table_full[["run"]]
    metrics_list_copy[["ps"]][["fold"]] <- metrics_table_full[["fold"]]
    metrics_list_copy[["ps"]][["task"]] <- metrics_table_full[["task"]]
    metrics_list_copy[["ps"]][["stage"]] <- metrics_table_full[["stage"]] # NULL
    metrics_list_copy[["ps"]][["name"]] <- metrics_table_full[["name"]] # NULL
    #ps_col_ind <- grep("^ps_|_parameter_search_", colnames(parameter_performance_table))
    metrics_list_copy[["test"]] <- metrics_table_full[,-ps_col_ind]
  } else {
    # Parameter search test set performance
    # (prefix, postfix, sign)
    performance_criteria <- list(list("patient_", "_reconstruction_loss", -0.25, TRUE), 
                                 list("cl_", "_reconstruction_loss", -0.25, TRUE),
                                 list("patient_", "_surv_c", 1, TRUE),
                                 list("cl_", "_drug_response_mse", -1, TRUE),
                                 list("combined_", "_batch_dsc", -1, TRUE))
    metrics_list_copy <- metrics_list
  }
  
  # need to identify best parameter set (task) for each run and fold
  runs <- unique(metrics_list_copy[["ps"]][["run"]])
  folds <- unique(metrics_list_copy[["ps"]][["fold"]])
  task_perf <- function(
    x, 
    pc = performance_criteria, 
    stage = "final", 
    train_set = "train", 
    test_set = "test", 
    scale_metrics = FALSE
  ) {
    x_stage <- x[["stage"]]
    if (!is.null(x_stage)) x <- x[x_stage == stage,]
    obj_train_list <- lapply(
      pc, 
      function(y) {
        scale(
          x[[paste0(y[[1]], train_set[y[[4]]], y[[2]])]], 
          center = scale_metrics, 
          scale = scale_metrics
        ) * y[[3]]
      }
    )
    obj_test_list <- lapply(
      pc, 
      function(y) {
        scale(
          x[[paste0(y[[1]], test_set[y[[4]]], y[[2]])]], 
          center = scale_metrics, 
          scale = scale_metrics
        ) * y[[3]]
      }
    )
    obj_train_list <- obj_train_list[sapply(obj_train_list, length) > 0]
    #obj_train_list <- obj_train_list[sapply(obj_train_list, function(x) sum(!is.na(x))) > 0]
    obj_test_list <- obj_test_list[sapply(obj_test_list, length) > 0]
    #obj_test_list <- obj_test_list[sapply(obj_test_list, function(x) sum(!is.na(x))) > 0]
    x$train_objective <- Reduce(.Primitive("+"), obj_train_list)
    x$test_objective <- Reduce(.Primitive("+"), obj_test_list)
    return(x)
  }
  best_task_f <- function(x, ...) {
    x <- task_perf(x, ...)
    return(x[which.max(x$test_objective), c("task", "train_objective", "test_objective")])
  }
  
  best_by <- c("fold", "run", "stage")
  
  if (TRUE) {
    best_task_ps <- plyr::ddply(
      metrics_list_copy[["ps"]], 
      best_by[best_by %in% colnames(metrics_list_copy[["ps"]])], 
      best_task_f)
  } else { # Before random search
    task_number <- 1 # 1 = wasserstein, 2 = cross-entropy
    best_task_ps <- metrics_list_copy[["ps"]]
    best_task_ps <- best_task_ps[best_task_ps[["task"]] == task_number & 
                                   best_task_ps[["stage"]] == "final", 
                                 c("fold", "run", "task")]
  }
  if ("stage" %in% colnames(best_task_ps)) {
    best_task_ps <- best_task_ps[best_task_ps[["stage"]] == "final",]
  }
  best_by_ps <- best_by[best_by %in% colnames(best_task_ps)]
  best_task_test <- plyr::join(best_task_ps[, c(best_by_ps, "task"), drop = FALSE], 
                               metrics_list_copy[["test"]], #[metrics_list_copy[["test"]][["stage"]] == "final",], 
                               by = c(best_by_ps, "task"))
  best_task_ps <- plyr::join(best_task_ps[, c(best_by_ps, "task"), drop = FALSE], 
                             metrics_list_copy[["ps"]], #[metrics_list_copy[["ps"]][["stage"]] == "final",], 
                             by = c(best_by_ps, "task"))
  best_task_test[["stage"]] <- NULL
  best_task_ps[["stage"]] <- NULL
  best_by_ps <- best_by_ps[best_by_ps != "stage"]
  
  write.csv(best_task_ps, paste0(plot_path, "best_task_ps_performance.csv"))
  write.csv(best_task_test, paste0(plot_path, "best_task_test_performance.csv"))
  
  test_criteria_names <- sapply(performance_criteria, function(x) paste0(x[[1]], "test"[x[[4]]], x[[2]]))
  train_criteria_names <- gsub("_test_", "_train_", test_criteria_names)
  
  test_perf <- sapply(test_criteria_names, function(x) mean(best_task_test[[x]], na.rm = TRUE))
  retrain_perf <- sapply(train_criteria_names, function(x) mean(best_task_test[[x]], na.rm = TRUE))
  valid_perf <- sapply(test_criteria_names, function(x) mean(best_task_ps[[x]], na.rm = TRUE))
  train_perf <- sapply(train_criteria_names, function(x) mean(best_task_ps[[x]], na.rm = TRUE))
  
  metric_names <- c("P_NLOGMSE", "CL_NLOGMSE", "SURV_C", "DR_R2", "CL_T_BACC", "B_RFC_BACC")#"B_RFC_AUC")#"B_SVM_AUC")
  metric_names_f <- factor(metric_names, levels = metric_names[c(1:2,5:6,3:4)])
  
  temp <- Reduce("rbind", 
                 list(data.frame(set = "train", metric_name = metric_names_f, metric = train_perf),
                      data.frame(set = "valid", metric_name = metric_names_f, metric = valid_perf),
                      data.frame(set = "retrain",  metric_name = metric_names_f, metric = retrain_perf),
                      data.frame(set = "test",  metric_name = metric_names_f, metric = test_perf)))
  temp[["set"]] <- factor(as.character(temp[["set"]]), levels = c("train", "valid", "retrain", "test"))
  
  save_figure_safe(
    ggplot(temp, aes(x = metric_name, y = metric, fill = set)) + 
      geom_bar(stat = "identity", position = "dodge") + theme_bw() + 
      scale_fill_brewer(palette = "Dark2") +
      ggh4x::facet_wrap2(metric_name ~ ., scales = "free", ncol = 6) + 
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank()) + 
      ggtitle(model_name), 
    png, 
    paste0(plot_path, "metrics_plot.png"), 
    width = plot_width, 
    height = plot_height / 2, 
    res = plot_res, 
    units = plot_units
  )
  
  # Select best settings for external validation
  best_by_external <- best_by[!best_by %in% c("fold", "run")]
  best_external <- plyr::ddply(
    metrics_list_copy[["test"]], 
    best_by_external[best_by_external %in% colnames(metrics_list_copy[["test"]])], 
    task_perf)
  best_external[["survival_overfit"]] <- best_external[["patient_train_surv_c"]] - best_external[["patient_test_surv_c"]]
  best_external_mean <- plyr::ddply(
    best_external, 
    c("task", best_by_external), 
    function(x) data.frame(
      train_objective = mean(x[["train_objective"]], na.rm = TRUE), 
      test_objective = mean(x[["test_objective"]], na.rm = TRUE), 
      survival_overfit = mean(x[["survival_overfit"]], na.rm = TRUE))
  )
  external_perf_mean <- best_external_mean
  best_external_mean <- best_external_mean[which.max(best_external_mean[["test_objective"]]),]
  # Delete stage to properly join with diagnostics and embeddings
  best_external_mean[["stage"]] <- NULL
  best_by_external <- best_by_external[best_by_external != "stage"]
  # metrics_list_copy[["test"]][metrics_list_copy[["test"]]$task == best_external_mean$task,]
  writeLines(
    paste(best_external_mean[1,"task"]), # format as string
    paste0(save_path, "best_task.txt"))
  
  # Check correlatedness of drug predictions
  drug_correlation <- NULL
  drug_correlation_mean <- NULL
  drug_correlation_files <- all_files[grep("drug_response_cor\\.csv\\.gz$", all_files)]
  if (length(drug_correlation_files) > 0) {
    drug_correlation <- list()
    for (i in 1:length(drug_correlation_files)) {
      drug_correlation_i <- read_result_files(
        drug_correlation_files[i], path = save_path, model_name = model_name)
      drug_correlation_i[["X"]] <- NULL
      
      if (TRUE) {
        hist(drug_correlation_i[["dr_cor_train"]])
        hist(drug_correlation_i[["dr_cor_test"]])
        
        drug_correlation_best <- plyr::join(
          best_external_mean, 
          drug_correlation_i, 
          by = "task", 
          type = "left")
        
        fold_inds <- split(1:nrow(drug_correlation_best), f = drug_correlation_best[,c("run", "fold")])
        fold_inds <- fold_inds[1]
        fold_correlations <- list()
        for (fold_ind in fold_inds) {
          one_ind <- fold_ind#with(drug_correlation_best, run == 0 & fold == 0)
          
          if (FALSE) {
            drug_correlation_mat <- reshape2::acast(
              drug_correlation_best[one_ind,], 
              drug_i ~ drug_j, 
              value.var = "dr_cor_train")
            
            full_cor <- matrix(
              NA, 
              nrow = nrow(drug_correlation_mat)+1, 
              ncol = ncol(drug_correlation_mat)+1)
            #drug_hc_order <- hclust(as.dist(drug_correlation_mat))$order
            #apply(is.na(drug_correlation_mat), 1, sum)
            cols <- paste0("dr_pred_", 1:ncol(drug_correlation_mat))
            rows <- paste0("dr_pred_", 0:(nrow(drug_correlation_mat)-1))
            
            ordered_cor_mat <- drug_correlation_mat[rows, cols]
            
            full_cor[upper.tri(full_cor, diag = FALSE)] <- 
              ordered_cor_mat[upper.tri(ordered_cor_mat, diag = TRUE)]
            
            full_cor[lower.tri(full_cor, diag = FALSE)] <- 
              t(ordered_cor_mat)[lower.tri(ordered_cor_mat, diag = TRUE)]
            diag(full_cor) <- 1
          } else {
            rows <- unique(drug_correlation_best[["drug_i"]])
            cols <- unique(drug_correlation_best[["drug_j"]])
            
            if (length(rows) != length(cols)) stop("Non-square correlation matrix.")
            full_cor <- matrix(
              NA, 
              nrow = length(rows)+1, 
              ncol = length(cols)+1)
            diag(full_cor) <- 1
            rownames(full_cor) <- colnames(full_cor) <- paste0("dr_pred_", 1:ncol(full_cor)-1)
            
            for (dr_cor_rowi in one_ind) {
              drugi <- drug_correlation_best[dr_cor_rowi, "drug_i"]
              drugj <- drug_correlation_best[dr_cor_rowi, "drug_j"]
              corij <- drug_correlation_best[dr_cor_rowi, "dr_cor_train"]
              full_cor[drugi, drugj] <- corij
              full_cor[drugj, drugi] <- corij
            }
          }
          
          missing <- apply(is.na(full_cor), 1, mean) > 0.8
          set.seed(0)
          full_cor_impute <- impute::impute.knn(full_cor[!missing, !missing], k = 10)$data
          
          fold_correlations <- c(fold_correlations, list(full_cor_impute))
        }
        
        if (FALSE) {
          # Check difference in correlation between folds
          cor_diff <- fold_correlations[[2]] - fold_correlations[[5]]
          
          cor_diff_order <- hclust(
            dist(
              cor_diff, 
              method = "euclidean"
            ), 
            method = "complete"
          )$order
          
          #cor_diff <- fold_correlations[[4]] - fold_correlations[[5]]
          
          cor_diff_hm <- ComplexHeatmap::Heatmap(
            cor_diff[cor_diff_order,cor_diff_order], 
            name = "AAC cor diff", 
            cluster_rows = FALSE, 
            cluster_columns = FALSE, 
            show_column_names = FALSE, 
            show_row_names = FALSE
          )
          print(cor_diff_hm)
        }
        
        full_cor_impute <- fold_correlations[[1]]
        #cor_dist <- as.dist(1 - full_cor_impute)
        cor_dist <- dist(full_cor_impute, method = "euclidean")
        hc_order <- hclust(cor_dist, method = "complete")$order
        drug_order <- c(which(missing), which(!missing)[hc_order])
        
        cl_pred_hm <- ComplexHeatmap::Heatmap(
          full_cor[drug_order, drug_order], 
          name = "cl", 
          show_row_names = FALSE, 
          show_column_names = FALSE, 
          cluster_rows = FALSE, 
          cluster_columns = FALSE)
        
        #pdf(paste0(plot_path, "drug_prediction_correlations_example.pdf"), width = plot_width / 25.4 * 1.25, height = plot_height / 25.4 * 2.25)
        save_figure_safe(
          cl_pred_hm, 
          pdf, 
          paste0(plot_path, "drug_prediction_correlations_example.pdf"), 
          width = plot_width / 25.4 * 1.25, 
          height = plot_height / 25.4 * 2.25
        )
      }
      
      drug_cor_sum <- plyr::ddply(
        drug_correlation_i, 
        c("task", "run", "fold", "name"), 
        function(x) data.frame(
          dr_abs_cor_train_mean = mean(abs(x[["dr_cor_train"]]), na.rm = TRUE), 
          dr_abs_cor_test_mean = mean(abs(x[["dr_cor_test"]]), na.rm = TRUE)))
      
      drug_correlation[[i]] <- drug_cor_sum
    }
    drug_correlation <- Reduce(outer_join, drug_correlation)
    drug_correlation_mean <- plyr::ddply(
      drug_correlation, 
      c("task", best_by_external), 
      function(x) data.frame(
        dr_abs_cor_train_mean = mean(x$dr_abs_cor_train_mean), 
        dr_abs_cor_test_mean = mean(x$dr_abs_cor_test_mean)))
  }
  
  # Re-plot objectives
  pareto_table_obj <- plyr::join(pareto_table, external_perf_mean)
  #pareto_table_obj[["highest sum"]] <- FALSE
  indicator_df <- best_external_mean[, "task", drop = FALSE]
  indicator_df[["Optimal"]] <- TRUE
  pareto_table_obj <- plyr::join(pareto_table_obj, indicator_df, type = "left")
  pareto_table_obj[is.na(pareto_table_obj[["Optimal"]]), "Optimal"] <- FALSE
  color_limits <- c(
    min(pareto_table_obj[["test_objective"]], na.rm = TRUE), 
    max(pareto_table_obj[["test_objective"]], na.rm = TRUE))
  #pdf(paste0(plot_path, "objectives_plot_colored.pdf"), width = plot_width / 25.4 * 1.25, height = plot_height / 25.4 * 2.25)
  pareto_table_obj[["Objective"]] <- pareto_table_obj[["test_objective"]]
  save_figure_safe(
    COPS::pareto_plot(
    #pareto_plot(
      pareto_table_obj, 
      metrics = objective_names[1:4], 
      color_var = Objective, 
      shape_var = Optimal, 
      size_var = NULL, 
      plot_pareto_front = TRUE, 
      front_color = "black", 
      metric_comparators = objective_comparators[1:4], 
      point_args = list(size = 2), 
      color_scale = scale_color_gradientn(
        limits = color_limits, 
        colours = pals::kovesi.rainbow(100)),
      #color_scale = scale_color_distiller(palette = "RdBu"), 
      #color_guide = guide_colourbar(title = "Sum of \nobjectives", ncol = 1, order = 1), 
      shape_scale = scale_shape_manual(values = c(3,8))), 
    pdf, 
    paste0(plot_path, "objectives_plot_colored.pdf"), 
    width = plot_width / 25.4 * 0.85, #1.25, 
    height = plot_width / 25.4 * 0.8 #plot_height / 25.4 * 2.25
  )
  #setEPS()
  save_figure_safe(
    COPS::pareto_plot(
    #pareto_plot(
      pareto_table_obj, 
      metrics = objective_names[1:4], 
      color_var = Objective, 
      shape_var = Optimal, 
      size_var = NULL, 
      plot_pareto_front = TRUE, 
      front_color = "black", 
      metric_comparators = objective_comparators[1:4], 
      point_args = list(size = 2), 
      color_scale = scale_color_gradientn(
        limits = color_limits, 
        colours = pals::kovesi.rainbow(100)),
      #color_scale = scale_color_distiller(palette = "Greens", direction = 1), 
      color_guide = guide_colourbar(title = "Sum of \nmetrics", ncol = 1, order = 1), 
      shape_scale = scale_shape_manual(values = c(3,8))), 
    #pdf, 
    #paste0(plot_path, "objectives_plot_colored_grad.pdf"), 
    #postscript, 
    #paste0(plot_path, "objectives_plot_colored_grad.eps"), 
    png, 
    paste0(plot_path, "objectives_plot_colored_grad.png"), 
    units = "in", 
    res = 300, 
    width = plot_width / 25.4 * 0.85, #1.25, 
    height = plot_width / 25.4 * 0.8 #plot_height / 25.4 * 2.25
  )
  #color_breaks <- c(
  #  min(pareto_table_obj[["test_objective"]], na.rm = TRUE), 
  #  max(pareto_table_obj[["test_objective"]], na.rm = TRUE))
  #colours = circlize::colorRamp2(color_breaks, c("blue", "red"))
  
  #png(paste0(plot_path, "objectives_plot_colored_surv_overfit.png"), width = plot_width / 25.4 * 1.25, height = plot_height / 25.4 * 2.25, units = "in", res = 300)
  #pareto_table_obj <- plyr::join(pareto_table, external_perf_mean)
  color_limits <- c(
    min(pareto_table_obj[["survival_overfit"]], na.rm = TRUE), 
    max(pareto_table_obj[["survival_overfit"]], na.rm = TRUE))
  save_figure_safe(
    COPS::pareto_plot(
      pareto_table_obj, 
      metrics = objective_names[1:4], 
      color_var = survival_overfit, 
      shape_var = Optimal, 
      size_var = NULL, 
      plot_pareto_front = TRUE, 
      front_color = "black", 
      metric_comparators = objective_comparators[1:4], 
      point_args = list(size = 2), 
      color_scale = scale_color_gradientn(
        limits = color_limits, 
        colours = pals::kovesi.rainbow(100)),
      shape_scale = scale_shape_manual(values = c(3,8))), 
    png, 
    paste0(plot_path, "objectives_plot_colored_surv_overfit.png"), 
    width = plot_width / 25.4 * 0.85, #1.25, 
    height = plot_width / 25.4 * 0.8, #plot_height / 25.4 * 2.25
    units = "in", 
    res = 300
  )
  
  if (!is.null(drug_correlation_mean)) {
    #png(paste0(plot_path, "objectives_plot_colored_drug_abscor.png"), width = plot_width / 25.4 * 1.25, height = plot_height / 25.4 * 2.25, units = "in", res = 300)
    pareto_table_obj <- plyr::join(pareto_table_obj, drug_correlation_mean)
    color_limits <- c(
      min(pareto_table_obj[["dr_abs_cor_train_mean"]], na.rm = TRUE), 
      max(pareto_table_obj[["dr_abs_cor_train_mean"]], na.rm = TRUE))
    pareto_table_obj[["dr_abs_cor_mean"]] <- pareto_table_obj[["dr_abs_cor_train_mean"]]
    save_figure_safe(
      COPS::pareto_plot(
        pareto_table_obj, 
        metrics = objective_names[1:3], 
        color_var = dr_abs_cor_mean, 
        shape_var = Optimal, 
        size_var = NULL, 
        plot_pareto_front = TRUE, 
        front_color = "black", 
        metric_comparators = objective_comparators[1:3], 
        point_args = list(size = 2), 
        color_scale = scale_color_gradientn(
          limits = color_limits, 
          colours = pals::kovesi.rainbow(100)),
        shape_scale = scale_shape_manual(values = c(3,8))), 
      png, 
      paste0(plot_path, "objectives_plot_colored_drug_abscor.png"), 
      width = plot_width / 25.4 * 0.85, #1.25, 
      height = plot_width / 25.4 * 0.8, #plot_height / 25.4 * 2.25
      units = "in", 
      res = 300
    )
  }
  
  if (!is.null(tissue_classifier_stability_processed)) {
    color_limits <- c(
      min(pareto_table_obj[["Objective"]], na.rm = TRUE), 
      max(pareto_table_obj[["Objective"]], na.rm = TRUE))
    pareto_table_obj[["primary_cl_dot"]] <- pareto_table_obj[["tissue_class_cl_train_primary_soft_max_mean_dot"]]
    pareto_table_obj[["metastatic_cl_dot"]] <- pareto_table_obj[["tissue_class_cl_train_metastatic_soft_max_mean_dot"]]
    tissue_objectives <- c("CL tissue BACC", "primary_cl_dot", "metastatic_cl_dot", "batch RFC BACC")
    save_figure_safe(
      COPS::pareto_plot(
        pareto_table_obj, 
        metrics = tissue_objectives, 
        color_var = Objective, 
        shape_var = Optimal, 
        size_var = NULL, 
        plot_pareto_front = TRUE, 
        front_color = "black", 
        metric_comparators = list(
          .Primitive(">"), 
          .Primitive(">"), 
          .Primitive(">"), 
          .Primitive("<")), 
        point_args = list(size = 2), 
        color_scale = scale_color_gradientn(
          limits = color_limits, 
          colours = pals::kovesi.rainbow(100)),
        shape_scale = scale_shape_manual(values = c(3,8))), 
      png, 
      paste0(plot_path, "objectives_plot_colored_tissue_stability_ari.png"), 
      width = plot_width / 25.4 * 0.85, #1.25, 
      height = plot_width / 25.4 * 0.8, #plot_height / 25.4 * 2.25
      units = "in", 
      res = 300
    )
  }
  
  
  # Baseline metrics
  if (FALSE) {
    eln_aac = pd.read_csv(ctrp_path + 'elasticnet_aac1.csv', header = 0, index_col = 0)
    eln_aac_old = pd.read_csv(ctrp_path + 'elasticnet_aac_old.csv', header = 0, index_col = 0)
    eln_aac_old_fixed = pd.read_csv(ctrp_path + 'elasticnet_aac_old_fixed.csv', header = 0, index_col = 0)
  }
  
  # Diagnostics
  diagnostics_list <- list()
  for (file_list in list(ps_files, t_files)) {
    diag_files <- file_list[grep("diagnostics", file_list)]
    
    diags <- lapply(diag_files, read_result_files, path = save_path, model_name = model_name)
    diags <- Reduce(COPS::rbind_fill, diags)
    
    #if (!"iteration" %in% colnames(x)){
    #  diags <- plyr::ddply(diags, c("fold", "run", "task", "stage", "name"),
    #                       function(x) {x$iteration <- 1:nrow(x);return(x)})
    #}
    diagnostics_list <- c(diagnostics_list, list(diags))
  }
  names(diagnostics_list) <- c("ps", "test")
  
  ## Check diagnostic plots of best results
  for (i in c("ps", "test")) {
    for (j in c("train", "valid")) {
      temp <- diagnostics_list[[i]]
      if (length(temp) > 0) {
        temp_best <- plyr::join(best_external_mean[, c(best_by_external, "task"), drop = FALSE], 
                                temp, by = c(best_by_external, "task"))
        temp_best_shaped <- reshape2::melt(temp_best, 
                                           id.vars = c("iteration", "fold", "run", "task", "stage", "name"),
                                           variable.name = "loss_type", value.name = "loss")
        
        temp_best_shaped_valid <- temp_best_shaped[grepl(paste0("^", j, "_"), temp_best_shaped[["stage"]]),]
        
        stage_map <- c("Auto_Encoder", "Batch_Detection", "Batch_Correction", "Survival_Risk", "Drug_Response", "Joint")
        names(stage_map) <- paste0(j, "_losses", c("_ae", "_bd", "_bc", "_sr", "_dr", ""))
        temp_best_shaped_valid[["Training"]] <- factor(stage_map[as.character(temp_best_shaped_valid[["stage"]])], 
                                                       levels = stage_map)
        
        loss_map <- c("P_MSE", "CL_MSE", "D_FNORM", "B_SCORE", "B_DSC", "SR_LL", "DR_MSE", "REG")
        names(loss_map) <- c("reconstruction_loss_dataset1", "reconstruction_loss_dataset2",
                             "confounder_alignment_norm" , "batch_cross_entropy", "batch_dsc", "survival_log_likelihood", 
                             "drug_response_mse", "regularization")
        temp_best_shaped_valid[["Loss_Type"]] <- factor(loss_map[as.character(temp_best_shaped_valid[["loss_type"]])], 
                                                        levels = loss_map)
        
        temp_best_shaped_valid[["instance"]] <- paste(temp_best_shaped_valid[["run"]], 
                                                      temp_best_shaped_valid[["fold"]], 
                                                      temp_best_shaped_valid[["task"]], 
                                                      sep = "_")
        
        if (FALSE) {
          temp_best_shaped_valid <- temp_best_shaped_valid[temp_best_shaped_valid[["run"]] == 0 & 
                                                             temp_best_shaped_valid[["fold"]] == 0 &
                                                             temp_best_shaped_valid[["iteration"]] %in% 1:1e6,]
        }
        
        #png(paste0(plot_path, "diagnostic_plot_", i, "_", j, ".png"), 
        #    width = plot_width * 2.5, height = plot_height * 2, res = plot_res, units = plot_units)
        save_figure_safe(
          ggplot(temp_best_shaped_valid, 
                 aes(iteration, loss, color = instance)) + 
            geom_line() + theme_bw() + scale_color_brewer(palette = "Dark2") + xlab("Epoch") + 
            #facet_wrap(loss_type ~., scale = "free_y") + 
            ggh4x::facet_grid2(Loss_Type ~ Training, scales = "free", independent = "none",
                               space = "free_x") + 
            theme(#axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank()) + 
            guides(color = "none"),# + ggtitle(paste("Model losses", i,"phase", j, "set")), 
          pdf, #png, 
          paste0(plot_path, "diagnostic_plot_", i, "_", j, ".pdf"), 
          #width = plot_width * 2.5, 
          #height = plot_height * 2, 
          width = plot_width / 25.4 * 1.5, 
          height = plot_width / 25.4 * 1#, 
          #res = plot_res, 
          #units = plot_units
        )
        #png(paste0(plot_path, "diagnostic_plot_", i, "_", j, "_batch_only.png"), 
        #    width = plot_width * 2.5, height = plot_height * 2, res = plot_res, units = plot_units)
        save_figure_safe(
          ggplot(temp_best_shaped_valid[temp_best_shaped_valid[["Training"]] %in% c("Batch_Detection", "Batch_Correction"),], 
                 aes(iteration, loss, color = instance)) + 
            geom_line() + theme_bw() + scale_color_brewer(palette = "Dark2") + xlab("Epoch") + 
            #facet_wrap(loss_type ~., scale = "free_y") + 
            ggh4x::facet_grid2(Loss_Type ~ Training, scales = "free", independent = "none",
                               space = "free_x") + 
            theme(#axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank()) + 
            guides(color = "none") + ggtitle(paste("Model losses", i,"phase", j, "set")), 
          png, 
          paste0(plot_path, "diagnostic_plot_", i, "_", j, "_batch_only.png"), 
          width = plot_width * 2.5, 
          height = plot_height * 2,
          res = plot_res, 
          units = plot_units
        )
      }
    }
  }
  
  # Finished iterations vs losses
  for (i in c("ps", "test")) {
    if (!is.null(diagnostics_list[[i]])) {
      temp <- plyr::ddply(
        diagnostics_list[[i]], 
        c("stage", "task", "run", "fold"), 
        function(x) data.frame(iterations = sum(!is.na(x[["reconstruction_loss_dataset1"]])))
      )
      # Only works properly without pre-train
      temp_mean <- plyr::ddply(
        temp, 
        c("task"), 
        function(x) data.frame(iterations = mean(x[["iterations"]]))
      )
      metricsi_mean <- plyr::ddply(
        metrics_list[[i]], 
        c("stage", "task"), 
        function(x) as.data.frame(lapply(x[1:26], mean, na.rm = TRUE))
      )
      metricsi_mean_melt <- reshape2::melt(
        metricsi_mean, 
        id.vars = c("stage", "task"), 
        variable.name = c("loss_name"), 
        value.name = c("loss_value")
      )
      temp_comb <- plyr::join(metricsi_mean_melt, temp_mean, by = c("task"))
      
      ggplot(temp_comb, aes(x = iterations, y = loss_value)) + 
        geom_point() + 
        theme_bw() + 
        facet_wrap(. ~ loss_name, scales = "free_y")
      
      temp_par <- reshape2::melt(
        pp_table[,c("task", pp_relevant[pp_relevant %in% colnames(pp_table)])], 
        id.vars = c("task"), 
        variable.name = c("par_name"), 
        value.name = c("par_value")
      )
      temp_par <- plyr::join(temp_par, temp_mean, by = c("task"))
      
      ggplot(temp_par, aes(x = par_value, y = iterations)) + 
        geom_point() + 
        theme_bw() + 
        facet_wrap(. ~ par_name, scales = "free_x")
    }
  }
  
  # TODO: update from above
  if (!exists("parameters")) {
    # Parameters
    param_files <- all_files[grepl("_parameters_task[0-9]+\\.csv", all_files)]
    parameters <- lapply(paste0(save_path, param_files), 
                         read.csv, header = TRUE, row.names = NULL)
    parameters <- Reduce(COPS::rbind_fill, parameters)
    parameters[["name"]] <- gsub(".*/", "", parameters[["file_name_prefix"]])
    parameters[["name"]] <- gsub("_$", "", parameters[["name"]])
    parameters[["task"]] <- gsub("_$", "", parameters[["task_id"]])
    enc_layers <- strsplit(gsub("\\[|\\]|\\s", "", parameters[["encoder_layers"]]), split = ",")
    for (layer_i in 1:(max(sapply(enc_layers, length)))) {
      parameters[[paste0("encoder_layer", layer_i)]] <- as.numeric(sapply(enc_layers, function(x) x[layer_i]))
    }
    parameters[["bottle_neck"]] <- as.numeric(sapply(enc_layers, function(x) x[length(x)]))
    dec_layers <- strsplit(gsub("\\[|\\]|\\s", "", parameters[["decoder_layers"]]), split = ",")
    for (layer_i in 1:(max(sapply(dec_layers, length)))) {
      parameters[[paste0("decoder_layer", layer_i)]] <- as.numeric(sapply(dec_layers, function(x) x[layer_i]))
    }
    cl_layers <- strsplit(gsub("\\[|\\]|\\s", "", parameters[["classifier_layers"]]), split = ",")
    for (layer_i in 1:(max(sapply(cl_layers, length)))) {
      parameters[[paste0("classifier_layer", layer_i)]] <- as.numeric(sapply(cl_layers, function(x) x[layer_i]))
    }
    sr_layers <- strsplit(gsub("\\[|\\]|\\s", "", parameters[["survival_model_layers"]]), split = ",")
    for (layer_i in 1:(max(sapply(sr_layers, length)))) {
      parameters[[paste0("survival_model_layer", layer_i)]] <- as.numeric(sapply(sr_layers, function(x) x[layer_i]))
    }
    #ow <- lapply(strsplit(parameters[["objective_weights"]], split = "\\n "), function(x) as.numeric(gsub("\\[|array|\\(|\\)|\\]", "", x)))
    ow <- lapply(strsplit(parameters[["objective_weights"]], split = " "), function(x) as.numeric(gsub("\\[|array|\\(|\\)|\\]", "", x)))
    parameters[["reconstruction_weight"]] <- sapply(ow, function(x) x[[1]])
    parameters[["classifier_weight"]] <- sapply(ow, function(x) x[[2]])
    parameters[["survival_weight"]] <- sapply(ow, function(x) x[[3]])
    parameters[["batch_correction_weight"]] <- sapply(ow, function(x) x[[4]])
    parameters[["drug_sensitivity_weight"]] <- sapply(ow, function(x) x[[5]])
    parameters[["confounder_norm_weight"]] <- sapply(ow, function(x) try(x[[6]]))
  }
  
  best_parameters <- plyr::join(best_external_mean[, "task", drop = FALSE], parameters)
  shared_embedding_names <- paste0("z", 1:best_parameters[["bottle_neck"]][1])
  # Embeddings
  embedding_list <- list()
  for (file_list in list(ps_files, t_files)) {
    embedding_files <- file_list[grep("embeddings", file_list)]
    embedding_files <- embedding_files[grep(paste0("task", best_external_mean[1,"task"]), embedding_files)]
    
    embeddings <- lapply(embedding_files, read_result_files, path = save_path, model_name = model_name)
    embeddings <- Reduce(COPS::rbind_fill, embeddings)
    
    embedding_list <- c(embedding_list, list(embeddings))
  }
  names(embedding_list) <- c("ps", "test")
  
  embeddings_best <- plyr::join(best_external_mean[, c(best_by_external, "task"), drop = FALSE], 
                                embedding_list[["test"]], 
                                by = c(best_by_external, "task"))
  single_embedding <- embeddings_best[embeddings_best[["run"]] == 0 &
                                        embeddings_best[["fold"]] == 0, ]
  
  single_embedding[["sample_type"]] <- NA
  single_embedding[["sample_type"]] <- factor(ifelse(grepl("TCGA", single_embedding[["X"]]), "patient", "cell-line"), levels = c("patient", "cell-line"))
  shape_scale <- scale_shape_manual(values = c(3,1))
  
  patient_expression_file_name <- 'mrna.csv.gz'
  tcga_datasets <- list.dirs(patient_expression_root_dir, full.names = FALSE, recursive = FALSE)
  tcga_samples <- list()
  tcga_cancer_type_map <- list()
  for (i in tcga_datasets) {
    raw_str <- readLines(paste0(patient_expression_root_dir, i, "/", patient_expression_file_name), n = 1)
    raw_str <- gsub("\"", "", raw_str)
    raw_str <- gsub("\\.", "-", raw_str)
    tcga_samples[[i]] <- strsplit(raw_str, ",")[[1]][-1]
    tcga_cancer_type_map[[i]] <- rep(i, length(tcga_samples[[i]]))
    names(tcga_cancer_type_map[[i]]) <- tcga_samples[[i]]
  }
  tcga_cancer_type_map <- Reduce("c", tcga_cancer_type_map)
  
  cell_line_oncotree_mapping_file <- 'Model_oncotree.csv'
  cell_line_oncotree_mappings <- read.csv(paste0(cell_line_expression_root_dir, cell_line_oncotree_mapping_file), header = TRUE, row.names = 1)
  cell_line_brca_map <- cell_line_oncotree_mappings[["level_2"]] == "BRCA"
  names(cell_line_brca_map) <- rownames(cell_line_oncotree_mappings)
  single_embedding[["brca_cell_line"]] <- single_embedding[["X"]] %in% names(which(cell_line_brca_map))
  
  single_embedding[["brca_class"]] <- NA
  if (grepl("scanb", save_dir)) {
    single_embedding[["brca_class"]][single_embedding[["brca_cell_line"]]] <- "BRCA CL"
    scanb_y <- read.csv(paste0(scanb_path, "scanb_pheno.csv.gz"), row.names = 1, header = TRUE)
    brca_match <- match(single_embedding[["X"]], scanb_y[["GEX.assay"]])
    single_embedding[["brca_class"]][!is.na(brca_match)] <- scanb_y[["NCN.PAM50"]][brca_match[!is.na(brca_match)]]
    single_embedding[["brca_class"]][is.na(single_embedding[["brca_class"]])] <- "NA"
  } else {
    # Assume TCGA
    single_embedding[["brca_class"]][!grepl("TCGA", single_embedding[["X"]]) & single_embedding[["brca_cell_line"]]] <- "BRCA CL"
    brca_subtypes <- TCGAbiolinks::TCGAquery_subtype("BRCA")
    brca_match <- match(substr(single_embedding[["X"]], 1, 12), brca_subtypes[["patient"]])
    single_embedding[["brca_class"]][!is.na(brca_match)] <- brca_subtypes[["BRCA_Subtype_PAM50"]][brca_match[!is.na(brca_match)]]
    single_embedding[["brca_class"]][is.na(single_embedding[["brca_class"]])] <- "NA"
  }
  
  table(single_embedding[["brca_class"]])
  
  brca_class_ind <- which(single_embedding[["brca_class"]] != "NA")
  if (length(brca_class_ind) > 0) {
    #png(paste0(plot_path, "embedding_umap_brca_classes.png"), width = plot_width, height = plot_height, res = plot_res, units = plot_units)
    set.seed(0)
    save_figure_safe(
      COPS::umap_viz(single_embedding[brca_class_ind, shared_embedding_names], 
                     category = single_embedding[brca_class_ind, "brca_class"], 
                     category_label = "BRCA class", pre_manifold_pca = FALSE, 
                     umap_args = list(init = "normlaplacian", min_dist = 0.5)), 
      png, 
      paste0(plot_path, "embedding_umap_brca_classes.png"), 
      width = plot_width, 
      height = plot_height, 
      res = plot_res, 
      units = plot_units
    )
    
    set.seed(0)
    save_figure_safe(
      COPS::tsne_viz(single_embedding[brca_class_ind, shared_embedding_names], 
                     category = single_embedding[brca_class_ind, "brca_class"], 
                     category_label = "BRCA class", pre_manifold_pca = FALSE), 
      png, 
      paste0(plot_path, "embedding_tsne_brca_classes.png"), 
      width = plot_width, 
      height = plot_height, 
      res = plot_res, 
      units = plot_units
    )
    #png(paste0(plot_path, "embedding_pca_brca_classes.png"), width = plot_width, height = plot_height, res = plot_res, units = plot_units)
    set.seed(0)
    save_figure_safe(
      COPS::pca_viz(single_embedding[brca_class_ind, shared_embedding_names], 
                    category = single_embedding[brca_class_ind, "brca_class"], 
                    category_label = "BRCA class"), 
      png, 
      paste0(plot_path, "embedding_pca_brca_classes.png"), 
      width = plot_width, 
      height = plot_height, 
      res = plot_res, 
      units = plot_units
    )
  }
  
  if (FALSE) {
    if (length(brca_class_ind) > 0) {
      if (!dir.exists(paste0(plot_path, "pairwise_embeddings"))) dir.create(paste0(plot_path, "pairwise_embeddings"))
      color_scale = scale_color_brewer(palette = "Dark2")
      z_ind <- grep("^z[0-9]+$", colnames(single_embedding))
      for (zi in rev(z_ind[-1])) {
        for (zj in z_ind[z_ind < zi]) {
          zi_name <- colnames(single_embedding)[zi]
          zj_name <- colnames(single_embedding)[zj]
          #png(paste0(plot_path, "pairwise_embeddings/", zj_name, "_", zi_name, "_brca_classes.png"), width = plot_width, height = plot_height, res = plot_res, units = plot_units)
          save_figure_safe(
            ggplot(single_embedding[brca_class_ind, c(zj_name, zi_name, "brca_class", "sample_type")], 
                   aes(!!ggplot2::ensym(zj_name), !!ggplot2::ensym(zi_name), color = brca_class, shape = sample_type)) + 
              geom_point(size = 3) + theme_bw() + 
              color_scale + shape_scale + labs(color = "BRCA class", shape = "Sample type"),
            png, 
            paste0(plot_path, "pairwise_embeddings/", zj_name, "_", zi_name, "_brca_classes.png"), 
            width = plot_width, 
            height = plot_height, 
            res = plot_res, 
            units = plot_units
          )
        }
      }
    }
  }
  
  ## Get cancer types from sample-dataset maps from above
  single_embedding[["cancer_type"]] <- NA
  #single_embedding[["cancer_type"]][!grepl("TCGA", single_embedding[["X"]])] <- "CL"
  cancer_map_ind <- single_embedding[["X"]] %in% names(tcga_cancer_type_map)
  single_embedding[["cancer_type"]][cancer_map_ind] <- tcga_cancer_type_map[single_embedding[["X"]][cancer_map_ind]]
  if (all(is.na(single_embedding[["cancer_type"]]))) {
    # Assume we are using SCANB + CCLE
    single_embedding[["cancer_type"]] <- "SCANB" # CCLE is identified from ID in next step
  }
  single_embedding[["cancer_type"]][grepl("ACH", single_embedding[["X"]])] <- "CL"
  single_embedding[["cancer_type"]][is.na(single_embedding[["cancer_type"]])] <- "NA"
  table(single_embedding[["cancer_type"]])
  
  n_labels <- length(unique(single_embedding[["cancer_type"]]))
  #png(paste0(plot_path, "embedding_umap_cancer_types.png"), width = plot_width, height = plot_height, res = plot_res, units = plot_units)
  set.seed(0)
  save_figure_safe(
    COPS::umap_viz(single_embedding[,shared_embedding_names], 
                   category = single_embedding[["cancer_type"]], 
                   category_label = "Cancer Type", pre_manifold_pca = FALSE, 
                   color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=n_labels)), 
                   umap_args = list(init = "normlaplacian", min_dist = 0.5)),
    png, 
    paste0(plot_path, "embedding_umap_cancer_types.png"), 
    width = plot_width, 
    height = plot_height, 
    res = plot_res, 
    units = plot_units
  )
  
  set.seed(0)
  save_figure_safe(
    COPS::tsne_viz(single_embedding[,shared_embedding_names], 
                   category = single_embedding[["cancer_type"]], 
                   category_label = "Cancer Type", pre_manifold_pca = FALSE, 
                   color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=n_labels))),
    png, 
    paste0(plot_path, "embedding_tsne_cancer_types.png"), 
    width = plot_width, 
    height = plot_height, 
    res = plot_res, 
    units = plot_units
  )
  
  #png(paste0(plot_path, "embedding_pca_cancer_types.png"), width = plot_width, height = plot_height, res = plot_res, units = plot_units)
  set.seed(0)
  save_figure_safe(
    COPS::pca_viz(single_embedding[,shared_embedding_names], 
                  category = single_embedding[["cancer_type"]], 
                  category_label = "Cancer Type", 
                  color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=n_labels))),
    png, 
    paste0(plot_path, "embedding_pca_cancer_types.png"), 
    width = plot_width,
    height = plot_height, 
    res = plot_res, 
    units = plot_units
  )
  
  #png(paste0(plot_path, "embedding_umap_sample_types.png"), width = plot_width, height = plot_height, res = plot_res, units = plot_units)
  set.seed(0)
  save_figure_safe(
    COPS::umap_viz(single_embedding[,shared_embedding_names], 
                   category = ifelse(single_embedding[["cancer_type"]] == "CL", "CL", "patient"), 
                   category_label = "Sample type", pre_manifold_pca = FALSE, 
                   umap_args = list(init = "normlaplacian", min_dist = 0.5)),
    png, 
    paste0(plot_path, "embedding_umap_sample_types.png"), 
    width = plot_width, 
    height = plot_height, 
    res = plot_res, 
    units = plot_units
  )
  
  set.seed(0)
  save_figure_safe(
    COPS::tsne_viz(single_embedding[,shared_embedding_names], 
                   category = ifelse(single_embedding[["cancer_type"]] == "CL", "CL", "patient"), 
                   category_label = "Sample type", pre_manifold_pca = FALSE),
    png, 
    paste0(plot_path, "embedding_tsne_sample_types.png"), 
    width = plot_width, 
    height = plot_height, 
    res = plot_res, 
    units = plot_units
  )
  
  #png(paste0(plot_path, "embedding_pca_sample_types.png"), width = plot_width, height = plot_height, res = plot_res, units = plot_units)
  set.seed(0)
  save_figure_safe(
    COPS::pca_viz(single_embedding[,shared_embedding_names], 
                  category = ifelse(single_embedding[["cancer_type"]] == "CL", "CL", "patient"), 
                  category_label = "Sample type"),
    png,
    paste0(plot_path, "embedding_pca_sample_types.png"), 
    width = plot_width, 
    height = plot_height, 
    res = plot_res, 
    units = plot_units
  )
  
  if (FALSE) {
    if (!dir.exists(paste0(plot_path, "pairwise_embeddings"))) dir.create(paste0(plot_path, "pairwise_embeddings"))
    color_scale <- scale_color_manual(values = pals::kovesi.rainbow(n=n_labels))
    z_ind <- grep("^z[0-9]+$", colnames(single_embedding))
    for (zi in rev(z_ind[-1])) {
      for (zj in z_ind[z_ind < zi]) {
        zi_name <- colnames(single_embedding)[zi]
        zj_name <- colnames(single_embedding)[zj]
        #png(paste0(plot_path, "pairwise_embeddings/", zj_name, "_", zi_name, "_cancer_types.png"), width = plot_width, height = plot_height, res = plot_res, units = plot_units)
        save_figure_safe(
          ggplot(single_embedding[, c(zj_name, zi_name, "cancer_type", "sample_type")], 
                 aes(!!ggplot2::ensym(zj_name), !!ggplot2::ensym(zi_name), color = cancer_type, shape = sample_type)) + 
            geom_point(size = 3) + theme_bw() + 
            color_scale + shape_scale + labs(color = "Cancer Type", shape = "Sample type"),
          png, 
          paste0(plot_path, "pairwise_embeddings/", zj_name, "_", zi_name, "_cancer_types.png"), 
          width = plot_width, 
          height = plot_height, 
          res = plot_res, 
          units = plot_units
        )
      }
    }
  }
  
  # Similarity matrix
  if (FALSE) {
    X <- as.matrix(single_embedding[, shared_embedding_names])
    K <- X %*% t(X)
    hist(K)
    K[K < -1000] <- -1000
    K[K > 5000] <- 5000
    hist(K)
    #COPS::plot_similarity_matrix(K, limits = c(-20,20))
    
    table(single_embedding[,"cancer_type"])
    n_ct <- length(unique(single_embedding[,"cancer_type"]))
    col <- pals::kovesi.rainbow(n_ct)
    names(col) <- unique(single_embedding[,"cancer_type"])
    col <- list(cancer_type = col) 
    
    top_annot <- ComplexHeatmap::HeatmapAnnotation(
      df = single_embedding[,"cancer_type", drop = FALSE], 
      col = col)
    left_annot <- ComplexHeatmap::rowAnnotation(
      df = single_embedding[,"cancer_type", drop = FALSE], 
      col = col)
    hm <- ComplexHeatmap::Heatmap(
      K, 
      top_annotation = top_annot, 
      left_annotation = left_annot, 
      show_column_dend = FALSE, 
      show_column_names = FALSE, 
      show_row_dend = FALSE, 
      show_row_names = FALSE)
    
    #png(paste0(plot_path, "linear_kernel_cancer_types.png"), width = plot_width, height = plot_height, res = plot_res, units = plot_units)
    save_figure_safe(
      hm, 
      png, 
      paste0(plot_path, "linear_kernel_cancer_types.png"), 
      width = plot_width, 
      height = plot_height, 
      res = plot_res, 
      units = plot_units
    )
  }
  
  # Prediction correlations, from parameter search data
  predictions_list <- list()
  for (file_list in list(ps_files, t_files)) {
    predictions_files <- file_list[grep("predictions", file_list)]
    predictions_files <- predictions_files[grep(paste0("task", best_external_mean[1,"task"]), predictions_files)]
    
    predictions <- lapply(predictions_files, read_result_files, path = save_path, model_name = model_name)
    predictions <- Reduce(COPS::rbind_fill, predictions)
    
    predictions_list <- c(predictions_list, list(predictions))
  }
  names(predictions_list) <- c("ps", "test")
  
  predictions_best <- plyr::join(
    best_external_mean[, c(best_by_external, "task"), drop = FALSE], 
    predictions_list[["test"]], 
    by = c(best_by_external, "task"))
  single_predictions <- predictions_best[
    predictions_best[["run"]] == 0 &
    predictions_best[["fold"]] == 0, ]
  drug_names <- readLines(paste0(ctrp_path, "drug_names.txt"))
  
  dr_pred_mat_cl <- single_predictions[
    with(single_predictions, dataset %in% c("cl_train")), #"cl_test")), 
    grep("^dr_pred_[0-9]+$", colnames(single_predictions))]
  
  dr_pred_cor_cl <- cor(dr_pred_mat_cl, method = "spearman")
  cl_pred_hm <- ComplexHeatmap::Heatmap(
    dr_pred_cor_cl[drug_order, drug_order], 
    name = "cl", 
    show_row_names = FALSE, 
    show_column_names = FALSE, 
    cluster_rows = FALSE, 
    cluster_columns = FALSE)
  
  #pdf(paste0(plot_path, "drug_prediction_correlations_example2_cl.pdf"), width = plot_width / 25.4 * 1.25, height = plot_height / 25.4 * 2.25)
  save_figure_safe(
    cl_pred_hm, 
    pdf, 
    paste0(plot_path, "drug_prediction_correlations_example2_cl.pdf"), 
    width = plot_width / 25.4 * 1.25, 
    height = plot_height / 25.4 * 2.25
  )
  
  dr_pred_mat_patient <- single_predictions[
    with(single_predictions, dataset %in% c("patient_train", "patient_test")), 
    grep("^dr_pred_[0-9]+$", colnames(single_predictions))]
  
  dr_pred_cor_patient <- cor(dr_pred_mat_patient, method = "spearman")
  patient_pred_hm <- ComplexHeatmap::Heatmap(
    dr_pred_cor_patient[drug_order, drug_order], 
    name = "patient", 
    show_row_names = FALSE, 
    show_column_names = FALSE, 
    cluster_rows = FALSE, 
    cluster_columns = FALSE)
  
  #pdf(paste0(plot_path, "drug_prediction_correlations_example2_patient.pdf"), width = plot_width / 25.4 * 1.25, height = plot_height / 25.4 * 2.25)
  save_figure_safe(
    patient_pred_hm, 
    pdf, 
    paste0(plot_path, "drug_prediction_correlations_example2_patient.pdf"), 
    width = plot_width / 25.4 * 1.25, 
    height = plot_height / 25.4 * 2.25
  )
  
  if (FALSE) {
    # Best drugs
    predictions_list <- list()
    for (file_list in list(ps_files, t_files)) {
      predictions_files <- file_list[grep("predictions", file_list)]
      predictions_files <- predictions_files[grep(paste0("task", best_external_mean[1,"task"]), predictions_files)]
      
      predictions <- lapply(predictions_files, read_result_files, path = save_path, model_name = model_name)
      predictions <- Reduce(COPS::rbind_fill, predictions)
      
      predictions_list <- c(predictions_list, list(predictions))
    }
    names(predictions_list) <- c("ps", "test")
    
    predictions_best <- plyr::join(
      best_external_mean[, c(best_by_external, "task"), drop = FALSE], 
      predictions_list[["test"]], 
      by = c(best_by_external, "task"))
    single_predictions <- predictions_best[predictions_best[["run"]] == 0 &
                                             predictions_best[["fold"]] == 0, ]
    drug_names <- readLines(paste0(ctrp_path, "drug_names.txt"))
    
    drug_prediction_ind <- grepl("dr_pred", colnames(single_predictions))
    drug_predictions <- single_predictions[
      , colnames(single_predictions) %in% c("fold","run", "task", "X") | 
        drug_prediction_ind]
    
    drug_predictions <- reshape2::melt(
      drug_predictions, 
      id.vars = c("fold","run", "task", "X"), 
      variable.name = "drug", 
      value.name = "response")
    write.csv(drug_predictions, gzfile(paste0(plot_path, "best_drug_predictions.csv.gz")))
    
    top10_drugs <- plyr::ddply(
      drug_predictions, 
      c("fold","run", "task", "X"), 
      function(x) x[nrow(x) - rank(x[["response"]]) < 10,])
    
    top10_drugs[["drug_name"]] <- drug_names[as.numeric(gsub("dr_pred_", "", top10_drugs[["drug"]]))+1]
    table(top10_drugs[["drug_name"]])
    
    top1_drugs <- plyr::ddply(
      drug_predictions, 
      c("fold","run", "task", "X"), 
      function(x) x[nrow(x) - rank(x[["response"]]) < 1,])
    
    top1_drugs[["drug_name"]] <- drug_names[as.numeric(gsub("dr_pred_", "", top1_drugs[["drug"]]))+1]
    
    n_labels <- length(table(top1_drugs[["drug_name"]]))
    single_embedding_drugged <- plyr::join(single_embedding, top1_drugs)
    brca_ind <- single_embedding_drugged[["cancer_type"]] == "BRCA"
    #png(paste0(plot_path, "brca_embedding_umap_best_drug.png"), width = plot_width, height = plot_height, res = plot_res, units = plot_units)
    set.seed(0)
    save_figure_safe(
      COPS::umap_viz(
        single_embedding_drugged[brca_ind, shared_embedding_names], 
        category = single_embedding_drugged[brca_ind, "drug_name"], 
        category_label = "Best Drug", pre_manifold_pca = FALSE, 
        color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=n_labels)), 
        umap_args = list(init = "normlaplacian", min_dist = 0.5)), 
      png, 
      paste0(plot_path, "brca_embedding_umap_best_drug.png"), 
      width = plot_width, 
      height = plot_height, 
      res = plot_res, 
      units = plot_units
    )
    
    #png(paste0(plot_path, "brca_embedding_pca_best_drug.png"), width = plot_width, height = plot_height, res = plot_res, units = plot_units)
    set.seed(0)
    save_figure_safe(
      COPS::pca_viz(
        single_embedding_drugged[brca_ind, shared_embedding_names], 
        category = single_embedding_drugged[brca_ind, "drug_name"], 
        category_label = "Best Drug", pre_manifold_pca = FALSE, 
        color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=n_labels))), 
      png,
      paste0(plot_path, "brca_embedding_pca_best_drug.png"), 
      width = plot_width, 
      height = plot_height, 
      res = plot_res, 
      units = plot_units
    )
  }
}

if (FALSE) {
  brca_exp <- read.csv(paste0(base_dir, "brca_data.csv"), row.names = 1, header = TRUE)
  ccle_exp <- read.csv(paste0(base_dir, "ccle_data.csv"), row.names = 1, header = TRUE)
  
  set.seed(0)
  #png(paste0(plot_path, "gex_data_umap.png"), width = plot_width * 1.15, height = plot_height, res = plot_res, units = plot_units)
  save_figure_safe(
    COPS::umap_viz(rbind(brca_exp, ccle_exp), 
                 category = rep(c("TCGA BRCA", "CCLE"), c(nrow(brca_exp), nrow(ccle_exp))), 
                 category_label = "dataset", 
                 umap_args = list(init = "normlaplacian", min_dist = 0.5)) + 
    theme_bw(base_size = 20) + geom_point(shape = "+", size = rel(0.75)), 
    png, 
    paste0(plot_path, "gex_data_umap.png"), 
    width = plot_width * 1.15, 
    height = plot_height, 
    res = plot_res, 
    units = plot_units
  )
  
  set.seed(0)
  #png(paste0(plot_path, "gex_data_pca.png"), width = plot_width * 1.15, height = plot_height, res = plot_res, units = plot_units)
  save_figure_safe(
    COPS::pca_viz(rbind(brca_exp, ccle_exp), 
                category = rep(c("TCGA BRCA", "CCLE"), c(nrow(brca_exp), nrow(ccle_exp))), 
                category_label = "dataset") + 
    theme_bw(base_size = 20) + geom_point(shape = "+", size = rel(0.75)), 
    png, 
    paste0(plot_path, "gex_data_pca.png"), 
    width = plot_width * 1.15, 
    height = plot_height, 
    res = plot_res, 
    units = plot_units
  )
  
  zero_var_genes <- lapply(list(brca_exp, ccle_exp), function(x) which(apply(x, 2, var) == 0))
  zero_var_genes <- Reduce(union, zero_var_genes)
  if (length(zero_var_genes) > 0) {
    brca_exp_filtered <- brca_exp[,-zero_var_genes]
    ccle_exp_filtered <- ccle_exp[,-zero_var_genes]
  } else {
    brca_exp_filtered <- brca_exp
    ccle_exp_filtered <- ccle_exp
  }
  set.seed(0)
  #png(paste0(plot_path, "gex_data_scaled_umap.png"), width = plot_width * 1.15, height = plot_height, res = plot_res, units = plot_units)
  save_figure_safe(
    COPS::umap_viz(rbind(scale(brca_exp_filtered), scale(ccle_exp_filtered)), 
                 category = rep(c("TCGA BRCA", "CCLE"), c(nrow(brca_exp), nrow(ccle_exp))), 
                 category_label = "dataset", 
                 umap_args = list(init = "normlaplacian", min_dist = 0.5)) + 
    theme_bw(base_size = 20) + geom_point(shape = "+", size = rel(0.75)), 
    png,
    paste0(plot_path, "gex_data_scaled_umap.png"), 
    width = plot_width * 1.15, 
    height = plot_height, 
    res = plot_res, 
    units = plot_units
  )
  
  set.seed(0)
  #png(paste0(plot_path, "gex_data_scaled_pca.png"), width = plot_width * 1.15, height = plot_height, res = plot_res, units = plot_units)
  save_figure_safe(
    COPS::pca_viz(rbind(scale(brca_exp_filtered), scale(ccle_exp_filtered)), 
                category = rep(c("TCGA BRCA", "CCLE"), c(nrow(brca_exp), nrow(ccle_exp))), 
                category_label = "dataset") + 
    theme_bw(base_size = 20) + geom_point(shape = "+", size = rel(0.75)), 
    png, 
    paste0(plot_path, "gex_data_scaled_pca.png"), 
    width = plot_width * 1.15, 
    height = plot_height, 
    res = plot_res, 
    units = plot_units
  )
  
  n_labels <- length(unique(single_embedding[["cancer_type"]]))
  #png(paste0(plot_path, "embedding_umap_cancer_types.png"), width = plot_width * 1.15, height = plot_height, res = plot_res, units = plot_units)
  set.seed(0)
  save_figure_safe(
    COPS::umap_viz(single_embedding[,shared_embedding_names], 
                   category = single_embedding[["cancer_type"]], 
                   category_label = "Cancer Type", pre_manifold_pca = FALSE, 
                   color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=n_labels)), 
                   umap_args = list(init = "normlaplacian", min_dist = 0.5)) + 
      theme_bw(base_size = 20) + geom_point(shape = "+", size = rel(0.75)), 
    png, 
    paste0(plot_path, "embedding_umap_cancer_types.png"), 
    width = plot_width * 1.15, 
    height = plot_height, 
    res = plot_res, 
    units = plot_units
  )
  
  #png(paste0(plot_path, "embedding_pca_cancer_types.png"), width = plot_width * 1.15, height = plot_height, res = plot_res, units = plot_units)
  set.seed(0)
  save_figure_safe(
    COPS::pca_viz(single_embedding[,shared_embedding_names], 
                  category = single_embedding[["cancer_type"]], 
                  category_label = "Cancer Type", 
                  color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=n_labels))) + 
      theme_bw(base_size = 20) + geom_point(shape = "+", size = rel(0.75)), 
    png, 
    paste0(plot_path, "embedding_pca_cancer_types.png"), 
    width = plot_width * 1.15, 
    height = plot_height, 
    res = plot_res, 
    units = plot_units
  )
}







