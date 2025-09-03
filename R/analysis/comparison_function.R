generate_comparisons <- function(
    base_dir, 
    save_path, 
    plot_path, 
    #model_name, 
    source_surv_dataset, 
    target_surv_dataset, 
    gene_preselection, 
    ext_survival_dir = NULL
) {
  # Best results
  best_task_ps <- read.csv(paste0(plot_path, "best_task_ps_performance.csv"))
  best_task_test <- read.csv(paste0(plot_path, "best_task_test_performance.csv"))
  
  external_path <- paste0(save_path, "external_evaluation/", ext_survival_dir, "/")
  
  ext_surv_fn <- paste0(external_path, "external_survival_validation_res.csv")
  if (file.exists(ext_surv_fn)) {
    best_task_surv_ext <- read.csv(ext_surv_fn, header = TRUE, row.names = 1)
    best_ext_survival_perf <- data.frame(
      Model = rep("MODAE", 2), 
      Set = c(
        toupper(source_surv_dataset),
        toupper(target_surv_dataset)
      ), 
      Scope = "external validation", 
      concordance = c(
        best_task_surv_ext["0", "survival_concordance"], 
        best_task_surv_ext["1", "survival_concordance"]), 
      concordance_sd = NA)
  } else {
    best_ext_survival_perf <- data.frame()
  }
  
  best_survival_perf <- data.frame(
    Model = rep("MODAE", 2), 
    Set = paste(toupper(source_surv_dataset), c("train", "test")),
    Scope = paste(toupper(source_surv_dataset), "cross-validation"), 
    concordance = c(
      mean(best_task_test[["patient_train_surv_c"]]), 
      mean(best_task_test[["patient_test_surv_c"]])), 
    concordance_sd = c(
      sd(best_task_test[["patient_train_surv_c"]]), 
      sd(best_task_test[["patient_test_surv_c"]])))
  
  best_survival_perf <- rbind(best_survival_perf, best_ext_survival_perf)
  
  best_dr_sens_perf <- data.frame(
    Model = rep("MODAE", 2), 
    Set = paste(toupper(source_surv_dataset), c("train", "test")), #c("TCGA train", "TCGA test"), 
    Scope = paste(toupper(source_surv_dataset), "cross-validation"), #rep("TCGA cross-validation", 2),
    dr_R2_Q90 = c(
      mean(best_task_test[["dr_cl_train_R2_Q90"]]), 
      mean(best_task_test[["dr_cl_test_R2_Q90"]])), 
    dr_R2_Q90_sd = c(
      sd(best_task_test[["dr_cl_train_R2_Q90"]]), 
      sd(best_task_test[["dr_cl_test_R2_Q90"]])))
  
  # Baseline survival comparison
  
  base_surv_dir <- paste0(output_dir, "baseline_results/survival/", ext_survival_dir, "/", source_surv_dataset, "/")
  if (gene_preselection) {
    base_surv_dir <- paste0(base_surv_dir, "gene_preselection/")
  }
  
  baseline_survival_files <- dir(
    base_surv_dir, 
    pattern = "*_validation.csv", 
    full.names = TRUE)
  age_model_ind <- grep("age_model", baseline_survival_files)
  if (length(age_model_ind) > 0) {
    # Age model is a sanity check with different formatting ... 
    # Easiest solution is to remove it. 
    baseline_survival_files <- baseline_survival_files[-age_model_ind]
  }
  baseline_survival_list <- list()
  for (i in baseline_survival_files) {
    fn <- gsub(".*/", "", i)
    model_name <- gsub("_.*", "", fn)
    dataset <- gsub(".*_tcga_", "", fn)
    dataset <- gsub("\\.csv", "", dataset)
    temp <- read.csv(i, header = TRUE)
    temp[["model_name"]] <- model_name
    #temp[["dataset"]] <- dataset
    temp[["X"]] <- NULL # invalid row index from pandas
    baseline_survival_list[[paste(model_name, dataset)]] <- temp
  }
  ext_ind <- grep("external_validation", names(baseline_survival_list))
  cv_ind <- grep("cross_validation", names(baseline_survival_list))
  baseline_survival_ext <- Reduce(COPS::rbind_fill, baseline_survival_list[ext_ind])
  baseline_survival_cv <- Reduce(COPS::rbind_fill, baseline_survival_list[cv_ind])
  
  if (!is.null(baseline_survival_ext)) {
    baseline_survival <- plyr::join(
      baseline_survival_cv, 
      baseline_survival_ext, 
      type = "full")
  } else {
    baseline_survival <- baseline_survival_cv
  }
  
  baseline_survival[["Model"]] <- ifelse(
    is.na(baseline_survival[["npc"]]), 
    "", 
    paste0(baseline_survival[["npc"]], " ")
  )
  model_names <- c(
    dr = "PCs + EN", 
    fs = "UFS + EN"
  )
  baseline_survival[["Model"]] <- paste0(
    baseline_survival[["Model"]], 
    model_names[baseline_survival[["model_name"]]])
  
  surv_res_cols <- paste0(
    c(
      paste0(source_surv_dataset, c("_train", "_test")), 
      source_surv_dataset, 
      target_surv_dataset
    ), 
    "_c"
  )
  
  surv_res_cols <- surv_res_cols[surv_res_cols %in% colnames(baseline_survival)]
  baseline_survival_df <- reshape2::melt(
    baseline_survival, 
    id.vars = c("Model"), 
    measure.vars = surv_res_cols, 
    variable.name = "dataset", 
    value.name = "concordance"
  )
  dataset_names <- c(
    paste(toupper(source_surv_dataset), c("train", "test")), 
    toupper(source_surv_dataset), 
    toupper(target_surv_dataset)
  )
  dataset_names <- dataset_names[1:length(dataset_names)]
  names(dataset_names) <- surv_res_cols
  dataset_scope <- c(
    rep(paste(toupper(source_surv_dataset), "cross-validation"), 2), 
    rep("external validation", 2)
  )
  names(dataset_scope) <- surv_res_cols
  
  baseline_survival_df[["Set"]] <- dataset_names[baseline_survival_df[["dataset"]]]
  baseline_survival_df[["Scope"]] <- dataset_scope[baseline_survival_df[["dataset"]]]
  
  baseline_survival_df <- plyr::ddply(
    baseline_survival_df, c("Model", "dataset", "Set", "Scope"), 
    function(x) data.frame(
      concordance = mean(x[["concordance"]]), 
      concordance_sd = sd(x[["concordance"]])))
  baseline_survival_df[["concordance_sd"]][baseline_survival_df[["concordance_sd"]] == 0] <- NA
  
  surv_perf_df <- rbind(best_survival_perf, baseline_survival_df[,colnames(best_survival_perf)])
  surv_perf_df[["Set"]] <- factor(
    as.character(surv_perf_df[["Set"]]), 
    levels = best_survival_perf[["Set"]])
  surv_perf_df[["Scope"]] <- factor(
    as.character(surv_perf_df[["Scope"]]), 
    levels = best_survival_perf[c(1,3), "Scope"])
    #levels = c("TCGA cross-validation", "external validation"))
  surv_perf_df[["Model"]] <- factor(
    as.character(surv_perf_df[["Model"]]), 
    levels = c("10 PCs + EN", "100 PCs + EN", "UFS + EN", "MODAE"))
  
  surv_comp_plot <- ggplot(
    surv_perf_df, 
    aes(x = Set, y = concordance, fill = Model)) + 
    geom_bar(position = "dodge", stat = "identity") + 
    geom_errorbar(
      aes(ymin = concordance - concordance_sd, 
          ymax = concordance + concordance_sd), 
      position = position_dodge(0.9), width = 0.2) +
    theme_bw() + scale_fill_brewer(palette = "Dark2") + 
    theme(axis.title.x = element_blank()) + ylim(0,1) +
    facet_wrap(Scope ~ ., scales = "free_x") 
  
  surv_comp_cv_plot <- ggplot(
    surv_perf_df[grepl("cross-validation", surv_perf_df[["Scope"]]), ], 
    aes(x = Set, y = concordance, fill = Model)) + 
    geom_bar(position = "dodge", stat = "identity") + 
    geom_errorbar(
      aes(ymin = concordance - concordance_sd, 
          ymax = concordance + concordance_sd), 
      position = position_dodge(0.9), width = 0.2) +
    theme_bw() + scale_fill_brewer(palette = "Dark2") + 
    theme(axis.title.x = element_blank()) + ylim(0,1)
  
  surv_comp_ext_plot <- ggplot(
    surv_perf_df[surv_perf_df[["Scope"]] == "external validation", ], 
    aes(x = Set, y = concordance, fill = Model)) + 
    geom_bar(position = "dodge", stat = "identity") + 
    geom_errorbar(
      aes(ymin = concordance - concordance_sd, 
          ymax = concordance + concordance_sd), 
      position = position_dodge(0.9), width = 0.2) +
    theme_bw() + scale_fill_brewer(palette = "Dark2") + 
    theme(axis.title.x = element_blank()) + ylim(0,1)
  
  # Baseline drug sensitivity prediction
  if (gene_preselection) {
    pca10_elnet_dr_r2 <- read.csv(paste0(xia_path, "pca10_elasticnet_aac_identical.csv"))
    pca100_elnet_dr_r2 <- read.csv(paste0(xia_path, "pca100_elasticnet_aac_identical.csv"))
    fs_elnet_dr_r2 <- read.csv(paste0(xia_path, "fs_elasticnet_aac_identical.csv"))
  } else {
    pca10_elnet_dr_r2 <- read.csv(paste0(xia_path, "pca10_elasticnet_aac_old.csv"))
    pca100_elnet_dr_r2 <- read.csv(paste0(xia_path, "pca100_elasticnet_aac_old.csv"))
    fs_elnet_dr_r2 <- read.csv(paste0(xia_path, "fs_elasticnet_aac_old.csv"))
  }
  
  fold_fixer <- function(x) plyr::ddply(x, c("drug"), function(y) {y[["fold"]] <- 1:5; return(y)})
  quantiler <- function(
    x, 
    cols = c("train_r2", "test_r2"), 
    col_labels = c("CTRP train", "CTRP test"), 
    q_probs = 0:10/10
  ) {
    qf <- function(y) {
      qs <- lapply(cols, function(col) quantile(y[[col]], probs = q_probs, na.rm = TRUE))
      qs <- Reduce(cbind, qs)
      qs[!is.finite(qs)] <- NA
      qs <- as.data.frame(qs)
      colnames(qs) <- col_labels
      qs[["quantile"]] <- paste0("Q", q_probs * 100)
      qs <- reshape2::melt(qs, id.vars = c("quantile"), variable.name = "Set", value.name = "R^2")
      return(qs)
    }
    out <- plyr::ddply(x, c("fold"), qf)
    out <- plyr::ddply(out, c("quantile", "Set"), 
                       function(y) data.frame(
                         R2 = mean(y[["R^2"]], na.rm = TRUE), 
                         R2_sd = sd(y[["R^2"]], na.rm = TRUE)
                       ))
    return(out)
  }
  pca10_elnet_dr_r2_qs <- quantiler(fold_fixer(pca10_elnet_dr_r2))
  pca100_elnet_dr_r2_qs <- quantiler(fold_fixer(pca100_elnet_dr_r2))
  fs_elnet_dr_r2_qs <- quantiler(fold_fixer(fs_elnet_dr_r2))
  pca10_elnet_dr_r2_qs[["Model"]] <- "10 PCs + EN"
  pca100_elnet_dr_r2_qs[["Model"]] <- "100 PCs + EN"
  fs_elnet_dr_r2_qs[["Model"]] <- "UFS + EN"
  
  baseline_dr_r2_qs <- Reduce(
    rbind, list(pca10_elnet_dr_r2_qs, pca100_elnet_dr_r2_qs, fs_elnet_dr_r2_qs))
  
  best_dr_perf <- data.frame(
    Model = rep("MODAE", 2), 
    Set = c("CTRP train", "CTRP test"), 
    R2 = c(
      mean(best_task_test[["dr_cl_train_R2_Q90"]]), 
      mean(best_task_test[["dr_cl_test_R2_Q90"]])), 
    R2_sd = c(
      sd(best_task_test[["dr_cl_train_R2_Q90"]]), 
      sd(best_task_test[["dr_cl_test_R2_Q90"]])))
  
  drug_perf_df <- rbind(
    best_dr_perf, 
    baseline_dr_r2_qs[baseline_dr_r2_qs[["quantile"]] == "Q90", c("Model", "Set", "R2", "R2_sd")])
  drug_perf_df[["Set"]] <- factor(
    as.character(drug_perf_df[["Set"]]), 
    levels = c("CTRP train", "CTRP test"))
  drug_perf_df[["Model"]] <- factor(
    as.character(drug_perf_df[["Model"]]), 
    levels = c("10 PCs + EN", "100 PCs + EN", "UFS + EN", "MODAE"))
  
  drug_comp_plot <- ggplot(drug_perf_df, aes(x = Set, y = R2, fill = Model)) + 
    geom_bar(position = "dodge", stat = "identity") + 
    geom_errorbar(
      aes(ymin = R2 - R2_sd, 
          ymax = R2 + R2_sd), 
      position = position_dodge(0.9), width = 0.2) +
    theme_bw() + scale_fill_brewer(palette = "Dark2") + 
    theme(axis.title.x = element_blank()) + #ylim(0,0.5) +
    ylab("R^2")# 90th percentile")
  
  model_legend <- cowplot::get_legend(surv_comp_cv_plot)
  add_theme <- theme(
    legend.position = "none", 
    plot.caption = element_text(hjust = 0.5, face = "bold"), 
    axis.text.x = element_text(size = 6))
  combined_comp_plot <- gridExtra::grid.arrange(
    surv_comp_cv_plot + add_theme + labs(caption = "A"), 
    surv_comp_ext_plot + add_theme + labs(caption = "B"), 
    drug_comp_plot + add_theme + labs(caption = "C"), 
    model_legend, 
    ncol = 4)
  
  if (dir.exists(paste0(save_path, "reverse_external_evaluation"))) {
    rev_external_path <- paste0(save_path, "reverse_external_evaluation/")
    best_task_surv_ext_rev <- read.csv(
      paste0(rev_external_path, "external_survival_validation_res.csv"), 
      header = TRUE, row.names = 1)
    
    best_survival_perf <- data.frame(
      Model = rep("MODAE", 4), 
      Set = c("TCGA train", "TCGA test", "TCGA", "SCANB"), #"TCGA (train)", "SCANB (test)"), 
      Scope = rep(c("TCGA cross-validation", "external validation"), c(2,2)), 
      concordance = c(
        mean(best_task_test[["patient_train_surv_c"]]), 
        mean(best_task_test[["patient_test_surv_c"]]), 
        best_task_surv_ext_rev["1", "survival_concordance"], 
        best_task_surv_ext_rev["0", "survival_concordance"]), 
      concordance_sd = c(
        sd(best_task_test[["patient_train_surv_c"]]), 
        sd(best_task_test[["patient_test_surv_c"]]), 
        NA, NA))
  }
  
  out <- list(
    surv_cv = surv_comp_cv_plot, 
    surv_ext = surv_comp_ext_plot, 
    surv_both = surv_comp_plot, 
    drug_cv = drug_comp_plot, 
    all = combined_comp_plot
  )
  return(out)
}