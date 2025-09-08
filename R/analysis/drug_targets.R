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
  
  high_confidence_drugs <- dr_r2_best[dr_r2_best[["cl_test"]] > 0.1, "drug"]
  
  # Final predictions
  int_pred <- get_final_internal_predictions()
  internal_predictions_patient <- int_pred[["patient"]]
  best_patient_predictions <- internal_predictions_patient
  
  #gex_fn <- paste0(save_path, "external_evaluation/internal/patient_mrna.csv.gz")
  #gene_expression_internal_filtered <- readr::read_csv(gex_fn)
  #gene_expression_internal_filtered <- as.data.frame(gene_expression_internal_filtered)
  #rownames(gene_expression_internal_filtered) <- gene_expression_internal_filtered[[1]]
  #gene_expression_internal_filtered[[1]] <- NULL
  ##gex_fn <- paste0(patient_expression_root_dir, "mrna_expression.csv.gz")
  ##gene_expression_internal <- read.csv(gex_fn, row.names = 1 , header = TRUE)
  
  gex_list <- lapply(
    pan_cancer_types, 
    function(x) readr::read_csv(
      paste0(patient_expression_root_dir, x, "/mrna.csv.gz"), 
      show_col_types = FALSE
    )
  )
  names(gex_list) <- pan_cancer_types
  gex_fun <- function(x) {
    ids <- colnames(x)[-1]
    rows <- x[[1]]
    x[[1]] <- NULL
    x <- t(as.matrix(x))
    colnames(x) <- rows
    rownames(x) <- ids
    return(x)
  }
  gex_list_processed <- lapply(gex_list, gex_fun)
  cancer_type_map <- c()
  for (i in names(gex_list_processed)) {
    cancer_type_map[rownames(gex_list_processed[[i]])] <- i
  }
  gene_expression_internal <- Reduce(rbind, gex_list_processed)
  rm(gex_list_processed)
  rownames(gene_expression_internal) <- gsub("-", ".", rownames(gene_expression_internal))
  #rownames(gene_expression_internal) <- gene_expression_internal[[1]]
  #gene_expression_internal[[1]] <- NULL
  
  if (FALSE) {
    # Check data consistency (some differences expected due to updated symbols)
    tcga_surv_data <- read.csv(
      paste0(patient_expression_root_dir, "survival.csv.gz"), 
      header = TRUE, row.names = 1)
    
    ids <- gsub("-", ".", rownames(gene_expression_internal_filtered))
    id_ind <- ids %in% rownames(gene_expression_internal)
    table(id_ind)
    
    type_from_survival <- tcga_surv_data[
      match(
        substr(rownames(gene_expression_internal_filtered), 1, 12), 
        tcga_surv_data[["bcr_patient_barcode"]]
      ), 
      "type"
    ]
    table(type_from_survival[!id_ind], useNA = "always")
    
    shared_genes <- intersect(
      colnames(gene_expression_internal), 
      colnames(gene_expression_internal_filtered)
    )
    
    dat1 <- gene_expression_internal[ids[id_ind], shared_genes]
    dat1 <- as.matrix(dat1)
    dat1 <- log2(dat1 + 1)
    dat2 <- as.matrix(gene_expression_internal_filtered[id_ind, shared_genes])
    hist(dat1 - dat2)
  }
  
  best_patient_dr_predictions <- best_patient_predictions[,grep("dr_pred_[0-9]+", colnames(best_patient_predictions))]
  rownames(best_patient_dr_predictions) <- best_patient_predictions[["X"]]
  ctrp_drugs <- get_ctrp_drugs(n_drugs = ncol(best_patient_dr_predictions))
  
  # Cancerwise mean sensitivity
  fn <- paste0(
    best_dr_res_path, 
    "drug_target_expression/cancer_sensitivity_table.csv"
  )
  if (file.exists(fn)) {
    mean_sens_df <- mean_sens_df <- readr::read_csv(fn)
  } else {
    mean_sens_df <- data.frame()
    for (drugi in ctrp_drugs) {
      sensi <- best_patient_dr_predictions[
        , 
        which(tolower(ctrp_drugs) == tolower(drugi)), 
        drop = FALSE
      ]
      sensi_mu <- mean(sensi[,1], na.rm = TRUE)
      sensi_sd <- sd(sensi[,1], na.rm = TRUE)
      for (canceri in pan_cancer_types) {
        pt_idx <- colnames(gex_list[[canceri]])[-1] # first col is gene id (rownames)
        pt_idx <- gsub("\\.", "-", pt_idx)
        pt_idx <- intersect(pt_idx, rownames(sensi))
        
        mean_sens_df <- rbind(
          mean_sens_df, 
          data.frame(
            drug = drugi, 
            cancer = canceri, 
            global_sensitivity_mean = sensi_mu, 
            global_sensitivity_sd = sensi_sd, 
            cancer_sensitivity_mean = mean(sensi[pt_idx,1]), 
            cancer_sensitivity_sd = sd(sensi[pt_idx,1])
          )
        )
      }
    }
    dir.create(
      paste0(best_dr_res_path, "drug_target_expression/"), 
      recursive = TRUE
    )
    readr::write_csv(mean_sens_df, file = fn)
  }
  
  # Drug targets
  drug_target_sources <- c(
    "drug_bank", 
    "open_targets", 
    "pharmgkb", 
    "ctd_expr", 
    "ctd_10int"
  )[-4]
  
  save_drug_target_figures <- TRUE
  # Regression based associations (just correlation plots for now)
  for (drug_target_source in drug_target_sources) {
    drug_target_genes_json <- readLines(paste0(ctrp_path, drug_target_source, "_drug_target_genes.json"))
    drug_target_genes <- rjson::fromJSON(drug_target_genes_json)
    names(drug_target_genes) <- only_alphanumericals(names(drug_target_genes))
    ctrp_drugs_adjusted <- only_alphanumericals(ctrp_drugs)
    ctrp_drugs_with_targets <- which(ctrp_drugs_adjusted %in% names(drug_target_genes))
    
    drug_target_expression_list <- list()
    dir.create(
      paste0(best_dr_res_path, "drug_target_expression_regression/", drug_target_source), 
      recursive = TRUE
    )
    
    for (drugi in ctrp_drugs_with_targets) {
      targetsi <- drug_target_genes[[ctrp_drugs_adjusted[drugi]]]
      targetsi <- targetsi[!is.na(targetsi)]
      targetsi <- targetsi[targetsi %in% colnames(gene_expression_internal)]
      if (length(targetsi)>0) {
        sensi <- best_patient_dr_predictions[, which(tolower(ctrp_drugs) == tolower(ctrp_drugs[drugi])), drop = FALSE]
        if (ncol(sensi) > 0) {
          # Should be true for all but Tipifarnib
          if (ncol(sensi) > 1) stop("More than one prediction!")
          sens_label_df <- data.frame(
            patient = gsub("\\.", "-", rownames(best_patient_dr_predictions)), 
            sensitivity = sensi[[1]])
          gex_mat <- as.matrix(gene_expression_internal[, targetsi, drop = FALSE])
          gex_df <- reshape2::melt(gex_mat, varnames = c("patient", "gene"), value.name = "TPM")
          gex_df[["cancer"]] <- cancer_type_map[gex_df[["patient"]]]
          gex_df[["patient"]] <- gsub("\\.", "-", gex_df[["patient"]])
          gex_df <- plyr::join(gex_df, sens_label_df, by = "patient", type = "inner")
          gex_df[["drug"]] <- ctrp_drugs[drugi]
          
          
          if (TRUE || tolower(ctrp_drugs[drugi]) %in% drugs_of_interest) {
            #png(paste0(best_dr_res_path, "drug_target_expression/", drug_target_source, "/", ctrp_drugs[drugi], ".png"), 
            #    width = plot_width, height = plot_height, res = plot_res, units = plot_units)
            if (FALSE) {
              # Beta regression for testing sensitivity ~ expression associations.
              # Test code using heavy computation to account for zero-inflation
              # in sensitity values. 
              for (targeti in targetsi) {
                mod1 <- betareg::betareg(sensitivity ~ log2(TPM+1):cancer + cancer, data = dplyr::filter(gex_df, gene == targeti))
                mod2 <- brms::brm(
                  brms::bf(
                    sensitivity ~ (1+log2(TPM+1) | cancer), 
                    zi ~ (1 | cancer)
                  ), 
                  data = dplyr::filter(gex_df, gene == targeti), 
                  family = brms::zero_inflated_beta(), 
                  chains = 4, 
                  iter = 2000, 
                  warmup = 1000, 
                  cores = 10, 
                  seed = 0, 
                  backend = "rstan"
                )
              }
            }
            
            if (save_drug_target_figures) {
              save_figure_safe(
                ggplot(
                  gex_df, 
                  aes(x = log2(TPM+1), y = sensitivity, color = cancer)
                ) + 
                  geom_point(shape = 3) + 
                  ggpubr::stat_cor(method = "pearson") + 
                  theme_bw() + 
                  theme(axis.text.x = element_blank()) + 
                  scale_color_manual(values = pals::kovesi.rainbow(length(pan_cancer_types))) + 
                  facet_wrap(gene ~ ., scales = "free_x"), 
                png, 
                paste0(
                  best_dr_res_path, "drug_target_expression_regression/", 
                  drug_target_source, "/", ctrp_drugs[drugi], ".png"), 
                width = plot_width, 
                height = plot_height, 
                res = plot_res, 
                units = plot_units
              )
            }
          }
        }
      }
    }
  }
  
  # Group based on predicted sensitivity
  sens_func <- function(x) {
    if (FALSE) {
      sens_q <- quantile(x, probs = c(0.05, 0.95))
      sens_label <- ifelse(x < sens_q[1], "resistant", "average")
      sens_label <- ifelse(x > sens_q[2], "sensitive", sens_label)
    } else if(FALSE) {
      sens_label <- ifelse(scale(x) < -1, "resistant", "average")
      sens_label <- ifelse(scale(x) > 1, "sensitive", sens_label)
    } else {
      sens_label <- ifelse(
        x <= mean(x), 
        "resistant", 
        "sensitive"
      )
    }
    return (sens_label)
  }
  
  for (drug_target_source in drug_target_sources) {
    drug_target_genes_json <- readLines(paste0(ctrp_path, drug_target_source, "_drug_target_genes.json"))
    drug_target_genes <- rjson::fromJSON(drug_target_genes_json)
    names(drug_target_genes) <- only_alphanumericals(names(drug_target_genes))
    ctrp_drugs_adjusted <- only_alphanumericals(ctrp_drugs)
    ctrp_drugs_with_targets <- which(ctrp_drugs_adjusted %in% names(drug_target_genes))
    
    drug_target_expression_list <- list()
    cancer_drug_target_expression_list <- list()
    dir.create(
      paste0(best_dr_res_path, "drug_target_expression/", drug_target_source), 
      recursive = TRUE
    )
    
    for (drugi in ctrp_drugs_with_targets) {
      targetsi <- drug_target_genes[[ctrp_drugs_adjusted[drugi]]]
      targetsi <- targetsi[!is.na(targetsi)]
      targetsi <- targetsi[targetsi %in% colnames(gene_expression_internal)]
      if (length(targetsi)>0) {
        sensi <- best_patient_dr_predictions[, which(tolower(ctrp_drugs) == tolower(ctrp_drugs[drugi])), drop = FALSE]
        if (ncol(sensi) > 0) {
          # Should be true for all but Tipifarnib
          if (ncol(sensi) > 1) stop("More than one prediction!")
          sens_label <- sens_func(sensi[[1]])
          
          sens_label_df <- data.frame(
            patient = gsub("\\.", "-", rownames(best_patient_dr_predictions)), 
            sensitivity = sensi[[1]], 
            sensitivity_class = factor(sens_label, levels = c("resistant", "average", "sensitive")))
          gex_mat <- as.matrix(gene_expression_internal[, targetsi, drop = FALSE])
          gex_df <- reshape2::melt(gex_mat, varnames = c("patient", "gene"), value.name = "TPM")
          gex_df[["patient"]] <- gsub("\\.", "-", gex_df[["patient"]])
          gex_df <- plyr::join(gex_df, sens_label_df, by = "patient", type = "inner")
          gex_df[["drug"]] <- ctrp_drugs[drugi]
          
          if (TRUE || tolower(ctrp_drugs[drugi]) %in% drugs_of_interest) {
            #png(paste0(best_dr_res_path, "drug_target_expression/", drug_target_source, "/", ctrp_drugs[drugi], ".png"), 
            #    width = plot_width, height = plot_height, res = plot_res, units = plot_units)
            if (save_drug_target_figures) {
              save_figure_safe(
                ggplot(gex_df, aes(y = log2(TPM+1), fill = sensitivity_class)) + 
                  geom_boxplot() + theme_bw() + theme(axis.text.x = element_blank()) + 
                  scale_fill_brewer(palette = "Dark2") + 
                  facet_wrap(gene ~ ., scales = "free_y"), 
                png, 
                paste0(best_dr_res_path, "drug_target_expression/", drug_target_source, "/", ctrp_drugs[drugi], ".png"), 
                width = plot_width, 
                height = plot_height, 
                res = plot_res, 
                units = plot_units
              )
            }
            cancer_drug_target_expression_list[[ctrp_drugs[drugi]]] <- list()
            for (canceri in pan_cancer_types) {
              pt_idx <- colnames(gex_list[[canceri]])[-1] # first col is gene id (rownames)
              pt_idx <- gsub("\\.", "-", pt_idx)
              pt_idx <- intersect(pt_idx, rownames(sensi))
              
              canceri_path <- paste0(best_dr_res_path, "drug_target_expression/", drug_target_source, "/", canceri, "/")
              dir.create(canceri_path, recursive = TRUE)
              
              c_sens <- sensi[pt_idx, 1]
              temp <- sens_func(c_sens)
              temp <- factor(temp, levels = c("resistant", "average", "sensitive"))
              #gex_df_ptr <- match(pt_idx, gex_df[["patient"]])
              temp <- data.frame(patient = pt_idx, sensitivity = c_sens, sensitivity_class = temp)
              temp <- plyr::join(
                temp, 
                gex_df[, !colnames(gex_df) %in% c("sensitivity", "sensitivity_class"), drop = FALSE], 
                by = "patient"
              )
              if (save_drug_target_figures) {
                save_figure_safe(
                  ggplot(temp, aes(y = log2(TPM+1), fill = sensitivity_class)) + 
                    geom_boxplot() + theme_bw() + theme(axis.text.x = element_blank()) + 
                    scale_fill_brewer(palette = "Dark2") + 
                    facet_wrap(gene ~ ., scales = "free_y"), 
                  png, 
                  paste0(canceri_path, ctrp_drugs[drugi], ".png"), 
                  width = plot_width, 
                  height = plot_height, 
                  res = plot_res, 
                  units = plot_units
                )
                save_figure_safe(
                  ggplot(temp, aes(x = sensitivity, y = log2(TPM+1))) + 
                    geom_point(shape = "+", size = 3) + 
                    theme_bw() + 
                    ggpubr::stat_cor(aes(color = "red")) +
                    scale_fill_brewer(palette = "Dark2") + 
                    facet_wrap(gene ~ ., scales = "free_y") + 
                    theme(legend.position = "none"), 
                  png, 
                  paste0(canceri_path, ctrp_drugs[drugi], "_cor.png"), 
                  width = plot_width, 
                  height = plot_height, 
                  res = plot_res, 
                  units = plot_units
                )
              }
              
              temp[["cancer"]] <- canceri
              
              cancer_drug_target_expression_list[[ctrp_drugs[drugi]]] <- c(
                cancer_drug_target_expression_list[[ctrp_drugs[drugi]]], 
                list(temp)
              )
            }
          }
        }
        
        drug_target_expression_list[[ctrp_drugs[drugi]]] <- gex_df
      } else {
        warning(paste0(
          "None of the targets (", 
          paste(targetsi, collapse = ", "), 
          ") were found in patient gene-expression data."))
      }
    }
    drug_target_expression_df <- do.call(
      dplyr::bind_rows, 
      args = drug_target_expression_list
    )
    drug_target_expression_df <- dplyr::filter(
      drug_target_expression_df,
      !is.na(sensitivity_class)
    )
    cancer_drug_target_expression_list <- lapply(
      cancer_drug_target_expression_list, 
      function(x) do.call(dplyr::bind_rows, args = x)
    )
    cancer_drug_target_expression_df <- do.call(
      dplyr::bind_rows, 
      args = cancer_drug_target_expression_list
    )
    cancer_drug_target_expression_df <- dplyr::filter(
      cancer_drug_target_expression_df,
      !is.na(sensitivity_class)
    )
    
    cancer_drug_target_expression_cor_df <- plyr::ddply(
      cancer_drug_target_expression_df, 
      c("drug", "gene", "cancer"), 
      function(x) {data.frame(cor = with(x, cor(TPM, sensitivity, method = "pearson")))}
    )
    ggplot(cancer_drug_target_expression_cor_df, aes(cor)) +
      geom_density() + 
      theme_bw() +
      facet_wrap(cancer ~.)
    
    drugs_of_interest <- high_confidence_drugs
    plot1 <- ggplot(drug_target_expression_df[tolower(drug_target_expression_df[["drug"]]) %in% drugs_of_interest, ], 
                    aes(x = drug, y = log2(TPM+1), fill = sensitivity_class, group = interaction(sensitivity_class, drug))) + 
      geom_boxplot() + theme_bw() + 
      scale_fill_brewer(palette = "Dark2") + 
      facet_wrap(gene ~ ., scales = "free")
    
    drug_target_expression_tested <- plyr::ddply(
      drug_target_expression_df, 
      c("drug", "gene"), 
      function(x) {
        sensitive <- log2(x[x[["sensitivity_class"]] == "sensitive", "TPM"]+1)
        resistant <- log2(x[x[["sensitivity_class"]] == "resistant", "TPM"]+1)
        sensitive_n = length(sensitive)
        resistant_n = length(resistant)
        if (sensitive_n > 1 & resistant_n > 1) {
          out <- data.frame(
            drug = x[1,"drug"], 
            gene = x[1,"gene"], 
            sensitive_mean = mean(sensitive), 
            resistant_mean = mean(resistant), 
            sensitive_n = sensitive_n, 
            resistant_n = resistant_n, 
            p_value = t.test(
              x = sensitive, 
              y = resistant, 
              alternative = "two.sided")[["p.value"]])
        } else {
          out <- data.frame()
        }
        return(out)
      })
    p_counts <- table(drug_target_expression_tested[["p_value"]] < 0.05)
    #p_counts / sum(p_counts)
    drug_target_expression_tested[["log2FC"]] <- with(drug_target_expression_tested, resistant_mean - sensitive_mean)
    
    hist(drug_target_expression_tested[["resistant_n"]])
    hist(drug_target_expression_tested[["sensitive_n"]])
    plot(drug_target_expression_tested[,c("sensitive_n", "resistant_n")])
    plot(
      drug_target_expression_tested[["log2FC"]], 
      -log10(drug_target_expression_tested[["p_value"]])
    )
    #drug_target_expression_tested[
    #  drug_target_expression_tested[["sensitive_n"]] > 400 &
    #    drug_target_expression_tested[["resistant_n"]] > 400,
    #]
    fn <- paste0(best_dr_res_path, "drug_target_expression/", drug_target_source, "/t_tested_differences.csv")
    write.csv(drug_target_expression_tested, file = fn)
    
    temp <- drug_target_expression_tested[
      , c("drug", "gene", "log2FC", "p_value")]
    
    temp[["p_value"]] <- ifelse(
      temp[["p_value"]] < 0.001, 
      "<0.001", 
      as.character(round(temp[["p_value"]], 3)))
    
    print(
      xtable::xtable(temp, digits = 3), 
      include.rownames=FALSE, 
      type = "latex"
    )
    
    # By cancer type
    cancer_drug_target_expression_tested <- plyr::ddply(
      cancer_drug_target_expression_df, 
      c("drug", "gene", "cancer"), 
      function(x) {
        sensitive <- log2(x[x[["sensitivity_class"]] == "sensitive", "TPM"]+1)
        resistant <- log2(x[x[["sensitivity_class"]] == "resistant", "TPM"]+1)
        sensitive_n = length(sensitive)
        resistant_n = length(resistant)
        if (sensitive_n > 1 & resistant_n > 1) {
          out <- data.frame(
            drug = x[1,"drug"], 
            gene = x[1,"gene"], 
            sensitive_mean = mean(sensitive), 
            resistant_mean = mean(resistant), 
            sensitive_n = sensitive_n, 
            resistant_n = resistant_n, 
            p_value = t.test(
              x = sensitive, 
              y = resistant, 
              alternative = "two.sided")[["p.value"]])
        } else {
          out <- data.frame()
        }
        return(out)
      })
    cancer_p_counts <- table(cancer_drug_target_expression_tested[["p_value"]] < 0.05)
    #cancer_p_counts / sum(cancer_p_counts)
    cancer_drug_target_expression_tested[["log2FC"]] <- with(
      cancer_drug_target_expression_tested, 
      resistant_mean - sensitive_mean
    )
    
    fn <- paste0(best_dr_res_path, "drug_target_expression/", drug_target_source, "/t_tested_differences_by_cancer.csv")
    write.csv(cancer_drug_target_expression_tested, file = fn)
    
    ggplot(cancer_drug_target_expression_tested, aes(x = log2FC, y = -log10(p_value))) + 
      geom_point(shape = "+", size = 3) + 
      theme_bw() + 
      facet_wrap(cancer ~ .)
    
    hcd_ind <- cancer_drug_target_expression_tested[["drug"]] %in% high_confidence_drugs
    ggplot(cancer_drug_target_expression_tested[hcd_ind,], aes(x = log2FC, y = -log10(p_value))) + 
      geom_point(shape = "+", size = 3) + 
      theme_bw() + 
      facet_wrap(cancer ~ .)
  }
  
  
}