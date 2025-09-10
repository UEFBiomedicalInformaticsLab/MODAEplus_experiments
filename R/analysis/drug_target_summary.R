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
  
  fn <- paste0(ctrp_path, "drug_bank_indications.csv")
  db_indications_df <- readr::read_csv(fn, show_col_types = FALSE)
  fn <- paste0(ctrp_path, "open_targets_indications.csv")
  ot_indications_df <- readr::read_csv(fn, show_col_types = FALSE)
  
  # Cancer sensitivity filter
  fn <- paste0(best_dr_res_path, "drug_target_expression/cancer_sensitivity_table.csv")
  mean_sens_df <- readr::read_csv(fn)
  
  with(
    mean_sens_df, 
    table(
      global_diff = abs(cancer_sensitivity_mean - global_sensitivity_mean) / global_sensitivity_sd > 1, 
      cancer_diff = cancer_sensitivity_sd / global_sensitivity_sd > 0.75
    )
  )
  
  mean_sens_filter_df <- dplyr::filter(
    mean_sens_df, 
    abs(cancer_sensitivity_mean - global_sensitivity_mean) / global_sensitivity_sd > 1 |
      cancer_sensitivity_sd / global_sensitivity_sd > 0.75
  )
  
  # CV drug performance filter
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
  table(unique(mean_sens_df$drug) %in% unique(dr_r2_best$drug))
  setdiff(unique(mean_sens_df$drug), unique(dr_r2_best$drug))
  table(only_alphanumericals(unique(mean_sens_df$drug)) %in% only_alphanumericals(unique(dr_r2_best$drug)))
  
  mean_sens_filter_df[["drug"]] <- only_alphanumericals(mean_sens_filter_df[["drug"]])
  dr_r2_best[["drug"]] <- only_alphanumericals(dr_r2_best[["drug"]])
  
  dr_r2_best_filter <- dplyr::filter(
    dr_r2_best, 
    cl_test > 0.1
  )
  
  # CV cancer specific drug performance
  drug_response_by_cancer_r2 <- read_result_files("drug_response_by_cancer_r2.csv.gz", path = save_path, model_name = model_name)
  drug_response_by_cancer_r2_best <- drug_response_by_cancer_r2[drug_response_by_cancer_r2[["task"]] == best_task,]
  drug_response_by_cancer_r2_best[["drug"]] <- process_ctrp_drug_names(
    tolower(
      fix_ctrp_drug_names(
        drug_response_by_cancer_r2_best[["drug"]], 
        n_drugs = n_drugs
      )
    )
  )
  dr_by_cancer_r2_best_long <- plyr::ddply(
    drug_response_by_cancer_r2_best, 
    c("cancer", "drug", "dataset"), 
    function(x) data.frame(r2 = mean(x$dr_r2))
  )
  dr_by_cancer_r2_best <- reshape2::dcast(dr_by_cancer_r2_best_long, drug ~ cancer + dataset, value.var = "r2")
  dr_by_cancer_r2_best <- reshape2::acast(dr_by_cancer_r2_best_long, drug ~ cancer ~ dataset, value.var = "r2")
  drug_cancer_mean_r2 <- apply(dr_by_cancer_r2_best[,,"cl_test"], 1, mean, na.rm = TRUE)
  drug_cancer_mean_r2[drug_cancer_mean_r2 > 0]
  
  cancer_drug_r2_cor <- cor(dr_by_cancer_r2_best[,,"cl_test"], use = "pairwise.complete.obs", method = "spearman")
  nna_tissue <- !apply(is.na(cancer_drug_r2_cor), 1, all)
  cancer_drug_r2_cor <- cancer_drug_r2_cor[nna_tissue, nna_tissue]
  cancer_hc_order <- hclust(as.dist(1-cancer_drug_r2_cor), method = "average")$order
  cancer_hc_order <- colnames(cancer_drug_r2_cor)[cancer_hc_order]
  
  drug_cancer_r2_cor <- cor(t(dr_by_cancer_r2_best[,,"cl_test"]), use = "pairwise.complete.obs", method = "spearman")
  nna_drugs <- !apply(is.na(drug_cancer_r2_cor), 1, all)
  drug_cancer_r2_cor <- drug_cancer_r2_cor[nna_drugs, nna_drugs]
  drug_cancer_r2_cor_imputed <- impute::impute.knn(drug_cancer_r2_cor)$data
  
  drug_hc_order <- hclust(as.dist(1-drug_cancer_r2_cor_imputed), method = "average")$order
  drug_hc_order <- colnames(drug_cancer_r2_cor)[drug_hc_order]
  
  accurate_cancer_specific_drugs <- apply(dr_by_cancer_r2_best[,,"cl_test"] > 0.5, 1, any, na.rm = TRUE)
  accurate_cancer_specific_drugs <- names(which(accurate_cancer_specific_drugs))
  
  ComplexHeatmap::Heatmap(
    pmax(dr_by_cancer_r2_best[intersect(drug_hc_order, accurate_cancer_specific_drugs), cancer_hc_order, "cl_test"], 0), 
    name = "R²"
  )
  
  # Drug targets
  drug_target_sources <- c(
    "drug_bank", 
    "open_targets", 
    "pharmgkb", 
    "ctd_expr", 
    "ctd_10int"
  )[-4]
  
  dt_path <- paste0(best_dr_res_path, "drug_target_expression/")
  target_res_list <- list()
  
  sensitivity_filter <- FALSE
  plot_drug_target_summaries <- FALSE
  for (filter_results in c(FALSE, TRUE)) {
    for (filter_indications in c(FALSE, TRUE)) {
      drug_target_res_list <- list()
      drug_min_p_list <- list()
      cancer_drug_target_res_list <- list()
      cancer_p_significance_list <-list()
      cancer_drug_min_p_list <- list()
      cancer_drug_min_p_sig_list <- list()
      
      for (i in drug_target_sources) {
        fn <- paste0(dt_path, i, "/t_tested_differences.csv")
        drug_target_res_list[[i]] <- readr::read_csv(fn, show_col_types = FALSE)[,-1]
        drug_target_res_list[[i]] <- dplyr::mutate(
          drug_target_res_list[[i]], 
          drug = only_alphanumericals(drug)
        )
        drug_target_res_list[[i]][["target_source"]] <- i
        drug_target_res_list[[i]] <- na.omit(drug_target_res_list[[i]])
        if (filter_results) {
          drug_target_res_list[[i]] <- dplyr::inner_join(
            drug_target_res_list[[i]], 
            dr_r2_best_filter, 
            by = c("drug")
          )
        }
        
        drug_min_p_list[[i]] <- plyr::ddply(
          drug_target_res_list[[i]], 
          c("drug"), 
          function(x) tibble::tibble(
            min_p = min(x[["p_value"]]), 
            log2FC = x[["log2FC"]][which.min(x[["p_value"]])]
          )
        )
        drug_min_p_list[[i]][["target_source"]] <- i
        
        fn <- paste0(dt_path, i, "/t_tested_differences_by_cancer.csv")
        cancer_drug_target_res_list[[i]] <- readr::read_csv(fn, show_col_types = FALSE)[,-1]
        cancer_drug_target_res_list[[i]] <- dplyr::mutate(
          cancer_drug_target_res_list[[i]], 
          drug = only_alphanumericals(drug)
        )
        cancer_drug_target_res_list[[i]][["target_source"]] <- i
        cancer_drug_target_res_list[[i]] <- na.omit(cancer_drug_target_res_list[[i]])
        if (filter_results) {
          if (sensitivity_filter) {
            cancer_drug_target_res_list[[i]] <- dplyr::inner_join(
              cancer_drug_target_res_list[[i]], 
              mean_sens_filter_df, 
              by = c("drug", "cancer")
            )
          }
          cancer_drug_target_res_list[[i]] <- dplyr::inner_join(
            cancer_drug_target_res_list[[i]], 
            dr_r2_best_filter, 
            by = c("drug")
          )
        }
        if (filter_indications) {
          cancer_drug_target_res_list[[i]] <- dplyr::inner_join(
            cancer_drug_target_res_list[[i]], 
            dplyr::distinct(ot_indications_df[,c("tcga_type", "drug")]), 
            by = dplyr::join_by(drug == drug, cancer == tcga_type)
          )
        }
        
        cancer_p_significance_list[[i]] <- plyr::ddply(
          cancer_drug_target_res_list[[i]], 
          c("cancer", "drug"), 
          function(x) tibble::tibble(
            "p<0.05 fraction" = mean(p.adjust(x[["p_value"]], method = "BH") < 0.05, na.rm = TRUE), 
            "p<0.01 fraction" = mean(p.adjust(x[["p_value"]], method = "BH") < 0.01, na.rm = TRUE), 
            "p<0.001 fraction" = mean(p.adjust(x[["p_value"]], method = "BH") < 0.001, na.rm = TRUE)
          )
        )
        cancer_p_significance_list[[i]][["target_source"]] <- i
        
        cancer_drug_min_p_list[[i]] <- plyr::ddply(
          cancer_drug_target_res_list[[i]], 
          c("cancer", "drug"), 
          function(x) tibble::tibble(
            min_p = min(x[["p_value"]]), 
            log2FC = x[["log2FC"]][which.min(x[["p_value"]])]
          )
        )
        cancer_drug_min_p_list[[i]][["target_source"]] <- i
        
        cancer_drug_min_p_sig_list[[i]] <- plyr::ddply(
          cancer_drug_min_p_list[[i]], 
          "cancer", 
          function(x) tibble::tibble(
            "p<0.05 fraction" = mean(p.adjust(x[["min_p"]], method = "BH") < 0.05, na.rm = TRUE), 
            "p<0.01 fraction" = mean(p.adjust(x[["min_p"]], method = "BH") < 0.01, na.rm = TRUE), 
            "p<0.001 fraction" = mean(p.adjust(x[["min_p"]], method = "BH") < 0.001, na.rm = TRUE)
          )
        )
        cancer_drug_min_p_sig_list[[i]][["target_source"]] <- i
      }
      
      drug_target_res <- dplyr::bind_rows(drug_target_res_list)
      drug_min_p <- dplyr::bind_rows(drug_min_p_list)
      cancer_drug_target_res <- dplyr::bind_rows(cancer_drug_target_res_list)
      cancer_p_significance <- dplyr::bind_rows(cancer_p_significance_list)
      cancer_drug_min_p <- dplyr::bind_rows(cancer_drug_min_p_list)
      cancer_drug_min_p_sig <- dplyr::bind_rows(cancer_drug_min_p_sig_list)
      
      res_str <- paste0(
        ifelse(filter_results, "filtered", "unfiltered"), 
        "_",
        ifelse(filter_indications, "indicated", "unindicated")
      )
      target_res_list[[res_str]] <- list(
        drug_target_res = drug_target_res, 
        drug_min_p = drug_min_p, 
        cancer_drug_target_res = cancer_drug_target_res, 
        cancer_p_significance = cancer_p_significance, 
        cancer_drug_min_p = cancer_drug_min_p, 
        cancer_drug_min_p_sig = cancer_drug_min_p_sig
      )
      
      # Tested across all cancers
      if (plot_drug_target_summaries) {
        p1 <- ggplot(
          drug_target_res, 
          aes(y = pmin(-log10(p_value), 5), 
              x = target_source, 
              color = target_source)
        ) + 
          ggbeeswarm::geom_beeswarm() +
          scale_color_manual(values = pals::kovesi.rainbow(length(drug_target_sources))) + 
          geom_hline(yintercept = -log10(0.05), color = "red") +
          theme_bw() + 
          facet_wrap(drug~.) + 
          theme(
            axis.text.x = element_blank(), #element_text(angle = 90, vjust = 0.5, hjust = 1), 
            axis.title.x = element_blank()
          ) + 
          ggtitle("t-test p-value between pan-cancer sensitive and resistant groups' drug target expression")
        
        save_figure_safe(
          p1, 
          png, 
          paste0(dt_path, res_str, "_drug_target_expression_pan_cancer_t_test.png"), 
          width = plot_width*1.4*ifelse(filter_results | filter_indications, 1, 1.5), 
          height = plot_width*ifelse(filter_results | filter_indications, 1, 1.5), 
          res = plot_res, 
          units = plot_units
        )
        
        p1 <- ggplot(
          drug_target_res, 
          aes(y = abs(log2FC), 
              x = target_source, 
              color = target_source)
        ) + 
          ggbeeswarm::geom_beeswarm() +
          scale_color_manual(values = pals::kovesi.rainbow(length(drug_target_sources))) + 
          geom_hline(yintercept = log2(1.5), color = "red") +
          theme_bw() + 
          facet_wrap(drug~.) + 
          theme(
            axis.text.x = element_blank(), #element_text(angle = 90, vjust = 0.5, hjust = 1), 
            axis.title.x = element_blank()
          ) + 
          ggtitle("effect-size between pan-cancer sensitive and resistant groups' drug target expression")
        
        save_figure_safe(
          p1, 
          png, 
          paste0(dt_path, res_str, "_drug_target_expression_pan_cancer_effect_size.png"), 
          width = plot_width*1.4*ifelse(filter_results | filter_indications, 1, 1.5), 
          height = plot_width*ifelse(filter_results | filter_indications, 1, 1.5), 
          res = plot_res, 
          units = plot_units
        )
        
        p1 <- ggplot(
          drug_min_p, 
          aes(y = pmin(-log10(min_p), 5), 
              x = target_source, 
              color = target_source)
        ) + 
          ggbeeswarm::geom_beeswarm() +
          scale_color_manual(values = pals::kovesi.rainbow(length(drug_target_sources))) + 
          geom_hline(yintercept = -log10(0.05), color = "red") +
          theme_bw() + 
          facet_wrap(drug~.) + 
          theme(
            axis.text.x = element_blank(), #element_text(angle = 90, vjust = 0.5, hjust = 1), 
            axis.title.x = element_blank()
          ) + 
          ggtitle("minimum t-test p-value between pan-cancer sensitive and resistant groups' drug target expression")
        
        save_figure_safe(
          p1, 
          png, 
          paste0(dt_path, res_str, "_minimum_drug_target_expression_pan_cancer_t_test.png"), 
          width = plot_width*1.4*ifelse(filter_results | filter_indications, 1, 1.5), 
          height = plot_width*ifelse(filter_results | filter_indications, 1, 1.5), 
          res = plot_res, 
          units = plot_units
        )
        
        # Tested within cancers
        p1 <- ggplot(
          cancer_drug_target_res, 
          aes(y = pmin(-log10(p_value), 5), 
              x = target_source, 
              color = target_source)
        ) + 
          ggbeeswarm::geom_beeswarm() +
          scale_color_manual(values = pals::kovesi.rainbow(length(drug_target_sources))) + 
          geom_hline(yintercept = -log10(0.05), color = "red") +
          theme_bw() + 
          facet_wrap(drug~.) + 
          theme(
            axis.text.x = element_blank(), #element_text(angle = 90, vjust = 0.5, hjust = 1), 
            axis.title.x = element_blank()
          ) + 
          ggtitle("aggregated t-test p-values between cancer specific sensitive and resistant groups' drug target expression")
        
        save_figure_safe(
          p1, 
          png, 
          paste0(dt_path, res_str, "_drug_target_expression_cancer_specific_t_test_drug_aggregated.png"), 
          width = plot_width*1.4*ifelse(filter_results | filter_indications, 1, 1.5), 
          height = plot_width*ifelse(filter_results | filter_indications, 1, 1.5), 
          res = plot_res, 
          units = plot_units
        )
        
        p1 <- ggplot(
          cancer_p_significance, 
          aes(y = `p<0.05 fraction`, 
              x = target_source, 
              color = target_source)
        ) + 
          ggbeeswarm::geom_beeswarm() +
          scale_color_manual(values = pals::kovesi.rainbow(length(drug_target_sources))) + 
          theme_bw() + 
          facet_wrap(cancer~.) + 
          theme(
            axis.text.x = element_blank(), #element_text(angle = 90, vjust = 0.5, hjust = 1), 
            axis.title.x = element_blank()
          ) + 
          ggtitle("aggregated t-test significance fraction between cancer specific sensitive and resistant groups' drug target expression")
        
        save_figure_safe(
          p1, 
          png, 
          paste0(dt_path, res_str, "_drug_target_expression_cancer_specific_t_test_significance_fraction_cancer_aggregated.png"), 
          width = plot_width*1.4, 
          height = plot_width, 
          res = plot_res, 
          units = plot_units
        )
        
        p1 <- ggplot(
          cancer_drug_target_res, 
          aes(y = pmin(-log10(p_value), 5), 
              x = target_source, 
              color = target_source)
        ) + 
          ggbeeswarm::geom_beeswarm() +
          scale_color_manual(values = pals::kovesi.rainbow(length(drug_target_sources))) + 
          geom_hline(yintercept = -log10(0.05), color = "red") +
          theme_bw() + 
          facet_wrap(cancer~.) + 
          theme(
            axis.text.x = element_blank(), #element_text(angle = 90, vjust = 0.5, hjust = 1), 
            axis.title.x = element_blank()
          ) + 
          ggtitle("aggregated t-test p-values between cancer specific sensitive and resistant groups' drug target expression")
        
        save_figure_safe(
          p1, 
          png, 
          paste0(dt_path, res_str, "_drug_target_expression_cancer_specific_t_test_cancer_aggregated.png"), 
          width = plot_width*1.4, 
          height = plot_width, 
          res = plot_res, 
          units = plot_units
        )
        
        p1 <- ggplot(
          cancer_drug_min_p, 
          aes(y = pmin(-log10(min_p), 5), 
              x = target_source, 
              color = target_source)
        ) + 
          ggbeeswarm::geom_beeswarm() +
          scale_color_manual(values = pals::kovesi.rainbow(length(drug_target_sources))) + 
          geom_hline(yintercept = -log10(0.05), color = "red") +
          theme_bw() + 
          facet_wrap(cancer~.) + 
          theme(
            axis.text.x = element_blank(), #element_text(angle = 90, vjust = 0.5, hjust = 1), 
            axis.title.x = element_blank()
          ) + 
          ggtitle("aggregated minimum t-test p-values between cancer specific sensitive and resistant groups' drug target expression")
        
        save_figure_safe(
          p1, 
          png, 
          paste0(dt_path, res_str, "_minimum_drug_target_expression_cancer_specific_t_test_cancer_aggregated.png"), 
          width = plot_width*1.4, 
          height = plot_width, 
          res = plot_res, 
          units = plot_units
        )
        
        p1 <- ggplot(
          cancer_drug_min_p, 
          aes(y = pmin(-log10(min_p), 5), 
              x = target_source, 
              color = target_source)
        ) + 
          ggbeeswarm::geom_beeswarm() +
          scale_color_manual(values = pals::kovesi.rainbow(length(drug_target_sources))) + 
          geom_hline(yintercept = -log10(0.05), color = "red") +
          theme_bw() + 
          facet_wrap(drug~.) + 
          theme(
            axis.text.x = element_blank(), #element_text(angle = 90, vjust = 0.5, hjust = 1), 
            axis.title.x = element_blank()
          ) + 
          ggtitle("aggregated minimum t-test p-values between cancer specific sensitive and resistant groups' drug target expression")
        
        save_figure_safe(
          p1, 
          png, 
          paste0(dt_path, res_str, "_minimum_drug_target_expression_cancer_specific_t_test_drug_aggregated.png"), 
          width = plot_width*1.4*ifelse(filter_results | filter_indications, 1, 1.5), 
          height = plot_width*ifelse(filter_results | filter_indications, 1, 1.5), 
          res = plot_res, 
          units = plot_units
        )
      }
    }
  }
  
  dt_sig_frac_table_name_map <- c(
    cancer_drug_min_p = "minimum\ncancer-specific\np-value", 
    drug_min_p = "minimum\npan-cancer\np-value", 
    cancer_drug_target_res = "cancer-specific\np-values", 
    drug_target_res = "pan-cancer\np-values"
  )
  dt_sig_frac_subset_name_map <- c(
    filtered_indicated = "R^2, sensitivity, and indication filter", 
    filtered_unindicated = "R^2 and sensitivity filter", 
    unfiltered_indicated = "indication filterg", 
    unfiltered_unindicated = "no filtering"
  )
  target_res_df_list <- list()
  dt_sig_frac_summary_df <- data.frame()
  for (filter_results in c(FALSE, TRUE)) {
    for (filter_indications in c(FALSE, TRUE)) {
      res_str <- paste0(
        ifelse(filter_results, "filtered", "unfiltered"), 
        "_",
        ifelse(filter_indications, "indicated", "unindicated")
      )
      for (j in list(
        list(table = c("cancer_drug_min_p", "drug_min_p"), var = "min_p"),
        list(table = c("cancer_drug_target_res", "drug_target_res"), var = "p_value")
      )) {
        for (k in j[["table"]]) {
          dt_sig_frac_summary_df <- rbind(
            dt_sig_frac_summary_df, 
            data.frame(
              filtering_subset = dt_sig_frac_subset_name_map[res_str], 
              p_value_table = dt_sig_frac_table_name_map[k], 
              plyr::ddply(
                target_res_list[[res_str]][[k]], 
                "target_source", 
                function(x) {
                  data.frame(
                    p_significant_fraction = mean(x[j[["var"]]] < 0.05), 
                    p_and_logfc_significant_fraction = mean(x[j[["var"]]] < 0.05 & abs(x[["log2FC"]]) > log2(1.5))
                  )
                }
              )
            )
          )
        }
      }
      target_res_df_list[[res_str]] <- dplyr::bind_rows(
        target_res_list[[res_str]][["cancer_drug_target_res"]], 
        data.frame(cancer = "pan_cancer", target_res_list[[res_str]][["drug_target_res"]])
      )
      target_res_df_list[[res_str]][["subset"]] <- res_str
    }
  }
  target_res_df <- dplyr::bind_rows(target_res_df_list)
  fn <- paste0(dt_path, "combined_t_tests.csv.gz")
  if (!file.exists(fn)) readr::write_csv(target_res_df, file = gzfile(fn))
  #target_res_df <- readr::read_csv(fn)
  
  dt_sig_frac_summary_df |>
    #dplyr::filter(filtering_subset == "indication filter") |> 
    tidyr::pivot_wider(
      id_cols = c("filtering_subset", "p_value_table"), 
      names_from = "target_source", 
      values_from = "p_significant_fraction", 
    )
  dt_sig_frac_summary_df |>
    #dplyr::filter(filtering_subset == "indication filter") |> 
    tidyr::pivot_wider(
      id_cols = c("filtering_subset", "p_value_table"), 
      names_from = "target_source", 
      values_from = "p_and_logfc_significant_fraction", 
    )
  
  dt_sig_frac_summary_df |>
    #dplyr::filter(filtering_subset == "indication filter") |> 
    tidyr::pivot_wider(
      id_cols = c("target_source", "p_value_table"), 
      names_from = "filtering_subset", 
      values_from = "p_significant_fraction", 
    )
  dt_sig_frac_summary_df |>
    #dplyr::filter(filtering_subset == "indication filter") |> 
    tidyr::pivot_wider(
      id_cols = c("target_source", "p_value_table"), 
      names_from = "filtering_subset", 
      values_from = "p_and_logfc_significant_fraction", 
    )
  
  p1 <- ggplot(
    dt_sig_frac_summary_df, 
    aes(x = target_source, 
        y = p_significant_fraction, 
        fill = target_source)
  ) + 
    geom_bar(stat = "identity") +
    scale_fill_manual(values = pals::kovesi.rainbow(length(drug_target_sources))) + 
    theme_bw() + 
    facet_grid(table~subset) + 
    theme(
      axis.text.x = element_blank(), #element_text(angle = 90, vjust = 0.5, hjust = 1), 
      axis.title.x = element_blank()
    ) + 
    ggtitle("fraction of significant p-values in drug-target expression")
  
  save_figure_safe(
    p1, 
    png, 
    paste0(dt_path, "drug_target_expression_t_test_significance_summary.png"), 
    width = plot_width, 
    height = plot_width, 
    res = plot_res, 
    units = plot_units
  )
  
  p1 <- ggplot(
    dt_sig_frac_summary_df, 
    aes(x = target_source, 
        y = p_and_logfc_significant_fraction, 
        fill = target_source)
  ) + 
    geom_bar(stat = "identity") +
    scale_fill_manual(values = pals::kovesi.rainbow(length(drug_target_sources))) + 
    theme_bw() + 
    facet_grid(table~subset) + 
    theme(
      axis.text.x = element_blank(), #element_text(angle = 90, vjust = 0.5, hjust = 1), 
      axis.title.x = element_blank()
    ) + 
    ggtitle("fraction of significant p-values and effect in drug-target expression")
  
  save_figure_safe(
    p1, 
    png, 
    paste0(dt_path, "drug_target_expression_t_test_significance_summary_with_lfc.png"), 
    width = plot_width, 
    height = plot_width, 
    res = plot_res, 
    units = plot_units
  )
  
  # Overall significance rate of minimum p-value
  with(target_res_list[["filtered_indicated"]][["cancer_drug_min_p"]], tapply(min_p < 0.05, target_source, mean))
  with(target_res_list[["unfiltered_indicated"]][["cancer_drug_min_p"]], tapply(min_p < 0.05, target_source, mean))
  with(target_res_list[["filtered_indicated"]][["drug_min_p"]], tapply(min_p < 0.05, target_source, mean))
  with(target_res_list[["unfiltered_indicated"]][["drug_min_p"]], tapply(min_p < 0.05, target_source, mean))
  
  with(target_res_list[["filtered_indicated"]][["cancer_drug_min_p"]], tapply(min_p < 0.05 & abs(log2FC) > log2(1.5), target_source, mean))
  with(target_res_list[["unfiltered_indicated"]][["cancer_drug_min_p"]], tapply(min_p < 0.05 & abs(log2FC) > log2(1.5), target_source, mean))
  with(target_res_list[["filtered_indicated"]][["drug_min_p"]], tapply(min_p < 0.05 & abs(log2FC) > log2(1.5), target_source, mean))
  with(target_res_list[["unfiltered_indicated"]][["drug_min_p"]], tapply(min_p < 0.05 & abs(log2FC) > log2(1.5), target_source, mean))
  
  # Overall significance rate of all p-values
  with(target_res_list[["filtered_indicated"]][["cancer_drug_target_res"]], tapply(p_value < 0.05, target_source, mean))
  with(target_res_list[["unfiltered_indicated"]][["cancer_drug_target_res"]], tapply(p_value < 0.05, target_source, mean))
  with(target_res_list[["filtered_indicated"]][["drug_target_res"]], tapply(p_value < 0.05, target_source, mean))
  with(target_res_list[["unfiltered_indicated"]][["drug_target_res"]], tapply(p_value < 0.05, target_source, mean))
  
  # Overall significance and effect-size exceeding rate of all p-values
  with(target_res_list[["filtered_indicated"]][["cancer_drug_target_res"]], tapply(p_value < 0.05 & abs(log2FC) > log2(1.5), target_source, mean))
  with(target_res_list[["unfiltered_indicated"]][["cancer_drug_target_res"]], tapply(p_value < 0.05 & abs(log2FC) > log2(1.5), target_source, mean))
  with(target_res_list[["filtered_indicated"]][["drug_target_res"]], tapply(p_value < 0.05 & abs(log2FC) > log2(1.5), target_source, mean))
  with(target_res_list[["unfiltered_indicated"]][["drug_target_res"]], tapply(p_value < 0.05 & abs(log2FC) > log2(1.5), target_source, mean))
}