source("../setup.R")

for (save_dir in save_dirs) {
  if(exists("var_list")) try(detach(var_list))
  var_list <- get_paths_and_parameters(save_dir, remove_incomplete = FALSE)
  attach(var_list)
  best_task <- readLines(paste0(save_path, "best_task.txt"))
  param_best_task_ind <- match(best_task, parameters[["task"]])
  shared_embedding_names <- paste0("z", 1:parameters[param_best_task_ind, "bottle_neck"])
  
  best_dr_res_path <- paste0(save_path, "../plots/", model_name, "_drug_responses/")
  dir.create(best_dr_res_path, recursive = TRUE)
  
  pancan_subtypes <- TCGAbiolinks::PanCancerAtlas_subtypes()
  fn <- paste0(best_dr_res_path, "drug_target_expression/t_tested_differences_by_cancer.csv")
  cancer_drug_target_expression_tested <- read.csv(fn, header = TRUE, row.names = 1)
  
  fn <- paste0(ctrp_path, "drug_bank_indications.csv")
  indications_df <- read.csv(fn, row.names = 1, header = TRUE)
  
  # Final predictions
  prediction_paths <- paste0(save_path, "external_evaluation/")
  if (dir.exists(paste0(save_path, "reverse_external_evaluation"))) {
    prediction_paths[2] <- paste0(save_path, "reverse_external_evaluation/")
  }
  for (prediction_path in prediction_paths) {
    fn <- paste0(prediction_path, "internal/unified_patient_output_table.csv.gz")
    tcga_final_output_table <- readr::read_csv(fn)
    fn <- paste0(plot_path, "best_model_drugwise_performance.csv")
    dr_perf_table <- readr::read_csv(fn)
    
    sensitivity_cols <- grep("_sensitivity$", colnames(tcga_final_output_table))
    sensitivity_cols <- colnames(tcga_final_output_table)[sensitivity_cols]
    # should be identical order to names saved in file
    sensitivity_drug_names <- get_ctrp_drugs(n_drugs = length(sensitivity_cols))
    sensitivity_drug_names <- process_ctrp_drug_names(sensitivity_drug_names)
    names(sensitivity_drug_names) <- sensitivity_cols
    names(sensitivity_cols) <- sensitivity_drug_names
    
    fn <- paste0(prediction_path, "internal/NRF2_global_sens_cor.csv.gz")
    try(global_correlation_df <- readr::read_csv(fn))
    fn <- paste0(prediction_path, "internal/NRF2_cancer_sens_cor.csv.gz")
    cancer_correlation_df <- readr::read_csv(fn)
    table(unique(cancer_correlation_df[["treatment"]]) %in% sensitivity_cols)
    cancer_correlation_df[["drug"]] <- sensitivity_drug_names[cancer_correlation_df[["treatment"]]]
    fn <- paste0(prediction_path, "internal/NRF2_cancer_sens_betareg.csv.gz")
    try(cancer_betareg_df <- readr::read_csv(fn))
    fn <- paste0(prediction_path, "internal/NRF2_cancer_sens_group_ttest.csv.gz")
    cancer_t_test_df <- readr::read_csv(fn)
    table(unique(cancer_t_test_df[["treatment"]]) %in% sensitivity_cols)
    cancer_t_test_df[["drug"]] <- sensitivity_drug_names[cancer_t_test_df[["treatment"]]]
    
    # Heatmaps
    heatmap_path <- paste0(plot_path, "nrf2_heatmaps/")
    dir.create(heatmap_path, recursive = TRUE)
    for (canceri in c("LUSC", "LUAD")) {
      cancer_i_cor_df <- cancer_correlation_df[cancer_correlation_df[["cancer"]] == canceri,]
      ggplot(cancer_i_cor_df, aes(rho, color = method)) + 
        geom_density() +
        theme_bw()
      
      cancer_i_cor_df <- cancer_i_cor_df[cancer_i_cor_df[["method"]] == "pearson",]
      cancer_i_cor_df[["method"]] <- NULL
      cancer_i_cor_df[["NRF2_cor"]] <- cancer_i_cor_df[["rho"]]
      cancer_i_cor_df[["rho"]] <- NULL
      cancer_i_cor_df[["NRF2_cor_p"]] <- cancer_i_cor_df[["p.value"]]
      cancer_i_cor_df[["p.value"]] <- NULL
      
      cancer_i_combined_df <- plyr::join(dr_perf_table, cancer_i_cor_df)
      cancer_i_combined_df <- plyr::join(cancer_i_combined_df, cancer_t_test_df[cancer_t_test_df[["cancer"]] == canceri,])
      cancer_i_combined_df <- cancer_i_combined_df[!is.na(cancer_i_combined_df[["cancer"]]), ]
      cancer_i_combined_df[["NRF2_t_test_p"]] <- cancer_i_combined_df[["t.test.p.val"]]
      cancer_i_combined_df[["t.test.p.val"]] <- NULL
      cancer_i_combined_df[["NRF2_effect_size"]] <- cancer_i_combined_df[["effect"]]
      cancer_i_combined_df[["effect"]] <- NULL
      cancer_i_combined_df[["mean_NRF2_in_resistant"]] <- cancer_i_combined_df[["mean.in.group.resistant"]]
      cancer_i_combined_df[["mean.in.group.resistant"]] <- NULL
      cancer_i_combined_df[["mean_NRF2_in_sensitive"]] <- cancer_i_combined_df[["mean.in.group.sensitive"]]
      cancer_i_combined_df[["mean.in.group.sensitive"]] <- NULL
      
      cancer_i_combined_df[["drug_bank_indication"]] <- cancer_i_combined_df[["drug"]] %in% with(indications_df, drug[type == canceri])
      
      ggplot(cancer_i_combined_df, aes(NRF2_effect_size, NRF2_cor)) + 
        geom_point(shape = "+", size = 3) + 
        theme_bw()
      ggplot(
        cancer_i_combined_df, 
        aes(
          x = -log10(p.adjust(NRF2_t_test_p, "BH")), 
          y = -log10(p.adjust(NRF2_cor_p, "BH"))
        )
      ) + 
        geom_point(shape = "+", size = 3) + 
        theme_bw()
      ggplot(cancer_i_combined_df, aes(pmax(0, r2_test_mean), NRF2_cor)) + 
        geom_point(shape = "+", size = 3) + 
        theme_bw()
      
      with(cancer_i_combined_df, table(
        t = p.adjust(NRF2_t_test_p, "BH") < 0.05, 
        cor = p.adjust(NRF2_cor_p, "BH") < 0.05
      ))
      with(cancer_i_combined_df, table(
        significant = p.adjust(NRF2_t_test_p, "BH") < 0.05 & p.adjust(NRF2_cor_p, "BH") < 0.05, 
        confidence = r2_test_mean > 0.1
      ))
      with(cancer_i_combined_df, table(
        significant = p.adjust(NRF2_t_test_p, "BH") < 0.05 & 
          p.adjust(NRF2_cor_p, "BH") < 0.05 & 
          r2_test_mean > 0.1, 
        drug_bank_indication
      ))
      with(cancer_i_combined_df, table(
        significant = p.adjust(NRF2_t_test_p, "BH") < 0.05 & 
          p.adjust(NRF2_cor_p, "BH") < 0.05 & 
          r2_test_mean > 0.1, 
        abs(NRF2_effect_size) >= 1
      ))
      cancer_i_drugs_of_interest_high_effect <- with(cancer_i_combined_df, 
        drug[
          p.adjust(NRF2_t_test_p, "BH") < 0.05 & 
          p.adjust(NRF2_cor_p, "BH") < 0.05 & 
          r2_test_mean > 0.1 &
          abs(NRF2_effect_size) >= 1
        ]
      )
      cancer_i_drugs_of_interest <- with(
        cancer_i_combined_df, 
        drug[
          p.adjust(NRF2_t_test_p, "BH") < 0.05 & 
            p.adjust(NRF2_cor_p, "BH") < 0.05 & 
            r2_test_mean > 0.1 
        ]
      )
      with(cancer_i_combined_df, hist(mean_sensitivity_in_sensitive - mean_sensitivity_in_resistant))
      cancer_i_drugs_of_interest_high_sens <- with(
        cancer_i_combined_df, 
        drug[
          p.adjust(NRF2_t_test_p, "BH") < 0.05 & 
            p.adjust(NRF2_cor_p, "BH") < 0.05 & 
            r2_test_mean > 0.1 &
            mean_sensitivity_in_sensitive - mean_sensitivity_in_resistant > 0.1
        ]
      )
      
      cancer_i_combined_df_best <- cancer_i_combined_df[
        cancer_i_combined_df[["drug"]] %in% cancer_i_drugs_of_interest, 
      ]
      
      hist(with(cancer_i_combined_df, mean_sensitivity_in_sensitive[drug %in% cancer_i_drugs_of_interest]))
      cancer_i_drugs_of_interest_high_sens %in% names(drug_target_genes)
      
      cancer_i_gex <- read.csv(paste0(patient_expression_root_dir, canceri, "/mrna.csv.gz"), row.names = 1, header = TRUE)
      
      drug_filter <- cancer_i_drugs_of_interest_high_effect
      #drug_filter <- cancer_i_drugs_of_interest_high_sens
      
      nrf2_genes <- c("CBR1", "SRXN1", "GCLC", "GCLM", "AKR1C3", "ME1")
      cancer_i_genes <- union(nrf2_genes, Reduce(union, drug_target_genes[drug_filter]))
      cancer_i_genes <- intersect(cancer_i_genes, rownames(cancer_i_gex))
      
      canceri_ind <- with(tcga_final_output_table, type == canceri & !is.na(type))
      annotation_df <- tcga_final_output_table[canceri_ind,]
      subtype_idx <- match(substr(annotation_df[["id"]], 1, 12), pancan_subtypes[["pan.samplesID"]])
      annotation_df[["PanCanAtlas_subtype"]] <- pancan_subtypes[subtype_idx, "Subtype_Selected", drop = TRUE]
      annotation_df[["age"]] <- annotation_df[["age_at_initial_pathologic_diagnosis"]]
      annotation_df[["ajcc_stage"]] <- annotation_df[["ajcc_pathologic_tumor_stage"]]
      annot_cols <- c("gender", "PanCanAtlas_subtype", "age", "ajcc_stage", "clinical_stage", "histological_grade", "NRF2_score", "NRF2_mutations")
      annot_nna <- sapply(annot_cols, function(x) !all(is.na(annotation_df[[x]])))
      annot_continuous_map <- c(
        gender = FALSE, 
        PanCanAtlas_subtype = FALSE, 
        age = TRUE, 
        ajcc_stage = FALSE, 
        clinical_stage = FALSE, 
        histological_grade = FALSE, 
        NRF2_score = TRUE, 
        NRF2_mutations = FALSE
      )
      annot_size <- sapply(annot_cols[!annot_continuous_map[annot_cols]], function(x) length(table(annotation_df[[x]])))
      annot_mintwo <- annot_size > 1
      annot_colors <- list(
        gender = c("#67A9CF", "#EF8A62"),
        PanCanAtlas_subtype = RColorBrewer::brewer.pal(annot_size["PanCanAtlas_subtype"], "Set1"), 
        ajcc_stage = pals::kovesi.rainbow(annot_size["ajcc_stage"]), 
        clinical_stage = RColorBrewer::brewer.pal(annot_size["clinical_stage"], "Dark2"),
        histological_grade = RColorBrewer::brewer.pal(annot_size["histological_grade"], "Dark2"), 
        NRF2_mutations = RColorBrewer::brewer.pal(annot_size["NRF2_mutations"], "PiYG")
      )
      annot_colors <- annot_colors[annot_mintwo]
      annot_final <- annot_cols[annot_nna & (annot_continuous_map[annot_cols] | annot_mintwo[annot_cols])]
      
      cancer_heatmap_annots <- COPS::heatmap_annotations(
        annotation_df[, annot_final, drop = FALSE], 
        factor_color_sets = annot_colors, 
        annotation_legend_param = list(nrow = 6)
      )
      #betareg_filter <- with(cancer_betareg_df, cancer == canceri)
      #drug_filter <- with(
      #  cancer_betareg_df[betareg_filter, , drop = FALSE], 
      #  treatment[abs(effect) > 0.05 & p.value < 0.05]
      #)
      #drug_filter <- drug_filter[!is.na(drug_filter)]
      #drug_filter <- paste0(known_cancer_drugs[["drug"]], "_sensitivity")
      drug_filter <- sensitivity_cols[drug_filter]
      sens_limits <- max(annotation_df[,drug_filter], na.rm = TRUE)
      #sens_limits <- c(0, sens_limits * 0.4, sens_limits * 0.8)
      sens_limits <- c(0, sens_limits * 0.5, sens_limits * 1.0)
      
      sens_cor <- cor(
        annotation_df[,drug_filter], 
        method = "pearson", 
        use = "pairwise.complete.obs"
      )
      sens_cor[is.na(sens_cor)] <- 0
      sens_order <- hclust(as.dist(1-sens_cor), method = "average")$order
      
      sample_cor <- cor(
        t(annotation_df[,drug_filter]), 
        method = "pearson", 
        use = "pairwise.complete.obs"
      )
      sample_cor[is.na(sample_cor)] <- 0
      sample_order <- hclust(as.dist(1-sample_cor), method = "average")$order
      
      cancer_heatmap <- ComplexHeatmap::Heatmap(
        #t(prediction_mat_patient)[brca_drug_correlation_order, brca_ind][minmax_drug_ind,], 
        t(annotation_df[sample_order, drug_filter[sens_order]]),
        name = "predicted DSS", 
        cluster_rows = FALSE, 
        cluster_columns = FALSE, 
        clustering_distance_columns = "pearson", 
        clustering_method_columns = "average", 
        show_column_names = FALSE, 
        show_column_dend = FALSE, 
        col = circlize::colorRamp2(sens_limits, c("blue", "white", "red"))
      )
      
      gex_mat <- log2(as.matrix(cancer_i_gex[cancer_i_genes, gsub("-", ".", annotation_df[["id"]][sample_order])])+1)
      #gex_lims <- c(0, median(gex_mat), max(gex_mat))
      gex_lims <- c(-2, 0, 2)
      gex_heatmap <- ComplexHeatmap::Heatmap(
        t(scale(t(gex_mat))), 
        name = "log2(TPM+1) z-score", 
        cluster_rows = TRUE, 
        cluster_columns = FALSE, 
        clustering_distance_columns = "pearson", 
        clustering_method_columns = "average", 
        show_column_names = FALSE, 
        show_column_dend = FALSE, 
        col = circlize::colorRamp2(gex_lims, c("blue", "white", "red"))
      )
      
      cancer_heatmap_list <- list(cancer_heatmap_annots, cancer_heatmap, gex_heatmap)
      
      save_figure_safe(
        Reduce(ComplexHeatmap::`%v%`, cancer_heatmap_list), 
        png, 
        paste0(heatmap_path, canceri, "_sensitivity_heatmap.png"), 
        width = plot_width * 1.5, 
        height = plot_width * 1.5, 
        res = plot_res, 
        units = plot_units
      )
    }
  }
}