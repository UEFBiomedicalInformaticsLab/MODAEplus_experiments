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
  
  #pancan_subtypes <- TCGAbiolinks::PanCancerAtlas_subtypes()
  #fn <- paste0(best_dr_res_path, "drug_target_expression/t_tested_differences_by_cancer.csv")
  #cancer_drug_target_expression_tested <- read.csv(fn, header = TRUE, row.names = 1)
  
  fn <- paste0(ctrp_path, "drug_bank_indications.csv")
  db_indications_df <- readr::read_csv(fn, show_col_types = FALSE)
  fn <- paste0(ctrp_path, "open_targets_indications.csv")
  ot_indications_df <- readr::read_csv(fn, show_col_types = FALSE)
  
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
    names(sensitivity_cols) <- only_alphanumericals(sensitivity_drug_names)
    
    heatmap_path <- paste0(plot_path, "sensitivity_heatmaps/")
    dir.create(heatmap_path, recursive = TRUE)
    
    for (canceri in pan_cancer_types) {
      canceri_ind <- with(tcga_final_output_table, type == canceri & !is.na(type))
      annotation_df <- tcga_final_output_table[canceri_ind,]
      #subtype_idx <- match(substr(annotation_df[["id"]], 1, 12), pancan_subtypes[["pan.samplesID"]])
      #annotation_df[["PanCanAtlas_subtype"]] <- pancan_subtypes[subtype_idx, "Subtype_Selected", drop = TRUE]
      annotation_df[["age"]] <- annotation_df[["age_at_initial_pathologic_diagnosis"]]
      annotation_df[["ajcc_stage"]] <- simplify_cancer_stages(annotation_df[["ajcc_pathologic_tumor_stage"]])
      annotation_df[["clinical_stage"]] <- simplify_cancer_stages(annotation_df[["clinical_stage"]])
      annot_cols <- c("gender", "Subtype_Selected", "age", "ajcc_stage", "clinical_stage", "histological_grade")#, "NRF2_score", "NRF2_mutations")
      annot_nna <- sapply(annot_cols, function(x) !all(is.na(annotation_df[[x]])))
      annot_continuous_map <- c(
        gender = FALSE, 
        Subtype_Selected = FALSE, 
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
        Subtype_Selected = RColorBrewer::brewer.pal(annot_size["Subtype_Selected"], "Set1"), 
        ajcc_stage = pals::kovesi.rainbow(annot_size["ajcc_stage"]), 
        clinical_stage = RColorBrewer::brewer.pal(annot_size["clinical_stage"], "Dark2"),
        histological_grade = RColorBrewer::brewer.pal(annot_size["histological_grade"], "Dark2")#, 
        #NRF2_mutations = RColorBrewer::brewer.pal(annot_size["NRF2_mutations"], "PiYG")
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
      #drug_filter <- sensitivity_cols[cancer_i_drugs_of_interest_high_sens]
      library(dplyr)
      drug_filter <- ot_indications_df |> filter(tcga_type == canceri) |> select(drug) |> distinct()
      #drug_filter <- db_indications_df |> filter(type == canceri) |> select(drug) |> distinct()
      drug_filter <- drug_filter[[1]]
      drug_perf_ind <- match(drug_filter, only_alphanumericals(dr_perf_table[["drug"]]))
      drug_perf <- dr_perf_table[["r2_test_mean"]][drug_perf_ind]
      drug_filter <- drug_filter[drug_perf > 0.1]
      drug_filter <- sensitivity_cols[drug_filter]
      drug_filter <- drug_filter[!is.na(drug_filter)]
      
      if(length(drug_filter) == 0) next
      
      sens_limits <- max(annotation_df[,drug_filter], na.rm = TRUE)
      
      sens_limits <- c(0, sens_limits * 0.4, sens_limits * 0.8)
      #sens_limits <- c(0, sens_limits * 0.5, sens_limits * 1.0)
      
      if (length(drug_filter) > 1) {
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
        
      } else {
        sens_order <- 1
        sample_order <- order(annotation_df[[drug_filter]])
      }
      
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
      
      cancer_heatmap_list <- list(cancer_heatmap_annots, cancer_heatmap)
      
      save_figure_safe(
        Reduce(ComplexHeatmap::`%v%`, cancer_heatmap_list), 
        png, 
        #paste0(heatmap_path, canceri, "_sensitivity_heatmap.png"),
        paste0(heatmap_path, canceri, "_sensitivity_heatmap_squashed.png"), 
        width = plot_width * 1.5, 
        height = plot_width * 1.5, 
        res = plot_res, 
        units = plot_units
      )
    }
  }
}