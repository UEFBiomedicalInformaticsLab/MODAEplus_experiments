script_root <- Sys.getenv("MODAE_SCRIPT_PATH")
if (script_root == "") {
  stop("Please define the MODAE_SCRIPT_PATH environment variable.")
}
source(paste0(script_root, "R/setup.R"))
library(dplyr)

for (save_dir in save_dirs) {
  if(exists("var_list")) try(detach(var_list))
  var_list <- get_paths_and_parameters(save_dir, remove_incomplete = FALSE)
  attach(var_list)
  best_task <- readLines(paste0(save_path, "best_task.txt"))
  param_best_task_ind <- match(best_task, parameters[["task"]])
  shared_embedding_names <- paste0("z", 1:parameters[param_best_task_ind, "bottle_neck"])
  
  best_dr_res_path <- paste0(save_path, "../plots/", model_name, "_drug_responses/")
  dir.create(best_dr_res_path, recursive = TRUE)
  
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
    tcga_final_output_table <- readr::read_csv(fn, show_col_types = FALSE)
    
    sensitivity_cols <- grep("_sensitivity$", colnames(tcga_final_output_table))
    sensitivity_cols <- colnames(tcga_final_output_table)[sensitivity_cols]
    sensitivity_drug_names <- get_ctrp_drugs(n_drugs = length(sensitivity_cols))
    sensitivity_drug_names <- process_ctrp_drug_names(sensitivity_drug_names)
    names(sensitivity_drug_names) <- sensitivity_cols
    names(sensitivity_cols) <- only_alphanumericals(sensitivity_drug_names)
    
    if (pan_cancer) {
      cancer_drug_pairs <- ot_indications_df |> 
        filter(only_alphanumericals(drug) %in% names(sensitivity_cols)) |>
        select(tcga_type, drug) |>
        distinct()
      for (cancer in unique(cancer_drug_pairs[["tcga_type"]])) {
        e <- tcga_final_output_table |>
          filter(type == cancer)
        e <- e[,shared_embedding_names]
        set.seed(0)
        res_umap <- uwot::umap(
          e, 
          n_components = 2, 
          n_neighbors = 20, 
          pca = NULL, 
          init = "normlaplacian", 
          min_dist = 0.5, 
          nn_method = "annoy", 
          metric = "correlation",
          verbose = FALSE
        )
        cancer_path <- paste0(
          plot_path, 
          "sensitivity_umaps/", 
          cancer, 
          "/"
        )
        dir.create(cancer_path, recursive = TRUE)
        
        cancer_drugs <- cancer_drug_pairs |> 
          filter(tcga_type == cancer)
        for (drug in unique(cancer_drugs[["drug"]])) {
          p <- tcga_final_output_table |>
            filter(type == cancer)
          p <- p[[sensitivity_cols[drug]]]
          umap_df <- data.frame(
            UMAP1 = res_umap[,1], 
            UMAP2 = res_umap[,2], 
            sensitivity = p
          )
          
          sens_limits <- c(0, max(p))
          
          p1 <- ggplot(umap_df, aes(UMAP1, UMAP2, color = sensitivity)) + 
            geom_point(shape = 3) + 
            theme_bw() +
            scale_color_distiller(
              palette = "RdBu", 
              limits = sens_limits
            ) + 
            theme(
              axis.title.x = element_blank(), 
              axis.title.y = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.ticks.y = element_blank(), 
              axis.text.x = element_blank(), 
              axis.text.y = element_blank()
            )
          save_figure_safe(
            p1, 
            png, 
            paste0(cancer_path, drug, ".png"), 
            width = plot_width * 0.65, 
            height = plot_width * 0.5, 
            res = plot_res, 
            units = plot_units
          )
        }
      }
    }
    patient_output <- tcga_final_output_table |>
      filter(!grepl("^ACH", id))
    e <- patient_output[,shared_embedding_names]
    set.seed(0)
    res_umap <- uwot::umap(
      e, 
      n_components = 2, 
      n_neighbors = 20, 
      pca = NULL, 
      init = "normlaplacian", 
      min_dist = 0.5, 
      nn_method = "annoy", 
      metric = "correlation",
      verbose = FALSE
    )
    umap_path <- paste0(
      plot_path, 
      "sensitivity_umaps/"
    )
    for (drug in names(sensitivity_cols)) {
      p <- patient_output[[sensitivity_cols[drug]]]
      umap_df <- data.frame(
        UMAP1 = res_umap[,1], 
        UMAP2 = res_umap[,2], 
        sensitivity = p
      )
      sens_limits <- c(0, max(p))
      p1 <- ggplot(umap_df, aes(UMAP1, UMAP2, color = sensitivity)) + 
        geom_point(shape = 3) + 
        theme_bw() +
        scale_color_distiller(
          palette = "RdBu", 
          limits = sens_limits
        ) + 
        theme(
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.text.x = element_blank(), 
          axis.text.y = element_blank()
        )
      save_figure_safe(
        p1, 
        png, 
        paste0(umap_path, drug, ".png"), 
        width = plot_width * 1.2, 
        height = plot_width * 0.9, 
        res = plot_res, 
        units = plot_units
      )
    }
  }
}