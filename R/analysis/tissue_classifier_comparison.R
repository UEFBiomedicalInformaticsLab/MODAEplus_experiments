source("../setup.R")
source("comparison_function.R")

for (save_dir in save_dirs) {
  if(exists("var_list")) try(detach(var_list))
  var_list <- get_paths_and_parameters(save_dir, remove_incomplete = FALSE)
  attach(var_list)
  
  baseline_correction_path <- paste0(
    output_dir, 
    "/baseline_results/batch_correction/pan_cancer/"
  )
  fn <- paste0(baseline_correction_path, "metrics_list.rds")
  baseline_tissue_metrics <- readRDS(fn)
  for(i in names(baseline_tissue_metrics)) {
    baseline_tissue_metrics[[i]][["model"]] <- i
  }
  
  fn <- paste0(plot_path, "tissue_metrics.csv")
  modae_tissue_metrics <- readr::read_csv(fn, show_col_types = FALSE)
  modae_tissue_metrics[["model"]] <- "MODAE"
  
  super_dirs <- dir(paste0(save_path, "../"))
  ablation_dir <- super_dirs[grepl("ablation", super_dirs)]
  ablation_tissue_metrics <- list()
  if (length(ablation_dir) > 0) {
    ablation_strings <- dir(paste0(save_path, "../", ablation_dir))
    file_extensions <- paste(
      paste0("\\.", c("tfrecords", "json", "csv", "csv\\.gz"), "$"), 
      collapse = "|"
    )
    ablation_strings <- ablation_strings[!grepl(file_extensions, ablation_strings)]
    ablation_strings <- ablation_strings[ablation_strings != "plots"]
    for (abl_str in ablation_strings) {
      fn <- paste0(
        save_path, "../", ablation_dir, 
        "/plots/", model_name, "_", abl_str, "_",
        "tissue_metrics.csv"
      )
      ablation_tissue_metrics[[abl_str]] <- readr::read_csv(
        fn, 
        show_col_types = FALSE
      )
      ablation_tissue_metrics[[abl_str]][["model"]] <- paste0("MODAE_", abl_str)
    }
  }
  tissue_metrics_df <- dplyr::bind_rows(
    dplyr::bind_rows(baseline_tissue_metrics), 
    dplyr::bind_rows(ablation_tissue_metrics), 
    modae_tissue_metrics
  )
  # Fix bad names
  tissue_metrics_df[["model"]] <- ifelse(
    tissue_metrics_df[["model"]] == "MODAE_no_private_without_deconfounding", 
    "MODAE_private_embedding_no_deconfounding", 
    tissue_metrics_df[["model"]]
  )
  tissue_metrics_df[["model"]] <- ifelse(
    tissue_metrics_df[["model"]] == "MODAE_no_deconfounding", 
    "MODAE_no_private_embedding_no_deconfounding", 
    tissue_metrics_df[["model"]]
  )
  col_order <- c(
    "model", 
    "accuracy", 
    "balanced_accuracy", 
    "f1_macro", 
    colnames(tissue_metrics_df)[-(ncol(tissue_metrics_df) - 0:3)]
  )
  fn <- paste0(plot_path, "combined_tissue_metrics.html")
  print(xtable::xtable(tissue_metrics_df[,col_order]), type = "html", file = fn)
  fn <- paste0(plot_path, "combined_tissue_metrics.csv")
  readr::write_csv(tissue_metrics_df[,col_order], file = fn)
}