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
  dr_r2_best[["rank"]] <- rank(-dr_r2_best[["cl_test"]])
  
  plot1 <- ggplot(dr_r2_best_long, aes(pmax(0, r2), fill = dataset)) + 
    geom_density(alpha = 0.4, adjust = 0.5) + 
    theme_bw() + scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
  save_figure_safe(
    plot1, 
    png, 
    paste0(best_dr_res_path, "drug_r2_density_plot.png"), 
    width = plot_width, 
    height = plot_height, 
    res = plot_res, 
    units = plot_units
  )
  
  # External drug performance
  ext_dr_results <- read.csv(
    paste0(save_path, "external_evaluation/external_drug_response_validation_performance.csv"), 
    row.names = NULL, header = TRUE)
  ext_dr_res <- plyr::join(ext_dr_results, dr_r2_best, by = c("drug"), type = "left")
  
  # One table
  drug_performance <- plyr::join(dr_r2_best, ext_dr_results, type = "left")
  drug_performance <- drug_performance[,c("drug", "cl_train", "cl_test", "PearsonR", "SpearmanR", "dataset")]
  colnames(drug_performance) <- c("drug", "CTRP_train_R2", "CTRP_test_R2", "ext_PearsonCor", "ext_SpearmanCor", "dataset")
  
  write.csv(drug_performance, paste0(plot_path, "external_drug_response_performance.csv"), row.names = FALSE)
  
  temp <- drug_performance#[drug_performance[["drug"]] %in% drugs_of_interest,]
  plot1 <- ggplot(temp, aes(x = CTRP_test_R2, y = ext_SpearmanCor, label = drug, color = dataset)) + 
    geom_point(shape = "+", size = 3) + 
    ggrepel::geom_text_repel() +
    theme_bw() + 
    xlim(c(-0.5, 0.4))
  save_figure_safe(
    plot1, 
    png, 
    paste0(plot_path, "external_validation_drug_performance_text.png"), 
    width = plot_width, 
    height = 0.7 * plot_width, 
    res = plot_res, 
    units = plot_units
  )
}