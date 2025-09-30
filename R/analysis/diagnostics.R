script_root <- Sys.getenv("MODAE_SCRIPT_PATH")
if (script_root == "") {
  stop("Please define the MODAE_SCRIPT_PATH environment variable.")
}
source(paste0(script_root, "R/setup.R"))

for (save_dir in save_dirs) {
  if(exists("var_list")) try(detach(var_list))
  var_list <- get_paths_and_parameters(save_dir, remove_incomplete = FALSE)
  attach(var_list)
  
  dn_file <- paste0(save_path, "external_evaluation/diagnostics.csv.gz")
  diag_ext <- read.csv(dn_file)
  
  ## Check diagnostic plots of best results
  for (j in c("train", "valid")) {
    temp <- diag_ext[grepl(paste0("^", j, "_"), diag_ext[["stage"]]),]
    loss_transform_map <- c(
      reconstruction_loss_dataset1 = log10, 
      reconstruction_loss_dataset2 = log10, 
      confounder_alignment_norm = log10, 
      batch_cross_entropy = identity, 
      batch_dsc = log10, 
      cross_entropy = log10, 
      survival_log_likelihood = function(x) log10(-x), 
      drug_response_mse = log10, 
      regularization = log10
    )
    for (trans_col in names(loss_transform_map)[names(loss_transform_map) %in% colnames(temp)]) {
      temp[[trans_col]] <- loss_transform_map[[trans_col]](temp[[trans_col]])
    }
    
    temp_shaped <- reshape2::melt(
      temp, 
      id.vars = c("iteration", "stage"),
      variable.name = "loss_type", 
      value.name = "loss"
    )
    
    stage_map <- c(
      "Auto_Encoder", 
      "Batch_Detection", 
      "Batch_Correction", 
      "Survival_Risk", 
      "Drug_Response", 
      "Joint"
    )
    stage_cols <- c("_ae", "_bd", "_bc", "_sr", "_dr", "")
    names(stage_map) <- paste0(j, "_losses", stage_cols)
    
    temp_shaped[["Training"]] <- factor(
      stage_map[as.character(temp_shaped[["stage"]])], 
      levels = stage_map
    )
    
    loss_map <- c(
      reconstruction_loss_dataset1 = "P_MSE", 
      reconstruction_loss_dataset2 = "CL_MSE", 
      confounder_alignment_norm = "D_FNORM", 
      batch_cross_entropy = "B_SCORE", 
      batch_dsc = "B_DSC", 
      cross_entropy = "CLASS_CE", 
      survival_log_likelihood = "SR_LL", 
      drug_response_mse = "DR_MSE", 
      regularization = "REG"
    )
    
    temp_shaped[["Loss_Type"]] <- factor(
      loss_map[as.character(temp_shaped[["loss_type"]])], 
      levels = loss_map
    )
    
    plot1 <- ggplot(
      temp_shaped, 
      aes(iteration, loss)
    ) + 
      geom_line() + 
      theme_bw() + 
      xlab("Epoch") + 
      #facet_wrap(loss_type ~., scale = "free_y") + 
      ggh4x::facet_grid2(
        Loss_Type ~ Training, 
        scales = "free", 
        independent = "none", 
        space = "free_x"
      ) + 
      theme(#axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      ) + 
      guides(color = "none")
    save_figure_safe(
      plot1,# + ggtitle(paste("Model losses", i,"phase", j, "set")), 
      pdf, #png, 
      paste0(plot_path, "external_diagnostic_plot_", j, ".pdf"), 
      #width = plot_width * 2.5, 
      #height = plot_height * 2, 
      width = plot_width / 25.4 * 1.5, 
      height = plot_width / 25.4 * 1#, 
      #res = plot_res, 
      #units = plot_units
    )
  }
}