script_root <- Sys.getenv("MODAE_SCRIPT_PATH")
if (script_root == "") {
  stop("Please define the MODAE_SCRIPT_PATH environment variable.")
}
source(paste0(script_root, "R/setup.R"))
source(paste0(script_root, "R/analysis/comparison_function.R"))

for (save_dir in save_dirs) {
  if(exists("var_list")) try(detach(var_list))
  var_list <- get_paths_and_parameters(save_dir, remove_incomplete = FALSE)
  attach(var_list)
  
  modae_file <- paste0(plot_path, "final_embedding_OT_level1.rds")
  modae_plot <- readRDS(modae_file)
  
  baseline_files <- paste0(
    output_dir, 
    "baseline_results/batch_correction/pan_cancer/", 
    c(
      "Celligner_OT_level1.rds", 
      "CODE_AE_OT_level1.rds", 
      "ComBat_OT_level1.rds", 
      "Unadjusted_OT_level1.rds"
    )
  )
  baseline_plots <- lapply(baseline_files, readRDS)
  
  dummy_plot <- ggplot(
    data.frame(
      dataset = factor(c("TCGA", "CCLE"), levels = c("TCGA", "CCLE")), 
      x = 1:2, 
      y = 1:2
    ), 
    aes(x, y, shape = dataset)
  ) + 
    geom_point() +
    theme_bw() + 
    scale_shape_manual(values = c(3,21))
  
  grob_list <- lapply(
    c(list(modae_plot), baseline_plots), 
    function(x) with(x, tissue_plot)
  )
  
  add_theme <-  theme(
    plot.caption = element_text(
      hjust = 0., 
      face = "bold", 
      size = 10
    ), 
    plot.tag = element_text(
      hjust = 0., 
      face = "bold", 
      size = 10
    ), 
    plot.tag.position = "topleft", 
    axis.title = element_blank()
  )
  for (i in 1:length(grob_list)) {
    grob_list[[i]] <- grob_list[[i]] + 
      labs(tag = LETTERS[i]) + 
      add_theme
  }
  
  grob_list <- c(grob_list, list(cowplot::get_legend(dummy_plot)))
  
  plot_grid_wide <- gridExtra::grid.arrange(
    grobs = grob_list, 
    layout_matrix = matrix(c(1,1,1,1,2,3,4,5,6,6), 2, 5),
    nrow = 2,
    ncol = 5, 
    widths = c(1,1,1,1,0.2), 
    heights = c(1,1)
  )
  
  save_figure_safe(
    grid::grid.draw(plot_grid_wide), 
    png, 
    paste0(plot_path, "embedding_comparison.png"), 
    width = plot_width * 2.4, 
    height = plot_width * 1.1, 
    res = plot_res, 
    units = plot_units
  )
  
  plot_grid_tall <- gridExtra::grid.arrange(
    grobs = grob_list, 
    layout_matrix = matrix(c(1,1,2,3,1,1,4,5,6,6,NA,NA), 4, 3),
    nrow = 4,
    ncol = 3, 
    widths = c(1,1,0.2), 
    heights = c(1,1,1,1)
  )
  
  save_figure_safe(
    grid::grid.draw(plot_grid_tall), 
    png, 
    paste0(plot_path, "embedding_comparison_tall.png"), 
    width = plot_width * 1.2, 
    height = plot_width * 1.8, 
    res = plot_res, 
    units = plot_units
  )
  
  plot_grid_square <- gridExtra::grid.arrange(
    grobs = grob_list, 
    layout_matrix = matrix(c(1,1,2,1,1,3,6,5,4), 3, 3),
    nrow = 3,
    ncol = 3, 
    widths = c(1,1,1), 
    heights = c(1,1,1)
  )
  
  save_figure_safe(
    grid::grid.draw(plot_grid_square), 
    png, 
    paste0(plot_path, "embedding_comparison_square.png"), 
    width = plot_width * 1.2, 
    height = plot_width * 1.2, 
    res = plot_res, 
    units = plot_units
  )
}