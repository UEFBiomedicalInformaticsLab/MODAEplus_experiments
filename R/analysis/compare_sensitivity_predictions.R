script_root <- Sys.getenv("MODAE_SCRIPT_PATH")
if (script_root == "") {
  stop("Please define the MODAE_SCRIPT_PATH environment variable.")
}
source(paste0(script_root, "R/setup.R"))

save_dirs <- c(
  "20250410_random_search/pancan_test/", 
  "20250415_random_search/pancan_test/" # primary labels only
)

pred_list <- list()
dr_perf_list <- list()

for (save_dir in save_dirs) {
  if(exists("var_list")) try(detach(var_list))
  var_list <- get_paths_and_parameters(save_dir, remove_incomplete = FALSE)
  attach(var_list)
  best_task <- readLines(paste0(save_path, "best_task.txt"))
  param_best_task_ind <- match(best_task, parameters[["task"]])
  shared_embedding_names <- paste0("z", 1:parameters[param_best_task_ind, "bottle_neck"])
  
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
    pred_list[[prediction_path]] <- tcga_final_output_table
    dr_perf_list[[prediction_path]] <- dr_perf_table
  }
}

colnames(dr_perf_list[[2]]) <- paste0(colnames(dr_perf_list[[2]]), "_primary_only")
colnames(dr_perf_list[[2]])[1] <- "drug"

combined_perf_df <- plyr::join(dr_perf_list[[1]], dr_perf_list[[2]], by = "drug", type = "full")

truncate_below <- function(x, a) ifelse(x < a, a, x)

p1 <- ggplot(
  combined_perf_df, 
  aes(
    truncate_below(r2_test_mean, 0), 
    truncate_below(r2_test_mean_primary_only, 0)
  )
) + 
  geom_point(shape = "+", size = 3) + 
  theme_bw()
fn <- paste0(plot_path, "best_model_drug_test_r2_cor_to_previous.png")
save_figure_safe(
  ggExtra::ggMarginal(p1, type = "densigram"), 
  png, 
  fn, 
  width = plot_width * 1.5, 
  height = plot_width * 1.3, 
  res = plot_res, 
  units = plot_units
)

p1 <- ggplot(
  combined_perf_df, 
  aes(
    truncate_below(r2_train_mean, 0), 
    truncate_below(r2_train_mean_primary_only, 0)
  )
) + 
  geom_point(shape = "+", size = 3) + 
  theme_bw()
fn <- paste0(plot_path, "best_model_drug_train_r2_cor_to_previous.png")
save_figure_safe(
  ggExtra::ggMarginal(p1, type = "densigram"), 
  png, 
  fn, 
  width = plot_width * 1.5, 
  height = plot_width * 1.3, 
  res = plot_res, 
  units = plot_units
)

sensitivity_cols <- grep("_sensitivity$", colnames(tcga_final_output_table))
sensitivity_cols <- colnames(tcga_final_output_table)[sensitivity_cols]
# should be identical order to names saved in file
sensitivity_drug_names <- get_ctrp_drugs(n_drugs = length(sensitivity_cols))
sensitivity_drug_names <- process_ctrp_drug_names(sensitivity_drug_names)

patient_pred_cor_df <- data.frame()
patient_sens_pred_cor_path <- paste0(plot_path, "patient_sensitivity_cor/")
dir.create(patient_sens_pred_cor_path, recursive = TRUE)
for (sensi in sensitivity_cols) {
  sens_df <- plyr::join(
    data.frame(
      id = pred_list[[1]][["id"]], 
      type = pred_list[[1]][["type"]], 
      sensitivity = pred_list[[1]][[sensi]]
    ),
    data.frame(
      id = pred_list[[2]][["id"]], 
      type = pred_list[[2]][["type"]], 
      sensitivity_primary_only = pred_list[[2]][[sensi]]
    ),
    by = "id", 
    type = "full"
  )
  n_labels <- length(unique(sens_df[["type"]]))
  p1 <- ggplot(sens_df, aes(x = sensitivity, y = sensitivity_primary_only, color = type)) + 
    geom_point(shape = "+", size = 3) + 
    theme_bw() + 
    ggtitle(sensi) + 
    scale_color_manual(values = pals::kovesi.rainbow(n=n_labels))
  fn <- paste0(patient_sens_pred_cor_path, sensi, ".png")
  save_figure_safe(
    ggExtra::ggMarginal(p1, type = "densigram"), 
    png, 
    fn, 
    width = plot_width * 1.5, 
    height = plot_width * 1.3, 
    res = plot_res, 
    units = plot_units
  )
  patient_pred_cor_df <- rbind(
    patient_pred_cor_df, 
    data.frame(
      col = sensi, 
      cor = with(sens_df, cor(sensitivity, sensitivity_primary_only, method = "pearson"))
    )
  )
}
patient_pred_cor_df[with(patient_pred_cor_df, !is.na(cor) & cor > 0.6), ]
patient_pred_cor_df[patient_pred_cor_df[["col"]] == "gemcitabine_sensitivity",] # 0.538

fnw <- paste0(patient_sens_pred_cor_path, "patient_pred_cor.csv")
readr::write_csv(patient_pred_cor_df, fnw)
#patient_pred_cor_df <- as.data.frame(readr::read_csv(fnw))

sensitivity_drug_names[sensitivity_drug_names != combined_perf_df[["drug"]]]
combined_perf_df[["drug"]][sensitivity_drug_names != combined_perf_df[["drug"]]]
combined_perf_df[["col"]] <- sensitivity_cols

combined_perf_df <- plyr::join(combined_perf_df, patient_pred_cor_df, by = "col")

p1 <- ggplot(combined_perf_df, aes(x = cor, y = truncate_below(r2_test_mean, 0))) + 
  geom_point(shape = "+", size = 3) + 
  theme_bw() + 
  ggtitle("performance vs. correlation")# + 
  #scale_color_manual(values = pals::kovesi.rainbow(n=n_labels))
ggExtra::ggMarginal(p1, type = "densigram")

p1 <- ggplot(combined_perf_df, aes(x = cor, y = truncate_below(r2_test_mean_primary_only, 0))) + 
  geom_point(shape = "+", size = 3) + 
  theme_bw() + 
  ggtitle("performance vs. correlation")# + 
#scale_color_manual(values = pals::kovesi.rainbow(n=n_labels))
ggExtra::ggMarginal(p1, type = "densigram")


with(combined_perf_df, hist(truncate_below(r2_test_mean_primary_only, 0)))
with(combined_perf_df, hist(truncate_below(r2_test_mean, 0)))

with(combined_perf_df, quantile(r2_test_mean_primary_only, probs = (1:9)/10))
with(combined_perf_df, quantile(r2_test_mean, probs = (1:9)/10))

with(
  combined_perf_df, 
  table(
    full = r2_test_mean > 0, 
    primary = r2_test_mean_primary_only > 0
  )
)

above_zero_r2_cols <- with(
  combined_perf_df, 
  col[r2_test_mean > 0 & r2_test_mean_primary_only > 0]
)

nonzero_pred_cols <- sapply(
  sensitivity_cols, 
  function(i) any(pred_list[[1]][[i]] != 0) & any(pred_list[[2]][[i]] != 0)
)
nonzero_pred_cols <- sensitivity_cols[nonzero_pred_cols]

id_match <- match(
  pred_list[[1]][["id"]], 
  pred_list[[2]][["id"]]
)
sens_diff <- (
  pred_list[[1]][id_match, intersect(above_zero_r2_cols, nonzero_pred_cols)] - 
  pred_list[[2]][, intersect(above_zero_r2_cols, nonzero_pred_cols)]
)
hist(as.matrix(sens_diff))

fn <- paste0(save_path, "external_evaluation/internal/patient_mrna.csv.gz")
gex <- readr::read_csv(fn)
id_match <- match(
  gex[[1]], 
  pred_list[[2]][["id"]]
)

sens_diff_cor <- cor(gex[id_match,-1], sens_diff)

p1 <- hist(as.matrix(sens_diff_cor))
fn <- paste0(patient_sens_pred_cor_path, "sensitivity_difference_gene_correlation_histogram.png")
save_figure_safe(
  plot(p1), 
  png, 
  fn, 
  width = plot_width * 0.5, 
  height = plot_width * 0.5, 
  res = plot_res, 
  units = plot_units
)

kegg_pws <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
tf_pws <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C3", subcollection = "TFT:TFT_LEGACY")
kegg_pw_lists <- with(kegg_pws, split(gene_symbol, f =  gs_name))
tf_pw_lists <- with(tf_pws, split(gene_symbol, f =  gs_name))
hist(sapply(kegg_pw_lists, length))
hist(sapply(tf_pw_lists, length))

kegg_enrichment_list <- list()
tft_enrichment_list <- list()
for (i in colnames(sens_diff_cor)) {
  kegg_enrichment_list[[i]] <- fgsea::fgsea(kegg_pw_lists, stats = sens_diff_cor[,i])
  tft_enrichment_list[[i]] <- fgsea::fgsea(tf_pw_lists, stats = sens_diff_cor[,i])
}






