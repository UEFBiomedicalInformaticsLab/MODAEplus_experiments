source("../setup.R")

codeae_data_path <- paste0(codeae_path, "data/preprocessed_dat/")
fn <- paste0(codeae_data_path, "xena_samples.txt")
tcga_samples <- readLines(fn)
fn <- paste0(codeae_data_path, "ccle_samples.txt")
ccle_samples <- readLines(fn)

codeae_pred_path <- paste0(
  codeae_path, 
  "XieResearchGroup-CODE-AE-6dc17a5/", 
  "intermediate_results/tcga_prediction/"
)
fn <- paste0(codeae_pred_path, "tcga_scores.csv")
codeae_predictions <- readr::read_csv(fn, show_col_types = FALSE)

table(codeae_predictions[["Sample"]] == tcga_samples) 

ctrp_drugs <- get_ctrp_drugs(n_drugs = 544)

codae_drugs <- only_alphanumericals(colnames(codeae_predictions)[-1:-2])

plot(codeae_predictions$Sorafenib...20, codeae_predictions$Sorafenib...51)

codae_drugs[tolower(codae_drugs) %in% ctrp_drugs]

