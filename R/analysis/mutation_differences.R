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
  #plot(dr_r2_best[["cl_train"]], dr_r2_best[["cl_test"]])
  dr_r2_best[["rank"]] <- rank(-dr_r2_best[["cl_test"]])
  #dr_r2_best[dr_r2_best[["rank"]] < 20,]
  
  high_confidence_drugs <- dr_r2_best[dr_r2_best[["cl_test"]] > 0.1, "drug"]
  
  # Final predictions
  int_pred <- get_final_internal_predictions()
  internal_predictions_patient <- int_pred[["patient"]]
  best_patient_predictions <- internal_predictions_patient
  
  # MAF from https://doi.org/10.1016/j.cels.2018.03.002
  fn <- paste0(patient_expression_root_dir, "maf_full.csv.gz")
  if (file.exists(fn)) {
    maf <- readr::read_csv(fn)
  }
}