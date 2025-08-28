source("../setup.R")

for (save_dir in save_dirs) {
  if(exists("var_list")) try(detach(var_list))
  var_list <- get_paths_and_parameters(save_dir, remove_incomplete = FALSE)
  attach(var_list)
  #best_task <- readLines(paste0(save_path, "best_task.txt"))
  #param_best_task_ind <- match(best_task, parameters[["task"]])
  #shared_embedding_names <- paste0("z", 1:parameters[param_best_task_ind, "bottle_neck"])
  
  dr_dir <- paste0(
    base_dir, 
    "../CODE-AE-v1.0/data/tcga/"
  )
  dr_files <- dir(path = dr_dir, pattern = "\\.csv$")
  
  dr_tables <- lapply(paste0(dr_dir, dr_files), readr::read_csv)
  names(dr_tables) <- gsub("\\.csv$", "", dr_files)
  
  tcga_first_treatments <- dr_tables[["tcga_drug_first_treatment"]]
  tcga_treats <- unique(tcga_first_treatments[["pharmaceutical_therapy_drug_name"]])
  
  ctrp_drugs <- get_ctrp_drugs()
  tcga_drugs <- tcga_treats[tolower(tcga_treats) %in% ctrp_drugs]
  
  
  
  
  brca_ind <- which(tcga_first_treatments[["tcga_project"]] == "BRCA")
  tcga_brca_treats <- unique(tcga_first_treatments[brca_ind, "pharmaceutical_therapy_drug_name"][[1]])
  tcga_drugs <- tcga_treats[tolower(tcga_treats) %in% ctrp_drugs]
  tcga_drugs_brca <- tcga_brca_treats[tolower(tcga_brca_treats) %in% ctrp_drugs]
  
  brca_drugs_ind <- which(tcga_first_treatments[["pharmaceutical_therapy_drug_name"]] %in% tcga_drugs_brca)
  with(
    tcga_first_treatments[
      intersect(brca_ind, brca_drugs_ind),
    ], 
    table(pharmaceutical_therapy_drug_name)
  )
  
  tcga_all_treatments <- dr_tables[["tcga_drug_treatment"]]
  tcga_treats <- unique(tcga_all_treatments[["pharmaceutical_therapy_drug_name"]])
  brca_ind <- which(tcga_all_treatments[["tcga_project"]] == "BRCA")
  tcga_brca_treats <- unique(tcga_all_treatments[brca_ind, "pharmaceutical_therapy_drug_name"][[1]])
  tcga_drugs <- tcga_treats[tolower(tcga_treats) %in% ctrp_drugs]
  tcga_drugs_brca <- tcga_brca_treats[tolower(tcga_brca_treats) %in% ctrp_drugs]
  
  brca_drugs_ind <- which(tcga_all_treatments[["pharmaceutical_therapy_drug_name"]] %in% tcga_drugs_brca)
  with(
    tcga_all_treatments[
      intersect(brca_ind, brca_drugs_ind),
    ], 
    table(pharmaceutical_therapy_drug_name)
  )
  ggplot(
    tcga_all_treatments[intersect(brca_ind, brca_drugs_ind),], 
    aes(x = pharmaceutical_tx_started_days_to)
  ) + 
    geom_histogram() + 
    theme_bw() + 
    facet_wrap(pharmaceutical_therapy_drug_name ~ ., scales = "free")
}
