source("../setup.R")

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
    if (FALSE) {
      # Old table for BRCA only
      best_dr_res_path <- paste0(save_path, "../plots/", model_name, "_drug_responses/")
      fn <- paste0(
        best_dr_res_path, 
        "annotated_patient_embeddings_with_predicted_aac_brca_only.csv.gz"
      )
      brca_heatmap_df <- readr::read_csv(file = fn)
    } else {
      fn <- paste0(prediction_path, "internal/unified_patient_output_table.csv.gz")
      tcga_final_output_table <- readr::read_csv(fn)
      tcga_final_output_table <- tcga_final_output_table[!is.na(tcga_final_output_table[["type"]]),]
    }
    
    cancers <- unique(tcga_final_output_table[["type"]])
    for (canceri in cancers) {
      
    }
    
    
    dir.create(paste0(best_dr_res_path, "target_clustering"), recursive = TRUE)
    
    library(randomForest)
    #library(ggdendro)
    #library(dendextend)
    
    pheno_cols <- c(
      "age_group", 
      "grade", 
      "subtype"
    )
    surv_cols <- c(
      "surv_risk_z", 
      "surv_risk", 
      "surv_risk_z_bin"
    )
    z_cols <- colnames(brca_heatmap_df)[grep("z[0-9]+", colnames(brca_heatmap_df))]
    drug_cols <- setdiff(
      colnames(brca_heatmap_df), 
      c("SCANB_id", pheno_cols, surv_cols, z_cols)
    )
    
    cluster_list <- list()
    for (drug_col in drug_cols) {
      data_cols <- c("age_group",drug_col,"subtype","surv_risk_z","grade")
      data <- brca_heatmap_df[, data_cols]
      data <- as.data.frame(data)
      data <- data[which(complete.cases(data)),]
      data$subtype <- factor(data$subtype)
      data$grade <- factor(data$grade)
      
      set.seed(0)
      rf.fit <- randomForest(
        x = data, 
        y = NULL, 
        ntree = 100, 
        proximity = TRUE, 
        oob.prox = TRUE
      )
      hclust.rf <- stats::hclust(
        as.dist(1-rf.fit$proximity), 
        method = "ward.D2"
      )
      rf.cluster = cutree(hclust.rf, k=20)
      table(rf.cluster)
      data$Cluster <- as.factor(rf.cluster)
      cluster_list[[drug_col]] <- as.factor(rf.cluster)
      
      # Function to plot either boxplot or barplot based on the type of variable
      plot_cluster_comparison <- function(variable_name) {
        if (is.numeric(data[[variable_name]])) {
          # If the variable is numeric
          ggplot(data, aes(x = Cluster, y = data[[variable_name]], fill = Cluster)) +
            geom_boxplot() +
            labs(title = paste("Boxplot of", variable_name, "by Cluster"),
                 x = "Cluster", y = variable_name) +
            theme_minimal() + 
            theme(legend.position = "none")
        } else if (is.factor(data[[variable_name]]) || is.character(data[[variable_name]])) {
          # If the variable is a factor or character
          ggplot(data, aes(x = Cluster, fill = data[[variable_name]])) +
            geom_bar(position = "dodge") +
            labs(title = paste("Barplot of", variable_name, "by Cluster"),
                 x = "Cluster", fill = variable_name) +
            theme_minimal()
        } else {
          stop("Variable is neither numeric nor factor.")
        }
      }
      
      plot_list <- list()
      for (col in setdiff(colnames(data), "Cluster")) {
        plot_list[[col]] <- plot_cluster_comparison(col)
      }
      g1 <- gridExtra::grid.arrange(grobs = plot_list)
      
      best_dr_res_path <- paste0(save_path, "../plots/", model_name, "_drug_responses/")
      fn <- paste0(
        best_dr_res_path, 
        "target_clustering/",
        drug_col, 
        "_clinical_clustering_res.png"
      )
      png(fn, width = 2*plot_width, height = 3*plot_height, res = plot_res, units = plot_units)
      grid::grid.draw(g1)
      dev.off()
    }
  }
}