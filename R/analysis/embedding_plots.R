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
  best_fn <- paste0(par_path, "best_task.txt")
  if (file.exists(best_fn)) {
    best_task <- readLines(best_fn)
    param_best_task_ind <- match(best_task, parameters[["task"]])
  } else {
    best_task <- NA
    param_best_task_ind <- 1
  }
  
  bn_size <- parameters[param_best_task_ind, "bottle_neck"]
  dcl_size <- parameters[param_best_task_ind, "deconfounder_layers_per_batch"]
  shared_embedding_names <- paste0("z", 1:bn_size)
  if (dcl_size > 0) {
    pt_private_embedding_names <- paste0("z", (bn_size + 1):(bn_size + dcl_size))
    cl_private_embedding_names <- paste0("z", (bn_size + dcl_size + 1):(bn_size + 2*dcl_size))
  } else {
    pt_private_embedding_names <- c()
    cl_private_embedding_names <- c()
  }
  
  embedding_file_names <- c(
    "internal_survival_validation_final_embeddings.csv.gz", 
    "internal_survival_validation_ae_pt_embeddings.csv.gz", 
    "internal_survival_validation_bc_pt_embeddings.csv.gz", 
    "internal_drug_response_validation_final_embeddings.csv.gz", 
    "internal_drug_response_validation_ae_pt_embeddings.csv.gz", 
    "internal_drug_response_validation_bc_pt_embeddings.csv.gz", 
    "external_survival_validation_embeddings.csv.gz"
  )
  embedding_names <- c(
    "internal_patient_final_embeddings", 
    "internal_patient_ae_pt_embeddings", 
    "internal_patient_bc_pt_embeddings", 
    "internal_cl_final_embeddings", 
    "internal_cl_ae_pt_embeddings", 
    "internal_cl_bc_pt_embeddings", 
    "scanb_embeddings"
  )
  embedding_list <- list()
  for (i in 1:length(embedding_file_names)) {
    fn <- paste0(save_path, "external_evaluation/", embedding_file_names[i])
    if (file.exists(fn)) {
      ei <- read.csv(fn, header = TRUE, row.names = 1)
      embedding_list[[embedding_names[i]]] <- ei
    }
  }
  
  embedding_combined_list <- list(
    final = rbind(
      embedding_list[["internal_patient_final_embeddings"]], 
      embedding_list[["internal_cl_final_embeddings"]]), 
    ae_pt = rbind(
      embedding_list[["internal_patient_ae_pt_embeddings"]], 
      embedding_list[["internal_cl_ae_pt_embeddings"]]), 
    bc_pt = rbind(
      embedding_list[["internal_patient_bc_pt_embeddings"]], 
      embedding_list[["internal_cl_bc_pt_embeddings"]])
  )
  
  if (FALSE) {
    row_ids <- rownames(embedding_combined_list[["final"]])
    tcga_inds <- grep("^TCGA", row_ids)
    ccle_inds <- grep("^ACH", row_ids)
    label <- c()
    label[tcga_inds] <- "TCGA"
    label[ccle_inds] <- "CCLE"
    
    embedding_mats <- lapply(
      embedding_combined_list, 
      function(x) if(!is.null(x)) as.matrix(x[row_ids, shared_embedding_names]))
    
    final_silh <- cluster::silhouette(as.integer(factor(label)), dist = dist(embedding_mats[[1]]))
    ae_pt_silh <- cluster::silhouette(as.integer(factor(label)), dist = dist(embedding_mats[[2]]))
    bc_pt_silh <- cluster::silhouette(as.integer(factor(label)), dist = dist(embedding_mats[[3]]))
    
    hist(final_silh[,"sil_width"])
    hist(ae_pt_silh[,"sil_width"])
    hist(bc_pt_silh[,"sil_width"])
    
    mean((as.vector(embedding_mats[[2]]) - as.vector(embedding_mats[[3]]))^2)
    mean((as.vector(embedding_mats[[1]]) - as.vector(embedding_mats[[3]]))^2)
    mean((as.vector(embedding_mats[[1]]) - as.vector(embedding_mats[[2]]))^2)
    plot(as.vector(embedding_mats[[2]]), as.vector(embedding_mats[[3]]))
    plot(as.vector(embedding_mats[[2]][tcga_inds,]), as.vector(embedding_mats[[3]][tcga_inds,]))
    plot(as.vector(embedding_mats[[2]][ccle_inds,]), as.vector(embedding_mats[[3]][ccle_inds,]))
    
    plot(as.vector(embedding_mats[[1]][ccle_inds,]), as.vector(embedding_mats[[2]][ccle_inds,]))
    plot(as.vector(embedding_mats[[1]][tcga_inds,]), as.vector(embedding_mats[[2]][tcga_inds,]))
  }
  
  if (pan_cancer) {
    e_final <- embedding_combined_list[["final"]]
    cell_line_info_file <- 'Model_augmented.csv'
    fn <- paste0(cell_line_expression_root_dir, cell_line_info_file)
    cell_line_info <- read.csv(fn, header = TRUE, row.names = 1)
    ccle_primary_ptr <- rownames(cell_line_info)[cell_line_info[["PrimaryOrMetastasis"]] == "Primary"]
    
    fn <- paste0(patient_expression_root_dir, "oncotree_level1.csv")
    patient_ot_level1_labels <- read.csv(fn, row.names = 1, header = TRUE)
    patient_ot_level1_labels[["level_1"]] <- tolower(gsub("_", " ", patient_ot_level1_labels[["level_1"]]))
    fn <- paste0(ccle_path, "oncotree_level1.csv")
    cl_ot_level1_labels <- read.csv(fn, row.names = 1, header = TRUE)
    cl_ot_level1_labels[["level_1"]] <- tolower(gsub("_", " ", cl_ot_level1_labels[["level_1"]]))
    n_labels <- length(union(
      unique(patient_ot_level1_labels[["level_1"]]), 
      unique(cl_ot_level1_labels[["level_1"]])
    ))
    tissue_label <- rep_len(NA, nrow(e_final))
    names(tissue_label) <- rownames(e_final)
    
    ccle_ptr <- intersect(rownames(e_final), rownames(cl_ot_level1_labels))
    tissue_label[ccle_ptr] <- cl_ot_level1_labels[ccle_ptr, "level_1"]
    tissue_label[rownames(patient_ot_level1_labels)] <- patient_ot_level1_labels[["level_1"]]
    dataset <- factor(toupper(e_final[["dataset"]]), levels = c("TCGA", "CCLE"))
    nna_ind <- !is.na(tissue_label)
    
    tissue_knn_metrics <- tissue_classifier_evaluation(
      data = e_final[,shared_embedding_names], 
      dataset_label = dataset, 
      reference_dataset_class = "TCGA", 
      predicted_dataset_class = "CCLE", 
      tissue_label = tissue_label, 
      labeled_ind = nna_ind, 
      excluded_labels = c("adrenal gland", "testis"),
      scale_datasets_separately = FALSE, 
      knn_method = correlation_knn
    )
    
    fn <- paste0(plot_path, "tissue_metrics.csv")
    readr::write_csv(tissue_knn_metrics, file = fn)
    
    set.seed(0)
    tissue_plot <- tissue_visualizer(
      e_final[,shared_embedding_names], 
      labeled_ind = nna_ind, 
      color_var = tissue_label, 
      color_name = "tissue",
      shape_var = dataset, 
      shape_name = "dataset",
      reference_shape_label = "TCGA", 
      umap_neighbors = 20, 
      pre_manifold_pca = FALSE, 
      scale_datasets_separately = FALSE, 
      umap_args = list(init = "normlaplacian", min_dist = 0.5, nn_method = "annoy", metric = "correlation"),
      primary_color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=n_labels)), 
      primary_fill_scale = scale_fill_manual(values = pals::kovesi.rainbow(n=n_labels)), 
      secondary_color_scale = scale_color_manual(values = c(NA, "#000000")), 
      primary_shape_scale = scale_shape_manual(values = c(3,NA)),
      secondary_shape_scale = scale_shape_manual(values = c(NA,21)), 
      point_size = 3, 
      point_alpha = 0.5, 
      annotation_size = 6, 
      annotation_force = 10, 
      knn_accuracy = FALSE,
      knn_predicted_dataset = "CCLE", 
      knn_title_method_name = "MODAE"
    )
    saveRDS(tissue_plot, file = paste0(plot_path, "final_embedding_OT_level1.rds"))
    save_figure_safe(
      with(tissue_plot, tissue_plot), 
      png, 
      paste0(plot_path, "final_embedding_OT_level1.png"), 
      width = plot_width * 1.2, 
      height = plot_width * 0.9, 
      res = plot_res, 
      units = plot_units
    )
    
    # Primary and non-NA tissues only
    ccle_ptr <- intersect(ccle_ptr, ccle_primary_ptr)
    tissue_label <- rep_len(NA, nrow(e_final))
    names(tissue_label) <- rownames(e_final)
    tissue_label[ccle_ptr] <- cl_ot_level1_labels[ccle_ptr, "level_1"]
    tissue_label[rownames(patient_ot_level1_labels)] <- patient_ot_level1_labels[["level_1"]]
    nna_ind <- !is.na(tissue_label)
    
    set.seed(0)
    tissue_plot <- tissue_visualizer(
      e_final[,shared_embedding_names], 
      labeled_ind = nna_ind, 
      color_var = tissue_label, 
      color_name = "tissue",
      shape_var = dataset, 
      shape_name = "dataset",
      reference_shape_label = "TCGA", 
      umap_neighbors = 20, 
      pre_manifold_pca = FALSE, 
      scale_datasets_separately = FALSE, 
      umap_args = list(init = "normlaplacian", min_dist = 0.5, nn_method = "annoy", metric = "correlation"),
      primary_color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=n_labels)), 
      primary_fill_scale = scale_fill_manual(values = pals::kovesi.rainbow(n=n_labels)), 
      secondary_color_scale = scale_color_manual(values = c(NA, "#000000")), 
      primary_shape_scale = scale_shape_manual(values = c(3,NA)),
      secondary_shape_scale = scale_shape_manual(values = c(NA,21)), 
      point_size = 3, 
      point_alpha = 0.5, 
      annotation_size = 6, 
      annotation_force = 10, 
      knn_accuracy = FALSE,
      knn_predicted_dataset = "CCLE", 
      knn_title_method_name = "MODAE"
    )
    save_figure_safe(
      with(tissue_plot, tissue_plot), 
      png, 
      paste0(plot_path, "final_embedding_OT_level1_primary_cl_only.png"), 
      width = plot_width * 1.2, 
      height = plot_width * 0.9, 
      res = plot_res, 
      units = plot_units
    )
    
    # Metastatic cl only
    ccle_metas_ptr <- rownames(cell_line_info)[cell_line_info[["PrimaryOrMetastasis"]] == "Metastatic"]
    ccle_ptr <- intersect(rownames(e_final), rownames(cl_ot_level1_labels))
    ccle_ptr <- intersect(ccle_ptr, ccle_metas_ptr)
    tissue_label <- rep_len(NA, nrow(e_final))
    names(tissue_label) <- rownames(e_final)
    tissue_label[ccle_ptr] <- cl_ot_level1_labels[ccle_ptr, "level_1"]
    tissue_label[rownames(patient_ot_level1_labels)] <- patient_ot_level1_labels[["level_1"]]
    nna_ind <- !is.na(tissue_label)
    
    set.seed(0)
    tissue_plot <- tissue_visualizer(
      e_final[,shared_embedding_names], 
      labeled_ind = nna_ind, 
      color_var = tissue_label, 
      color_name = "tissue",
      shape_var = dataset, 
      shape_name = "dataset",
      reference_shape_label = "TCGA", 
      umap_neighbors = 20, 
      pre_manifold_pca = FALSE, 
      scale_datasets_separately = FALSE, 
      umap_args = list(init = "normlaplacian", min_dist = 0.5, nn_method = "annoy", metric = "correlation"),
      primary_color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=n_labels)), 
      primary_fill_scale = scale_fill_manual(values = pals::kovesi.rainbow(n=n_labels)), 
      secondary_color_scale = scale_color_manual(values = c(NA, "#000000")), 
      primary_shape_scale = scale_shape_manual(values = c(3,NA)),
      secondary_shape_scale = scale_shape_manual(values = c(NA,21)), 
      point_size = 3, 
      point_alpha = 0.5, 
      annotation_size = 6, 
      annotation_force = 10, 
      knn_accuracy = FALSE,
      knn_predicted_dataset = "CCLE", 
      knn_title_method_name = "MODAE"
    )
    save_figure_safe(
      with(tissue_plot, tissue_plot), 
      png, 
      paste0(plot_path, "final_embedding_OT_level1_metastatic_cl_only.png"), 
      width = plot_width * 1.2, 
      height = plot_width * 0.9, 
      res = plot_res, 
      units = plot_units
    )
    
    # Collection site
    
    sample_site_map <- c(
      abdomen = NA, 
      ascites = NA, 
      autonomic_ganglia = "peripheral nervous system", 
      biliary_tract = "biliary tract",
      bone = "bone", 
      bone_marrow = NA, #table(cell_line_info[cell_line_info[["SampleCollectionSite"]] == "bone_marrow","OncotreePrimaryDisease"])
      breast = "breast", 
      central_nervous_system = "brain", #table(cell_line_info[cell_line_info[["SampleCollectionSite"]] == "central_nervous_system","OncotreePrimaryDisease"])
      cervix = "cervix", 
      colon = "colon", 
      embryonal = NA, 
      endometrium = "uterus", 
      eye = "eye", 
      fibroblast = NA, 
      haematopoietic_and_lymphoid_tissue = NA, 
      kidney = "kidney", 
      large_intestine = "bowel", 
      liver = "liver", 
      lung = "lung", 
      lymph_node = "lymphoid", 
      matched_normal_tissue = NA, 
      oesophagus = "stomach", # gastrointestinal tract with bowel?
      ovary = "ovary", 
      pancreas = "pancreas", 
      pericardial_effusion = NA, 
      placenta = NA, 
      pleura = "pleura", 
      pleural_effusion = "pleura", 
      prostate = "prostate", 
      salivary_gland = NA, 
      sinonasal = NA, 
      skin = "skin", 
      small_intestine = "bowel", 
      soft_tissue = "soft tissue", 
      spleen = NA, 
      stomach = "stomach", 
      testes = "testis", 
      thyroid = "thyroid", 
      Unknown = NA, 
      upper_aerodigestive_tract = "head neck", 
      urinary_tract = "bladder", 
      uvea = "eye"
    )
    #table(tissue_label[rownames(patient_ot_level1_labels)])
    #table(cell_line_info[["SampleCollectionSite"]])
    
    ccle_metas_ptr <- rownames(cell_line_info)[cell_line_info[["PrimaryOrMetastasis"]] == "Metastatic"]
    ccle_ptr <- intersect(rownames(e_final), rownames(cl_ot_level1_labels))
    ccle_ptr <- intersect(ccle_ptr, ccle_metas_ptr)
    tissue_label <- rep_len(NA, nrow(e_final))
    names(tissue_label) <- rownames(e_final)
    tissue_label[ccle_ptr] <- sample_site_map[tolower(cell_line_info[ccle_ptr,"SampleCollectionSite"])]
    tissue_label[rownames(patient_ot_level1_labels)] <- patient_ot_level1_labels[["level_1"]]
    nna_ind <- !is.na(tissue_label)
    
    set.seed(0)
    tissue_plot <- tissue_visualizer(
      e_final[,shared_embedding_names], 
      labeled_ind = nna_ind, 
      color_var = tissue_label, 
      color_name = "tissue",
      shape_var = dataset, 
      shape_name = "dataset",
      reference_shape_label = "TCGA", 
      umap_neighbors = 20, 
      pre_manifold_pca = FALSE, 
      scale_datasets_separately = FALSE, 
      umap_args = list(init = "normlaplacian", min_dist = 0.5, nn_method = "annoy", metric = "correlation"),
      primary_color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=n_labels)), 
      primary_fill_scale = scale_fill_manual(values = pals::kovesi.rainbow(n=n_labels)), 
      secondary_color_scale = scale_color_manual(values = c(NA, "#000000")), 
      primary_shape_scale = scale_shape_manual(values = c(3,NA)),
      secondary_shape_scale = scale_shape_manual(values = c(NA,21)), 
      point_size = 3, 
      point_alpha = 0.5, 
      annotation_size = 6, 
      annotation_force = 10, 
      knn_accuracy = FALSE,
      knn_predicted_dataset = "CCLE", 
      knn_title_method_name = "MODAE"
    )
    save_figure_safe(
      with(tissue_plot, tissue_plot), 
      png, 
      paste0(plot_path, "final_embedding_OT_level1_collection_site_metastatic_cl_only.png"), 
      width = plot_width * 1.2, 
      height = plot_width * 0.9, 
      res = plot_res, 
      units = plot_units
    )
    
    # Private embeddings
    tissue_label <- rep_len(NA, nrow(e_final))
    names(tissue_label) <- rownames(e_final)
    
    ccle_ptr <- intersect(rownames(e_final), rownames(cl_ot_level1_labels))
    tissue_label[ccle_ptr] <- cl_ot_level1_labels[ccle_ptr, "level_1"]
    tissue_label[rownames(patient_ot_level1_labels)] <- patient_ot_level1_labels[["level_1"]]
    dataset <- factor(toupper(e_final[["dataset"]]), levels = c("TCGA", "CCLE"))
    nna_ind <- !is.na(tissue_label)
    metastasis_label <- rep_len("primary", nrow(e_final))
    names(metastasis_label) <- rownames(e_final)
    metastasis_label[ccle_metas_ptr[!is.na(ccle_metas_ptr)]] <- "metastasis"
    
    set.seed(0)
    meta_plot <- COPS::umap_viz(
      e_final[ccle_ptr, cl_private_embedding_names], 
      metastasis_label[ccle_ptr], 
      "type", 
      color_scale = scale_color_brewer(palette = "Set1"), 
      pre_manifold_pca = FALSE, 
      umap_args = list(init = "normlaplacian", min_dist = 0.5)
    )
    save_figure_safe(
      meta_plot, 
      png, 
      paste0(plot_path, "cl_private_embedding_metastasis_label.png"), 
      width = plot_width * 1.2, 
      height = plot_width * 0.9, 
      res = plot_res, 
      units = plot_units
    )
    
    set.seed(0)
    tissue_plot <- tissue_visualizer(
      e_final[ccle_ptr, cl_private_embedding_names], 
      labeled_ind = nna_ind[ccle_ptr], 
      color_var = tissue_label[ccle_ptr], 
      color_name = "tissue",
      shape_var = metastasis_label[ccle_ptr], 
      shape_name = "dataset",
      reference_shape_label = "TCGA", 
      umap_neighbors = 20, 
      pre_manifold_pca = FALSE, 
      scale_datasets_separately = FALSE, 
      umap_args = list(init = "normlaplacian", min_dist = 0.5, nn_method = "annoy", metric = "correlation"),
      primary_color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=n_labels)), 
      primary_fill_scale = scale_fill_manual(values = pals::kovesi.rainbow(n=n_labels)), 
      secondary_color_scale = scale_color_manual(values = c(NA, "#000000")), 
      primary_shape_scale = scale_shape_manual(values = c(3,NA)),
      secondary_shape_scale = scale_shape_manual(values = c(NA,21)), 
      point_size = 3, 
      point_alpha = 0.5, 
      annotation_size = 6, 
      annotation_force = 10
    )
    save_figure_safe(
      with(tissue_plot, tissue_plot), 
      png, 
      paste0(plot_path, "cl_private_embedding_OT_level1.png"), 
      width = plot_width * 1.2, 
      height = plot_width * 0.9, 
      res = plot_res, 
      units = plot_units
    )
  }
  
  if (FALSE) {
    embedding_plot_list <- list()
    private_embedding_plot_list <- list()
    for (i in 1:length(embedding_combined_list)) {
      if (!is.null(embedding_combined_list[[i]])) {
        n_labels <- length(table(embedding_combined_list[[i]][["dataset"]]))
        row_ids <- rownames(embedding_combined_list[[i]])
        ccle_inds <- grep("^ACH", row_ids)
        patient_inds <- row_ids[-ccle_inds]
        ccle_inds <-row_ids[ccle_inds]
        
        set.seed(0)
        umap_ploti <- COPS::umap_viz(
          embedding_combined_list[[i]][,shared_embedding_names], 
          category = toupper(embedding_combined_list[[i]][["dataset"]]), 
          category_label = "Dataset", 
          color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=n_labels)), 
          umap_neighbors = 20)
        
        embedding_plot_list[[names(embedding_combined_list)[i]]] <- umap_ploti
        
        png(paste0(
          plot_path, 
          names(embedding_combined_list)[i], 
          "_embedding_umap_sample_types.png"), 
          width = plot_width * 0.5, height = plot_width * 0.4, res = plot_res, units = plot_units)
        print(umap_ploti)
        dev.off()
        
        pca_ploti <- COPS::pca_viz(
          embedding_combined_list[[i]][,shared_embedding_names], 
          category = toupper(embedding_combined_list[[i]][["dataset"]]), 
          category_label = "Dataset", 
          color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=n_labels)))
        
        png(paste0(
          plot_path, 
          names(embedding_combined_list)[i], 
          "_embedding_pca_sample_types.png"), 
          width = plot_width * 0.5, height = plot_width * 0.4, res = plot_res, units = plot_units)
        print(pca_ploti)
        dev.off()
        
        # Cell-line private embeddings
        set.seed(0)
        umap_ploti <- COPS::umap_viz(
          embedding_combined_list[[i]][ccle_inds, cl_private_embedding_names], 
          category = toupper(embedding_combined_list[[i]][ccle_inds, "dataset"]), 
          category_label = "Dataset", 
          color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=n_labels)), 
          umap_neighbors = 20)
        
        if (FALSE) {
          cell_line_oncotree_mapping_file <- 'Model_oncotree.csv'
          cell_line_oncotree_mappings <- read.csv(
            paste0(cell_line_expression_root_dir, cell_line_oncotree_mapping_file), 
            header = TRUE, 
            row.names = 1
          )
          cl_tissue_type <- cell_line_oncotree_mappings[ccle_inds, "level_1"]
          n_labels <- length(unique(cl_tissue_type))
          set.seed(0)
          umap_ploti <- COPS::umap_viz(
            embedding_combined_list[[i]][ccle_inds, cl_private_embedding_names], 
            category = toupper(cl_tissue_type), 
            category_label = "Tissue type", 
            color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=n_labels)), 
            umap_neighbors = 20)
        }
        
        png(paste0(
          plot_path, 
          names(embedding_combined_list)[i], 
          "_embedding_umap_cl_private.png"), 
          width = plot_width * 0.5, height = plot_width * 0.4, res = plot_res, units = plot_units)
        print(umap_ploti)
        dev.off()
        
        # Patient private embeddings
        set.seed(0)
        umap_ploti <- COPS::umap_viz(
          embedding_combined_list[[i]][patient_inds, pt_private_embedding_names], 
          category = toupper(embedding_combined_list[[i]][patient_inds, "dataset"]), 
          category_label = "Dataset", 
          color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=n_labels)[2]), 
          umap_neighbors = 20)
        
        png(paste0(
          plot_path, 
          names(embedding_combined_list)[i], 
          "_embedding_umap_pt_private.png"), 
          width = plot_width * 0.5, height = plot_width * 0.4, res = plot_res, units = plot_units)
        print(umap_ploti)
        dev.off()
      }
    }
    
    final_embeds <- embedding_combined_list[["final"]]
    
    cell_line_oncotree_mapping_file <- 'Model_oncotree.csv'
    cell_line_oncotree_mappings <- read.csv(
      paste0(cell_line_expression_root_dir, cell_line_oncotree_mapping_file), 
      header = TRUE, row.names = 1)
    cell_line_brca_map <- cell_line_oncotree_mappings[["level_2"]] == "BRCA"
    names(cell_line_brca_map) <- rownames(cell_line_oncotree_mappings)
    
    final_embeds[["BRCA"]] <- ifelse(
      final_embeds[["dataset"]] %in% c("tcga", "scanb"), 
      TRUE, 
      rownames(final_embeds) %in% names(which(cell_line_brca_map)))
    
    final_embeds[["brca_subtype"]] <- NA
    final_embeds[["brca_subtype"]][
      final_embeds[["dataset"]] == "ccle" & 
        final_embeds[["BRCA"]]] <- "CCLE_BRCA"
    tcga_brca_subtypes <- TCGAbiolinks::TCGAquery_subtype("BRCA")
    tcga_brca_match <- match(
      substr(rownames(final_embeds), 1, 12), 
      tcga_brca_subtypes[["patient"]])
    if (any(!is.na(tcga_brca_match))) {
      final_embeds[["brca_subtype"]][!is.na(tcga_brca_match)] <- 
        as.data.frame(tcga_brca_subtypes)[
          tcga_brca_match[!is.na(tcga_brca_match)], 
          "BRCA_Subtype_PAM50"]
    }
    scanb_y_fn <- paste0(scanb_path, "scanb_pheno.csv.gz")
    scanb_y <- read.csv(scanb_y_fn, row.names = 1, header = TRUE)
    scanb_match <- match(rownames(final_embeds), scanb_y[["GEX.assay"]])
    if (any(!is.na(scanb_match))) {
      scanb_subtype_map <- c(
        Basal = "Basal", 
        Her2 = "Her2", 
        LumA = "LumA", 
        LumB = "LumB", 
        Normal = "Normal", 
        unclassified = NA)
      final_embeds[["brca_subtype"]][!is.na(scanb_match)] <- scanb_subtype_map[
        scanb_y[
          scanb_match[!is.na(scanb_match)], 
          "NCN.PAM50"
        ]
      ]
    }
    
    table(final_embeds[["brca_subtype"]], useNA = "always")
    
    #pals::kovesi.rainbow(7)
    type_color_scale <- scale_color_manual(
      values = c(
        Basal = "#4DA910", 
        Her2 = "#B3C120", 
        LumA = "#FCC228", 
        Normal = "#FF8410", 
        LumB = "#FD3000", 
        CCLE_BRCA = "#0034F5"))
    
    set.seed(0)
    umap_ploti <- COPS::umap_viz(
      final_embeds[final_embeds[["BRCA"]], shared_embedding_names], 
      category = final_embeds[final_embeds[["BRCA"]], "brca_subtype"], 
      category_label = "Subtype/CL", 
      color_scale = type_color_scale, 
      umap_neighbors = 20, 
      pre_manifold_pca = FALSE, 
      umap_args = list(init = "normlaplacian", min_dist = 0.5))
    
    png(paste0(
      plot_path, 
      "final_embedding_umap_brca_subtypes.png"), 
      width = plot_width * 0.5, height = plot_width * 0.4, res = plot_res, units = plot_units)
    print(umap_ploti)
    dev.off()
    
    set.seed(0)
    temp_umap <- uwot::umap(
      final_embeds[final_embeds[["BRCA"]], shared_embedding_names],
      n_neighbors = 20, 
      n_components = 2, 
      pca = NULL, 
      verbose = FALSE, 
      init = "normlaplacian", 
      min_dist = 0.5)
    temp_umap_df <- data.frame(Dim.1 = temp_umap[,1], Dim.2 = temp_umap[,2])
    temp_umap_df <- cbind(temp_umap_df, final_embeds[final_embeds[["BRCA"]], "brca_subtype"])
    colnames(temp_umap_df)[3] <- "category"
    
    temp_umap_df[["CCLEID"]] <- with(
      final_embeds[final_embeds[["BRCA"]],], 
      ifelse(
        dataset == "ccle", 
        rownames(temp_umap_df), 
        NA
      )
    )
    cell_line_model_info <- read.csv(
      paste0(cell_line_expression_root_dir, "Model_augmented.csv"), 
      header = TRUE, row.names = 1)
    temp_umap_df[["cell_line_name"]] <- cell_line_model_info[temp_umap_df[["CCLEID"]], "CellLineName"]
    
    umap_ploti <- ggplot(temp_umap_df, aes(Dim.1, Dim.2, color = category)) + 
      geom_point(shape = "+", size = 3) + 
      theme_bw() + 
      type_color_scale + 
      labs(x = "Z1", y = "Z2", color = "Subtype/CL")
    
    png(paste0(
      plot_path, 
      "final_embedding_umap_brca_subtypes_cl_names.png"), 
      width = plot_width * 2.2, height = plot_width * 2, res = plot_res, units = plot_units)
    print(umap_ploti + ggrepel::geom_text_repel(aes(label = cell_line_name)))
    dev.off()
    
    patient_expression_fn <- paste0(save_path, "external_evaluation/internal/patient_mrna.csv.gz")
    cl_expression_fn <- paste0(save_path, "external_evaluation/internal/cl_mrna.csv.gz")
    
    patient_exp <- read.csv(patient_expression_fn, row.names = 1, header = TRUE)
    cl_exp <- read.csv(cl_expression_fn, row.names = 1, header = TRUE)
    
    exp_df <- rbind(patient_exp, cl_exp)
    
    patient_dataset <- ifelse("scanb" %in% final_embeds[["dataset"]], "SCANB", "TCGA")
    dataset_label <- rep(c(patient_dataset, "CCLE"), c(nrow(patient_exp), nrow(cl_exp)))
    
    set.seed(0)
    umap_plot_gex <- COPS::umap_viz(
      scale(as.matrix(exp_df), center = TRUE, scale = TRUE), 
      category = dataset_label, 
      category_label = "Dataset", 
      color_scale = scale_color_manual(values = pals::kovesi.rainbow(n=2)), 
      umap_neighbors = 20, 
      umap_args = list(init = "normlaplacian", min_dist = 0.5))
    
    if (FALSE) {
      g_variance <- apply(exp_df, 2, var)
      top100_varg <- rank(-g_variance) <= 100
      
      brca_gex_hm <- ComplexHeatmap::Heatmap(
        scale(as.matrix(exp_df[dataset_label == "SCANB", top100_varg]), scale = FALSE), 
        name = "scanb", 
        show_row_dend = FALSE, 
        show_column_dend = FALSE, 
        show_row_names = FALSE, 
        show_column_names= FALSE
      )
      ccle_gex_hm <- ComplexHeatmap::Heatmap(
        scale(as.matrix(exp_df[dataset_label == "CCLE", top100_varg]), scale = FALSE), 
        name = "ccle", 
        show_row_dend = FALSE, 
        show_column_dend = FALSE, 
        show_row_names = FALSE, 
        show_column_names= FALSE
      )
      png(
        paste0(
          plot_path, 
          "scanb_top100_gene_heatmap.png"
        ), 
        width = plot_width / 2, 
        height = plot_width / 2, 
        unit = plot_units, 
        res = plot_res
      )
      print(brca_gex_hm)
      dev.off()
      png(
        paste0(
          plot_path, 
          "ccle_top100_gene_heatmap.png"
        ), 
        width = plot_width / 2, 
        height = plot_width / 2, 
        unit = plot_units, 
        res = plot_res
      )
      print(ccle_gex_hm)
      dev.off()
    }
    
    model_legend <- cowplot::get_legend(umap_plot_gex)
    add_theme <-  theme(
      legend.position = "none", 
      plot.caption = element_text(hjust = 0., face = "bold", size = 10))
    plot_captions <- LETTERS
    grid_plot_list <- list(umap_plot_gex + add_theme + labs(caption = plot_captions[1]))
    plot_captions <- plot_captions[-1]
    if (!is.null(embedding_plot_list[["bc_pt"]])) {
      temp <- embedding_plot_list[["bc_pt"]] + add_theme + labs(caption = plot_captions[1])
      grid_plot_list <- c(grid_plot_list, list(temp))
      plot_captions <- plot_captions[-1]
    } 
    temp <- embedding_plot_list[["final"]] + add_theme + labs(caption = plot_captions[1])
    grid_plot_list <- c(grid_plot_list, list(temp))
    plot_captions <- plot_captions[-1]
    grid_plot_list <- c(grid_plot_list, list(model_legend))
    
    combined_embedding_plot <- gridExtra::grid.arrange(
      grobs = grid_plot_list, 
      ncol = 2, 
      widths = c(4,4))
    pdf(paste0(plot_path, "combined_embedding_plot.pdf"), width = plot_width / 25.4 * 1, height = plot_width / 25.4 * 1)
    grid::grid.draw(combined_embedding_plot)
    dev.off()
    
    if (scanb || tcga_brca) {
      source_surv_dataset_map <- c(
        TCGA = "tcga", 
        SCANB = "scanb"
      )
      target_surv_dataset_map <- c(
        TCGA = "scanb", 
        SCANB = "tcga"
      )
      source_surv_dataset = source_surv_dataset_map[patient_dataset]
      target_surv_dataset = target_surv_dataset_map[patient_dataset]
    } else {
      source_surv_dataset = "tcga"
      target_surv_dataset = NULL
    }
    
    perf_plots <- generate_comparisons(
      base_dir = base_dir, 
      save_path = save_path, 
      plot_path = plot_path, 
      #model_name = model_name, 
      source_surv_dataset = source_surv_dataset, 
      target_surv_dataset = target_surv_dataset, 
      gene_preselection = FALSE, # feature selection applied before training models
      ext_survival_dir = ifelse(scanb || tcga_brca, "brca", "pancan")
    )
    
    perf_legend <- cowplot::get_legend(perf_plots[["surv_cv"]])
    #embed_legend <- cowplot::get_legend(umap_plot_gex)
    
    for (ploti in perf_plots[c("surv_cv", "drug_cv")]) {
      temp <- ploti + add_theme + labs(caption = plot_captions[1])
      grid_plot_list <- c(grid_plot_list, list(temp))
      plot_captions <- plot_captions[-1]
    }
    grid_plot_list <- c(grid_plot_list, list(perf_legend))
    
    if (length(grid_plot_list) == 7) {
      combined_embedding_plot <- gridExtra::grid.arrange(
        grobs = grid_plot_list, 
        ncol = 5, 
        nrow = 2, 
        widths = c(4,4,2,1,1), 
        heights = c(3,2),
        layout_matrix = rbind(c(1,2,3,3,3), c(5,6,4,7,7)))
    } else if (length(grid_plot_list) == 6) {
      combined_embedding_plot <- gridExtra::grid.arrange(
        grobs = grid_plot_list, 
        ncol = 3, 
        nrow = 2, 
        widths = c(4,4,2), 
        heights = c(3,2),
        layout_matrix = rbind(c(1,2,3), c(4,5,6)))
    } else {
      stop("Unexpected plot list length")
    }
    
    png(paste0(plot_path, "performance_composite_plot.png"), width = plot_width / 25.4 * 1, height = plot_width / 25.4 * 0.6, res = 300, unit = "in")
    grid::grid.draw(combined_embedding_plot)
    dev.off()
    
    png(paste0(plot_path, "png"), width = plot_width / 25.4 * 0.5, height = plot_width / 25.4 * 0.6 * 0.4, res = 300, unit = "in")
    perf_plots[["surv_ext"]]
    dev.off()
    
    png(paste0(plot_path, "survival_performance_plot.png"), width = plot_width / 25.4 * 0.5, height = plot_width / 25.4 * 0.5, res = 300, unit = "in")
    print(perf_plots$surv_cv)
    dev.off()
    
    png(paste0(plot_path, "drug_performance_plot.png"), width = plot_width / 25.4 * 0.5, height = plot_width / 25.4 * 0.5, res = 300, unit = "in")
    print(perf_plots$drug_cv)
    dev.off()
  }
}