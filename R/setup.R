data_root <- Sys.getenv("MODAE_DATA_PATH")
if (data_root == "") {
  stop("Please define the MODAE_DATA_PATH environment variable.")
}
output_dir <- Sys.getenv("MODAE_OUTPUT_PATH")
if (output_dir == "") {
  stop("Please define the MODAE_OUTPUT_PATH environment variable.")
}

patient_expression_root_dir <- paste0(data_root, "tcga/pan_cancer/data_full/")
cell_line_expression_root_dir <- paste0(data_root, "ccle/")
cell_line_drug_response_root_dir <- paste0(data_root, "ctrp/")
bruna_path <- paste0(data_root, "breast_pdtc_bruna/")
mundi_path <- paste0(data_root, "drug_response_mundi/")
xia_path <- paste0(data_root, "drug_response_dataset/")
gao_path <- paste0(data_root, "pdtc_gao/")
ccle_path <- paste0(data_root, "ccle/")
ctrp_path <- paste0(data_root, "ctrp/")
pgx_path <- paste0(data_root, "pharmacogx/")
ctrdb_path <- paste0(data_root, "ctrdb/")
drugbank_path <- paste0(data_root, "Common_databases/DrugBank/")
opentargets_path <- paste0(data_root, "opentargets/")
pharmgkb_path <- paste0(data_root, "pharmgkb/")
ctd_path <- paste0(data_root, "ctd/CTD_chem_gene_ixns.csv.gz")
scanb_path <- paste0(data_root, "scanb_preprocessed/")
biomart_path <- paste0(data_root, "biomart_results/")
celligner_path <- paste0(data_root, "Celligner/")
codeae_path <- paste0(data_root, "CODE-AE-v1.0/")
harkonen_path <- paste0(data_root, "lung_cancer/harkonen_2023/")

ctrp_drug_target_file <- ""
cell_line_oncotree_mapping_file <- "Model_oncotree.csv"

library(ggplot2)
poster_size <- FALSE
if (poster_size) {
  plot_units <- "in"
  plot_width <- 16 # 86 or 178
  plot_height <- 11.45 #50 #86
  plot_base_size <- 34
} else {
  plot_units <- "mm"
  # Bioinformatics (OUP)
  plot_width <- 178 #86 or 178 
  plot_height <- 86 #50 #86
  # Genome Biology (BMC)
  #plot_width <- 170
  #plot_height <- 125
  plot_base_size <- 8
}
plot_res <- 300
plot_ncol <- 3

save_dirs <- c(
  "20250410_random_search/pancan_test/" 
  #"20250410_random_search/pancan_ablation_test/no_classifier/" 
  #"20250410_random_search/pancan_ablation_test/no_famo/" 
  #"20250410_random_search/pancan_ablation_test/no_deconfounding/" 
  #"20250410_random_search/pancan_ablation_test/no_private_without_deconfounding/"
)

get_paths_and_parameters <- function(savedir, remove_incomplete = FALSE) {
  out <- list()
  save_path <- paste0(output_dir, savedir)
  plot_path <- paste0(save_path, "../plots")
  if (!dir.exists(plot_path)) {
    dir.create(plot_path)
  }
  if (grepl("ablation", save_dir)) {
    #par_path <- paste0(output_dir, gsub("_ablation|/.*?/$", "", save_dir), "/")
    par_path <- paste0(output_dir, save_dir, "../")
    additional_model_string <- gsub(".*?/", "", gsub("/$", "", save_dir))
    additional_model_string <- paste0(additional_model_string, "_")
  } else {
    par_path <- paste0(output_dir, save_dir)
    additional_model_string <- ""
  }
  all_files <- dir(save_path)
  par_files <- dir(par_path)
  parameter_string <- paste0(additional_model_string, "parameters_")
  model_name <- par_files[grep(parameter_string, par_files)][1]
  model_name <- gsub("_parameters_.*", "", model_name)
  
  plot_path <- paste0(plot_path, "/", model_name, "_")
  
  parameter_jsons <- par_files[grep(paste0(parameter_string, "task[0-9|None]+\\.json"), par_files)]
  data_files <- grep("\\.tfrecords$|\\.json$", all_files)
  if (length(data_files) > 0) {
    all_files <- all_files[-data_files]
  }
  incomplete <- grep("run[0-9]+_fold[0-9]", all_files)
  if (remove_incomplete & length(incomplete) > 0) {
    incomplete_files <- all_files[incomplete]
    all_files <- all_files[-incomplete]
  } else {
    incomplete_files <- c()
  }
  
  parameter_files <- par_files[grep("parameters", all_files)]
  if (length(parameter_files) == 0) {
    parameter_files <- parameter_jsons
  }
  parameters <- lapply(paste0(par_path, parameter_files), parameter_parser)
  parameters <- Reduce(COPS::rbind_fill, parameters)
  if (!any(grepl("\\.json$", parameter_files))) {
    parameters <- unpack_parameters(parameters)
  }
  # Check if analysing pan-cancer or breast cancer
  pan_cancer <- grepl("pancan", savedir)
  tcga_brca <- grepl("brca", savedir) && !pan_cancer
  scanb <- grepl("scanb", savedir) && !tcga_brca && !pan_cancer
  out <- list(
    model_name = model_name, 
    save_path = save_path, 
    par_path = par_path, 
    plot_path = plot_path, 
    all_files = all_files, 
    incomplete_files = incomplete_files, 
    parameters = parameters, 
    pan_cancer = pan_cancer, 
    tcga_brca = tcga_brca, 
    scanb = scanb
  )
  return(out)
}

read_result_files <- function(f, path, model_name) {
  out <- read.csv(paste0(path, f), header = TRUE)
  if (nrow(out)>0) {
    out[["name"]] <- model_name
  }
  return(out)
}

outer_join <- function(a,b) plyr::join(a, b, by = intersect(colnames(a), colnames(b)), type = "full")

parameter_parser <- function(fn) {
  if (grepl("\\.json$", fn)) {
    return(parameter_parser_json(fn))
  } else {
    return(parameter_parser_csv(fn))
  }
}

# New parameters & save format
parameter_parser_json <- function(fn) {
  parameters <- jsonlite::fromJSON(fn)
  parameters_to_remove <- grep(
    paste(c(
      "serialized_data", 
      "ps_test_sets", 
      "ps_validation_sets", 
      "survival_evaluation_brier_times", 
      "nthreads_interop", 
      "nthreads", 
      "gpu_memory", 
      "parallel", 
      "omics_layer", 
      "file_name_prefix", 
      "save", 
      "data_batch_args", 
      "data_dict"), 
      collapse = "|"), 
    names(parameters))
  if (length(parameters_to_remove) > 0) {
    parameters <- parameters[-parameters_to_remove]
  }
  for (i in names(parameters)) {
    if (!inherits(parameters[[i]], "list")) {
      temp <- list(parameters[[i]])
      names(temp) <- i
      parameters[[i]] <- temp
    }
  }
  parameters <- Reduce(c, parameters)
  
  enc_layers <- parameters[["encoder_layers"]]
  for (layer_i in 1:length(enc_layers)) {
    parameters[[paste0("encoder_layer", layer_i)]] <- enc_layers[layer_i]
  }
  parameters[["bottle_neck"]] <- enc_layers[length(enc_layers)]
  dec_layers <- parameters[["decoder_layers"]]
  for (layer_i in 1:length(dec_layers)) {
    parameters[[paste0("decoder_layer", layer_i)]] <- dec_layers[layer_i]
  }
  cl_layers <- parameters[["classifier_layers"]]
  if (length(cl_layers) > 0) for (layer_i in 1:length(cl_layers)) {
    parameters[[paste0("classifier_layer", layer_i)]] <- cl_layers[layer_i]
  }
  sr_layers <- parameters[["survival_model_layers"]]
  if (length(sr_layers) > 0) for (layer_i in 1:length(sr_layers)) {
    parameters[[paste0("survival_model_layer", layer_i)]] <- sr_layers[layer_i]
  }
  dr_layers <- parameters[["drug_response_model_layers"]]
  if (length(dr_layers) > 0) for (layer_i in 1:length(dr_layers)) {
    parameters[[paste0("drug_response_model_layer", layer_i)]] <- dr_layers[layer_i]
  }
  ba_layers <- parameters[["batch_adversarial_model_layers"]]
  if (length(ba_layers) > 0) for (layer_i in 1:length(ba_layers)) {
    parameters[[paste0("batch_adversarial_model_layer", layer_i)]] <- ba_layers[layer_i]
  }
  for (i in names(parameters)) {
    if (length(parameters[[i]]) > 1) {
      if (inherits(parameters[[i]], "list")) {
        temp <- as.data.frame(parameters[[i]])
        #colnames(temp) <- paste0(i, "_", colnames(temp))
        parameters[[i]] <- temp
      } else {
        parameters[[i]] <- paste(parameters[[i]], collapse = ",")
      }
    }
    if (length(parameters[[i]]) == 0) parameters[[i]] <- NA
  }
  parameters[["task"]] <- parameters[["task_id"]]
  parameters_df <- as.data.frame(parameters)
  return(parameters_df)
}

parameter_parser_csv <- function(fn) {
  parameter_content <- readLines(fn)
  parameter_cols <- strsplit(parameter_content[1], split = ",")[[1]]
  
  parameter_cols_to_remove <- grep(
    paste(c(
      "serialized_data", 
      "ps_test_sets", 
      "ps_validation_sets", 
      "survival_evaluation_brier_times", 
      "nthreads_interop", 
      "nthreads", 
      "gpu_memory", 
      "parallel", 
      "omics_layer", 
      "file_name_prefix", 
      "save", 
      "data_batch_args", 
      "data_dict"), 
      collapse = "|"), 
    parameter_cols)
  
  parameter_body <- paste(parameter_content[-1], collapse = "")
  parameter_quote_split <- strsplit(parameter_body, split = "\"")[[1]]
  
  parameters <- c()
  j <- 1
  for (i in 1:length(parameter_quote_split)) {
    if (nchar(parameter_quote_split[[i]]) > 0) {
      if (i %% 2 == 0) {
        # Inside quotes
        parameters[j] <- parameter_quote_split[[i]]
        j <- j + 1
      } else {
        cont <- gsub("^,|,$", "", parameter_quote_split[[i]])
        if (nchar(cont) == 0) next
        vals <- strsplit(cont, split = ",")[[1]]
        #if (any(nchar(vals) > 0)) stop()
        #vals <- vals[nchar(vals) > 0]
        #if (length(vals) > 0) {
        #  parameters[j:(j+length(vals)-1)] <- vals
        #  j <- j + length(vals)
        #}
        parameters[j:(j+length(vals)-1)] <- vals
        j <- j + length(vals)
      }
    }
  }
  parameters_relevant <- parameters[-parameter_cols_to_remove]
  out <- t(parameters_relevant)
  out <- as.data.frame(out)
  colnames(out) <- parameter_cols[-parameter_cols_to_remove]
  
  gym_args <- out[["gym_args"]]
  if (!is.null(gym_args)) {
    out[["gym_args"]] <- NULL
    gym_comma_split <- strsplit(gym_args, split = ",")[[1]]
    gym_key_value_pairs <- strsplit(gym_comma_split, split = ": ")
    gym_key_value_pairs <- lapply(gym_key_value_pairs, function(x) gsub("\\{|\\'|\\}| ", "", x))
    for (i in 1:length(gym_key_value_pairs)) {
      key <- gym_key_value_pairs[[i]][1]
      val <- gym_key_value_pairs[[i]][2]
      out[[key]] <- python_map_f(val)
    }
  }
  
  for (i in colnames(out)) {
    num <- suppressWarnings(as.numeric(out[[i]]))
    if (!any(is.na(num))) out[[i]] <- num
  }
  
  return(out)
}

python_map_f <- function(v) {
  if (v == "True") return(TRUE)
  if (v == "False") return(FALSE)
  return(as.numeric(v))
}

unpack_parameters <- function(parameters) {
  #parameters[["name"]] <- gsub(".*/", "", parameters[["file_name_prefix"]])
  #parameters[["name"]] <- gsub("_$", "", parameters[["name"]])
  parameters[["task"]] <- gsub("_$", "", parameters[["task_id"]])
  enc_layers <- strsplit(gsub("\\[|\\]|\\s", "", parameters[["encoder_layers"]]), split = ",")
  for (layer_i in 1:(max(sapply(enc_layers, length)))) {
    parameters[[paste0("encoder_layer", layer_i)]] <- as.numeric(sapply(enc_layers, function(x) x[layer_i]))
  }
  parameters[["bottle_neck"]] <- as.numeric(sapply(enc_layers, function(x) x[length(x)]))
  dec_layers <- strsplit(gsub("\\[|\\]|\\s", "", parameters[["decoder_layers"]]), split = ",")
  for (layer_i in 1:(max(sapply(dec_layers, length)))) {
    parameters[[paste0("decoder_layer", layer_i)]] <- as.numeric(sapply(dec_layers, function(x) x[layer_i]))
  }
  cl_layers <- strsplit(gsub("\\[|\\]|\\s", "", parameters[["classifier_layers"]]), split = ",")
  for (layer_i in 1:(max(sapply(cl_layers, length)))) {
    parameters[[paste0("classifier_layer", layer_i)]] <- as.numeric(sapply(cl_layers, function(x) x[layer_i]))
  }
  sr_layers <- strsplit(gsub("\\[|\\]|\\s", "", parameters[["survival_model_layers"]]), split = ",")
  for (layer_i in 1:(max(sapply(sr_layers, length)))) {
    parameters[[paste0("survival_model_layer", layer_i)]] <- as.numeric(sapply(sr_layers, function(x) x[layer_i]))
  }
  dr_layers <- strsplit(gsub("\\[|\\]|\\s", "", parameters[["drug_response_model_layers"]]), split = ",")
  for (layer_i in 1:(max(sapply(dr_layers, length)))) {
    parameters[[paste0("drug_response_model_layer", layer_i)]] <- as.numeric(sapply(dr_layers, function(x) x[layer_i]))
  }
  ba_layers <- strsplit(gsub("\\[|\\]|\\s", "", parameters[["batch_adversarial_model_layers"]]), split = ",")
  for (layer_i in 1:(max(sapply(ba_layers, length)))) {
    parameters[[paste0("batch_adversarial_model_layer", layer_i)]] <- as.numeric(sapply(ba_layers, function(x) x[layer_i]))
  }
  #ow <- lapply(strsplit(parameters[["objective_weights"]], split = "\\n "), function(x) as.numeric(gsub("\\[|array|\\(|\\)|\\]", "", x)))
  #ow <- lapply(strsplit(parameters[["objective_weights"]], split = ", "), function(x) as.numeric(gsub("\\[|array|\\(|\\)|\\]", "", x)))
  ow <- lapply(strsplit(parameters[["objective_weights"]], split = "\\s+"), function(x) as.numeric(gsub("\\[|array|\\(|\\)|\\]", "", x)))
  parameters[["reconstruction_weight"]] <- sapply(ow, function(x) x[[1]])
  parameters[["classifier_weight"]] <- sapply(ow, function(x) x[[2]])
  parameters[["survival_weight"]] <- sapply(ow, function(x) x[[3]])
  parameters[["batch_weight"]] <- sapply(ow, function(x) x[[4]])
  parameters[["drug_weight"]] <- sapply(ow, function(x) x[[5]])
  
  return(parameters)
}

save_figure_safe <- function(fig, gfx_dev = png, ...) {
  gfx_dev(...)
  tryCatch(
    print(fig), 
    error = function(e) warning(paste("Encountered error while saving figure", list(...)[[1]], ":", e)),
    finally = dev.off()
  )
}

ggforest_fixed <- function(
    model, 
    data = NULL, 
    main = "Hazard ratio", 
    cpositions = c(0.02, 0.22, 0.4), 
    fontsize = 0.7, 
    refLabel = "reference", 
    noDigits = 2, 
    range_width_factor = 0.7
) {
  conf.high <- conf.low <- estimate <- NULL
  stopifnot(inherits(model, "coxph"))
  data <- survminer:::.get_data(model, data = data)
  terms <- attr(model$terms, "dataClasses")[-1]
  coef <- as.data.frame(broom::tidy(model, conf.int = TRUE))
  gmodel <- broom::glance(model)
  allTerms <- lapply(seq_along(terms), function(i) {
    var <- names(terms)[i]
    if (terms[i] %in% c("factor", "character")) {
      adf <- as.data.frame(table(data[, var]))
      cbind(var = var, adf, pos = 1:nrow(adf))
    }
    else if (terms[i] == "numeric") {
      data.frame(var = var, Var1 = "", Freq = nrow(data), pos = 1)
    }
    else {
      vars = grep(paste0("^", var, "*."), coef$term, value = TRUE)
      data.frame(var = vars, Var1 = "", Freq = nrow(data), pos = seq_along(vars))
    }
  })
  allTermsDF <- do.call(rbind, allTerms)
  colnames(allTermsDF) <- c("var", "level", "N", "pos")
  inds <- apply(allTermsDF[, 1:2], 1, paste0, collapse = "")
  rownames(coef) <- gsub(coef$term, pattern = "`", replacement = "")
  toShow <- cbind(allTermsDF, coef[inds, ])[
    , 
    c("var", 
      "level", 
      "N", 
      "p.value", 
      "estimate", 
      "conf.low", 
      "conf.high", 
      "pos")
  ]
  toShowExp <- toShow[, 5:7]
  toShowExp[is.na(toShowExp)] <- 0
  toShowExp <- format(exp(toShowExp), digits = noDigits)
  toShowExpClean <- data.frame(toShow, pvalue = signif(toShow[, 4], noDigits + 1), toShowExp)
  toShowExpClean$stars <- paste0(
    round(toShowExpClean$p.value, noDigits + 1), 
    " ", ifelse(toShowExpClean$p.value < 0.05, "*", ""), 
    ifelse(toShowExpClean$p.value < 0.01, "*", ""), 
    ifelse(toShowExpClean$p.value < 0.001, "*", ""))
  toShowExpClean$ci <- paste0(
    "(", toShowExpClean[, "conf.low.1"], " - ", toShowExpClean[, "conf.high.1"], ")")
  toShowExpClean$estimate.1[is.na(toShowExpClean$estimate)] = refLabel
  toShowExpClean$stars[which(toShowExpClean$p.value < 0.001)] = "<0.001 ***"
  toShowExpClean$stars[is.na(toShowExpClean$estimate)] = ""
  toShowExpClean$ci[is.na(toShowExpClean$estimate)] = ""
  toShowExpClean$estimate[is.na(toShowExpClean$estimate)] = 0
  toShowExpClean$var = as.character(toShowExpClean$var)
  toShowExpClean$var[duplicated(toShowExpClean$var)] = ""
  toShowExpClean$N <- paste0("(N=", toShowExpClean$N, ")")
  toShowExpClean <- toShowExpClean[nrow(toShowExpClean):1, ]
  rangeb <- range(toShowExpClean$conf.low, toShowExpClean$conf.high, na.rm = TRUE)
  breaks <- axisTicks(rangeb/2, log = TRUE, nint = 7)
  rangeplot <- rangeb
  rangeplot[1] <- rangeplot[1] - diff(rangeb) / range_width_factor
  rangeplot[2] <- rangeplot[2] + 0.15 * diff(rangeb)
  width <- diff(rangeplot)
  y_variable <- rangeplot[1] + cpositions[1] * width
  y_nlevel <- rangeplot[1] + cpositions[2] * width
  y_cistring <- rangeplot[1] + cpositions[3] * width
  y_stars <- rangeb[2]
  x_annotate <- seq_len(nrow(toShowExpClean))
  annot_size_mm <- fontsize * as.numeric(grid::convertX(unit(theme_get()$text$size, "pt"), "mm"))
  p <- ggplot(toShowExpClean, aes(seq_along(var), exp(estimate))) + 
    geom_rect(
      aes(
        xmin = seq_along(var) - 0.5, 
        xmax = seq_along(var) + 0.5, 
        ymin = exp(rangeplot[1]), 
        ymax = exp(rangeplot[2]), 
        fill = ordered(seq_along(var)%%2 + 1)
      )
    ) + 
    scale_fill_manual(values = c("#FFFFFF33", "#00000033"), guide = "none") + 
    geom_point(pch = 15, size = 4) + 
    geom_errorbar(
      aes(
        ymin = exp(conf.low), 
        ymax = exp(conf.high)
      ), width = 0.15) + 
    geom_hline(yintercept = 1, linetype = 3) + 
    coord_flip(ylim = exp(rangeplot)) + 
    ggtitle(main) + 
    scale_y_log10(
      name = "", 
      labels = sprintf("%g", breaks), 
      expand = c(0.02, 0.02), 
      breaks = breaks) + 
    theme_light() + 
    theme(
      panel.grid.minor.y = element_blank(), 
      panel.grid.minor.x = element_blank(), 
      panel.grid.major.y = element_blank(), 
      legend.position = "none", 
      panel.border = element_blank(), 
      axis.title.y = element_blank(), 
      axis.text.y = element_blank(), 
      axis.ticks.y = element_blank(), 
      plot.title = element_text(hjust = 0.5)) + 
    xlab("") + 
    annotate(
      geom = "text", 
      x = x_annotate, 
      y = exp(y_variable), 
      label = toShowExpClean$var, 
      fontface = "bold", hjust = 0, 
      size = annot_size_mm
    ) + 
    annotate(
      geom = "text", 
      x = x_annotate, 
      y = exp(y_nlevel), 
      hjust = 0, 
      label = toShowExpClean$level, 
      vjust = -0.1, 
      size = annot_size_mm
    ) + 
    annotate(
      geom = "text", 
      x = x_annotate, 
      y = exp(y_nlevel), 
      label = toShowExpClean$N, 
      fontface = "italic", 
      hjust = 0, 
      vjust = ifelse(toShowExpClean$level == "", 0.5, 1.1), 
      size = annot_size_mm
    ) + 
    annotate(
      geom = "text", 
      x = x_annotate, 
      y = exp(y_cistring), 
      label = toShowExpClean$estimate.1, 
      size = annot_size_mm, 
      vjust = ifelse(toShowExpClean$estimate.1 == "reference", 0.5, -0.1)
    ) + 
    annotate(
      geom = "text", 
      x = x_annotate, 
      y = exp(y_cistring), 
      label = toShowExpClean$ci, 
      size = annot_size_mm, 
      vjust = 1.1, 
      fontface = "italic"
    ) + 
    annotate(
      geom = "text", 
      x = x_annotate, 
      y = exp(y_stars), 
      label = toShowExpClean$stars, 
      size = annot_size_mm, 
      hjust = -0.2, 
      fontface = "italic"
    ) + 
    annotate(
      geom = "text", 
      x = 0.5, 
      y = exp(y_variable), 
      label = paste0(
        "# Events: ", 
        gmodel$nevent, 
        "; Global p-value (Log-Rank): ", 
        format.pval(gmodel$p.value.log, eps = ".001"), 
        " \nAIC: ", 
        round(gmodel$AIC, 2), 
        "; Concordance Index: ", 
        round(gmodel$concordance, 2)
      ), 
      size = annot_size_mm, 
      hjust = 0, 
      vjust = 1.2, 
      fontface = "italic")
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  ggpubr::as_ggplot(gt)
}

get_ctrp_drugs <- function(n_drugs = 545) {
  if (n_drugs == 545) {
    # PharmacoGx
    ctrp_drugs <- readLines(paste0(ctrp_path, "drug_names.txt"))
    ctrp_drugs_fix <- ctrp_drugs
    ctrp_drugs_fix[grep("tipifarnib", ctrp_drugs_fix)] <- "tipifarnib-P2///tipifarnib-P1"
  } else if (n_drugs == 544) {
    # Xia
    ctrp_drugs_fix <- readLines(paste0(xia_path, "processed_drug_names.txt"))
    ctrp_drugs_fix <- tolower(ctrp_drugs_fix)
    ctrp_drugs_fix[ctrp_drugs_fix == "azd-2281"] <- "olaparib"
    ctrp_drugs_fix[ctrp_drugs_fix == "5-fu"] <- "fluorouracil"
    ctrp_drugs_fix[ctrp_drugs_fix == "plx-4032"] <- "vemurafenib"
    ctrp_drugs_fix[ctrp_drugs_fix == "platin"] <- "cisplatin"
    ctrp_drugs_fix[ctrp_drugs_fix == "kx2-391"] <- "tirbanibulin"
    #ctrp_drugs_fix[ctrp_drugs_fix == "gsk1120212"] <- "trametinib"
    
  } else {
    stop("Uknown number of drugs")
  }
  
  return(ctrp_drugs_fix)
}

process_ctrp_drug_names <- function(x) {
  out <- gsub("\\(.*?\\)$", "", x)
  out <- gsub("\\s", "", out)
  out <- gsub("\\:", "+", out)
  return(out)
}

fix_ctrp_drug_names <- function(x, n_drugs = 545) {
   if (n_drugs == 544) {
     # Xia
     ctrp_drug_cols <- readLines(paste0(xia_path, "processed_drug_columns.txt"))
     ctrp_drug_names <- readLines(paste0(xia_path, "processed_drug_names.txt"))
     names(ctrp_drug_names) <- ctrp_drug_cols
     ctrp_drugs_fix <- ctrp_drug_names[as.character(x)]
     ctrp_drugs_fix <- tolower(ctrp_drugs_fix)
     ctrp_drugs_fix[ctrp_drugs_fix == "azd-2281"] <- "olaparib"
     ctrp_drugs_fix[ctrp_drugs_fix == "5-fu"] <- "fluorouracil"
     ctrp_drugs_fix[ctrp_drugs_fix == "plx-4032"] <- "vemurafenib"
     ctrp_drugs_fix[ctrp_drugs_fix == "platin"] <- "cisplatin"
     ctrp_drugs_fix[ctrp_drugs_fix == "kx2-391"] <- "tirbanibulin"
     #ctrp_drugs_fix[ctrp_drugs_fix == "gsk1120212"] <- "trametinib"
     x <- ctrp_drugs_fix
   }
  # Otherwise do nothing
  return(x)
}

get_ctrp_dr_df <- function() {
  ctrp_dr_data <- read.csv(
    paste0(ctrp_path, "screen_aac.csv.gz"), 
    row.names = NULL, header = TRUE)
  ctrp_dr_rowinfo <- read.csv(
    paste0(ctrp_path, "screen_rowinfo.csv.gz"), 
    row.names = NULL, header = TRUE)
  ctrp_cl_id_map <- read.csv(
    paste0(ctrp_path, "ctrp_ccle_clid_map.csv.gz"), 
    row.names = 1, header = TRUE)
  
  ctrp_cl_ccle_id <- ctrp_cl_id_map[ctrp_dr_rowinfo[["sampleid"]], "CCLE_model_id"]
  
  cell_line_oncotree_mapping_file <- 'Model_oncotree.csv'
  cell_line_oncotree_mappings <- read.csv(paste0(cell_line_expression_root_dir, cell_line_oncotree_mapping_file), header = TRUE, row.names = 1)
  ctrp_cl_cancer_type <- cell_line_oncotree_mappings[ctrp_cl_ccle_id, "level_2"]
  
  ctrp_dr_df <- data.frame(
    cl_id = ctrp_cl_ccle_id, 
    cl_type = ctrp_cl_cancer_type, 
    drug = ctrp_dr_rowinfo[["treatmentid_fixed"]],
    aac = ctrp_dr_data[["aac_recomputed"]])
  
  ctrp_dr_df <- plyr::ddply(ctrp_dr_df, c("drug"), function(x) {x[["aac_drugwise_z_score"]] <- scale(x[["aac"]]);return(x)})
  return(ctrp_dr_df)
}

get_final_internal_predictions <- function() {
  fn_patient <- paste0(save_path, "external_evaluation/internal_survival_validation_predictions.csv.gz")
  if (!file.exists(fn_patient)) {
    fn_patient <- paste0(save_path, "external_evaluation/internal_survival_validation_final_predictions.csv.gz")
  }
  fn_cl <- paste0(save_path, "external_evaluation/internal_drug_response_validation_predictions.csv.gz")
  if (!file.exists(fn_cl)) {
    fn_cl <- paste0(save_path, "external_evaluation/internal_drug_response_validation_final_predictions.csv.gz")
  }
  internal_predictions_patient <- read.csv(fn_patient, row.names = NULL, header = TRUE)
  internal_predictions_cl <- read.csv(fn_cl, row.names = NULL, header = TRUE)
  
  return(list(patient = internal_predictions_patient, cl = internal_predictions_cl))
}

tissue_visualizer <- function(
    data, 
    labeled_ind = c(TRUE), 
    color_var = NULL, 
    color_name = "color_name",
    shape_var = NULL, 
    shape_name = "shape_name",
    reference_shape_label = NULL, 
    umap_neighbors = 20, 
    pre_manifold_pca = TRUE, 
    max_pcs = 50, 
    scale_datasets_separately = FALSE, 
    umap_args = list(init = "normlaplacian"), 
    primary_color_scale = scale_color_brewer(palette = "Dark2"), 
    primary_fill_scale = scale_fill_brewer(palette = "Dark2"), 
    secondary_color_scale = scale_color_manual(values = c(NA, "#000000")), 
    primary_shape_scale = scale_shape_manual(values = c(3,NA)),
    secondary_shape_scale = scale_shape_manual(values = c(NA,21)), 
    point_size = 3, 
    point_alpha = 0.5, 
    annotation_size = 5, 
    annotation_force = 10, 
    annotation_max_overlaps = 10, 
    remove_legend = TRUE, 
    knn_accuracy = FALSE,
    knn_predicted_dataset = NA, 
    knn_title_method_name = ""
) {
  if (scale_datasets_separately) {
    data_split <- split(data, f = shape_var)
    data_split <- lapply(data_split, scale, center = TRUE, scale = TRUE)
    data <- Reduce(rbind, data_split)[rownames(data),]
  }
  umap_args[["X"]] <- data
  umap_args[["n_components"]] <- 2
  umap_args[["n_neighbors"]] <- umap_neighbors
  umap_args[["pca"]] <- if(pre_manifold_pca) min(max_pcs, dim(data)) else NULL
  umap_args[["verbose"]] <- FALSE
  res_umap <- do.call(uwot::umap, umap_args)
  res_umap <- data.frame(Dim.1 = res_umap[,1], Dim.2 = res_umap[,2])
  res_umap[, color_name] <- color_var
  res_umap[, shape_name] <- shape_var
  res_umap[!labeled_ind, color_name] <- NA
  res_umap[!labeled_ind, shape_name] <- NA
  
  p1 <- ggplot(res_umap, aes(Dim.1, Dim.2)) + 
    geom_point(
      aes(
        color = !!ggplot2::ensym(color_name), 
        shape = !!ggplot2::ensym(shape_name)
      ), 
      size = point_size, 
      alpha = point_alpha
    ) + 
    theme_bw() + 
    primary_color_scale + 
    primary_shape_scale +
    labs(x = "UMAP 1", y = "UMAP 2")
  
  p1 <- p1 + 
    ggnewscale::new_scale("color") + 
    ggnewscale::new_scale("shape") + 
    geom_point(
    data = res_umap, 
    aes(
      color = !!ggplot2::ensym(shape_name), 
      fill = !!ggplot2::ensym(color_name), 
      shape = !!ggplot2::ensym(shape_name)
    ), 
    size = point_size, 
    alpha = point_alpha
  ) + 
    primary_fill_scale +
    secondary_color_scale + 
    secondary_shape_scale
  
  p1_legend <- cowplot::get_legend(p1)
  if (remove_legend) {
    p1 <- p1 + theme(legend.position = "none")
  }
  
  # Annotations
  annotations <- plyr::ddply(
    if (is.null(reference_shape_label)) res_umap else res_umap[res_umap[[shape_name]] == reference_shape_label,], 
    .variables = color_name, 
    .fun = function(x) data.frame(
      Dim.1 = mean(x[["Dim.1"]]), 
      Dim.2 = mean(x[["Dim.2"]]), 
      label = x[[color_name]][1]
    )
  )
  
  p1 <- p1 + ggrepel::geom_text_repel(
    aes(Dim.1, Dim.2, label = label), 
    data = annotations, 
    size = annotation_size, 
    force = annotation_force, 
    max.overlaps = annotation_max_overlaps
  )
  
  if (knn_accuracy) {
    knn_preds <- correlation_knn(
      X = data[dataset == reference_shape_label & labeled_ind,], 
      Y = data[dataset == knn_predicted_dataset & labeled_ind,], 
      x_labels = tissue_label[dataset == reference_shape_label & labeled_ind], 
      k = 25
    )
    reference_labels <- unique(tissue_label[dataset == reference_shape_label & labeled_ind])
    predicted_labels <- unique(tissue_label[dataset == knn_predicted_dataset & labeled_ind])
    
    knn_bacc <- balanced_accuracy_score(
      pred_labels = knn_preds, 
      true_labels = tissue_label[dataset == knn_predicted_dataset & labeled_ind],
      available_labels = intersect(reference_labels, predicted_labels)
    )
    knn_acc <- accuracy_score(
      pred_labels = knn_preds, 
      true_labels = tissue_label[dataset == knn_predicted_dataset & labeled_ind],
      available_labels = intersect(reference_labels, predicted_labels)
    )
    p1 <- p1 + ggtitle(paste(
      knn_title_method_name, 
      knn_predicted_dataset, 
      "knn BACC:", 
      signif(knn_bacc, 3),
      "knn ACC:", 
      signif(knn_acc, 3)
    ))
  } else {
    knn_bacc <- NULL
    knn_acc <- NULL
  }
  
  out <- rlang::new_environment(
    data = list(
      tissue_plot = p1, 
      plot_legend = p1_legend, 
      res_umap = res_umap, 
      knn_bacc = knn_bacc, 
      knn_acc = knn_acc
    )
  )
  
  return(out)
}
tyrosine_kinase_inhibitors <- c(
  "lapatinib", 
  "sunitinib", 
  "imatinib", 
  "sorafenib", 
  "pazopanib", 
  "adavosertib"
)
pan_cancer_types <- c(
  "ACC",
  "BLCA",
  "BRCA",
  "CESC",
  "CHOL",
  "COAD",
  "DLBC",
  "ESCA",
  "GBM",
  "HNSC",
  "KICH",
  "KIRC",
  "KIRP",
  "LGG",
  "LIHC",
  "LUAD",
  "LUSC",
  "MESO",
  "OV",
  "PAAD",
  "PCPG",
  "PRAD",
  "READ",
  "SARC",
  "SKCM",
  "STAD",
  "TGCT",
  "THCA",
  "THYM",
  "UCEC",
  "UCS",
  "UVM"
)

only_alphanumericals <- function(x) gsub("[^[:alnum:]]", "", x)

ctrdb_tcga_types_map <- list(
  GSE20194 = "BRCA", 
  GSE20271 = "BRCA", 
  GSE14764 = "OVCA", 
  `E-MTAB-7083` = "OVCA", 
  GSE14209 = "STAD", 
  GSE196434 = "SKCM"
)

ctrdb_types_to_oncotree <- c(
  `Colorectal cancer` = "BOWEL",
  `Brain cancer` = "BRAIN",
  `Leukemia` = NA,
  `Kidney cancer` = "KIDNEY",
  `Melanoma` = "SKIN",
  `Lymphoma` = "LYMPH",
  `Ovarian cancer` = "OVARY",
  `Bone marrow cancer` = "BONE",
  `Breast cancer` = "BREAST",
  `Head and neck cancer` = "HEAD_NECK",
  `Esophageal cancer` = "STOMACH",
  `liver cancer` = "LIVER",
  `Urinary bladder cancer` = "BLADDER",
  `Lung cancer` = "LUNG",
  `Stomach cancer` = "STOMACH",
  `Gastrointestinal carcinoma` = "STOMACH",
  `Pancreatic cancer` = "PANCREAS",
  `Neuroendocrine carcinoma` = NA,
  `Oral cavity cancer` = "HEAD_NECK",
  `Transitional cell carcinoma` = NA,
  `Thymus cancer` = "THYMUS",
  `Gastroesophageal cancer` = "STOMACH",
  `Liver cancer` = "LIVER",
  `Anus cancer` = "BOWEL",
  `Tonsil cancer` = "HEAD_NECK",
  `Carcinoma` = NA,
  `Vulva cancer` = "VULVA",
  `Uterine cancer` = "UTERUS",
  `Muscle cancer` = "SOFT_TISSUE",
  `Soft tissue tumours` = "SOFT_TISSUE",
  `Sarcoma` = "SOFT_TISSUE",
  `Bone cancer` = "BONE",
  `Skin cancer` = "SKIN"
)


fn <- paste0(ctrdb_path, "non_tcga_datasets.json")
if (file.exists(fn)) {
  ctrdb_datasets <- rjson::fromJSON(readLines(fn))
} else {
  ctrdb_ids <- dir(ctrdb_path, pattern = "CTR_")
  ctrdb_response_list <- list()
  for (i in ctrdb_ids) {
    fn <- paste0(ctrdb_path, "/", i, "/cli.inf_.csv")
    ctrdb_response_list[[i]] <- readr::read_csv(fn)
    ctrdb_response_list[[i]][["CTRDB_id"]] <- i
    ctrdb_response_list[[i]] <- dplyr::distinct(
      ctrdb_response_list[[i]][,c("Source", "CTRDB_id", "Cancer_type_level1")]
    )
  }
  ctrdb_response_all <- dplyr::bind_rows(ctrdb_response_list)
  ctrdb_datasets <- plyr::dlply(
    ctrdb_response_all, 
    "Source", 
    function(x) unique(x[["CTRDB_id"]])
  )
  ctrdb_datasets <- ctrdb_datasets[!grepl("^TCGA", names(ctrdb_datasets))]
  
  ctrdb_datasets_json <- rjson::toJSON(ctrdb_datasets)
  fnw <- paste0(ctrdb_path, "non_tcga_datasets.json")
  writeLines(ctrdb_datasets_json, fnw)
}

if (FALSE) {
  # Manual selection
  ctrdb_datasets <- list(
    #BRCA
    GSE20194 = c(
      "CTR_Microarray_48-I",
      "CTR_Microarray_8-I",
      "CTR_Microarray_63-I",
      "CTR_Microarray_10-I",
      "CTR_Microarray_82-I",
      "CTR_Microarray_39-I",
      "CTR_Microarray_100-I"
    ), 
    GSE20271 = c(
      "CTR_Microarray_106-I", 
      "CTR_Microarray_89-I", 
      "CTR_Microarray_80-I", 
      "CTR_Microarray_45-I", 
      "CTR_Microarray_43-I", 
      "CTR_Microarray_3-I", 
      "CTR_Microarray_23-I", 
      "CTR_Microarray_24-I", 
      "CTR_Microarray_99-I", 
      "CTR_Microarray_68-I", 
      "CTR_Microarray_73-I", 
      "CTR_Microarray_69-I"
    ), 
    # OVCA
    GSE14764 = c(
      "CTR_Microarray_34-I", 
      "CTR_Microarray_12-I", 
      "CTR_Microarray_2-I", 
      "CTR_Microarray_33-I", 
      "CTR_Microarray_42-I", 
      "CTR_Microarray_87-I", 
      "CTR_Microarray_66-I"
    ), 
    `E-MTAB-7083` = c(
      "CTR_Microarray_109-I"
    ),
    # STAD (Gastric)
    GSE14209 = c(
      "CTR_Microarray_230-II"
    ),
    # SKCM (Cutaneous Melanoma)
    GSE196434 = c(
      "CTR_RNAseq_568-I",
      "CTR_RNAseq_563-I"
    )
  )
}

simplify_cancer_stages <- function(x) {
  x <- gsub("[A-C]*[0-9]*$", "", x)
  # Exclude Stage X
  x <- ifelse(x == "Stage X", NA, x)
  x <- ifelse(grepl("^Stage", x), x, NA)
  return(x)
}

tissue_classifier_evaluation <- function(
    data, 
    dataset_label, 
    reference_dataset_class, 
    predicted_dataset_class, 
    tissue_label, 
    labeled_ind = rep_len(TRUE, length(dataset_label)), 
    excluded_labels = c(), 
    scale_datasets_separately = FALSE, 
    apply_pca = FALSE, 
    npc = 100, 
    knn_method = correlation_knn
) {
  if (scale_datasets_separately) {
    data_split <- split(data, f = dataset_label)
    data_split <- lapply(data_split, scale, center = TRUE, scale = TRUE)
    data <- Reduce(rbind, data_split)[rownames(data),]
  }
  if (apply_pca) {
    pca_res <- FactoMineR::PCA(data, scale.unit = FALSE, ncp = npc, graph = FALSE)
    data <- pca_res$ind$coord
  }
  knn_preds <- knn_method(
    data[dataset_label == reference_dataset_class & labeled_ind,], 
    data[dataset_label == predicted_dataset_class & labeled_ind,], 
    tissue_label[dataset_label == reference_dataset_class & labeled_ind], 
    k = 25
  )
  reference_labels <- unique(tissue_label[dataset_label == reference_dataset_class & labeled_ind])
  predicted_labels <- unique(tissue_label[dataset_label == predicted_dataset_class & labeled_ind])
  available_labels <- intersect(reference_labels, predicted_labels)
  available_labels <- setdiff(
    intersect(reference_labels, predicted_labels), 
    excluded_labels
  )
  knn_cm <- wide_confusion_matrix(
    pred_labels = knn_preds, 
    true_labels = tissue_label[dataset_label == predicted_dataset_class & labeled_ind],
    available_labels = available_labels
  )
  knn_metrics <- metrics_from_cm(knn_cm)
  out <- data.frame(
    t(knn_metrics[,"f1"]), 
    t(knn_metrics[,"rec"]), 
    t(knn_metrics[,"pre"])
  )
  colnames(out) <- paste0(
    rep(rownames(knn_metrics), 3), 
    rep(c("_f1", "_recall", "_precision"), each = nrow(knn_metrics))
  )
  out[["accuracy"]] <- mean(knn_preds == tissue_label[dataset_label == predicted_dataset_class & labeled_ind])
  out[["balanced_accuracy"]] <- mean(knn_metrics[,"rec"])
  out[["f1_macro"]] <- mean(knn_metrics[,"f1"])
  return(out)
}

correlation_knn <- function(X, Y, x_labels, k = 25) {
  XY_cor <- cor(t(X), t(Y), method = "pearson")
  XY_cor_rank <- apply(-XY_cor, 2, rank)
  Y_knn_in_X <- which(XY_cor_rank <= k, arr.ind = TRUE)
  knn_x_ptrs <- split(Y_knn_in_X[,1], f = Y_knn_in_X[,2])
  knn_x_labels <- lapply(knn_x_ptrs, function(ptr) x_labels[ptr])
  y_pred <- sapply(knn_x_labels, function(x) names(which.max(table(x))))
  return(y_pred)
}

euclidean_knn <- function(X, Y, x_labels, k = 25) {
  XY_d2 <- matrix(NA, nrow = nrow(X), ncol = nrow(Y))
  for (j in 1:nrow(Y)) {
    d_j <- sweep(X, 2, Y[j,], "-")
    d_j <- apply(d_j^2, 1, sum)
    XY_d2[,j] <- d_j
  }
  XY_d2_rank <- apply(XY_d2, 2, rank)
  Y_knn_in_X <- which(XY_d2_rank <= k, arr.ind = TRUE)
  knn_x_ptrs <- split(Y_knn_in_X[,1], f = Y_knn_in_X[,2])
  knn_x_labels <- lapply(knn_x_ptrs, function(ptr) x_labels[ptr])
  y_pred <- sapply(knn_x_labels, function(x) names(which.max(table(x))))
  return(y_pred)
}

wide_confusion_matrix <- function(
    pred_labels, 
    true_labels, 
    available_labels = unique(true_labels)
) {
  lcm <- array(
    NA, 
    dim = c(length(available_labels), 4), 
    dimnames = list(label = available_labels, c("tp", "fn", "tn", "fp"))
  )
  for (labi in available_labels) {
    lcm[labi, "tp"] <- sum((pred_labels == labi) * (true_labels == labi))
    lcm[labi, "fn"] <- sum((pred_labels != labi) * (true_labels == labi))
    lcm[labi, "tn"] <- sum((pred_labels != labi) * (true_labels != labi))
    lcm[labi, "fp"] <- sum((pred_labels == labi) * (true_labels != labi))
  }
  return(lcm)
}

metrics_from_cm <- function(cm) {
  f1 <- 2*cm[,"tp"] / (2*cm[,"tp"] + cm[,"fp"] + cm[,"fn"])
  rec <- cm[,"tp"] / (cm[,"tp"] + cm[,"fn"])
  pre <- cm[,"tp"] / (cm[,"tp"] + cm[,"fp"])
  f1 <- 2 * (pre * rec) / (pre + rec)
  out <- cbind(f1, rec, pre)
  rownames(out) <- rownames(cm)
  colnames(out) <- c("f1", "rec", "pre")
  return(out)
}

balanced_accuracy_score <- function(pred_labels, true_labels, available_labels = unique(true_labels)) {
  recall_scores <- c()
  for (labi in available_labels) {
    tp <- sum((pred_labels == labi) * (true_labels == labi))
    fn <- sum((pred_labels != labi) * (true_labels == labi))
    recall_scores[labi] <- tp / (tp + fn)
  }
  out <- mean(recall_scores)
  attributes(out)[["recall_scores"]] <- recall_scores
  return(out)
}

accuracy_score <- function(pred_labels, true_labels, available_labels = unique(true_labels)) {
  available_ind <- true_labels %in% available_labels
  return(mean(pred_labels[available_ind] == true_labels[available_ind]))
}
