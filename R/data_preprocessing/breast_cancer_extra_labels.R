source("setup.R")

brca_mrna <- suppressMessages(
  curatedTCGAData::curatedTCGAData(
    diseaseCode = "BRCA", 
    assays = "RNASeq2GeneNorm", 
    version = "2.0.1", 
    dry.run = FALSE))

brca_coldata <- brca_mrna@colData

brca_labels <- TCGAbiolinks::TCGAquery_subtype(tumor = "BRCA")
brca_labels <- as.data.frame(brca_labels)
sample_id <- colnames(brca_mrna)[[1]]
patient_id <- substr(sample_id, 1, 12)
temp <- data.frame(sample = sample_id, patient = patient_id)
brca_labels <- plyr::join(temp, brca_labels, by = "patient", type = "left")
#colnames(brca_labels)[1] <- "sample"
#brca_labels$sample
rownames(brca_labels) <- brca_labels$sample
#brca_labels <- brca_labels[sample_id,]

brca_labels$ER.Status <- brca_coldata$ER.Status[
  match(brca_labels$patient, brca_coldata$patientID)]
ER_map <- c(1, -1)
names(ER_map) <- c("Positive", "Negative")
brca_labels$ER.Status <- factor(ER_map[brca_labels$ER.Status])

brca_labels$HER2.Final.Status <- brca_coldata$HER2.Final.Status[
  match(brca_labels$patient, brca_coldata$patientID)]
HER2_map <- c(1, -1)
names(HER2_map) <- c("Positive", "Negative")
brca_labels$HER2.Final.Status <- factor(HER2_map[
  brca_labels$HER2.Final.Status])

brca_labels$PR.Status <- brca_coldata$PR.Status[
  match(brca_labels$patient, brca_coldata$patientID)]
PR_map <- c(1, -1)
names(PR_map) <- c("Positive", "Negative")
brca_labels$PR.Status <- factor(PR_map[brca_labels$PR.Status])

write.csv(
  brca_labels, 
  paste0(
    patient_expression_root_dir, 
    "../../../tcga_extra/brca_labels_extended.csv"))

if (FALSE) {
  status_map <- c("+", "-")
  names(status_map) <- c("1", "-1")
  out <- as.data.frame(table(x$cluster))
  colnames(out) <- c("cluster", "total_count")
  temp_er <- as.data.frame(table(x$cluster, x$ER.Status, useNA = "always"))
  colnames(temp_er) <- c("cluster", "receptor_status", "count")
  temp_er <- temp_er[!is.na(temp_er$cluster),]
  temp_er$receptor_status <- factor(paste0("ER ", status_map[as.character(temp_er$receptor_status)]))
  temp_pr <- as.data.frame(table(x$cluster, x$PR.Status, useNA = "always"))
  colnames(temp_pr) <- c("cluster", "receptor_status", "count")
  temp_pr <- temp_pr[!is.na(temp_pr$cluster),]
  temp_pr$receptor_status <- factor(paste0("PR ", status_map[as.character(temp_pr$receptor_status)]))
  temp_her2 <- as.data.frame(table(x$cluster, x$HER2.Final.Status, useNA = "always"))
  colnames(temp_her2) <- c("cluster", "receptor_status", "count")
  temp_her2 <- temp_her2[!is.na(temp_her2$cluster),]
  temp_her2$receptor_status <- factor(paste0("HER2 ", status_map[as.character(temp_her2$receptor_status)]))
  out <- lapply(list(temp_er, temp_pr, temp_her2), function(x) plyr::join(out, x, by = "cluster"))
  out <- Reduce("rbind", out)
}


if (FALSE) {
  brca_exp <- MultiAssayExperiment::assay(brca_mrna)
  brca_sample_id <- colnames(brca_exp)
  brca_sample_type <- sapply(strsplit(brca_sample_id, split = "-"), function(x) x[[4]])
  brca_primary_tumor_ind <- grep("^01", brca_sample_type)
  table(brca_sample_type[brca_primary_tumor_ind])
  brca_subtype <- brca_labels[
    brca_sample_id[brca_primary_tumor_ind], 
    "BRCA_Subtype_PAM50"]
  table(brca_subtype, useNA = "always")
  
  #brca_subtype[is.na(brca_subtype)] <-  "NA"
  
  limma_model <- function(
    expr, 
    group, 
    max_cut = 3, 
    raw_counts = TRUE
  ) {
    colnames(expr) <- gsub("-", ".", colnames(expr))
    
    model <- model.matrix(~ 0 + group)
    if (raw_counts) {
      edger_obj <- edgeR::DGEList(counts = expr)
      edger_obj <- edgeR::calcNormFactors(edger_obj, method = "TMM")
      
      low_counts <- which(apply(edgeR::cpm(edger_obj), 1, max) < max_cut)
      if (length(low_counts) > 0) {
        expression_edger_filtered <- edger_obj[-low_counts,]
      } else {
        expression_edger_filtered <- edger_obj
      }
      expression_voom <- limma::voom(expression_edger_filtered, model, plot = FALSE)
      
      fit <- limma::lmFit(expression_voom, model)
    } else {
      low_expr <- which(apply(expr, 1, max) < log2(max_cut+1))
      if (length(low_expr) > 0) {
        expr_filtered <- expr[-low_expr,]
      } else {
        expr_filtered <- expr
      }
      fit <- limma::lmFit(expr_filtered, model)
    }
    
    return(fit)
  }
  
  limma_dea <- function(
    model,
    contrast
  ) {
    tmp <- limma::contrasts.fit(model, contrast)
    tmp <- limma::eBayes(tmp)
    degs <- limma::topTable(tmp, sort.by = "P", n = Inf)
    
    return(degs)
  }
  
  nna_subtype_ind <- !is.na(brca_subtype)
  limma_mod <- limma_model(
    expr = log2(brca_exp[,brca_primary_tumor_ind][,nna_subtype_ind]+1), 
    #expr = brca_exp[,brca_primary_tumor_ind][,nna_subtype_ind], 
    group = factor(brca_subtype[nna_subtype_ind]), 
    max_cut = 3, 
    raw_counts = FALSE
  )
  
  degs_list <- list()
  subtype_space <- names(table(brca_subtype))
  for (i in subtype_space) {
    contr_str <- paste0(
      "group", i, " - (", 
      paste0("group", setdiff(subtype_space, i), 
             collapse = " + "), 
      ") / ", length(subtype_space)-1)
    contr <- limma::makeContrasts(contr_str, levels = colnames(coef(limma_mod)))
    degs_list[[i]] <- limma_dea(limma_mod, contrast = contr)
    degs_list[[i]][["contrast"]] <- gsub("group", "", contr_str)
    degs_list[[i]][["gene"]] <- rownames(degs_list[[i]])
  }
  degs <- Reduce("rbind", degs_list)
  rownames(degs) <- 1:nrow(degs)
  
  write.csv(
    degs, 
    gzfile(paste0(
      patient_expression_root_dir, 
      "../../../tcga_extra/brca_subtype_degs.csv.gz")))
  table(FC = abs(degs[["logFC"]]) > log2(2), FDR = degs[["adj.P.Val"]] < 0.01)
  
  ggplot(degs, aes(logFC, -log10(adj.P.Val))) + geom_point(shape = "+") + theme_bw()
  
  deg_genes <- unique(degs[abs(degs[["logFC"]]) > log2(2) & degs[["adj.P.Val"]] < 0.01, "gene"])
  length(deg_genes)
  
  deg_only_df <- degs[abs(degs[["logFC"]]) > log2(2) & degs[["adj.P.Val"]] < 0.01, ]
  deg_splits <- split(deg_only_df, f = deg_only_df[["contrast"]])
  sapply(deg_splits, nrow)
  
  brca_nmf <- NMF::nmf(
    log2(brca_exp[deg_genes,brca_primary_tumor_ind]+1), 
    rank = 5, 
    method = "brunet", 
    seed = "nndsvd")
  
  w <- basis(brca_nmf)
  h <- coef(brca_nmf)
  
  clust <- apply(h, 2, which.max)
  
  table(clust, brca_subtype)
}