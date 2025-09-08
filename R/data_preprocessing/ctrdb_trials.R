script_root <- Sys.getenv("MODAE_SCRIPT_PATH")
if (script_root == "") {
  stop("Please define the MODAE_SCRIPT_PATH environment variable.")
}
source(paste0(script_root, "R/setup.R"))

# Loop through datasets and mark for log-transform if log-quantiles are 
# far apart (presumed RNA-Seq).
library(dplyr)
ctrdb_gex_list <- list()
ctrdb_clininfo_list <- list()
quants <- list()
na_gene_rows_count <- c()
qlogdiff_threshold <- 5
for (study_id in names(ctrdb_datasets)) {
  study_gex_list <- list()
  study_clininfo_list <- list()
  for (ctrdb_id in ctrdb_datasets[[study_id]]) {
    fn <- paste0(ctrdb_path, ctrdb_id, "/matrix_.csv")
    study_gex_list[[ctrdb_id]] <- readr::read_csv(fn)
    colnames(study_gex_list[[ctrdb_id]])[1] <- "gene_id"
    fn <- paste0(ctrdb_path, ctrdb_id, "/cli.inf_.csv")
    study_clininfo_list[[ctrdb_id]] <- readr::read_csv(fn)
  }
  if (length(study_gex_list) > 1) {
    ctrdb_gex <- Reduce(full_join, study_gex_list)
    ctrdb_clininfo <- bind_rows(study_clininfo_list)
  } else {
    ctrdb_gex <- study_gex_list[[1]]
    ctrdb_clininfo <- study_clininfo_list[[1]]
  }
  na_gex_rows <- which(apply(is.na(as.matrix(ctrdb_gex[,-1])), 1, any))
  na_gene_rows_count[study_id] <- length(na_gex_rows)
  if (na_gene_rows_count[study_id] > 0) {
    ctrdb_gex <- ctrdb_gex[-na_gex_rows,]
  }
  quants[[study_id]] <- quantile(
    as.vector(as.matrix(ctrdb_gex[,-1])), 
    probs = c(0.01, 0.1, 0.5, 0.9, 0.99), 
  )
  qlogdiff <- log2(quants[[study_id]][5]+1) - log2(quants[[study_id]][1]+1)
  if (qlogdiff > qlogdiff_threshold) {
    # No negative numbers allowed
    ctrdb_gex[,-1] <- pmax(as.matrix(ctrdb_gex[,-1]), 0)
  }
  fnw <- paste0(ctrdb_path, study_id, "_gex.csv.gz")
  readr::write_csv(ctrdb_gex, fnw)
  fnw <- paste0(ctrdb_path, study_id, "_response.csv")
  readr::write_csv(ctrdb_clininfo, fnw)
  fnw <- paste0(ctrdb_path, study_id, "_gene_mapping.csv")
  study_genes <- ctrdb_gex[["gene_id"]]
  readr::write_csv(
    data.frame(
      old_id = study_genes, 
      new_id = limma::alias2SymbolTable(study_genes)
    ), 
    fnw
  )
}

q_mat <- bind_rows(quants)
qdiff <- q_mat[[5]] - q_mat[[1]]
qlogdiff <- log2(q_mat[[5]]+1) - log2(q_mat[[1]]+1)

log2_transform_map <- qlogdiff > qlogdiff_threshold
names(log2_transform_map) <- names(ctrdb_datasets)

fn <- paste0(ctrdb_path, "log2_transform_map.json")
writeLines(rjson::toJSON(log2_transform_map), fn)