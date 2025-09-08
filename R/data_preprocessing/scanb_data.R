script_root <- Sys.getenv("MODAE_SCRIPT_PATH")
if (script_root == "") {
  stop("Please define the MODAE_SCRIPT_PATH environment variable.")
}
source(paste0(script_root, "R/setup.R"))

path <- paste0(data_root, "BRCA/SCANB/")
int_path <- paste0(data_root, "scanb_preprocessed/")
path <- paste0(
  path, 
  "RNA Sequencing-Based Single Sample Predictors of Molecular Subtype ", 
  "and Risk of Recurrence for Clinical Assessment of Early-Stage Breast Cancer/"
)

unadjusted_data <- FALSE
if (unadjusted_data) {
  path <- paste0(path, "StringTie FPKM Gene Data unadjusted/")
  load(paste0(path, "SCANB.9206.mymatrix.Rdata"))
  scanb_gex_mat <- SCANB.9206.mymatrix
} else {
  path <- paste0(path, "StringTie FPKM Gene Data LibProtocol adjusted/")
  load(paste0(path, "SCANB.9206.genematrix_noNeg.Rdata"))
  scanb_gex_mat <- SCANB.9206.genematrix_noNeg
}

if (!exists("gtf_genes")) {
  gtf <- rtracklayer::import(
    paste0(
      path, 
      "../GENCODE v27 GTF SCANB Hisat StringTie pipeline/", 
      "gencode.v27.primary_assembly.annotation.gtf"
    )
  )
  gtf_genes <- unique(gtf@elementMetadata[,c("gene_id", "gene_name")])
}
gene_symbols <- gtf_genes[["gene_name"]][
  match(rownames(scanb_gex_mat), gtf_genes[["gene_id"]])
]

# Update gene symbols and sum expression
gene_symbols_updated <- limma::alias2SymbolTable(gene_symbols, species = "Hs")

save_updated_gene_expression <- FALSE # or do it later in Python code
if (save_updated_gene_expression) {
  nna_ind <- !is.na(gene_symbols_updated)
  temp_mat <- scanb_gex_mat[nna_ind,]
  nna_gene_symbols_updated <- gene_symbols_updated[nna_ind]
  ndup_ind <- !duplicated(nna_gene_symbols_updated)
  fixed_matrix <- temp_mat[ndup_ind,]
  rownames(fixed_matrix) <- nna_gene_symbols_updated[ndup_ind]
  for (i in which(!ndup_ind)) {
    a <- fixed_matrix[nna_gene_symbols_updated[i],]
    a <- a + temp_mat[i,]
    fixed_matrix[nna_gene_symbols_updated[i],] <- a
  }
  dge <- edgeR::DGEList(counts=fixed_matrix)
} else {
  dge <- edgeR::DGEList(counts=scanb_gex_mat)
}
dge <- edgeR::calcNormFactors(dge, method = "upperquartile")
scanb_normalized <- edgeR::cpm(dge, normalized.lib.sizes = TRUE)

scanb_y <- readxl::read_xlsx(
  paste0(
    path, 
    "../Supplemental Data Table/Supplementary Data Table 1 - 2023-01-13.xlsx"
  ), 
  sheet = 1
)
table(colnames(scanb_normalized) == scanb_y$GEX.assay) # all
scanb_age_col <- "Age (5-year range, e.g., 35(31-35), 40(36-40), 45(41-45) etc.)"
scanb_y[["Age_group"]] <- scanb_y[["scanb_age_col"]]

write.csv(
  data.frame(
    old_id = rownames(scanb_gex_mat), 
    new_id = gene_symbols_updated
  ), 
  file = paste0(int_path, "scanb_gene_mapping.csv")
)
write.csv(scanb_normalized, gzfile(paste0(int_path, "scanb_tpm_uq.csv.gz")))
write.csv(scanb_y, gzfile(paste0(int_path, "scanb_pheno.csv.gz")))

