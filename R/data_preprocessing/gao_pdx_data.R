script_root <- Sys.getenv("MODAE_SCRIPT_PATH")
if (script_root == "") {
  stop("Please define the MODAE_SCRIPT_PATH environment variable.")
}
source(paste0(script_root, "R/setup.R"))

fn <- paste0(gao_path, "41591_2015_BFnm3954_MOESM10_ESM.xlsx")

sheets <- readxl::excel_sheets(fn)

gao_data <- lapply(sheets, function(x) readxl::read_xlsx(fn, sheet = x))

names(gao_data) <- sheets

gao_data[["PCT curve metrics"]][["Treatment"]] <- gsub("\"", "", gao_data[["PCT curve metrics"]][["Treatment"]])
gao_drugs <- unique(gao_data[["PCT curve metrics"]][["Treatment"]])
gao_drugs_split <- strsplit(gao_drugs, split = " \\+ ")
gao_drugs_unique <- Reduce(union, gao_drugs_split)

ctrp_drugs <- readLines(paste0(ctrp_path, "drug_names.txt"))

gao_drugs_unique[tolower(gao_drugs_unique) %in% tolower(ctrp_drugs)]
gao_drugs[tolower(gao_drugs) %in% tolower(ctrp_drugs)]

gao_rnaseq <- as.matrix(gao_data[["RNAseq_fpkm"]][,-1])
rownames(gao_rnaseq) <- gao_data[["RNAseq_fpkm"]][[1]]

gao_drug_response <- as.data.frame(gao_data[["PCT curve metrics"]])
gao_drug_response_raw <- as.data.frame(gao_data[["PCT raw data"]])

hist(log2(gao_rnaseq+1))

hist(log10(apply(gao_rnaseq, 2, sum)))

pcs <- FactoMineR::PCA(
  t(log2(gao_rnaseq+1)), 
  scale.unit = FALSE, 
  ncp = 2, 
  graph = FALSE)

print(plot(pcs, label = "none"))

gao_rnaseq_tpm <- sweep(
  gao_rnaseq, 
  2, 
  apply(gao_rnaseq, 2, sum), 
  "/") * 1e6

pcs <- FactoMineR::PCA(
  t(log2(gao_rnaseq_tpm+1)), 
  scale.unit = FALSE, 
  ncp = 2, 
  graph = FALSE)

print(plot(pcs, label = "none"))

edger_normalize <- function(x, group = NULL, tpm = FALSE, gene_lengths = NULL) {
  out <- edgeR::DGEList(counts = x, group = group)
  out <- edgeR::calcNormFactors(out, method = "TMM")
  if (tpm) {
    if (is.null(gene_lengths)) stop("Please provide gene lengths for TPM.")
    out <- edgeR::rpkm(out, gene.length = gene_lengths[rownames(x)])
  } else {
    out <- edgeR::cpm(out)
  }
  return(out)
}

gao_rnaseq_tmm <- edger_normalize(
  gao_rnaseq)

pcs <- FactoMineR::PCA(
  t(log2(gao_rnaseq_tmm+1)), 
  scale.unit = FALSE, 
  ncp = 2, 
  graph = FALSE)

print(plot(pcs, label = "none"))

dge <- edgeR::DGEList(counts=gao_rnaseq)
dge <- edgeR::calcNormFactors(dge, method = "upperquartile")
gao_rnaseq_tpm_uq <- edgeR::cpm(dge, normalized.lib.sizes = TRUE)

pcs <- FactoMineR::PCA(
  t(log2(gao_rnaseq_tpm_uq+1)), 
  scale.unit = FALSE, 
  ncp = 2, 
  graph = FALSE)

print(plot(pcs, label = "none"))

fn <- paste0(gao_path, "exp_tpm_uq.csv.gz")
write.csv(gao_rnaseq_tpm_uq, gzfile(fn))

gene_symbols_updated <- limma::alias2SymbolTable(
  rownames(gao_rnaseq_tpm_uq), species = "Hs")
write.csv(
  data.frame(
    old_id = rownames(gao_rnaseq_tpm_uq), 
    new_id = gene_symbols_updated), 
  file = paste0(gao_path, "gao_gene_mapping.csv"))

