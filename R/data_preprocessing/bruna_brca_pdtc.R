source("../setup.R")

readLines(paste0(bruna_path, "ExpressionModels.txt"), n = 2)

model_drug_response <- read.table(
  paste0(bruna_path, "DrugResponsesAUCModels.txt"), 
  header = TRUE, sep = "\t")

sample_drug_response <- read.table(
  paste0(bruna_path, "DrugResponsesAUCSamples.txt"), 
  header = TRUE, sep = "\t")

model_expression <- read.table(
  paste0(bruna_path, "ExpressionModels.txt"), 
  row.names = 1, header = TRUE, sep = "\t")

sample_expression <- read.table(
  paste0(bruna_path, "ExpressionSamples.txt"), 
  row.names = 1, header = TRUE, sep = "\t")

passing_ids <- unique(sample_drug_response[["ID"]])
table(passing_ids %in% gsub("\\.", "-", colnames(sample_expression)))
passing_ids[!passing_ids %in% gsub("\\.", "-", colnames(sample_expression))]

model_genes <- rownames(model_expression)
sample_genes <- rownames(sample_expression)

sample_expression <- sample_expression[model_genes, ]

table(model_genes %in% sample_genes)

model_genes_updated <- limma::alias2SymbolTable(model_genes)

ctrp_drugs <- readLines(paste0(ctrp_path, "drug_names.txt"))
ctrp_drugs_fix <- ctrp_drugs
ctrp_drugs_fix[grep("tipifarnib", ctrp_drugs_fix)] <- 
  "tipifarnib-P2///tipifarnib-P1"
ctrp_drug_info <- read.csv(
  paste0(ctrp_path, "drug_info.csv"), 
  row.names = 1, header = TRUE)
ctrp_drugs2 <- rownames(ctrp_drug_info)[
  match(ctrp_drugs_fix, ctrp_drug_info[["cpd_name"]])]
ctrp_drugs_cid <- ctrp_drug_info[
  match(ctrp_drugs_fix, ctrp_drug_info[["cpd_name"]]), "cid"]
bruna_drugs1 <- unique(model_drug_response[["Drug"]])
bruna_drugs2 <- unique(sample_drug_response[["Drug"]])
table(tolower(bruna_drugs1) %in% tolower(ctrp_drugs2))

matching_drugs <- bruna_drugs1[tolower(bruna_drugs1) %in% tolower(ctrp_drugs2)]

if (FALSE) {
  # Check mismatched drugs
  mismatching_drugs <- bruna_drugs[!tolower(bruna_drugs1) %in% tolower(ctrp_drugs)]
  drug_bank <- readRDS(paste0(drugbank_path, "parsed_DrugBank_data.rds"))
  
  # BRCA drugs
  table(drug_bank$drugs$pharmacology[["indication"]])
  
  bruna_drug_dbid <- dplyr::pull(drug_bank$drugs$synonyms[
    match(bruna_drugs1, drug_bank$drugs$synonyms$synonym), 
    "drugbank-id"])
  bruna_drug_cid <- dplyr::pull(drug_bank$drugs$general_information[
    match(bruna_drug_dbid, drug_bank$drugs$general_information$primary_key),
    "cas_number"])
  
  bruna_drug_cid[as.numeric(gsub("-", "", bruna_drug_cid)) %in% ctrp_drugs_cid]
  
  
  db_additionals <- mismatching_drugs[
    tolower(mismatching_drugs) %in% tolower(drug_bank$drugs$synonyms$synonym)]
  db_add_ids <- lapply(
    db_additionals, 
    function(x) drug_bank$drugs$synonyms[
      drug_bank$drugs$synonyms$synonym == x, 
      "drugbank-id"])
  
  table(sapply(db_add_ids, length)) # Only one name
  table(table(drug_bank$drugs$synonyms[["drugbank-id"]])) # Most drugs have only one name
}

# Missing values
model_missing_genes <- which(apply(is.na(model_expression), 1, sum)>0)
sample_missing_genes <- which(apply(is.na(sample_expression), 1, sum)>0)

apply(as.matrix(model_expression)[model_missing_genes, ], 1, mean, na.rm = TRUE)
apply(as.matrix(sample_expression)[sample_missing_genes, ], 1, mean, na.rm = TRUE)

hist(as.matrix(model_expression)[model_missing_genes[2],])
hist(as.matrix(sample_expression)[model_missing_genes[2],])

write.csv(model_expression, gzfile(paste0(bruna_path, "pdtc_exp.csv.gz")))
write.csv(sample_expression, gzfile(paste0(bruna_path, "pdtx_exp.csv.gz")))

model_expression_imputed <- impute::impute.knn(as.matrix(model_expression), k = 5, rng.seed = 0)
sample_expression_imputed <- impute::impute.knn(as.matrix(sample_expression), k = 10, rng.seed = 1)
write.csv(model_expression_imputed$data, gzfile(paste0(bruna_path, "pdtc_exp_imputed.csv.gz")))
write.csv(sample_expression_imputed$data, gzfile(paste0(bruna_path, "pdtx_exp_imputed.csv.gz")))
# No need to process further, not used as input for DL model
#write.csv(model_drug_response, gzfile(paste0(bruna_path, "pdtc_drug_response.csv.gz")))
#write.csv(sample_drug_response, gzfile(paste0(bruna_path, "pdtx_drug_response.csv.gz")))

write.csv(data.frame(old_id = model_genes, new_id = model_genes_updated), 
          file = paste0(bruna_path, "bruna_gene_mapping.csv"))


