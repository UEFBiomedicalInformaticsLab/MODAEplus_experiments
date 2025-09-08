script_root <- Sys.getenv("MODAE_SCRIPT_PATH")
if (script_root == "") {
  stop("Please define the MODAE_SCRIPT_PATH environment variable.")
}
source(paste0(script_root, "R/setup.R"))
library(dplyr)

# DrugBank
get_drug_bank_drug_targets <- function(
    drug_bank, 
    biomart_path, 
    ctrp_drugs, 
    cett_names = c("targets"), 
    debug = FALSE
) {
  db_lists <- lapply(
    cett_names, 
    function(i) drug_bank$cett[[i]]$general_information
  )
  
  db_targets <- do.call(dplyr::bind_rows, args = db_lists)
  
  db_targets <- db_targets[db_targets[["organism"]] == "Humans",]
  db_targets[["drug_id"]] <- db_targets[["parent_key"]]
  db_target_id <- tidyr::pivot_wider(
    drug_bank$cett$targets$polypeptides$external_identy, 
    id_cols = parent_key, 
    names_from = resource, 
    values_from = identifier
  )
  db_target_id[["id"]] <- db_target_id[["parent_key"]]
  db_targets <- plyr::join(
    db_targets[,c("id", "drug_id")], 
    db_target_id[,c("id", "UniProtKB")], 
    type = "left")
  
  prot_map_fn <- paste0(biomart_path, "protein_gene_map.csv")
  if (!file.exists(prot_map_fn)) {
    mart <- biomaRt::useEnsembl(
      biomart = "genes", 
      dataset = "hsapiens_gene_ensembl"
    )
    protein_map <- biomaRt::getBM(
      attributes = c("ensembl_gene_id", "hgnc_symbol", "uniprotswissprot"), 
      mart = mart
    )
    
    write.csv(protein_map, prot_map_fn)
  } else {
    protein_map <- read.csv(prot_map_fn, row.names = 1, header = TRUE)
  }
  protein_map <- unique(protein_map[,c("hgnc_symbol", "uniprotswissprot")])
  # 14k one-to-one matches
  if (debug) {
    print(
      table(
        duplicated(protein_map$hgnc_symbol), 
        duplicated(protein_map$uniprotswissprot)
      )
    )
  }
  
  # Calculate number of matches for uniprot
  uniprot_mapping_number <- table(protein_map$uniprotswissprot)
  # Most targets are one-to-one matching
  db_targets[["n_mapped_genes"]] <- uniprot_mapping_number[db_targets$UniProtKB]
  db_drug_names <- drug_bank$drugs$synonyms
  db_drug_names[["drug_id"]] <- db_drug_names[["drugbank-id"]]
  db_drug_names[["drug_name"]] <- db_drug_names[["synonym"]]
  db_targets <- plyr::join(
    db_targets, 
    db_drug_names[,c("drug_id", "drug_name")], 
    type = "left")
  
  ctrp_drugs_expanded <- gsub(" \\(.*?\\)", "", ctrp_drugs)
  ctrp_drugs_expanded <- strsplit(ctrp_drugs_expanded, split = "\\+") 
  ctrp_drugs_expanded <- Reduce(union, ctrp_drugs_expanded)
  
  db_targets_ctrp <- db_targets[
    tolower(db_targets[["drug_name"]]) %in% tolower(ctrp_drugs_expanded),]
  db_targets_ctrp <- db_targets_ctrp[!is.na(db_targets_ctrp[["UniProtKB"]]),]
  db_targets_ctrp[["uniprotswissprot"]] <- db_targets_ctrp[["UniProtKB"]]
  db_targets_ctrp <- plyr::join(
    db_targets_ctrp, 
    protein_map, 
    by = "uniprotswissprot"
  )
  drug_target_genes <- plyr::dlply(
    db_targets_ctrp, 
    "drug_name", 
    function(x) x[["hgnc_symbol"]]
  )
  
  return(drug_target_genes)
}

ctrp_drugs1 <- tolower(process_ctrp_drug_names(get_ctrp_drugs(n_drugs = 545)))
ctrp_drugs2 <- tolower(process_ctrp_drug_names(get_ctrp_drugs(n_drugs = 544)))
ctrp_drugs <- union(ctrp_drugs1, ctrp_drugs2)
drug_bank <- readRDS(paste0(drugbank_path, "parsed_DrugBank_data.rds"))
drug_target_genes <- get_drug_bank_drug_targets(
  drug_bank, 
  biomart_path, 
  ctrp_drugs
)
drug_target_genes_extended <- get_drug_bank_drug_targets(
  drug_bank, 
  biomart_path, 
  ctrp_drugs, 
  cett_names = c("targets", "carriers", "enzymes", "transporters")
)
names(drug_target_genes) <- tolower(names(drug_target_genes))
names(drug_target_genes_extended) <- tolower(names(drug_target_genes_extended))
drug_target_genes[["dinaciclib"]] <- c("CDK1", "CDK2", "CDK9")
drug_target_genes[["tanespimycin"]] <- c("HSP90AA1", "HSP90AB1")
drug_target_genes[["bi-2536"]] <- c("PLK1")

drugbank_id_map <- plyr::join(
  data.frame(synonym = tolower(ctrp_drugs)), 
  data.frame(
    synonym = tolower(drug_bank$drugs$synonyms[["synonym"]]), 
    id = drug_bank$drugs$synonyms[["drugbank-id"]]), 
  type = "left"
)

# BRCA drug-indications
indication_ind <- grepl(
  "breast cancer", 
  drug_bank$drugs$pharmacology[["indication"]], 
  ignore.case = TRUE
)
drugbank_brca_drug_id <- dplyr::pull(
  drug_bank$drugs$pharmacology[indication_ind, "drugbank_id"]
)

known_brca_drugs <- drugbank_id_map[
  drugbank_id_map[["id"]] %in% drugbank_brca_drug_id,
]
colnames(known_brca_drugs) <- c("drug", "drugbank_id")

drug_pairs <- strsplit(ctrp_drugs[grep("\\+", ctrp_drugs)], split = "\\+")
names(drug_pairs) <- ctrp_drugs[grep("\\+", ctrp_drugs)]
for (i in names(drug_pairs)) {
  u <- drug_pairs[[i]][1]
  v <- drug_pairs[[i]][2]
  a <- drug_target_genes[[u]]
  b <- drug_target_genes[[v]]
  drug_target_genes[[i]] <- union(a, b)
}

writeLines(
  known_brca_drugs[["drug"]], 
  paste0(ctrp_path, "drug_bank_brca_drugs.txt")
)

# All CTRP drug targets
drug_target_genes_json <- rjson::toJSON(drug_target_genes)
writeLines(
  drug_target_genes_json, 
  paste0(ctrp_path, "drug_bank_drug_target_genes.json")
)

# drug indications for all cancers
fn <- paste0(patient_expression_root_dir, "oncotree_level1.csv")
patient_ot_level1_labels <- read.csv(fn, row.names = 1, header = TRUE)
patient_ot_level1_labels[["level_1"]] <- tolower(
  gsub("_", " ", patient_ot_level1_labels[["level_1"]])
)
cancer_types <- unique(patient_ot_level1_labels[["type"]])
cancer_type_ot1_map <- match(cancer_types, patient_ot_level1_labels[["type"]])
cancer_type_ot1_map <- patient_ot_level1_labels[cancer_type_ot1_map, "level_1"]
names(cancer_type_ot1_map) <- cancer_types

indications <- list()
for (i in names(cancer_type_ot1_map)) {
  indication_ind <- grepl(
    "cancer", 
    drug_bank$drugs$pharmacology[["indication"]], 
    ignore.case = TRUE
  )
  indication_ind <- (
    indication_ind & 
      grepl(
        cancer_type_ot1_map[i], 
        drug_bank$drugs$pharmacology[["indication"]], 
        ignore.case = TRUE
      )
  )
  
  drugbank_cancer_drug_id <- dplyr::pull(
    drug_bank$drugs$pharmacology[indication_ind, "drugbank_id"]
  )
  
  known_cancer_drugs <- drugbank_id_map[
    drugbank_id_map[["id"]] %in% drugbank_cancer_drug_id,
  ]
  colnames(known_cancer_drugs) <- c("drug", "drugbank_id")
  if (nrow(known_cancer_drugs) > 0) {
    known_cancer_drugs[["type"]] <- i
    indications[[i]] <- known_cancer_drugs
  }
}
indications_df <- Reduce("rbind", indications)

fnw <- paste0(ctrp_path, "drug_bank_indications.csv")
write.csv(indications_df, fnw)


# Open Targets

ot_drug_mechanisms <- arrow::open_dataset(
  sources = paste0(
    opentargets_path, 
    "drug_mechanisms.snappy.parquet"
  ), 
  format = "parquet"
)
ot_drug_mechanisms <- as.data.frame(ot_drug_mechanisms)

ot_drug_info <- arrow::open_dataset(
  sources = paste0(
    opentargets_path, 
    "known_drug_part", 
    1:2, 
    ".snappy.parquet"
  ), 
  format = "parquet"
)
ot_drug_info <- as.data.frame(ot_drug_info)


test <- ot_drug_info[,c("drugId", "approvedSymbol")] |>
  group_by(drugId) |>
  summarize(targets = n_distinct(approvedSymbol)) |>
  collect()

table(test[["targets"]])

test2 <- ot_drug_info[,c("drugId", "targetId")] |>
  group_by(drugId) |>
  summarize(targets = n_distinct(targetId)) |>
  collect()

table(test2[["targets"]])

# drug synonyms
drug_id_table <- ot_drug_info[,c("drugId", "synonyms")] |>
  collect() |>
  tidyr::unnest(synonyms) |>
  distinct()
drug_id_table[["synonyms"]] <- tolower(drug_id_table[["synonyms"]])
drug_id_table <- distinct(drug_id_table)

table(ctrp_drugs %in% drug_id_table[["synonyms"]])
table(
  only_alphanumericals(ctrp_drugs) %in% 
    only_alphanumericals(drug_id_table[["synonyms"]])
)

hit_list <- ctrp_drugs %in% drug_id_table[["synonyms"]]
hit_list2 <- only_alphanumericals(ctrp_drugs) %in% 
  only_alphanumericals(drug_id_table[["synonyms"]])

ctrp_drugs[hit_list2 & ! hit_list]

test <- drug_id_table |> 
  group_by(synonyms) |>
  count()
table(test[["n"]])

syn_pointer <- match(
  only_alphanumericals(ctrp_drugs), 
  only_alphanumericals(test[["synonyms"]])
)
test_ind <- test[syn_pointer, "n"] > 1
test_ind[is.na(test_ind)] <- FALSE
duplicate_drugs <- test[syn_pointer[test_ind], "synonyms"]

duplicate_drug_ids <- drug_id_table[
  drug_id_table[["synonyms"]] %in% duplicate_drugs[["synonyms"]],
]

syn_ind <- only_alphanumericals(drug_id_table[["synonyms"]]) %in% 
  only_alphanumericals(ctrp_drugs)

drug_ids <- drug_id_table[syn_ind,]

# Targets
test <- drug_ids["drugId"] |> 
  left_join(ot_drug_info[,c("drugId", "approvedSymbol")], join_by(drugId))

target_table <- ot_drug_info[,c("drugId", "approvedSymbol")] |> 
  right_join(drug_ids["drugId"], join_by(drugId)) |> 
  distinct() |>
  collect()

target_list <- target_table |>
  group_by(drugId) |>
  group_split()
target_drug_ids <- lapply(target_list, pull, var = drugId)
target_drug_ids <- sapply(target_drug_ids, function(x) x[1])
target_list <- lapply(target_list, pull, var = approvedSymbol)
names(target_list) <- target_drug_ids

# Check duplicate drugs' targets
setequal(target_list[["CHEMBL2107358"]], target_list[["CHEMBL408194"]])
setequal(target_list[["CHEMBL84"]], target_list[["CHEMBL1607"]])
setequal(target_list[["CHEMBL1201179"]], target_list[["CHEMBL554"]])
setequal(target_list[["CHEMBL314854"]], target_list[["CHEMBL544665"]])
# all true

named_target_list <- target_list
name_ptr <- match(names(named_target_list), drug_ids[["drugId"]])
names(named_target_list) <- only_alphanumericals(drug_ids[["synonyms"]][name_ptr])

named_target_list_json <- rjson::toJSON(named_target_list)
writeLines(
  named_target_list_json, 
  paste0(ctrp_path, "open_targets_drug_target_genes.json")
)

# Comparison
fn <- paste0(ctrp_path, "drug_bank_drug_target_genes.json")
drug_target_genes_json <- readLines(fn)
drug_target_genes <- rjson::fromJSON(drug_target_genes_json)

fn <- paste0(ctrp_path, "open_targets_drug_target_genes.json")
named_target_list_json <- readLines(fn)
named_target_list <- rjson::fromJSON(named_target_list_json)

target_comp_list <- c()
target_diff1_list <- list()
target_diff2_list <- list()
target_intersect_list <- list()
for(i in names(drug_target_genes)) {
  j <- only_alphanumericals(i)
  if(j %in% names(named_target_list)) {
    target_comp_list[i] <- setequal(
      drug_target_genes[[i]], 
      named_target_list[[j]]
    )
    target_diff1_list[[i]] <- setdiff(
      drug_target_genes[[i]], 
      named_target_list[[j]]
    )
    target_diff2_list[[i]] <- setdiff(
      named_target_list[[j]], 
      drug_target_genes[[i]]
    )
    target_intersect_list[[i]] <- intersect(
      drug_target_genes[[i]], 
      named_target_list[[j]]
    )
  }
}

# PharmGKB
fn <- paste0(pharmgkb_path, "relationships/relationships.tsv")
pharmgkb_rels <- readr::read_tsv(fn)

with(pharmgkb_rels, table(Entity1_type, Entity2_type))

chem_to_gene <- pharmgkb_rels[
  pharmgkb_rels[["Entity1_type"]] == "Chemical" & 
    pharmgkb_rels[["Entity2_type"]] == "Gene" &
    pharmgkb_rels[["Association"]] == "associated",
]
gene_to_chem <- pharmgkb_rels[
  pharmgkb_rels[["Entity1_type"]] == "Gene" & 
    pharmgkb_rels[["Entity2_type"]] == "Chemical" &
    pharmgkb_rels[["Association"]] == "associated",
]

chem_to_gene_pairs <- with(
  chem_to_gene, 
  data.frame(chem = Entity1_name, gene = Entity2_name)
)
gene_to_chem_pairs <- with(
  gene_to_chem, 
  data.frame(chem = Entity2_name, gene = Entity1_name)
)

dim(chem_to_gene_pairs)
dim(gene_to_chem_pairs)
pgkb_drug_target_table <- distinct(
  bind_rows(
    chem_to_gene_pairs, 
    gene_to_chem_pairs
  )
)
dim(pgkb_drug_target_table)
pgkb_drug_target_list <- pgkb_drug_target_table |>
  group_by(chem) |>
  group_split()

pgkb_target_drug_ids <- lapply(pgkb_drug_target_list, pull, var = chem)
pgkb_target_drug_ids <- sapply(pgkb_target_drug_ids, function(x) x[1])
pgkb_drug_target_list <- lapply(pgkb_drug_target_list, pull, var = gene)
names(pgkb_drug_target_list) <- pgkb_target_drug_ids

pgkb_drug_target_list_json <- rjson::toJSON(pgkb_drug_target_list)
writeLines(
  pgkb_drug_target_list_json, 
  paste0(ctrp_path, "pharmgkb_drug_target_genes.json")
)

# PGKB omparison to open targets
target_comp_list <- c()
target_diff1_list <- list()
target_diff2_list <- list()
target_intersect_list <- list()
for(i in names(pgkb_drug_target_list)) {
  j <- only_alphanumericals(i)
  if(j %in% names(named_target_list)) {
    target_comp_list[i] <- setequal(
      pgkb_drug_target_list[[i]], 
      named_target_list[[j]]
    )
    target_diff1_list[[i]] <- setdiff(
      pgkb_drug_target_list[[i]], 
      named_target_list[[j]]
    )
    target_diff2_list[[i]] <- setdiff(
      named_target_list[[j]], 
      pgkb_drug_target_list[[i]]
    )
    target_intersect_list[[i]] <- intersect(
      pgkb_drug_target_list[[i]], 
      named_target_list[[j]]
    )
  }
}

# PGKB omparison to drug_bank
target_comp_list <- c()
target_diff1_list <- list()
target_diff2_list <- list()
target_intersect_list <- list()
for(i in names(pgkb_drug_target_list)) {
  if(i %in% names(drug_target_genes)) {
    target_comp_list[i] <- setequal(
      pgkb_drug_target_list[[i]], 
      drug_target_genes[[i]]
    )
    target_diff1_list[[i]] <- setdiff(
      pgkb_drug_target_list[[i]], 
      drug_target_genes[[i]]
    )
    target_diff2_list[[i]] <- setdiff(
      drug_target_genes[[i]], 
      pgkb_drug_target_list[[i]]
    )
    target_intersect_list[[i]] <- intersect(
      pgkb_drug_target_list[[i]], 
      drug_target_genes[[i]]
    )
  }
}

# CTD gene interactions
fn <- paste0(ctd_path, "CTD_chem_gene_ixns.csv.gz")
header <- readLines(fn, n = 28)[28]
cols <- strsplit(gsub("\\# ", "", header), split = ",")[[1]]
ctd_gene_interactions <- readr::read_csv(fn, skip = 29, col_names = cols)

ctd_chems <- unique(ctd_gene_interactions[["ChemicalName"]])

table(only_alphanumericals(ctrp_drugs) %in% 
        only_alphanumericals(tolower(ctd_chems)))

ctd_actions <- unique(ctd_gene_interactions[["InteractionActions"]])
ctd_actions_split <- strsplit(ctd_actions, split = "\\|")
ctd_actions_expanded <- Reduce(union, ctd_actions_split)
ctd_actions_pairs <- strsplit(ctd_actions_expanded, split = "\\^")
ctd_actions_expanded[
  sapply(ctd_actions_pairs, function(x) x[[2]]) == "expression"
]

ctd_sorafenib <- ctd_gene_interactions[
  ctd_gene_interactions[["ChemicalName"]] == "Sorafenib",
]

ctd_sorafenib_counts <- ctd_sorafenib |>
  group_by(GeneSymbol) |>
  count()

ctd_sorafenib_counts[ctd_sorafenib_counts[["n"]] > 9,]

ctd_sorafenib |>
  filter(GeneSymbol == "MAPK3") |>
  select(Interaction, InteractionActions, PubMedIDs)

ctd_sorafenib |>
  filter(grepl("expression", InteractionActions)) |>
  select(GeneSymbol, InteractionActions) |>
  distinct()

ctd_sorafenib |>
  filter(GeneSymbol %in% drug_target_genes[["sorafenib"]]) |>
  select(GeneSymbol, InteractionActions) |>
  distinct()

table(grepl("expression", ctd_sorafenib[["InteractionActions"]]))

table(
  grepl(
    "expression", 
    filter(
      ctd_gene_interactions, 
      ChemicalName == "Docetaxel"
    )[["InteractionActions"]]
  )
)

ctrp_ctd_ixs <- ctd_gene_interactions |>
  filter(
    only_alphanumericals(tolower(ChemicalName)) %in% 
      only_alphanumericals(ctrp_drugs)
  )

ctrp_ctd_ixs[["ChemicalName"]] <- only_alphanumericals(
  tolower(ctrp_ctd_ixs[["ChemicalName"]])
)
ctrp_ctd_ixs <- ctrp_ctd_ixs |>
  mutate(n_interactions = stringr::str_count(InteractionActions, "\\|") + 1)

ctrp_ctd_ixs_counts <- ctrp_ctd_ixs |>
  summarize(
    total_interactions = sum(n_interactions), 
    .by = c(ChemicalName, GeneSymbol)
  )

ctrp_ctd_gene_list <- ctrp_ctd_ixs_counts |>
  filter(total_interactions > 9) |>
  group_by(ChemicalName) |>
  group_split()

ctd_drug_names <- lapply(ctrp_ctd_gene_list, pull, var = "ChemicalName")
ctd_drug_names <- sapply(ctd_drug_names, function(x) x[1])
ctrp_ctd_gene_list <- lapply(ctrp_ctd_gene_list, pull, var = "GeneSymbol")
names(ctrp_ctd_gene_list) <- ctd_drug_names

ctrp_ctd_expr_gene_list <- ctrp_ctd_ixs |>
  filter(grepl("expression", InteractionActions)) |>
  select(ChemicalName, GeneSymbol) |>
  distinct() |>
  group_by(ChemicalName) |>
  group_split()
ctd_drug_names <- lapply(ctrp_ctd_expr_gene_list, pull, var = "ChemicalName")
ctd_drug_names <- sapply(ctd_drug_names, function(x) x[1])
ctrp_ctd_expr_gene_list <- lapply(
  ctrp_ctd_expr_gene_list, 
  pull, 
  var = "GeneSymbol"
)
names(ctrp_ctd_expr_gene_list) <- ctd_drug_names

names(ctrp_ctd_gene_list) <- only_alphanumericals(
  names(ctrp_ctd_gene_list)
)
names(ctrp_ctd_expr_gene_list) <- only_alphanumericals(
  names(ctrp_ctd_expr_gene_list)
)

ctrp_ctd_gene_list_json <- rjson::toJSON(ctrp_ctd_gene_list)
writeLines(
  ctrp_ctd_gene_list_json, 
  paste0(ctrp_path, "ctd_10int_drug_target_genes.json")
)

ctrp_ctd_expr_gene_list_json <- rjson::toJSON(ctrp_ctd_expr_gene_list)
writeLines(
  ctrp_ctd_expr_gene_list_json, 
  paste0(ctrp_path, "ctd_expr_drug_target_genes.json")
)

# Venn diagrams
names(drug_target_genes) <- only_alphanumericals(
  names(drug_target_genes)
)
names(drug_target_genes_extended) <- only_alphanumericals(
  names(drug_target_genes_extended)
)
names(named_target_list) <- only_alphanumericals(
  names(named_target_list)
)
names(pgkb_drug_target_list) <- only_alphanumericals(
  names(pgkb_drug_target_list)
)
names(ctrp_ctd_gene_list) <- only_alphanumericals(
  names(ctrp_ctd_gene_list)
)
names(ctrp_ctd_expr_gene_list) <- only_alphanumericals(
  names(ctrp_ctd_expr_gene_list)
)

common_drugs <- Reduce(
  intersect, 
  list(
    names(drug_target_genes), 
    #names(drug_target_genes_extended),
    names(named_target_list), 
    names(pgkb_drug_target_list), 
    names(ctrp_ctd_gene_list),
    names(ctrp_ctd_expr_gene_list)
  )
)

venn_diag_path <- paste0(ctrp_path, "target_venn_diagrams/")
dir.create(venn_diag_path, recursive = TRUE)
for (i in common_drugs) {
  dat <- list(
    DrugBank = drug_target_genes[[i]], 
    OpenTargets = named_target_list[[i]], 
    PharmGKB = pgkb_drug_target_list[[i]], 
    CTD_10int = ctrp_ctd_gene_list[[i]]#,
    #CTD_expr = ctrp_ctd_expr_gene_list[[i]]
  )
  p1 <- gplots::venn(dat, show.plot = FALSE)
  #par(mar = c(0,0,0,0)+0.1)
  save_figure_safe(
    plot(p1), 
    png, 
    paste0(venn_diag_path, i, "_target_venn.png"), 
    width = plot_width * 1, 
    height = plot_width * 1.1, 
    res = plot_res, 
    units = plot_units
  )
}




