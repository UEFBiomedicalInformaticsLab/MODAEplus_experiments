script_root <- Sys.getenv("MODAE_SCRIPT_PATH")
if (script_root == "") {
  stop("Please define the MODAE_SCRIPT_PATH environment variable.")
}
source(paste0(script_root, "R/setup.R"))

library(dplyr)
# Get drug ids from target table
ctrp_drugs1 <- tolower(process_ctrp_drug_names(get_ctrp_drugs(n_drugs = 545)))
ctrp_drugs2 <- tolower(process_ctrp_drug_names(get_ctrp_drugs(n_drugs = 544)))
ctrp_drugs <- union(ctrp_drugs1, ctrp_drugs2)

ot_drug_info <- arrow::open_dataset(
  sources = paste0(opentargets_path, "known_drug_part", 1:2, ".snappy.parquet"), 
  format = "parquet"
)
ot_drug_info <- as.data.frame(ot_drug_info)

# drug synonyms = names
drug_id_table <- ot_drug_info[,c("drugId", "synonyms")] |>
  collect() |>
  tidyr::unnest(synonyms) |>
  distinct()
drug_id_table[["synonyms"]] <- tolower(drug_id_table[["synonyms"]])
drug_id_table <- distinct(drug_id_table)

syn_ind <- only_alphanumericals(drug_id_table[["synonyms"]]) %in% 
  only_alphanumericals(ctrp_drugs)

drug_ids <- drug_id_table[syn_ind,] 

# Indications
ot_drug_indications <- arrow::open_dataset(
  sources = paste0(opentargets_path, "drug_indications.snappy.parquet"), 
  format = "parquet"
)
ot_drug_indications <- as.data.frame(ot_drug_indications)

ot_ctrp_drug_indications <- drug_ids |> 
  left_join(ot_drug_indications, by = join_by(drugId == id))

ot_di_list <- lapply(
  ot_ctrp_drug_indications[["indications"]], 
  function(x) x[["efoName"]]
)

names(ot_di_list) <- ot_ctrp_drug_indications[["synonyms"]]

ot_indications_subset <- Reduce(union, ot_di_list)

ot_indications_subset[grepl("adren", ot_indications_subset)]

# What to look for for each TCGA project
tcga_cancer_string_map <- list(
  ACC = "adrenal cortex carcinoma",
  BLCA = c(
    "urinary bladder carcinoma",
    "urinary bladder cancer",
    "bladder tumor"
  ), 
  BRCA = c(
    "breast cancer", 
    "breast neoplasm", 
    "breast carcinoma", 
    "triple-negative breast cancer", 
    "breast ductal carcinoma in situ", 
    "invasive lobular carcinoma",
    "her2 positive breast carcinoma"
  ), 
  CESC = c(
    "cervical cancer", 
    "cervical carcinoma", 
    "cervical adenocarcinoma"
  ),
  CHOL = "cholangiocarcinoma", 
  COAD = c(
    "colon carcinoma", 
    "colon adenocarcinoma", 
    "colonic neoplasm", 
    "malignant colonic neoplasm"
  ), 
  DLBC = c(
    "lymphoma", 
    "non-hodgkins lymphoma", 
    "diffuse large b-cell lymphoma", 
    "neoplasm of mature b-cells", 
    "b-cell neoplasm"
  ), 
  ESCA = c(
      "esophageal cancer", 
      "esophageal carcinoma", 
      "esophageal squamous cell carcinoma", 
      "neoplasm of esophagus"
  ),
  GBM = c(
    "glioblastoma multiforme", 
    "brain cancer"
  ), 
  HNSC = paste(
    "head and neck",
    c(
      "malignant neoplasia", 
      "squamous cell carcinoma", 
      "carcinoma"
    )
  ), 
  KICH = c(
    "kidney cancer", 
    "kidney neoplasm", 
    "renal carcinoma", 
    "renal cell carcinoma", 
    "chromophobe renal cell carcinoma"
  ), 
  KIRC = c(
    "kidney cancer", 
    "kidney neoplasm", 
    "renal carcinoma", 
    "renal cell carcinoma", 
    "clear cell renal carcinoma"
  ), 
  KIRP = c(
    "kidney cancer", 
    "kidney neoplasm", 
    "renal carcinoma", 
    "renal cell carcinoma", 
    "papillary renal cell carcinoma"
  ), 
  LGG = c(
    "brain cancer", 
    "glioma", 
    "astrocytoma", 
    "oligodendroglioma"
  ),
  LIHC = c(
    "liver cancer", 
    "hepatocellular carcinoma", 
    "liver neoplasm"
  ),
  LUAD = c(
    "lung cancer", 
    "lung neoplasm", 
    "lung carcinoma", 
    "lung adenocarcinoma"
  ),
  LUSC = c(
    "lung cancer", 
    "lung neoplasm", 
    "lung carcinoma", 
    "squamous cell lung carcinoma"
  ),
  MESO = "mesothelioma",
  OV = paste(
    "ovarian",
    c(
      "cancer", 
      "neoplasm", 
      "carcinoma"
    )
  ), 
  PAAD = c(
    "pancreatic neoplasm", 
    "malignant pancreatic neoplasm", 
    "pancreatic carcinoma", 
    "pancreatic ductal adenocarcinoma"
  ),
  PCPG = "adrenal gland pheochromocytoma", # no match for paraganglioma
  PRAD = paste(
    "prostate",
    c(
      "cancer", 
      "carcinoma", 
      "adenocarcinoma"
    )
  ), 
  READ = c(
    "colorectal cancer", 
    "colorectal neoplasm", 
    "colorectal carcinoma",
    "colorectal adenocarcinoma", 
    "rectum cancer",
    "rectal carcinoma"
  ),
  SARC = c( # no match for "myxofibrosarcoma"
    "sarcoma", 
    "liposarcoma", 
    "dedifferentiated liposarcoma", 
    "leiomyosarcoma",
    "histiocytoma", # partially matches undifferentiated pleomorphic sarcoma in oncotree
    "malignant peripheral nerve sheath tumor", 
    "synovial sarcoma"
  ),
  SKCM = c(
    "skin cancer", 
    "skin neoplasm", 
    "melanoma",
    "cutaneous melanoma"
  ),
  STAD = c(
    "gastric cancer",
    "gastric carcinoma",
    "gastric adenocarcinoma",
    "stomach neoplasm"
  ),
  TGCT = c(
    "testicular cancer",
    "testicular neoplasm",
    "testicular carcinoma",
    "testicular yolk sac tumor", 
    "seminoma"
  ),
  THCA = c(
    "thyroid cancer",
    "thyroid neoplasm",
    "thyroid carcinoma",
    "papillary thyroid carcinoma"
  ),
  THYM = c(
    "thymus cancer",
    "thymus neoplasm",
    "thymoma"
  ),
  UCEC = c(
    "endometrial cancer",
    "endometrial neoplasm",
    "endometrial carcinoma"
  ),
  UCS = c(
    "uterine cancer",
    "uterine neoplasm",
    "uterine carcinosarcoma"
  ),
  UVM = c(
    NA
  )
)

tcga_cancer_string_union <- Reduce(union, tcga_cancer_string_map)
tcga_di_ind <- sapply(
  ot_di_list, 
  function(x) any(x %in% tcga_cancer_string_union)
)

ot_di_df <- bind_rows(lapply(
  names(ot_di_list),
  function(i) data.frame(drug = i, indication = ot_di_list[[i]])
))

tcga_ot_di_df <- bind_rows(lapply(
  names(tcga_cancer_string_map), 
  function(i) left_join(
    data.frame(tcga_type = i, indication = tcga_cancer_string_map[[i]]), 
    ot_di_df, 
    join_by(indication)
  )
))
tcga_ot_di_df[["drug"]] <- only_alphanumericals(tcga_ot_di_df[["drug"]])
tcga_ot_di_df <- distinct(tcga_ot_di_df)

tcga_ot_di_df <- tcga_ot_di_df |> 
  filter(!is.na(indication))

tcga_ot_di_df |> 
  count(tcga_type)

fnw <- paste0(ctrp_path, "open_targets_indications.csv")
readr::write_csv(tcga_ot_di_df, fnw)

if (FALSE) {
  # String based matching to find cancer drugs
  table(sapply(ot_di_list, function(x) any(grepl("cancer", x))))
  table(sapply(ot_di_list, function(x) any(grepl("carcinoma", x))))
  table(
    cancer = sapply(ot_di_list, function(x) any(grepl("cancer", x))), 
    carcinoma = sapply(ot_di_list, function(x) any(grepl("carcinoma", x)))
  )
  table(
    cancer = sapply(ot_di_list, function(x) any(grepl("cancer", x))), 
    lymphoma = sapply(ot_di_list, function(x) any(grepl("lymphoma", x)))
  )
  table(
    carcinoma = sapply(ot_di_list, function(x) any(grepl("carcinoma", x))), 
    lymphoma = sapply(ot_di_list, function(x) any(grepl("lymphoma", x)))
  )
  table(
    carcinoma = sapply(ot_di_list, function(x) any(grepl("carcinoma", x))), 
    sarcoma = sapply(ot_di_list, function(x) any(grepl("sarcoma", x)))
  )
  table(
    carcinoma = sapply(ot_di_list, function(x) any(grepl("carcinoma", x))), 
    glioma = sapply(ot_di_list, function(x) any(grepl("glioma", x)))
  )
  cancer_strings <- paste(
    c(
      "cancer", 
      "carcinoma", 
      "lymphoma", 
      "sarcoma", 
      "glioma", 
      "leukemia", 
      "myeloma",
      "myelofibrosis", 
      "mesothelioma", 
      "blastoma", 
      "melanoma",
      "neoplasm", 
      "neoplasia", 
      "tumor"
    ), 
    collapse = "|"
  )
  ot_di_cancer_list <- lapply(
    ot_di_list, 
    function(x) x[grep(cancer_strings, x)]
  )
  cancer_drug_ind <- sapply(ot_di_cancer_list, length) > 0
  ot_di_cancer_list <- ot_di_cancer_list[cancer_drug_ind]
}

