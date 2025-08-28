ccle_model <- read.csv(paste0(ccle_path, "Model.csv"), row.names = NULL, header = TRUE)

ccle_model[["solid"]] <- "solid"
ccle_model[
  grep(
    "leukemia", 
    ccle_model[["OncotreePrimaryDisease"]], 
    ignore.case = TRUE), 
  "solid"] <- "liquid"
ccle_model[
  grep(
    "myelo", 
    ccle_model[["OncotreePrimaryDisease"]], 
    ignore.case = TRUE), 
  "solid"] <- "liquid"

ccle_model[["solid_primary"]] <- with(
  ccle_model, 
  ifelse(
    solid == "solid" & PrimaryOrMetastasis == "Primary", 
    "solid_primary", 
    "other")
)

ccle_model[["PrimaryOrMetastasis"]] <- with(ccle_model, ifelse(PrimaryOrMetastasis == "", NA, PrimaryOrMetastasis))
write.csv(ccle_model, paste0(ccle_path, "Model_augmented.csv"), row.names = FALSE)

oncotree <- jsonlite::fromJSON(
  "http://oncotree.mskcc.org/api/tumorTypes", 
  flatten = TRUE
)
oncotree_parent_map <- oncotree$parent
names(oncotree_parent_map) <- oncotree$code
oncotree_level_map <- oncotree$level
names(oncotree_level_map) <- oncotree$code

ccle_model_oncotree <- data.frame(ModelID = ccle_model$ModelID)
codesi <- ccle_model$OncotreeCode
for (i in rev(sort(unique(oncotree$level)))) {
  levelsi <- oncotree_level_map[codesi]
  ccle_model_oncotree[[paste0("level_", i)]] <- ifelse(levelsi == i, codesi, NA)
  codesi <- ifelse(levelsi == i, oncotree_parent_map[codesi], codesi)
}

rownames(ccle_model_oncotree) <- ccle_model_oncotree[["ModelID"]]
ccle_model_oncotree[["ModelID"]] <- NULL

write.csv(ccle_model_oncotree, file = paste0(ccle_path, "Model_oncotree.csv"))

# Code for interactive checks
if (FALSE) {
  ccle_model_oncotree[ccle_model_oncotree$level_2 == "BRCA" & !is.na(ccle_model_oncotree$level_2),]
  ccle_model_oncotree[ccle_model_oncotree$level_2 == "RCC" & !is.na(ccle_model_oncotree$level_2),]
  ccle_model_oncotree[ccle_model_oncotree$level_2 == "COADREAD" & !is.na(ccle_model_oncotree$level_2),]
  ccle_model_oncotree[ccle_model_oncotree$level_2 == "NSCLC" & !is.na(ccle_model_oncotree$level_2),]
  ccle_model_oncotree[ccle_model_oncotree$level_3 %in% c("LUAD", "LUSC") & !is.na(ccle_model_oncotree$level_3),]
  
  pan_cancer_list <- c(
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
    "UVM")
  
  tcga_oncotree <- data.frame(TCGA_histology = pan_cancer_list)
  codesi <- pan_cancer_list
  for (i in rev(sort(unique(oncotree$level)))) {
    levelsi <- oncotree_level_map[codesi]
    tcga_oncotree[[paste0("level_", i)]] <- ifelse(levelsi == i, codesi, NA)
    codesi <- ifelse(levelsi == i, oncotree_parent_map[codesi], codesi)
  }
  
}
