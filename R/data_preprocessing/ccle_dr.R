source("../setup.R")

ccle_info <- read.csv(paste0(ccle_path, "sample_info.csv"), 
                      header = TRUE, row.names = 1)
ccle_data0 <- read.csv(paste0(ccle_path, "prism/primary-screen-mfi-normalized.csv"), 
                       header = TRUE, row.names = 1)

real_colnames <- readLines(paste0(ccle_path, "prism/secondary-screen-replicate-collapsed-logfold-change.csv"), 1)
real_colnames <- strsplit(real_colnames, split = ",")[[1]]
real_colnames <- strsplit(real_colnames, split = "::")
ncol <- sapply(real_colnames, length)
extend_vector <- function(x, len) if(length(x) < len) replace(x, (length(x)+1):len, NA) else x
real_colnames_extended <- lapply(real_colnames, extend_vector, len = max(ncol))
drug_ids <- Reduce("rbind", real_colnames_extended[2:length(real_colnames_extended)])
rownames(drug_ids) <- 1:nrow(drug_ids)
colnames(drug_ids) <- c("broad_id", "dose", "screen_id", "compound_plate")

ccle_data <- read.csv(paste0(ccle_path, "prism/secondary-screen-replicate-collapsed-logfold-change.csv"), 
                      header = TRUE, row.names = 1)
ccle_dr_details <- read.csv(paste0(ccle_path, "prism/secondary-screen-replicate-treatment-info.csv"), 
                            header = TRUE, row.names = 1)

#table(drug_ids[!is.na(drug_ids[,"unknown"]), "unknown"])
#sapply(ccle_dr_details, function(x) any(grepl("PR500", x)))
#table(ccle_dr_details$compound_plate, useNA = "always") # no NA, so data IDs are incomplete

ccle_dr_details[,!colnames(ccle_dr_details) %in% c("moa", "target", "smiles")]


ccle_data[,1]

# dose-response curve
ccle_dr_data <- read.csv(paste0(ccle_path, "prism/secondary-screen-dose-response-curve-parameters.csv"), 
                         header = TRUE, row.names = NULL)

hist(ccle_dr_data[ccle_dr_data[,"broad_id"] == "BRD-A00077618-236-07-6","ic50"])

table(over_limit = ccle_dr_data$ec50 > ccle_dr_data$upper_limit, 
      under_limit = ccle_dr_data$ec50 < ccle_dr_data$lower_limit, 
      useNA = "always")

table(limit_inversion = ccle_dr_data$upper_limit < ccle_dr_data$lower_limit,
      ic50_na = is.na(ccle_dr_data$ic50))

table(limit_inversion = ccle_dr_data$upper_limit < ccle_dr_data$lower_limit,
      auc_over_1 = ccle_dr_data$auc > 1)

# dr curve summaries
test <- ccle_dr_data
test <- data.table::setDT(test)
test2 <- test[,mean(log10(ic50), na.rm = TRUE), broad_id]

# Expression data (RNA-Seq)
ccle_expression <- read.csv(paste0(ccle_path, "CCLE_expression.csv"), 
                            header = TRUE, row.names = 1)


COPS::pca_viz(log2(ccle_expression+1), 1, "none")


# GDSC
gdsc1 <- readxl::read_xlsx(paste0(gdsc_path, "GDSC1_fitted_dose_response_24Jul22.xlsx"))

table(over_limit = exp(gdsc1$LN_IC50) > gdsc1$MAX_CONC, 
      under_limit = exp(gdsc1$LN_IC50) < gdsc1$MIN_CONC, 
      useNA = "always")

gdsc2 <- readxl::read_xlsx(paste0(gdsc_path, "GDSC2_fitted_dose_response_24Jul22.xlsx"))

table(over_limit = exp(gdsc2$LN_IC50) > gdsc2$MAX_CONC, 
      under_limit = exp(gdsc2$LN_IC50) < gdsc2$MIN_CONC, 
      useNA = "always")

hist(gdsc2$AUC)

# PharmacoGx drug response data and old multi-omic data
if (TRUE) {
  # "Download" from disk
  psets <- PharmacoGx::availablePSets()
  
  GDSC1 <- PharmacoGx::downloadPSet("GDSC_2020(v1-8.2)", saveDir = gdsc_path)
  GDSC2 <- PharmacoGx::downloadPSet("GDSC_2020(v2-8.2)", saveDir = gdsc_path)
  
  CCLE <- PharmacoGx::downloadPSet("CCLE_2015", saveDir = ccle_path)
  PRISM <- PharmacoGx::downloadPSet("PRISM_2020", saveDir = ccle_path)
  
  CTRP <- PharmacoGx::downloadPSet("CTRPv2_2015", saveDir = ctrp_path)
} else {
  # readRDS does not work
  GDSC1 <- readRDS(paste0(gdsc_path, "GDSC_2020(v1-8.2).rds"))
  GDSC2 <- readRDS(paste0(gdsc_path, "GDSC_2020(v2-8.2).rds"))
  
  CCLE <- readRDS(paste0(ccle_path, "CCLE_2015.rds"))
  PRISM <- readRDS(paste0(ccle_path, "PRISM_2020.rds"))
}

PharmacoGx::sensitivityInfo(GDSC1)
PharmacoGx::sensitivityInfo(GDSC2)
PharmacoGx::sensitivityInfo(PRISM)

# CCLE
PharmacoGx::sensitivityInfo(CCLE)
PharmacoGx::sensitivityMeasures(CCLE)
x <- PharmacoGx::sensitivityProfiles(CCLE)[,"ic50_published"]
y <- PharmacoGx::sensitivityProfiles(CCLE)[,"ic50_recomputed"]
plot(x,y) # x is corrupted?

write.csv(
  PharmacoGx::sensitivityProfiles(CCLE)[,"aac_recomputed", drop = FALSE],
  gzfile(paste0(ccle_path, "drug_sensitivity/screen_aac.csv.gz")))

write.csv(
  PharmacoGx::drugInfo(CCLE), 
  paste0(ccle_path, "drug_sensitivity/drug_info.csv"))

dr_sens_info <- PharmacoGx::sensitivityInfo(CCLE)

table(
  rownames(PharmacoGx::sensitivityProfiles(CCLE)) == 
    rownames(dr_sens_info)
)
ccle_processed_info <- PharmacoGx::sampleInfo(CCLE)
intersect(colnames(dr_sens_info), colnames(ccle_processed_info))

table(
  matches = dr_sens_info[["sampleid"]] %in% ccle_processed_info[["sampleid"]], 
  na = is.na(dr_sens_info[["sampleid"]]), 
  useNA = "ifany")

ccle_model <- read.csv(paste0(ccle_path, "Model.csv"), row.names = NULL, header = TRUE)

table(
  matches = ccle_processed_info[["CCLE.name"]] %in% ccle_model[["CCLEName"]], 
  na = is.na(ccle_processed_info[["CCLE.name"]]), 
  useNA = "always")


matches = ccle_processed_info[["sampleid"]] %in% ccle_model[["CellLineName"]]
table(
  matches = matches, 
  na = is.na(ccle_processed_info[["sampleid"]]), 
  useNA = "always")

matches = ccle_processed_info[["CCLE_rnaseq.sampleid"]] %in% ccle_model[["CellLineName"]]
table(
  matches = matches, 
  na = is.na(ccle_processed_info[["CCLE_rnaseq.sampleid"]]), 
  useNA = "always")

ccle_model_ind <- match(
  ccle_processed_info[["CCLE.name"]], 
  ccle_model[["CCLEName"]])

ccle_processed_info[["CCLE_model_id"]] <- ccle_model[ccle_model_ind, "ModelID"]
write.csv(
  ccle_processed_info, 
  gzfile(paste0(ccle_path, "drug_sensitivity/ccle_clid_map.csv.gz")))

# CCLE screen info
new_dr_sens_info <- PharmacoGx::sensitivityInfo(CCLE)
table(rownames(new_dr_sens_info) == rownames(PharmacoGx::sensitivityProfiles(CCLE)))

table(table(new_dr_sens_info[c("sampleid", "treatmentid")])) 
# Up to 2 replicates

write.csv(
  PharmacoGx::sensitivityInfo(CCLE),
  gzfile(paste0(ccle_path, "drug_sensitivity/screen_rowinfo.csv.gz")))

# PRISM
PharmacoGx::sensitivityInfo(PRISM)
PharmacoGx::sensitivityMeasures(PRISM)
x <- PharmacoGx::sensitivityProfiles(PRISM)[,"ic50_published"]
y <- PharmacoGx::sensitivityProfiles(PRISM)[,"ic50_recomputed"]
plot(log10(x),log10(y))

# Compare to manual download
ccle_pgx_id <- paste(ccle_dr_data$broad_id, 
                     ccle_dr_data$screen_id,
                     ccle_dr_data$depmap_id, 
                     ccle_dr_data$name, sep = "::")
x <- ccle_dr_data$ic50
y <- PharmacoGx::sensitivityProfiles(PRISM)[ccle_pgx_id,"ic50_published"]
table(x == y, useNA = "always")
table(is.na(x))
table(is.na(y)) 
table(is.na(PharmacoGx::sensitivityProfiles(PRISM)[,"ic50_published"]))
table(is.na(x))

# Area above/under curve 
x <- PharmacoGx::sensitivityProfiles(PRISM)[,"aac_recomputed"]
y <- PharmacoGx::sensitivityProfiles(PRISM)[,"auc_published"]
y[is.na(y)] <- 0
plot(x, y)

ccle_rows <- rownames(PharmacoGx::sensitivityProfiles(CCLE))
prism_rows <- rownames(PharmacoGx::sensitivityProfiles(PRISM))
ctrp_rows <- rownames(PharmacoGx::sensitivityProfiles(CTRP))

colnames(PharmacoGx::sensitivityInfo(CCLE))
colnames(PharmacoGx::sensitivityInfo(PRISM))
colnames(PharmacoGx::sensitivityInfo(CTRP))

ccle_cell_line_ids <- unique(PharmacoGx::sensitivityInfo(CCLE)$sampleid)
prism_cell_line_ids <- unique(PharmacoGx::sensitivityInfo(PRISM)$sampleid)
ccle_prism_union <- union(ccle_cell_line_ids, prism_cell_line_ids)
table(ccle = ccle_prism_union %in% ccle_cell_line_ids, 
      prism = ccle_prism_union %in% prism_cell_line_ids,
      useNA = "always")

ctrp_cell_line_ids <- unique(PharmacoGx::sensitivityInfo(CTRP)$sampleid)
ccle_ctrp_union <- union(ccle_cell_line_ids, ctrp_cell_line_ids)
table(ccle = ccle_ctrp_union %in% ccle_cell_line_ids, 
      ctrp = ccle_ctrp_union %in% ctrp_cell_line_ids,
      useNA = "always")

# Check cell line overlap with omics data
ccle_rna <- read.csv(paste0(ccle_path, "CCLE_expression.csv"), row.names = 1, header = TRUE)
ccle_model <- read.csv(paste0(ccle_path, "Model.csv"), row.names = NULL, header = TRUE)

ccle_cell_line_names <- PharmacoGx::sampleInfo(CCLE)[ccle_cell_line_ids, "CCLE.name"]
table(id = ccle_cell_line_ids %in% ccle_model$CellLineName,
      name = ccle_cell_line_names %in% ccle_model$CCLEName)
ccle_cell_line_names[!ccle_cell_line_names %in% ccle_model$CCLEName]
# Need to use both sampleid and CCLE.name
ccle_model_id_ind1 <- match(ccle_cell_line_names, ccle_model$CCLEName)
ccle_model_id_ind2 <- match(ccle_cell_line_ids, ccle_model$CellLineName)
table(ccle_model_id_ind1 == ccle_model_id_ind2, useNA = "always")
# Not always true
test1 <- ccle_model_id_ind1[ccle_model_id_ind1 != ccle_model_id_ind2]
ccle_model[test1[!is.na(test1)],]
test2 <- ccle_model_id_ind2[ccle_model_id_ind2 != ccle_model_id_ind1]
ccle_model[test2[!is.na(test2)],]
test3 <- ccle_model_id_ind2 != ccle_model_id_ind1
test3[is.na(test3)] <- FALSE
PharmacoGx::sampleInfo(CCLE)[ccle_cell_line_ids[test3],]
# Some discrepancy between sampleid and CCLE model ids
ccle_cl_id <- PharmacoGx::sampleInfo(CCLE)[ccle_cell_line_ids, "Cell.line.primary.name"]
ccle_model_id_ind3 <- match(ccle_cl_id, ccle_model$CellLineName)
table(name = is.na(ccle_model_id_ind1), 
      id = is.na(ccle_model_id_ind3)) # 15 not matching
table(ccle_model_id_ind1 == ccle_model_id_ind3, useNA = "always")
# 472 true, 31 NA
ccle_model_id_ind <- ifelse(is.na(ccle_model_id_ind1), 
                            ccle_model_id_ind3, 
                            ccle_model_id_ind1)
ccle_model_id <- ccle_model$ModelID[ccle_model_id_ind]
ccle_model[ccle_model$ModelID %in% ccle_model_id[!ccle_model_id %in% rownames(ccle_rna)],]
ccle_model_id_ind4 <- match(ccle_cl_id, ccle_model$StrippedCellLineName)
table(comb = is.na(ccle_model_id_ind), 
      strip_id = is.na(ccle_model_id_ind4)) # gain one more match
table(ccle_model_id_ind == ccle_model_id_ind4, useNA = "always") # 92 true, 411 NA
# Probably not worth it to use this ID for 1 more sample
table(ccle_model_id %in% rownames(ccle_rna)) # 471 true, 32 false

prism_model_id <- PharmacoGx::sampleInfo(PRISM)[prism_cell_line_ids, "depmap_id"]
table(is.na(prism_model_id)) # 0
table(prism_model_id  %in% rownames(ccle_rna)) # 476 true, 4 false

# CTRP cell-line id map
ctrp_cell_line_ids <- unique(PharmacoGx::sensitivityInfo(CTRP)$sampleid)
ctrp_cl_id <- PharmacoGx::sampleInfo(CTRP)[ctrp_cell_line_ids, "ccl_name"]
table(ctrp_cl_id == ctrp_cell_line_ids, useNA = "always")
data.frame(id1 = ctrp_cell_line_ids, id2 = ctrp_cl_id)[ctrp_cl_id != ctrp_cell_line_ids,]
table(ctrp_cell_line_ids %in% ccle_model$CellLineName) # 183 false
table(ctrp_cl_id %in% ccle_model$StrippedCellLineName) # 47 false
ctrp_model_ind1 <- match(ctrp_cell_line_ids, ccle_model$CellLineName)
ctrp_model_ind2 <- match(ctrp_cl_id, ccle_model$StrippedCellLineName)
table(is.na(ctrp_model_ind1), is.na(ctrp_model_ind2))
ctrp_model_id <- ccle_model$ModelID[ctrp_model_ind2]
table(ctrp_model_id  %in% rownames(ccle_rna)) # 822 true, 65 false

ctrp_processed_info <- PharmacoGx::sampleInfo(CTRP)[ctrp_cell_line_ids,]
ctrp_processed_info[["CCLE_model_id"]] <- ctrp_model_id
prev <- read.csv(paste0(ctrp_path, "ctrp_ccle_clid_map.csv.gz"), row.names = 1, header = TRUE)
write.csv(ctrp_processed_info, gzfile(paste0(ctrp_path, "ctrp_ccle_clid_map.csv.gz")))


model_union <- Reduce(union, list(ccle_model_id, prism_model_id, ctrp_model_id))
table(ccle = model_union %in% ccle_model_id, 
      prism = model_union %in% prism_model_id)
table(ccle = model_union %in% ccle_model_id, 
      ctrp = model_union %in% ctrp_model_id)
table(prism = model_union %in% prism_model_id, 
      ctrp = model_union %in% ctrp_model_id)

# Drug overlap
drug_union <- Reduce(union, lapply(list(CCLE, PRISM, CTRP), PharmacoGx::drugNames))

table(ccle = drug_union %in% PharmacoGx::drugNames(CCLE),
      prism = drug_union %in% PharmacoGx::drugNames(PRISM))
table(ccle = drug_union %in% PharmacoGx::drugNames(CCLE),
      ctrp = drug_union %in% PharmacoGx::drugNames(CTRP))
table(prism = drug_union %in% PharmacoGx::drugNames(PRISM),
      ctrp = drug_union %in% PharmacoGx::drugNames(CTRP))

# CTRP
x <- PharmacoGx::sensitivityProfiles(CTRP)[,"aac_published"]
y <- PharmacoGx::sensitivityProfiles(CTRP)[,"aac_recomputed"]
plot(x, y)
table(is.na(PharmacoGx::sensitivityProfiles(CTRP)[,"ic50_recomputed"]))

write.csv(PharmacoGx::sensitivityProfiles(CTRP)[,"ic50_recomputed", drop = FALSE],
          gzfile(paste0(ctrp_path, "screen_ic50.csv.gz")))
write.csv(PharmacoGx::sensitivityProfiles(CTRP)[,"aac_recomputed", drop = FALSE],
          gzfile(paste0(ctrp_path, "screen_aac.csv.gz")))
#write.csv(PharmacoGx::sensitivityInfo(CTRP)[,c("sampleid", "treatmentid", "culture_media"), drop = FALSE],
#          gzfile(paste0(ctrp_path, "screen_rowinfo.csv.gz")))

write.csv(PharmacoGx::drugInfo(CTRP), paste0(ctrp_path, "drug_info.csv"))






# check PharmacoGx molecular profiles
CCLE_rna <- PharmacoGx::molecularProfiles(CCLE, mDataType = "rna")
colnames(CCLE_rna)

# Check duplicated experiments
dr_sens_info <- PharmacoGx::sensitivityInfo(CTRP)[,c("sampleid", "treatmentid", "culture_media"), drop = FALSE]
dr_index_dupl <- duplicated(dr_sens_info)
dr_dupl <- plyr::join(unique(dr_sens_info[dr_index_dupl,]), 
                      PharmacoGx::sensitivityInfo(CTRP), 
                      by = c("sampleid", "treatmentid", "culture_media"))
# experiment_id is sometimes different (replicates?)
dr_sens_info <- PharmacoGx::sensitivityInfo(CTRP)[,c("sampleid", "treatmentid", "culture_media", "experiment_id"), drop = FALSE]
dr_index_dupl <- duplicated(dr_sens_info)
dr_dupl <- plyr::join(unique(dr_sens_info[dr_index_dupl,]), 
                      PharmacoGx::sensitivityInfo(CTRP), 
                      by = c("sampleid", "treatmentid", "culture_media", "experiment_id"))

table(apply(is.na(dr_dupl), 1, sum))
dr_dupl[apply(is.na(dr_dupl), 1, any),]
table(dr_dupl$treatmentid) # all Tipifarnib

test <- PharmacoGx::sensitivityProfiles(CTRP)[dr_dupl$experimentIds,]
test$Phase <- gsub(".*?tipifarnib-", "", rownames(test))
test$Phase <- gsub("_.*", "", test$Phase)
test$EID_new <- gsub("-P[1-2]_", "_", rownames(test))

test <- split(test, f = test$Phase)
test[[2]] <- test[[2]][match(test[[2]]$EID_new, test[[1]]$EID_new),]
plot(test[[1]]$aac_recomputed, test[[2]]$aac_recomputed)

#PharmacoGx::sensitivityProfiles(CTRP)[,"aac_recomputed", drop = FALSE]
new_dr_sens_info <- PharmacoGx::sensitivityInfo(CTRP)
table(sapply(strsplit(new_dr_sens_info$experimentIds, split = "_"), length))
new_dr_sens_info$treatmentid_fixed <- sapply(strsplit(new_dr_sens_info$experimentIds, split = "_"), 
                                             function(x) x[[2]])

table(table(new_dr_sens_info[c("sampleid", "treatmentid_fixed", "culture_media")])) 
sum(duplicated(new_dr_sens_info[,c("sampleid", "treatmentid_fixed", "culture_media", "experiment_id")]))
# Up to 2 replicates

dr_sens_repl_dupl <- duplicated(new_dr_sens_info[,c("sampleid", "treatmentid_fixed", "culture_media")])
dr_sens_replicates <- plyr::join(unique(new_dr_sens_info[dr_sens_repl_dupl, c("sampleid", "treatmentid_fixed", "culture_media")]), 
                                 new_dr_sens_info, 
                                 by = c("sampleid", "treatmentid_fixed", "culture_media"))
table(tapply(dr_sens_replicates$chosen.min.range, 
             paste(dr_sens_replicates$sampleid, 
                   dr_sens_replicates$treatmentid_fixed, 
                   dr_sens_replicates$culture_media), var) > 0) # None
table(tapply(dr_sens_replicates$chosen.max.range, 
             paste(dr_sens_replicates$sampleid, 
                   dr_sens_replicates$treatmentid_fixed, 
                   dr_sens_replicates$culture_media), var) > 0) # None

test <- PharmacoGx::sensitivityProfiles(CTRP)[dr_sens_replicates$experimentIds,]
test$EID_new <- gsub("_[0-9]+", "", rownames(test))
table(table(test$EID_new)) # all 2
test$experiment_id <- sapply(strsplit(rownames(test), split = "_"), function(x) x[[4]])

# Check correlation of replicated experimental results
experiment_pairs <- tapply(test$experiment_id, test$EID_new, list)
experiment_pairs <- Reduce("rbind", experiment_pairs)
experiment_pairs <- unique(experiment_pairs)

experiment_pair_map <- rep(c(1,2), times = rep(nrow(experiment_pairs), 2))
names(experiment_pair_map) <- as.vector(experiment_pairs)
test$experiment_pair <- experiment_pair_map[test$experiment_id]

test_split <- split(test, f = test$experiment_pair)
test_ind <- match(test_split[[2]]$EID_new, test_split[[1]]$EID_new)
plot(test_split[[1]][,"aac_recomputed"], test_split[[2]][test_ind,"aac_recomputed"])

write.csv(new_dr_sens_info[,c("sampleid", "treatmentid_fixed", "culture_media", "experiment_id"), drop = FALSE],
          gzfile(paste0(ctrp_path, "screen_rowinfo.csv.gz")))

table(table(new_dr_sens_info[c("treatmentid_fixed", "culture_media")])) 
table(table(new_dr_sens_info[c("sampleid", "treatmentid_fixed", "experiment_id")])) # unique ID
table(table(new_dr_sens_info[c("sampleid", "treatmentid_fixed")]))

test <- new_dr_sens_info
test$EID_new <- gsub("_[0-9]+", "", rownames(test))
table(table(test$EID_new)) # all 2
test$experiment_id <- sapply(strsplit(rownames(test), split = "_"), function(x) x[[4]])

experiment_pairs <- tapply(test$experiment_id, test$EID_new, list)
experiment_pairs <- lapply(experiment_pairs, function(x) if (length(x)==1) return(c(1, NA)) else return(c(1,2)))
experiment_pairs <- Reduce("rbind", experiment_pairs)
experiment_pairs <- unique(experiment_pairs)

experiment_pair_map <- rep(c(1,2), times = rep(nrow(experiment_pairs), 2))
names(experiment_pair_map) <- as.vector(experiment_pairs)
test$experiment_pair <- experiment_pair_map[test$experiment_id]

table(table(new_dr_sens_info[c("sampleid", "treatmentid_fixed", "experiment_id")]))

# Try to decouple experiment_id from other IDs to check medium repeats
dr_sens_key <- c("sampleid", "treatmentid_fixed", "culture_media")
dr_sens_repl_dupl <- duplicated(new_dr_sens_info[,dr_sens_key])
dr_sens_replicates <- plyr::join(unique(new_dr_sens_info[dr_sens_repl_dupl, dr_sens_key]), 
                                 new_dr_sens_info, 
                                 by = dr_sens_key)

dr_sens_noreplicates <- plyr::join(unique(new_dr_sens_info[!dr_sens_repl_dupl, dr_sens_key]), 
                                   new_dr_sens_info, 
                                   by = dr_sens_key)
dr_sens_noreplicates$experiment_id_fixed <- 1
dr_sens_replicates$EID_new <- gsub("_[0-9]+", "", dr_sens_replicates$experimentIds)

experiment_pairs <- tapply(dr_sens_replicates$experiment_id, dr_sens_replicates$EID_new, list)
experiment_pairs <- Reduce("rbind", experiment_pairs)
experiment_pairs <- unique(experiment_pairs)

experiment_pair_map <- rep(c(1,2), times = rep(nrow(experiment_pairs), 2))
names(experiment_pair_map) <- as.vector(experiment_pairs)
dr_sens_replicates$experiment_id_fixed <- experiment_pair_map[dr_sens_replicates$experiment_id]

dr_sens_combined <- rbind(dr_sens_noreplicates, dr_sens_replicates[,colnames(dr_sens_noreplicates)])
table(table(dr_sens_combined[c("sampleid", "treatmentid_fixed", "experiment_id_fixed")]))

# Check sensitivity correlation
test <- PharmacoGx::sensitivityProfiles(CTRP)[dr_sens_replicates$experimentIds,]
test$EID_new <- gsub("_[0-9]+", "", rownames(test))
table(table(test$EID_new)) # all 2
test$experiment_id <- sapply(strsplit(rownames(test), split = "_"), function(x) x[[4]])

experiment_pairs <- tapply(test$experiment_id, test$EID_new, list)
experiment_pairs <- Reduce("rbind", experiment_pairs)
experiment_pairs <- unique(experiment_pairs)

experiment_pair_map <- rep(c(1,2), times = rep(nrow(experiment_pairs), 2))
names(experiment_pair_map) <- as.vector(experiment_pairs)
test$experiment_pair <- experiment_pair_map[test$experiment_id]



# Wrong key
dr_sens_key <- c("sampleid", "treatmentid", "culture_media", "experiment_id")
dr_sens_info <- new_dr_sens_info[,dr_sens_key, drop = FALSE]
dr_index_dupl <- duplicated(dr_sens_info)
dr_dupl <- plyr::join(unique(dr_sens_info[dr_index_dupl,]), 
                      new_dr_sens_info, 
                      by = dr_sens_key)

table(apply(is.na(dr_dupl), 1, sum))
dr_dupl[apply(is.na(dr_dupl), 1, any),]
table(dr_dupl$treatmentid) # all Tipifarnib

# Reverse selection
dr_nodupl <- plyr::join(unique(dr_sens_info[!dr_index_dupl,]), 
                        PharmacoGx::sensitivityInfo(CTRP), 
                        by = dr_sens_key)

# Between media correlation
test <- PharmacoGx::sensitivityProfiles(CTRP)[dr_dupl$experimentIds,]
test$Phase <- gsub(".*?tipifarnib-", "", rownames(test))
test$Phase <- gsub("_.*", "", test$Phase)
test$EID_new <- gsub("-P[1-2]_", "_", rownames(test))
