source("../setup.R")

CCLE <- PharmacoGx::downloadPSet("CCLE_2015", saveDir = ccle_path)
PRISM <- PharmacoGx::downloadPSet("PRISM_2020", saveDir = ccle_path)
CTRP <- PharmacoGx::downloadPSet("CTRPv2_2015", saveDir = ctrp_path)

# Cell-line info straight from DepMap
ccle_model <- read.csv(paste0(ccle_path, "Model.csv"), row.names = NULL, header = TRUE)

# Sensitivity tables can be saved as is
write.csv(
  PharmacoGx::sensitivityProfiles(CCLE),
  gzfile(paste0(ccle_path, "drug_sensitivity/screen_all.csv.gz"))
)
write.csv(
  PharmacoGx::sensitivityProfiles(PRISM),
  gzfile(paste0(ccle_path, "prism/screen_all.csv.gz"))
)
write.csv(
  PharmacoGx::sensitivityProfiles(CTRP),
  gzfile(paste0(ctrp_path, "screen_all.csv.gz"))
)

# Sensitivity info tables contains treatment ids
write.csv(
  PharmacoGx::sensitivityInfo(CCLE),
  gzfile(paste0(ccle_path, "drug_sensitivity/screen_rowinfo.csv.gz")))

# CTRP sensitivity info table must be fixed
new_dr_sens_info <- PharmacoGx::sensitivityInfo(CTRP)
# all row ids are structured the same
if (FALSE) table(sapply(strsplit(new_dr_sens_info$experimentIds, split = "_"), length))
new_dr_sens_info$treatmentid_fixed <- sapply(
  strsplit(new_dr_sens_info$experimentIds, split = "_"), 
  function(x) x[[2]])
# Sometimes CTRP cell-lines were cultured in two different mediums
if (FALSE) table(table(new_dr_sens_info[c("sampleid", "treatmentid_fixed", "culture_media")])) 
if (FALSE) sum(duplicated(new_dr_sens_info[,c("sampleid", "treatmentid_fixed", "culture_media", "experiment_id")]))

dr_sens_repl_dupl <- duplicated(new_dr_sens_info[,c("sampleid", "treatmentid_fixed", "culture_media")])
dr_sens_replicates <- plyr::join(
  unique(new_dr_sens_info[dr_sens_repl_dupl, c("sampleid", "treatmentid_fixed", "culture_media")]), 
  new_dr_sens_info, 
  by = c("sampleid", "treatmentid_fixed", "culture_media"))
# Range of concentrations is the same for all replicates
if (FALSE) {
  table(tapply(dr_sens_replicates$chosen.min.range, 
               paste(dr_sens_replicates$sampleid, 
                     dr_sens_replicates$treatmentid_fixed, 
                     dr_sens_replicates$culture_media), var) > 0) # None
  table(tapply(dr_sens_replicates$chosen.max.range, 
               paste(dr_sens_replicates$sampleid, 
                     dr_sens_replicates$treatmentid_fixed, 
                     dr_sens_replicates$culture_media), var) > 0) # None
}

test <- PharmacoGx::sensitivityProfiles(CTRP)[dr_sens_replicates$experimentIds,]
test$EID_new <- gsub("_[0-9]+", "", rownames(test))
if (FALSE) table(table(test$EID_new)) # all 2
test$experiment_id <- sapply(strsplit(rownames(test), split = "_"), function(x) x[[4]])

# Check correlation of replicated experimental results
if (FALSE) {
  experiment_pairs <- tapply(test$experiment_id, test$EID_new, list)
  experiment_pairs <- Reduce("rbind", experiment_pairs)
  experiment_pairs <- unique(experiment_pairs)
  
  experiment_pair_map <- rep(c(1,2), times = rep(nrow(experiment_pairs), 2))
  names(experiment_pair_map) <- as.vector(experiment_pairs)
  test$experiment_pair <- experiment_pair_map[test$experiment_id]
  
  test_split <- split(test, f = test$experiment_pair)
  test_ind <- match(test_split[[2]]$EID_new, test_split[[1]]$EID_new)
  plot(test_split[[1]][,"aac_recomputed"], test_split[[2]][test_ind,"aac_recomputed"])
  cor(
    test_split[[1]][,"aac_recomputed"], 
    test_split[[2]][test_ind,"aac_recomputed"], 
    method = "pearson", 
    use = "complete.obs"
  ) # 0.794
}

write.csv(
  new_dr_sens_info[
    ,c("sampleid", "treatmentid_fixed", "culture_media", "experiment_id"), 
    drop = FALSE
  ],
  gzfile(paste0(ctrp_path, "screen_rowinfo.csv.gz"))
)

## RNA-Seq and sensitivity data have different IDs that must be matched
# Cell line id map for CCLE
ccle_processed_info <- PharmacoGx::sampleInfo(CCLE)
# Have to use sample id to match sensitivity to gene-expression
if (FALSE) intersect(
  colnames(PharmacoGx::sensitivityInfo(CCLE)), 
  colnames(ccle_processed_info))
# PharmacoGx ids match
if (FALSE) table(
  matches = PharmacoGx::sensitivityInfo(CCLE)[["sampleid"]] %in% ccle_processed_info[["sampleid"]], 
  na = is.na(PharmacoGx::sensitivityInfo(CCLE)[["sampleid"]]), 
  useNA = "ifany")
# CCLE.name has 989 matches and a few NAs
if (FALSE) table(
  matches = ccle_processed_info[["CCLE.name"]] %in% ccle_model[["CCLEName"]], 
  na = is.na(ccle_processed_info[["CCLE.name"]]), 
  useNA = "always")
# sampleid has 848 matches and zero NAs
if (FALSE) table(
  matches = ccle_processed_info[["sampleid"]] %in% ccle_model[["CellLineName"]], 
  na = is.na(ccle_processed_info[["sampleid"]]), 
  useNA = "always")
# CCLE_rnaseq.sampleid has 821 matches and many NAs
if (FALSE) table(
  matches = ccle_processed_info[["CCLE_rnaseq.sampleid"]] %in% ccle_model[["CellLineName"]], 
  na = is.na(ccle_processed_info[["CCLE_rnaseq.sampleid"]]), 
  useNA = "always")
# CCLE.name maximizes matches
ccle_model_ind <- match(
  ccle_processed_info[["CCLE.name"]], 
  ccle_model[["CCLEName"]])

ccle_processed_info[["CCLE_model_id"]] <- ccle_model[ccle_model_ind, "ModelID"]
write.csv(
  ccle_processed_info, 
  gzfile(paste0(ccle_path, "drug_sensitivity/ccle_clid_map.csv.gz")))

# Cell line id map for CTRP
ctrp_cell_line_ids <- unique(PharmacoGx::sensitivityInfo(CTRP)$sampleid)
ctrp_cl_id <- PharmacoGx::sampleInfo(CTRP)[ctrp_cell_line_ids, "ccl_name"]
# Some CL ids do not match between PharmacoGx 
if (FALSE) table(ctrp_cl_id == ctrp_cell_line_ids, useNA = "always")
# Compare mis-matched ids, difference in capitalization and dashes
if (FALSE) data.frame(id1 = ctrp_cell_line_ids, id2 = ctrp_cl_id)[ctrp_cl_id != ctrp_cell_line_ids,]
# 183 ctrp ids do not match
if (FALSE) table(ctrp_cell_line_ids %in% ccle_model$CellLineName)
# Using stripped ids reduces mismatch to 47
if (FALSE) table(ctrp_cl_id %in% ccle_model$StrippedCellLineName)
ctrp_model_ind1 <- match(ctrp_cell_line_ids, ccle_model$CellLineName)
ctrp_model_ind2 <- match(ctrp_cl_id, ccle_model$StrippedCellLineName)
# Two additional matches can be made using non-stripped ids
if (FALSE) table(is.na(ctrp_model_ind1), is.na(ctrp_model_ind2))
# strip all non alpha-numeric characters
strip_nan <- function(x) gsub("^a-z0-9", "", x)
ctrp_model_ind3 <- match(
  strip_nan(tolower(ctrp_cell_line_ids)), 
  strip_nan(tolower(ccle_model$CellLineName))
)
ctrp_model_ind4 <- match(
  strip_nan(tolower(ctrp_cl_id)), 
  strip_nan(tolower(ccle_model$StrippedCellLineName))
)
# Manual stripping does not increase matches
if (FALSE) table(is.na(ctrp_model_ind3), is.na(ctrp_model_ind4))
# Compare mis-matched ids
if (FALSE) {
  data.frame(
    id1 = ctrp_cell_line_ids[is.na(ctrp_model_ind2) & !is.na(ctrp_model_ind1)], 
    id2 = ctrp_cl_id[is.na(ctrp_model_ind2) & !is.na(ctrp_model_ind1)])
  ccle_model[grep("TT", ccle_model$StrippedCellLineName),]
  # Stripping does not resolve this match
  data.frame(
    id3 = ctrp_cell_line_ids[is.na(ctrp_model_ind4) & !is.na(ctrp_model_ind3)], 
    id4 = ctrp_cl_id[is.na(ctrp_model_ind4) & !is.na(ctrp_model_ind3)])
  # No hit CLs
  data.frame(
    id3 = ctrp_cell_line_ids[is.na(ctrp_model_ind4) & is.na(ctrp_model_ind3)], 
    id4 = ctrp_cl_id[is.na(ctrp_model_ind4) & is.na(ctrp_model_ind3)])
  # Potential matches (1-11 not checked yet)
  ccle_model[grep("TM87", ccle_model$StrippedCellLineName),]
  ccle_model[grep("578", ccle_model$StrippedCellLineName),]
}
ctrp_model_id <- ccle_model$ModelID[ctrp_model_ind4]

ctrp_processed_info <- PharmacoGx::sampleInfo(CTRP)[ctrp_cell_line_ids,]
ctrp_processed_info[["CCLE_model_id"]] <- ctrp_model_id
#previous_version <- read.csv(paste0(ctrp_path, "ctrp_ccle_clid_map.csv.gz"), row.names = 1, header = TRUE)
write.csv(ctrp_processed_info, gzfile(paste0(ctrp_path, "ctrp_ccle_clid_map.csv.gz")))

