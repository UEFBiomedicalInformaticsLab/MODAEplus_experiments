script_root <- Sys.getenv("MODAE_SCRIPT_PATH")
if (script_root == "") {
  stop("Please define the MODAE_SCRIPT_PATH environment variable.")
}
source(paste0(script_root, "R/setup.R"))

fn <- paste0(xia_path, "combined_cl_metadata")
cl_meta <- read.table(fn, sep = "\t", header = TRUE)

cl_meta[cl_meta[["dataset"]] == "CTRP", "core_str"]

ccle_model <- read.csv(paste0(ccle_path, "Model.csv"), row.names = NULL, header = TRUE)

model_idx1 <- match(cl_meta[["core_str"]], ccle_model[["CellLineName"]])
model_idx2 <- match(cl_meta[["core_str"]], ccle_model[["StrippedCellLineName"]])

table(name = is.na(model_idx1), stripped = is.na(model_idx2))

# strip all non alpha-numeric characters
strip_nan <- function(x) gsub("^a-z0-9", "", x)
model_idx3 <- match(
  strip_nan(tolower(cl_meta[["core_str"]])), 
  strip_nan(tolower(ccle_model[["CellLineName"]]))
)
model_idx4 <- match(
  strip_nan(tolower(cl_meta[["core_str"]])), 
  strip_nan(tolower(ccle_model[["StrippedCellLineName"]]))
)

table(name = is.na(model_idx3), stripped = is.na(model_idx4))

with(cl_meta[!is.na(model_idx3),], table(dataset))

id_split <- strsplit(cl_meta[["sample_name"]], split = "\\.")

table(sapply(id_split, length))
table(sapply(id_split, function(x) x[1]))
id_exdat <- sapply(id_split, function(x) paste(x[-1], collapse = "."))

table(id_exdat == cl_meta[["core_str"]])

model_idx5 <- match(
  strip_nan(tolower(id_exdat)), 
  strip_nan(tolower(ccle_model[["CellLineName"]]))
)
model_idx6 <- match(
  strip_nan(tolower(id_exdat)), 
  strip_nan(tolower(ccle_model[["StrippedCellLineName"]]))
)

table(name = is.na(model_idx5), stripped = is.na(model_idx6))

with(cl_meta[!is.na(model_idx5),], table(dataset))
with(cl_meta[!is.na(model_idx6),], table(dataset))
with(cl_meta[!is.na(model_idx5)|!is.na(model_idx6),], table(dataset))

ccle_xia_ind <- !is.na(model_idx6) & cl_meta[,"dataset"] == "CCLE"
ccle_xia_ids <- cl_meta[ccle_xia_ind, "sample_name"]

xia_depmap_map <- ccle_model[model_idx6, "ModelID"]
names(xia_depmap_map) <- cl_meta[["sample_name"]]

ctrp_xia_ind <- !is.na(model_idx6) & cl_meta[,"dataset"] == "CTRP"
ctrp_xia_ids <- cl_meta[ctrp_xia_ind, "sample_name"]

fn <- paste0(xia_path, "cl_mapping")
cl_mapping <- read.table(fn, sep = "\t", header = TRUE)

mapping_id_split <- strsplit(cl_mapping[[2]], split = "\\.")
table(sapply(mapping_id_split, function(x) x[1]))

fn <- paste0(xia_path, "combined_rnaseq_data")
cl_rnaseq <- read.table(fn, sep = "\t", header = TRUE, row.names = 1)

cl_rnaseq_project <- sapply(
  strsplit(rownames(cl_rnaseq), split = "\\."), 
  function(x) x[1]
)
names(cl_rnaseq_project) <- rownames(cl_rnaseq)
table(cl_rnaseq_project)

fn <- paste0(xia_path, "combined_single_response_agg")
dr_data <- readr::read_tsv(fn)

dr_cell_ids <- unique(dr_data[["CELL"]])
table(cl_rnaseq_project[dr_cell_ids])

table(strip_nan(tolower(dr_cell_ids)) %in% strip_nan(tolower(rownames(cl_rnaseq))))

ctrp_xia_dr_cell_ids <- unique(dr_data[["CELL"]][dr_data[["SOURCE"]] == "CTRP"])
ctrp_xia_dr_cell_ccle_idx <- match(
  ctrp_xia_dr_cell_ids, 
  cl_mapping[[2]]
)
ctrp_xia_dr_cell_ccle_id <- cl_mapping[ctrp_xia_dr_cell_ccle_idx,1]
table(ctrp_xia_dr_cell_ccle_id %in% ccle_xia_ids)

ctrp_ccle_xia_id_map <- data.frame(
  depmap_id = xia_depmap_map[ctrp_xia_dr_cell_ccle_id]
)
rownames(ctrp_ccle_xia_id_map) <- ctrp_xia_dr_cell_ids

fnw <- paste0(xia_path, "xia_id_ccle_clid_map.csv.gz")
write.csv(ctrp_ccle_xia_id_map, gzfile(fnw))

fnw <- paste0(xia_path, "xia_ctrp_screen_dss.csv.gz")
write.csv(
  as.data.frame(dr_data)[
    dr_data[["CELL"]] %in% ctrp_xia_dr_cell_ids, 
    c("DSS1"), 
    drop = FALSE
  ], 
  gzfile(fnw)
)

fnw <- paste0(xia_path, "xia_ctrp_screen_info.csv.gz")
write.csv(
  as.data.frame(dr_data)[
    dr_data[["CELL"]] %in% ctrp_xia_dr_cell_ids, 
    c("CELL", "DRUG"), 
    drop = FALSE
  ], 
  gzfile(fnw)
)

fn <- paste0(xia_path, "drug_info")
xia_drug_info <- readr::read_tsv(fn)


