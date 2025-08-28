source("../setup.R")

ctrdb_dirs <- dir(ctrdb_path)
ctrdb_dirs <- ctrdb_dirs[-grep("data.tar", ctrdb_dirs)]
ctrdb_files <- grep("\\.csv|\\.txt|\\.png|\\.json", ctrdb_dirs)
if (length(ctrdb_files) > 0) {
  ctrdb_dirs <- ctrdb_dirs[-ctrdb_files]
}

ctrdb_table_list <- list()
for (i in ctrdb_dirs) {
  fni <- paste0(ctrdb_path, i, "/cli.inf_.csv")
  temp <- tryCatch(
    readr::read_csv(
      fni, 
      col_names = TRUE, 
      show_col_types = FALSE
    ), 
    error = function(e) {
      warning(paste(i , e));return(NULL)
    }
  )
  ctrdb_table_list[[i]] <- temp
}
#lapply(ctrdb_table_list, readr::problems)

lapply(ctrdb_table_list, function(x) table(x[["Resource"]], useNA = "always"))

tcga_ind <- sapply(ctrdb_table_list, function(x) "TCGA" %in% x[["Resource"]])
tcga_response <- Reduce("rbind", ctrdb_table_list[tcga_ind])

with(tcga_response, table(Source))
with(tcga_response, table(Source, Response))
with(tcga_response, table(Drug_list))

folfiri_idx <- grep("folfiri", tcga_response[["Drug_list"]], ignore.case = TRUE)
tcga_response[folfiri_idx, "Drug_list"] <- "Fluorouracil+Irinotecan+Leucovorin"
folfox_idx <- grep("folfox", tcga_response[["Drug_list"]], ignore.case = TRUE)
tcga_response[folfox_idx, "Drug_list"] <- "Fluorouracil+Leucovorin+Oxaliplatin"
tcga_drug_lists <- strsplit(tcga_response[["Drug_list"]], split = "\\+")
tcga_drug_lists <- lapply(tcga_drug_lists, tolower)
fnw <- paste0(ctrdb_path, "tcga_responses.csv")
readr::write_csv(tcga_response, fnw)

tcga_unique_drugs <- unique(Reduce("c", tcga_drug_lists))

tcga_drug_idxs <- lapply(tcga_drug_lists, match, table = tcga_unique_drugs)

n_drugs <- length(tcga_unique_drugs)
one_hot_from_idx <- function(idx) {
  out <- rep_len(0, n_drugs)
  out[idx] <- 1
  return(out)
}
tcga_drug_ind_mat <- Reduce("rbind", lapply(tcga_drug_idxs, one_hot_from_idx))

tcga_drug_ind_df <- data.frame(tcga_drug_ind_mat)
colnames(tcga_drug_ind_df) <- tcga_unique_drugs
rownames(tcga_drug_ind_df) <- tcga_response[["Sample_id"]]
tcga_drug_ind_df[["Source"]] <- tcga_response[["Source"]]
tcga_drug_ind_df[["Response"]] <- tcga_response[["Response"]]

UpSetR::upset(tcga_drug_ind_df, sets = tcga_unique_drugs, nsets = 50, cutoff = 5)

cancers <- unique(tcga_response[["Source"]])
for (i in cancers) {
  temp <- tcga_drug_ind_df[tcga_drug_ind_df[["Source"]] == i,]
  drugis <- names(which(apply(temp[,tcga_unique_drugs], 2, sum) > 0))
  tryCatch(
    {
      p1 <- UpSetR::upset(
        temp, 
        sets = drugis, 
        nsets = 50, 
        nintersects = NA, 
        queries = list(list(query = function(x) x[["Response"]] == "Non_response", color = "red", active = TRUE))
      )
      fn <- paste0(ctrdb_path, i, "_response_upset_plot.png")
      save_figure_safe(
        p1, 
        png, 
        fn, 
        width = plot_width * 1.5, 
        height = plot_height * 1.3, 
        res = plot_res, 
        units = plot_units
      )
    }, 
    error = function(e) {
      warning(e)
      p1 <- UpSetR::upset(
        temp, 
        sets = drugis, 
        nsets = 50, 
        nintersects = NA
      )
      fn <- paste0(ctrdb_path, i, "_response_upset_plot.png")
      save_figure_safe(
        p1, 
        png, 
        fn, 
        width = plot_width * 1.5, 
        height = plot_height * 1.3, 
        res = plot_res, 
        units = plot_units
      )
    }
  )
}
