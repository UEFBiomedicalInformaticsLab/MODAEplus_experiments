source("../setup.R")

fn <- paste0(data_root, "sc_dr/GSE114459_Polyak.csv.gz")
melanoma_mcf7_scdr <- readr::read_csv(fn)