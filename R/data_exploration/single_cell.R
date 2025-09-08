script_root <- Sys.getenv("MODAE_SCRIPT_PATH")
if (script_root == "") {
  stop("Please define the MODAE_SCRIPT_PATH environment variable.")
}
source(paste0(script_root, "R/setup.R"))

fn <- paste0(data_root, "sc_dr/GSE114459_Polyak.csv.gz")
melanoma_mcf7_scdr <- readr::read_csv(fn)