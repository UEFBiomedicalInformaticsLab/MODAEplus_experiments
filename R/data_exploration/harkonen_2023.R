script_root <- Sys.getenv("MODAE_SCRIPT_PATH")
if (script_root == "") {
  stop("Please define the MODAE_SCRIPT_PATH environment variable.")
}
source(paste0(script_root, "R/setup.R"))

fn <- paste0(harkonen_path, "1-s2.0-S2213231723000459-mmc3.xlsx")
nrf2_muts <- readxl::read_xlsx(fn, sheet = 1)
fn <- paste0(harkonen_path, "1-s2.0-S2213231723000459-mmc2.xlsx")
nrf2_targets <- readxl::read_xlsx(fn, sheet = 1)
fn <- paste0(harkonen_path, "1-s2.0-S2213231723000459-mmc1.xlsx")
nrf2_dea <- readxl::read_xlsx(fn, sheet = 1)

table(nrf2_muts[["Mutation"]])
table(LFC = abs(nrf2_dea$logFC) > log2(2), p = nrf2_dea$adj.P.Val < 0.01)