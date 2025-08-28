source("../setup.R")

# MAF from https://doi.org/10.1016/j.cels.2018.03.002
fn <- paste0(patient_expression_root_dir, "maf_full.csv.gz")
if (file.exists(fn)) {
  maf <- readr::read_csv(fn)
} else {
  maf <- TCGAbiolinks::getMC3MAF()
  readr::write_tsv(maf, gzfile(fn))
}

maftools_cols <- c(
  "Hugo_Symbol", 
  "Chromosome", 
  "Start_Position", 
  "End_Position", 
  "Reference_Allele", 
  "Tumor_Seq_Allele2", 
  "Variant_Classification", 
  "Variant_Type", 
  "Tumor_Sample_Barcode"
)

if (any(!maftools_cols %in% colnames(maf))) {
  stop(
    paste(
      "Missing mandatory MAF columns:", 
      paste(
        maftools_cols[!maftools_cols %in% colnames(maf)], 
        collapse = ", "
      )
    )
  )
}

maf <- maftools::read.maf(
  maf = maf, 
  
)



table(maf[, "Consequence"])
table(maf[, "CLIN_SIG"])
table(maf[, "IMPACT"])
table(maf[, "VARIANT_CLASS"])
maf[["all_effects"]][1:10]
maf[["SIFT"]][1:10]
maf[["PolyPhen"]][1:10]
table(maf[, "PHENO"])
table(maf[, "PolyPhen"])


clin_sig_relevant <- grep("pathogenic|drug_response", maf[["CLIN_SIG"]])
impact_high <- which(maf[["IMPACT"]] %in% c("HIGH", "MODERATE"))

f_union <- union(clin_sig_relevant, impact_high)

table(
  clin_sig = f_union %in% clin_sig_relevant, 
  impact = f_union %in% impact_high
)
filter_idx <- intersect(clin_sig_relevant, impact_high)

table(maf[filter_idx, "project_id"])

mut_genes <- c("KEAP1", "NFE2L2", "CUL3", "STK11", "SMARCA4")

gene_presel_idx <- which(maf[["SYMBOL"]] %in% mut_genes)

table(maf[gene_presel_idx, "SYMBOL"])


