script_root <- Sys.getenv("MODAE_SCRIPT_PATH")
if (script_root == "") {
  stop("Please define the MODAE_SCRIPT_PATH environment variable.")
}
source(paste0(script_root, "R/setup.R"))

# Updated multi-omics data (directly from CCLE)
# main omics
ccle_cnv <- read.csv(paste0(ccle_path, "CCLE_gene_cn.csv"), row.names = 1, header = TRUE)
ccle_wes_cnv <- read.csv(paste0(ccle_path, "CCLE_wes_gene_cn.csv"), row.names = 1, header = TRUE)
ccle_rna <- read.csv(paste0(ccle_path, "CCLE_expression.csv"), row.names = 1, header = TRUE)
#ccle_mut <- read.csv(paste0(ccle_path, "CCLE_mutations.csv"), row.names = NULL, header = TRUE)
#ccle_mut_genes <- plyr::ddply(ccle_mut, c("Hugo_Symbol", "DepMap_ID"), function(x) data.frame(damaging = sum(x[["Variant_annotation"]] == "damaging")))

ccle_rna_genes <- readLines(paste0(ccle_path, "CCLE_expression.csv"), n = 1L)
tcga_rna_example <- curatedTCGAData::curatedTCGAData("ACC", assays = "RNASeq2GeneNorm",
                                                   version = "2.0.1", dry.run = FALSE)
tcga_rna_genes <- rownames(SummarizedExperiment::rowData(MultiAssayExperiment::experiments(tcga_rna_example)[[1]]))
ccle_rna_symbols <- strsplit(gsub(" \\([0-9]+\\)", "", ccle_rna_genes), split = ",")[[1]][-1]
ccle_rna_entrez <- gsub(".*?\\(", "", strsplit(ccle_rna_genes, split = ",")[[1]][-1])
ccle_rna_entrez <- as.integer(gsub("\\)", "", ccle_rna_entrez))

symbols_intersect <- intersect(ccle_rna_symbols, tcga_rna_genes)
symbols_union <- union(ccle_rna_symbols, tcga_rna_genes)

table(CCLE = symbols_union %in% ccle_rna_symbols, TCGA = symbols_union %in% tcga_rna_genes)
# 16570 overlapping

ccle_rna_symbols[!ccle_rna_symbols %in% tcga_rna_genes][1:10]
tcga_rna_genes[!tcga_rna_genes %in% ccle_rna_symbols][1:10]

# Some of the missing genes seem to match an  older symbol
# So let's try matching through the biomaRt
mart <- biomaRt::useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
mart_attributes <- biomaRt::listAttributes(mart)
mart_results <- biomaRt::getBM(attributes = c("entrezgene_id", "hgnc_id", "hgnc_symbol"), mart = mart)

table(ccle_rna_entrez %in% mart_results$entrezgene_id) # Most match
table(tcga_rna_genes %in% mart_results$hgnc_symbol) # Most do not match

tcga_rna_symbols <- limma::alias2SymbolTable(tcga_rna_genes, species = "Hs")
tcga_rna_symbols_multi <- lapply(tcga_rna_genes, limma::alias2Symbol)
table(sapply(tcga_rna_symbols_multi, length))

ccle_rna_symbols_multi <- lapply(ccle_rna_symbols, limma::alias2Symbol)
table(sapply(ccle_rna_symbols_multi, length))

ccle_rna_symbols_updated <- limma::alias2SymbolTable(ccle_rna_symbols, species = "Hs")

symbols_union_updated <- union(ccle_rna_symbols_updated, tcga_rna_symbols)
symbols_union_updated <- symbols_union_updated[!is.na(symbols_union_updated)]

table(CCLE = symbols_union_updated %in% ccle_rna_symbols_updated, 
      TCGA = symbols_union_updated %in% tcga_rna_symbols)
# 18402 overlapping

# Save updated gene IDs in CCLE folder for matching
write.csv(data.frame(old_id = tcga_rna_genes, new_id = tcga_rna_symbols), 
          file = paste0(ccle_path, "TCGA_gene_mapping.csv"))
write.csv(data.frame(old_id = strsplit(ccle_rna_genes, split = ",")[[1]][-1], 
                     new_id = ccle_rna_symbols_updated), 
          file = paste0(ccle_path, "CCLE_gene_mapping.csv"))

# methylation
#readLines(paste0(ccle_path, "CCLE_RRBS_TSS_1kb_20180614.txt"), 1)
ccle_rrbs_tss <- read.table(paste0(ccle_path, "CCLE_RRBS_TSS_1kb_20180614.txt"), 
                            header = TRUE, row.names = NULL, sep = "\t")
ccle_rrbs_tss_arr <- lapply(ccle_rrbs_tss[,-(1:7)], as.numeric)
ccle_rrbs_tss_arr <- Reduce(cbind, ccle_rrbs_tss_arr)
rownames(ccle_rrbs_tss_arr) <- ccle_rrbs_tss[["TSS_id"]]
#hist(apply(is.na(ccle_rrbs_tss_arr), 2, sum))
#ccle_rrbs_tss[,-(1:7)] <- lapply(ccle_rrbs_tss[,-(1:7)], as.numeric)

#readLines(paste0(ccle_path, "CCLE_RRBS_TSS_CpG_clusters_20180614.txt"), 1)
ccle_rrbs_tss_clusters <- read.table(paste0(ccle_path, "CCLE_RRBS_TSS_CpG_clusters_20180614.txt"), 
                                     header = TRUE, row.names = NULL, sep = "\t")
ccle_rrbs_tss_clusters_arr <- lapply(ccle_rrbs_tss_clusters[,-(1:5)], as.numeric)
ccle_rrbs_tss_clusters_arr <- Reduce(cbind, ccle_rrbs_tss_clusters_arr)
#hist(apply(is.na(ccle_rrbs_tss_clusters_arr), 2, sum)) # NAs

# mirna
#temp <- readLines(paste0(ccle_path, "CCLE_miRNA_20180525.gct"), 5)
ccle_mirna <- read.table(paste0(ccle_path, "CCLE_miRNA_20180525.gct"), 
                         header = TRUE, row.names = NULL, sep = "\t", skip = 2)
ccle_mirna_arr <- as.matrix(ccle_mirna[,-(1:2)])
ccle_mirna_arr <- log2(ccle_mirna_arr + 1)
rownames(ccle_mirna_arr) <- ccle_mirna[["Name"]]

# RPPA
#readLines(paste0(ccle_path, "CCLE_RPPA_20180123.csv"), 2)
ccle_rppa <- read.csv(paste0(ccle_path, "CCLE_RPPA_20180123.csv"), row.names = 1, header = TRUE)
ccle_rppa <- as.matrix(ccle_rppa)
#readLines(paste0(ccle_path, "CCLE_RPPA_Ab_info_20180123.csv"), 2)
ccle_rppa_abs <- read.csv(paste0(ccle_path, "CCLE_RPPA_Ab_info_20180123.csv"), row.names = NULL, header = TRUE)

# sample info
#readLines(paste0(ccle_path, "CCLE_sample_info_file_2012-10-18.txt"), 5)
ccle_info <- read.table(paste0(ccle_path, "CCLE_sample_info_file_2012-10-18.txt"), 
                        header = TRUE, row.names = NULL, sep = "\t")
#readLines(paste0(ccle_path, "OmicsProfiles.csv"), 2)
ccle_profiles <- read.csv(paste0(ccle_path, "OmicsProfiles.csv"), row.names = NULL, header = TRUE)
#readLines(paste0(ccle_path, "OmicsDefaultModelProfiles.csv"), 2)
ccle_model_profiles <- read.csv(paste0(ccle_path, "OmicsDefaultModelProfiles.csv"), row.names = NULL, header = TRUE)
#readLines(paste0(ccle_path, "Model.csv"), 2)
ccle_model <- read.csv(paste0(ccle_path, "Model.csv"), row.names = NULL, header = TRUE)
#readLines(paste0(ccle_path, "ModelCondition.csv"), 2)
ccle_modelcon <- read.csv(paste0(ccle_path, "ModelCondition.csv"), row.names = NULL, header = TRUE)
#readLines(paste0(ccle_path, "Media.csv"), 2)
ccle_media <- read.csv(paste0(ccle_path, "Media.csv"), row.names = NULL, header = TRUE)

# Map all ids into ModelIDs
table(colnames(ccle_rrbs_tss)[-(1:7)] %in% ccle_model$CCLEName)
fixed_colnames <- gsub("^X", "", colnames(ccle_rrbs_tss)[-(1:7)])
ccle_rrbs_tss_matching_cols <- fixed_colnames %in% ccle_model$CCLEName
ccle_rrbs_tss_ids <- fixed_colnames[ccle_rrbs_tss_matching_cols]
fixed_colnames[!ccle_rrbs_tss_matching_cols] # print non-matching cell lines
ccle_rrbs_tss_arr <- ccle_rrbs_tss_arr[,ccle_rrbs_tss_matching_cols]
ccle_rrbs_tss_model_ids <- ccle_model$ModelID[match(ccle_rrbs_tss_ids, ccle_model$CCLEName)]

table(colnames(ccle_mirna_arr) %in% ccle_model$CCLEName)
fixed_colnames <- gsub("^X", "", colnames(ccle_mirna_arr))
ccle_mirna_matching_cols <- fixed_colnames %in% ccle_model$CCLEName
ccle_mirna_ids <- fixed_colnames[ccle_mirna_matching_cols]
fixed_colnames[!ccle_mirna_matching_cols]
ccle_mirna_arr <- ccle_mirna_arr[,ccle_mirna_matching_cols]
ccle_mirna_model_ids <- ccle_model$ModelID[match(ccle_mirna_ids, ccle_model$CCLEName)]

table(rownames(ccle_rppa) %in% ccle_model$CCLEName)
ccle_rppa_matching_cols <- rownames(ccle_rppa) %in% ccle_model$CCLEName
ccle_rppa_ids <- rownames(ccle_rppa)[ccle_rppa_matching_cols]
rownames(ccle_rppa)[!ccle_rppa_matching_cols]
ccle_rppa <- ccle_rppa[ccle_rppa_matching_cols,]
ccle_rppa_model_ids <- ccle_model$ModelID[match(ccle_rppa_ids, ccle_model$CCLEName)]

intersect_ids <- Reduce(intersect, list(rownames(ccle_cnv), rownames(ccle_rna), 
                                        ccle_rrbs_tss_model_ids, 
                                        ccle_mirna_model_ids, 
                                        ccle_rppa_model_ids))

type_counts <- table(ccle_model[match(intersect_ids, ccle_model$ModelID), "OncotreePrimaryDisease"])
breast_counts <- type_counts[grep("breast", names(type_counts), ignore.case = TRUE)]
breast_cancer_ids <- ccle_model$ModelID[grep("breast", ccle_model$OncotreePrimaryDisease, ignore.case = TRUE)]

# Create subset of only breast cancer cell line data
omic_names <- c("cnv", "rna", "rrbs", "mirna", "rppa")
data_sets <- list(as.matrix(ccle_cnv), 
                  as.matrix(ccle_rna), 
                  t(ccle_rrbs_tss_arr), 
                  t(ccle_mirna_arr), 
                  ccle_rppa)
names(data_sets) <- omic_names
sample_ids <- list(rownames(ccle_cnv), rownames(ccle_rna), ccle_rrbs_tss_model_ids, 
                   ccle_mirna_model_ids, ccle_rppa_model_ids)
names(sample_ids) <- omic_names
for (i in omic_names) {
  rownames(data_sets[[i]]) <- sample_ids[[i]]
}

for (i in names(data_sets)) {
  int_ids <- intersect(breast_cancer_ids, rownames(data_sets[[i]]))
  data_sets[[i]] <- data_sets[[i]][int_ids,]
}

rrbs_imputed <- impute::impute.knn(t(data_sets$rrbs), k = 5)

data_sets$rrbs <- t(rrbs_imputed$data)

saveRDS(data_sets, file = paste0(ccle_path, "processed_brca_cl_data.rds"))
## nearZeroVariance features

ccle_cnv_sd <- apply(ccle_cnv, 2, sd)
ccle_rna_sd <- apply(ccle_rna, 2, sd)
ccle_rrbs_tss_sd <- apply(ccle_rrbs_tss_arr, 1, sd)
ccle_mirna_sd <- apply(ccle_mirna_arr, 1, sd)
ccle_rppa_sd <- apply(ccle_rppa, 1, sd)

omic_names <- c("cnv", "rna", "rrbs", "mirna", "rrpa")
top_var_numbers <- c(200, 2000, 300, 100, 150)
names(top_var_numbers) <- omic_names

data_sets <- list(ccle_cnv, ccle_rna, t(ccle_rrbs_tss_arr), t(ccle_mirna_arr), 
                  ccle_rppa)
names(data_sets) <- omic_names

sample_ids <- list(rownames(ccle_cnv), rownames(ccle_rna), ccle_rrbs_tss_model_ids, 
                   ccle_mirna_model_ids, ccle_rppa_model_ids)
names(sample_ids) <- omic_names
for (i in omic_names) {
  rownames(data_sets[[i]]) <- sample_ids[[i]]
}

ccle_pairs_spls_tuned <- list()
ccle_pairs_spls_final <- list()
ccle_pairs_spls_final_perf <- list()

require(mixOmics)
omic_names_iterate <- c("cnv", "rna", "mirna", "rrpa")
for (i in omic_names_iterate[-length(omic_names_iterate)]) {
  remaining_omics <- omic_names_iterate[(which(omic_names_iterate == i)+1):length(omic_names_iterate)]
  for (j in remaining_omics) {
    intersection_ij <- intersect(sample_ids[[i]], sample_ids[[j]])
    intersection_ij <- intersect(intersection_ij, breast_cancer_ids)
    X <- data_sets[[i]][intersection_ij,]
    Y <- data_sets[[j]][intersection_ij,]
    list.keepX <- seq(top_var_numbers[i]/20, top_var_numbers[i]/2, top_var_numbers[i]/20)
    list.keepY <- seq(top_var_numbers[j]/20, top_var_numbers[j]/2, top_var_numbers[j]/20)
    X_var <- apply(X, 2, var)
    Y_var <- apply(Y, 2, var)
    X_top_var <- (length(X_var) - rank(X_var, ties = "first")) < (top_var_numbers[i])
    Y_top_var <- (length(Y_var) - rank(Y_var, ties = "first")) < (top_var_numbers[j])
    
    ccle_spls_tuned <- tune.spls(X = X[,X_top_var], 
                                 Y = Y[,Y_top_var], 
                                 ncomp = 2,
                                 test.keepX = list.keepX,
                                 test.keepY = list.keepY,
                                 nrepeat = 1, folds = 10, # use 10 folds
                                 mode = 'canonical', measure = 'cor') 
    optimal.keepX <- ccle_spls_tuned$choice.keepX
    optimal.keepY <- ccle_spls_tuned$choice.keepY
    optimal.ncomp <-  length(optimal.keepX)
    ccle_spls_final <- spls(X[,X_top_var], 
                            Y[,Y_top_var], 
                            ncomp = optimal.ncomp, 
                            keepX = optimal.keepX,
                            keepY = optimal.keepY,
                            mode = "canonical")
    ccle_spls_final_perf <- perf(ccle_spls_final, 
                                 folds = 5, nrepeat = 10, # use repeated cross-validation
                                 validation = "Mfold", 
                                 dist = "max.dist",  # use max.dist measure
                                 progressBar = FALSE)
    
    ccle_pairs_spls_tuned[[paste0(i, "_", j)]] <- ccle_spls_tuned
    ccle_pairs_spls_final[[paste0(i, "_", j)]] <- ccle_spls_final
    ccle_pairs_spls_final_perf[[paste0(i, "_", j)]] <- ccle_spls_final_perf
  }
}

# plot the stability of each feature for the first two components, 
# 'h' type refers to histogram
par(mfrow=c(1,2)) 
plot(ccle_pairs_spls_final_perf[["cnv_rna"]]$features$stability.X$comp1, 
     type = 'h',
     ylab = 'Stability',
     xlab = 'Features',
     main = '(a) Comp 1', las =2)
plot(ccle_pairs_spls_final_perf[["cnv_rna"]]$features$stability.X$comp2, 
     type = 'h',
     ylab = 'Stability',
     xlab = 'Features',
     main = '(b) Comp 2', las =2)

stability_df_x_list <- lapply(ccle_pairs_spls_final_perf, function(x) x$features$stability.X)
stability_df_y_list <- lapply(ccle_pairs_spls_final_perf, function(x) x$features$stability.Y)

stability_df_x <- list()
stability_df_y <- list()
for (i in names(ccle_pairs_spls_final_perf)) {
  stability_df_x[[i]] <- list()
  for (j in names(stability_df_x_list[[i]])) {
    stability_df_x[[i]][[j]] <- data.frame(omic_pair = i, variate = "x", comp = j, 
                                      feature = names(stability_df_x_list[[i]][[j]]), 
                                      freq = stability_df_x_list[[i]][[j]])
  }
  stability_df_x[[i]] <- Reduce(rbind, stability_df_x[[i]])
  stability_df_y[[i]] <- list()
  for (j in names(stability_df_y_list[[i]])) {
    stability_df_y[[i]][[j]] <- data.frame(omic_pair = i, variate = "y", comp = j, 
                                      feature = names(stability_df_y_list[[i]][[j]]), 
                                      freq = stability_df_y_list[[i]][[j]])
  }
  stability_df_y[[i]] <- Reduce(rbind, stability_df_y[[i]])
}
stability_df_x <- Reduce(rbind, stability_df_x)
stability_df_y <- Reduce(rbind, stability_df_y)
stability_df <- rbind(stability_df_x, stability_df_y)

require(ggplot2)
ggplot(stability_df, aes(x = reorder(feature, -freq), y = freq)) + 
  geom_bar(stat = "identity") + theme_bw() + theme(axis.text.x=element_blank()) + 
  ggh4x::facet_grid2(variate + comp ~ omic_pair, scales = "free_x", independent = "x")

# Plotting
require(mixOmics)
# mRNA miRNA
plotIndiv(ccle_pairs_spls_final[["rna_mirna"]], 
          ind.names = FALSE, 
          rep.space = "X-variate")

plotIndiv(ccle_pairs_spls_final[["rna_mirna"]], 
          ind.names = FALSE,
          rep.space = "Y-variate")

plotIndiv(ccle_pairs_spls_final[["rna_mirna"]], 
          ind.names = FALSE,
          rep.space = "XY-variate")

plotVar(ccle_pairs_spls_final[["rna_mirna"]], 
        var.names = c(FALSE, TRUE),
        cex = c(4, 4), cutoff = 0.5,
        title = "sPLS CCLE breast cancer comp 1-2")

# mRNA RPPA
plotIndiv(ccle_pairs_spls_final[["rna_rrpa"]], 
          ind.names = FALSE, 
          rep.space = "X-variate")

plotIndiv(ccle_pairs_spls_final[["rna_rrpa"]], 
          ind.names = FALSE,
          rep.space = "Y-variate")

plotIndiv(ccle_pairs_spls_final[["rna_rrpa"]], 
          ind.names = FALSE,
          rep.space = "XY-variate")

plotVar(ccle_pairs_spls_final[["rna_rrpa"]], 
        var.names = c(FALSE, TRUE),
        cex = c(4, 4), cutoff = 0.5,
        title = "sPLS CCLE breast cancer comp 1-2")

# mmiRNA RPPA
plotIndiv(ccle_pairs_spls_final[["mirna_rrpa"]], 
          ind.names = FALSE, 
          rep.space = "X-variate")

plotIndiv(ccle_pairs_spls_final[["mirna_rrpa"]], 
          ind.names = FALSE,
          rep.space = "Y-variate")

plotIndiv(ccle_pairs_spls_final[["mirna_rrpa"]], 
          ind.names = FALSE,
          rep.space = "XY-variate")

plotVar(ccle_pairs_spls_final[["mirna_rrpa"]], 
        var.names = c(TRUE, TRUE),
        cex = c(4, 4), cutoff = 0.5,
        title = "sPLS CCLE breast cancer comp 1-2")
ccle_mirna[ccle_mirna$Name == "nmir00344.2", 1:2]
ccle_rppa_abs[ccle_rppa_abs$Antibody_Name == "PKC.alpha",]

var_cor <- cor(X, CV.rcc.nutrimouse$variates$X, use = "pairwise")
## mixOmics
# CNV - RNA
# loop this 
#ccle_cnv_nzv <- caret::nearZeroVar(ccle_cnv[cnv_rna_intersection,])
#ccle_rna_nzv <- caret::nearZeroVar(ccle_rna[cnv_rna_intersection,])
cnv_rna_intersection <- Reduce(intersect, list(breast_cancer_ids, 
                                               rownames(ccle_cnv), 
                                               rownames(ccle_rna)))
X <- ccle_cnv[cnv_rna_intersection,]
Y <- ccle_rna[cnv_rna_intersection,]

if (length(ccle_cnv_nzv) > 0) X <- X[,-ccle_cnv_nzv]
if (length(ccle_rna_nzv) > 0) Y <- Y[,-ccle_rna_nzv]
ccle_cnv_rna_spls <- mixOmics::spls(X = X,
                                    Y = Y, 
                                    ncomp = 2, mode = 'canonical')

mixOmics::plotIndiv(ccle_spls, ind.names = FALSE, 
                    rep.space = "X-variate", # plot in X-variate subspace
                    #group = sapply(strsplit(total_intersection, split = "_"), function(x) x[2]), # colour by group
                    #col.per.group = mixOmics::color.mixo(1:2), 
                    #legend = TRUE, legend.title = 'Group',mix
                    title = "CNV-RNA canonical variates")

# set range of test values for number of variables to use from X dataframe
list.keepX <- seq(10, 100, 10)
# set range of test values for number of variables to use from Y dataframe
list.keepY <- seq(100, 1000, 100)

X_var <- apply(X, 2, var)
Y_var <- apply(Y, 2, var)

X_top_var <- (length(X_var) - rank(X_var)) < 200
Y_top_var <- (length(Y_var) - rank(Y_var)) < 2000

ccle_cnv_rna_spls_tuned <- mixOmics::tune.spls(X = X[,X_top_var], 
                                               Y = Y[,Y_top_var], 
                                               ncomp = 2,
                                               test.keepX = list.keepX,
                                               test.keepY = list.keepY,
                                               nrepeat = 1, folds = 10, # use 10 folds
                                               mode = 'canonical', measure = 'cor') 
plot(ccle_cnv_rna_spls_tuned)         # use the correlation measure for tuning

# extract optimal number of variables for X dataframe
optimal.keepX <- ccle_cnv_rna_spls_tuned$choice.keepX 

# extract optimal number of variables for Y datafram
optimal.keepY <- ccle_cnv_rna_spls_tuned$choice.keepY

optimal.ncomp <-  length(optimal.keepX) # extract optimal number of components

# use all tuned values from above
ccle_cnv_rna_spls_final <- mixOmics::spls(X, Y, ncomp = optimal.ncomp, 
                                          keepX = optimal.keepX,
                                          keepY = optimal.keepY,
                                          mode = "canonical")

plotIndiv(ccle_cnv_rna_spls_final, ind.names = FALSE, 
          rep.space = "X-variate")

plotIndiv(ccle_cnv_rna_spls_final, ind.names = FALSE,
          rep.space = "Y-variate")

plotIndiv(ccle_cnv_rna_spls_final, ind.names = FALSE,
          rep.space = "XY-variate")

# PharmacoGx drug response data and old multi-omic data
psets <- PharmacoGx::availablePSets()

GDSC1 <- PharmacoGx::downloadPSet("GDSC_2020(v1-8.2)", saveDir = gdsc_path)
GDSC2 <- PharmacoGx::downloadPSet("GDSC_2020(v2-8.2)", saveDir = gdsc_path)

CCLE <- PharmacoGx::downloadPSet("CCLE_2015", saveDir = ccle_path)
PRISM <- PharmacoGx::downloadPSet("PRISM_2020", saveDir = ccle_path)

# PharmacoGx multi-omics data (out of date)
CCLE <- readRDS("~/pharmacogx/CCLE_2015.rds")

PharmacoGx::mDataNames(GDSC1)
PharmacoGx::mDataNames(GDSC1)
PharmacoGx::mDataNames(CCLE)
PharmacoGx::mDataNames(PRISM)

CCLE_mut <- PharmacoGx::molecularProfiles(CCLE, mDataType = "mutation")
CCLE_cnv <- PharmacoGx::molecularProfiles(CCLE, mDataType = "cnv")
CCLE_rna <- PharmacoGx::molecularProfiles(CCLE, mDataType = "rna")

CCLE_cell_info <- PharmacoGx::cellInfo(CCLE)

#cnv_ids <- colnames(CCLE_cnv)[colnames(CCLE_cnv) %in% CCLE_cell_info[["SNP.arrays"]]]
#rna_ids <- colnames(CCLE_rna)[colnames(CCLE_rna) %in% CCLE_cell_info[["Expression.arrays"]]]

#cell_line_intersection <- CCLE_cell_info[["cellid"]][!is.na(CCLE_cell_info[["Expression.arrays"]]) & # RNA id
#                                                     is.na(CCLE_cell_info[["SNP.arrays"]] %in% cnv_ids)] # CNV id

cnv_pheno <- PharmacoGx::phenoInfo(CCLE, "cnv")
rna_pheno <- PharmacoGx::phenoInfo(CCLE, "rna")

table(rownames(cnv_pheno) == colnames(CCLE_cnv))
table(rownames(rna_pheno) == colnames(CCLE_rna))

cell_line_intersection <- intersect(cnv_pheno[["Cell.line.primary.name"]], 
                                    rna_pheno[["Cell.line.primary.name"]])

colnames(CCLE_cnv) <- cnv_pheno[["Cell.line.primary.name"]]
colnames(CCLE_rna) <- rna_pheno[["Cell.line.primary.name"]]

ccle_spls <- mixOmics::spls(X = t(CCLE_cnv[,cell_line_intersection]),
                            Y = t(CCLE_rna[,cell_line_intersection]), 
                            ncomp = 2, mode = 'regression')

ccle_tissues <- CCLE_cell_info[cell_line_intersection,"tissueid"]

mixOmics::plotIndiv(ccle_spls, ind.names = FALSE, 
                    rep.space = "X-variate", # plot in X-variate subspace
                    #group = sapply(strsplit(total_intersection, split = "_"), function(x) x[2]), # colour by group
                    #col.per.group = mixOmics::color.mixo(1:2), 
                    #legend = TRUE, legend.title = 'Group',
                    title = "CNV effect on RNA expression")

# Better filters
ccle_brca <- readRDS(paste0(ccle_path, "processed_brca_cl_data.rds"))
lapply(ccle_brca, dim)

sample_ids <- lapply(ccle_brca, rownames)
ccle_brca_var <- lapply(ccle_brca, function(x) apply(x, 2, var))
ccle_zv <- lapply(ccle_brca_var, function(x) x == 0)

ccle_brca_cross_cor_filter <- list()
for (i in 1:(length(ccle_brca)-1)) {
  for (j in (i+1):length(ccle_brca)) {
    pairname <- paste0(names(ccle_brca)[i], "_", names(ccle_brca)[j])
    ccle_brca_cross_cor_filter[[pairname]] <- list()
    intersection_ij <- intersect(sample_ids[[i]], sample_ids[[j]])
    X <- ccle_brca[[i]][intersection_ij,]
    Y <- ccle_brca[[j]][intersection_ij,]
    
    X_zv <- apply(X, 2, var) == 0
    Y_zv <- apply(Y, 2, var) == 0
    
    cross_cor <- cor(X[,!X_zv], Y[,!Y_zv], method = "spearman")
    #cross_cor <- ccle_brca_cross_cor[[pairname]]
    ccle_brca_cross_cor_filter[[pairname]]$X <- names(which(apply(abs(cross_cor) > 0.80, 1, any)))
    ccle_brca_cross_cor_filter[[pairname]]$Y <- names(which(apply(abs(cross_cor) > 0.80, 2, any)))
  }
}
saveRDS(ccle_brca_cross_cor_filter, file = paste0(ccle_path, "brca_cl_data_cc_filters.rds"))

