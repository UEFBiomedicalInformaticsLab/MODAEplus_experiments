int_path <- paste0(data_root, "scanb_preprocessed/")

scanb_y_fn <- paste0(int_path, "scanb_pheno.csv.gz")
scanb_y <- read.csv(scanb_y_fn, row.names = 1, header = TRUE)

scanb_normalized <- read.csv(
  paste0(int_path, "scanb_tpm_uq.csv.gz"), 
  row.names = 1, 
  header = TRUE
)

scanb_pcs <- FactoMineR::PCA(
  t(log2(scanb_normalized+1)),
  scale.unit = FALSE,
  ncp = 50,
  graph = FALSE
)

write.csv(
  scanb_pcs$ind$coord, 
  file = gzfile(paste0(int_path, "scanb_pca50.csv.gz"))
)

write.csv(
  scanb_pcs$var$contrib, 
  file = gzfile(paste0(int_path, "scanb_pca50_loadings.csv.gz"))
)

scanb_nmf <- NMF::nmf(
  t(log2(scanb_normalized+1)),
  rank = 2:10, 
  nrun = 30, 
  method = "brunet", 
  seed = "nndsvd"
)

set.seed(0)
scanb_umap <- uwot::umap(
  scanb_pcs$ind$coord, 
  n_neighbors = 50, 
  n_components = 2, 
  verbose = FALSE, 
  init = "normlaplacian"
)
colnames(scanb_umap) <- paste0("Z", 1:2)

write.csv(
  scanb_umap, 
  file = gzfile(paste0(int_path, "scanb_umap.csv.gz"))
)


COPS::pca_viz(
  t(log2(scanb_normalized+1)), 
  category = scanb_y$NCN.PAM50, 
  category_label = "PAM50"
)

n_unique <- sapply(scanb_y, function(x) length(levels(factor(x))))


scanb_assoc <- embedding_associations(
  log2(scanb_normalized+1), 
  class = scanb_y[,n_unique <= 20 & n_unique > 1],
  n_pc_max = 10
)

# TODO: check if there is anything useful in the code below
scanb_corrected <- sva::ComBat_seq(
  fixed_matrix, 
  batch = scanb_y[["LibraryProtocol"]]
)

dge <- edgeR::DGEList(counts=scanb_corrected)
dge <- edgeR::calcNormFactors(dge, method = "upperquartile")
scanb_normalized <- edgeR::cpm(dge, normalized.lib.sizes = TRUE)

COPS::pca_viz(t(log2(scanb_normalized+1)), category = scanb_y$NCN.PAM50, "PAM50")
plots <- COPS::triple_viz(t(log2(scanb_normalized+1)), category = scanb_y$NCN.PAM50, "PAM50", tsne = FALSE)


dat <- log2(scanb_normalized+1)
class <- scanb_y[,n_unique <= 600 & n_unique > 1]
n_pc_max <- 20

if (is.null(dim(class))) class <- data.frame(class = as.character(class))
class <- lapply(class, as.factor)
class <- class[sapply(class, function(x) length(levels(x)))>1]

pca_silh <- c()
pca_reg <- list()
DSC_res <- c()

dat_pca <- FactoMineR::PCA(t(dat),
                           scale.unit = FALSE,
                           ncp = min(c(n_pc_max, dim(dat))),
                           graph = FALSE)

for (i in 1:length(class)) {
  # PCA based
  if (FALSE) {
    pca_silh[i] <- kBET::batch_sil(list(x = dat_pca$ind$coord[!is.na(class[[i]]),]),
                                   class[[i]][!is.na(class[[i]])],
                                   nPCs = min(c(n_pc_max, dim(dat))))
    pca_reg[[i]] <- kBET::pcRegression(list(x = dat_pca$ind$coord[!is.na(class[[i]]),], 
                                            sdev = sqrt(dat_pca$eig[,"eigenvalue"])),
                                       class[[i]][!is.na(class[[i]])],
                                       n_top = min(c(n_pc_max, dim(dat))))
  }
  # Other
  DSC_res[i] <- COPS::DSC(dat[,!is.na(class[[i]])], class[[i]][!is.na(class[[i]])])
}

dat_comb <- cbind(data.frame(dat_pca$ind$coord), class)
swamp::confounding(dat_comb, note = FALSE)

PC_R2_max = sapply(pca_reg, function(a) a[["maxR2"]])

data.frame(var = names(class)[order(DSC_res)], DSC = sort(DSC_res))
data.frame(var = names(class)[order(pca_silh)], pca_silh = sort(pca_silh))
data.frame(var = names(class)[order(PC_R2_max)], PC_R2_max = sort(PC_R2_max))

set.seed(0)
res_umap <- uwot::umap(dat_pca$ind$coord, 
                       n_neighbors = 20, 
                       n_components = 2, 
                       pca = NULL, 
                       verbose = FALSE, 
                       init = "normlaplacian")

res_umap <- data.frame(Dim.1 = res_umap[,1], Dim.2 = res_umap[,2])

class_name <- "NCN.PAM50"
res_umap$category <- class[[class_name]]
ggplot(res_umap, aes(Dim.1, Dim.2, color = category)) + geom_point(shape = "+", size = 3) + 
  theme_bw() + scale_color_brewer(palette = "Dark2") + 
  labs(x = "Z1", y = "Z2", color = class_name)

class_name <- "SequencerSerial"
res_umap$category <- class[[class_name]]
ggplot(res_umap, aes(Dim.1, Dim.2, color = category)) + geom_point(shape = "+", size = 3) + 
  theme_bw() + scale_color_brewer(palette = "Dark2") + 
  labs(x = "Z1", y = "Z2", color = class_name)
ggplot(res_umap, aes(Dim.1, Dim.2, color = category == "SN587")) + geom_point(shape = "+", size = 3) + 
  theme_bw() + scale_color_brewer(palette = "Dark2") + 
  labs(x = "Z1", y = "Z2", color = class_name)

res_umap$category <- scanb_y[["PoolName"]]
ggplot(res_umap, aes(Dim.1, Dim.2, color = category)) + geom_point(shape = "+", size = 3) + 
  theme_bw() + #scale_color_brewer(palette = "Dark2") + 
  labs(x = "Z1", y = "Z2") + theme(legend.position="none")

res_umap[["outlier"]] <- res_umap$Dim.1 < 0 & res_umap$Dim.2 > 2.5
scanb_y[["outlier"]] <- res_umap[["outlier"]]
outlier_pools <- table(scanb_y[["outlier"]], scanb_y[["PoolName"]])
outlier_pool_f <- sweep(outlier_pools, 2, apply(outlier_pools, 2, sum), .Primitive("/"))
sum(outlier_pools[,outlier_pool_f[2,] > 0.])
sum(outlier_pools[,outlier_pool_f[2,] == 0.])

outlier_seqs <- table(scanb_y[["outlier"]], scanb_y[["SequencerSerial"]])
outlier_seq_f <- sweep(outlier_seqs, 2, apply(outlier_seqs, 2, sum), .Primitive("/"))
outlier_seqs[,outlier_seq_f[2,] > 0.1]
outlier_seqs[,outlier_seq_f[2,] == 0.]

table(apply(table(scanb_y[["PoolName"]], scanb_y[["SequencerSerial"]]) > 0, 1, sum))

dat_pca_df <- data.frame(dat_pca$ind$coord)
dat_pca_df[["ID"]] <- rownames(dat_pca_df)
#dat_pca_df[["outlier"]] <- scanb_y[["outlier"]]
dat_pca_df_long <- reshape2::melt(dat_pca_df, id.vars = "ID", variable.name = "PC", value.name = "value")

table(scanb_y[["GEX.assay"]] == rownames(dat_pca_df))
scanb_y <- as.data.frame(scanb_y)
rownames(scanb_y) <- scanb_y[["GEX.assay"]]
dat_pca_df_long[["outlier"]] <- scanb_y[dat_pca_df_long[["ID"]], "outlier"]

ggplot(dat_pca_df_long, aes(x = value, fill = outlier)) + geom_density(alpha = 0.4) + 
  theme_bw() + scale_fill_brewer(palette = "Dark2") + facet_wrap(PC ~.)

dat_pca_df_wide <- dat_pca_df
dat_pca_df_wide[["outlier"]] <- scanb_y[dat_pca_df_wide[["ID"]], "outlier"]
GGally::ggpairs(dat_pca_df_wide, aes(color = outlier, shape = "+"), columns = 1:20,
                diag = list(continuous = GGally::wrap("densityDiag", alpha = 0.5)),
                upper = list(continuous = "blank"), 
                lower = list(continuous = GGally::wrap("points", shape = "+", size = 3))) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + theme_bw()




dat_pca_df[["outlier"]] <- scanb_y[dat_pca_df[["ID"]], "outlier"]
dat_pca_df[["outlier2"]] <- dat_pca_df[["Dim.3"]] < -25
table(umap = dat_pca_df[["outlier"]], pc3 = dat_pca_df[["outlier2"]])
outlier2_pools <- table(dat_pca_df[["outlier2"]], scanb_y[["PoolName"]])
outlier2_pool_f <- sweep(outlier2_pools, 2, apply(outlier2_pools, 2, sum), .Primitive("/"))
sum(outlier2_pools[,outlier2_pool_f[2,] > 0.])
sum(outlier2_pools[,outlier2_pool_f[2,] == 0.])

class_filtered <- scanb_y[scanb_y[["PoolName"]] %in% colnames(outlier_pools[,outlier_pool_f[2,] == 0.]),]
dat_pca_filtered <- FactoMineR::PCA(
  t(dat)[scanb_y[["PoolName"]] %in% colnames(outlier_pools[,outlier_pool_f[2,] == 0.]),],
  scale.unit = FALSE,
  ncp = 20,
  graph = FALSE)
set.seed(0)
res_umap_filtered <- uwot::umap(
  dat_pca_filtered$ind$coord, 
  n_neighbors = 20, 
  n_components = 2, 
  pca = NULL, 
  verbose = FALSE, 
  init = "normlaplacian")
res_umap_filtered <- data.frame(Dim.1 = res_umap_filtered[,1], Dim.2 = res_umap_filtered[,2])

class_name <- "SequencerSerial"
res_umap_filtered$category <- class_filtered[[class_name]]
ggplot(res_umap_filtered, aes(Dim.1, Dim.2, color = category)) + geom_point(shape = "+", size = 3) + 
  theme_bw() + scale_color_brewer(palette = "Dark2") + 
  labs(x = "Z1", y = "Z2", color = class_name)

res_umap_filtered$category <- class_filtered[["PoolName"]]
ggplot(res_umap_filtered, aes(Dim.1, Dim.2, color = category)) + geom_point(shape = "+", size = 3) + 
  theme_bw() + #scale_color_brewer(palette = "Dark2") + 
  labs(x = "Z1", y = "Z2") + theme(legend.position="none")

res_umap_filtered$category <- class_filtered[["NCN.PAM50"]]
ggplot(res_umap_filtered, aes(Dim.1, Dim.2, color = category)) + geom_point(shape = "+", size = 3) + 
  theme_bw() + scale_color_brewer(palette = "Dark2") + 
  labs(x = "Z1", y = "Z2")



class_name <- "CClust"
res_umap$category <- class[[class_name]]
ggplot(res_umap, aes(Dim.1, Dim.2, color = category)) + geom_point(shape = "+", size = 3) + 
  theme_bw() + scale_color_brewer(palette = "Dark2") + 
  labs(x = "Z1", y = "Z2", color = class_name)

# Fix both LibraryProtocol and SequencerSerial
ss_filter <- scanb_y[["SequencerSerial"]] %in% c("NB501986", "NS500281", "SN587")
scanb_batch <- interaction(scanb_y[ss_filter,"LibraryProtocol"][[1]], 
                           scanb_y[ss_filter,"SequencerSerial"][[1]])
scanb_corrected <- sva::ComBat_seq(fixed_matrix[,ss_filter], batch = droplevels(scanb_batch))

dge <- edgeR::DGEList(counts=scanb_corrected)
dge <- edgeR::calcNormFactors(dge, method = "upperquartile")
scanb_normalized <- edgeR::cpm(dge, normalized.lib.sizes = TRUE)

# run batch-association analysis again
dat <- log2(scanb_normalized+1)
class <- scanb_y[ss_filter,n_unique <= 600 & n_unique > 1]
n_pc_max <- 20

if (is.null(dim(class))) class <- data.frame(class = as.character(class))
class <- lapply(class, as.factor)
class <- class[sapply(class, function(x) length(levels(x)))>1]

pca_silh <- c()
pca_reg <- list()
DSC_res <- c()

dat_pca <- FactoMineR::PCA(t(dat),
                           scale.unit = FALSE,
                           ncp = min(c(n_pc_max, dim(dat))),
                           graph = FALSE)

for (i in 1:length(class)) {
  # PCA based
  pca_silh[i] <- kBET::batch_sil(list(x = dat_pca$ind$coord[!is.na(class[[i]]),]),
                                 class[[i]][!is.na(class[[i]])],
                                 nPCs = min(c(n_pc_max, dim(dat))))
  pca_reg[[i]] <- try(kBET::pcRegression(list(x = dat_pca$ind$coord[!is.na(class[[i]]),], 
                                          sdev = sqrt(dat_pca$eig[1:n_pc_max,"eigenvalue"])),
                                     class[[i]][!is.na(class[[i]])],
                                     n_top = min(c(n_pc_max, dim(dat)))))
  
  # Other
  DSC_res[i] <- COPS::DSC(dat[,!is.na(class[[i]])], class[[i]][!is.na(class[[i]])])
}

PC_R2_max = sapply(pca_reg, function(a) ifelse(class(a) != "try-error", a[["maxR2"]], NA))

data.frame(var = names(class)[order(DSC_res)], DSC = sort(DSC_res))
data.frame(var = names(class)[order(pca_silh)], pca_silh = sort(pca_silh))
data.frame(var = names(class)[order(PC_R2_max)], PC_R2_max = sort(PC_R2_max, na.last = FALSE))

res_umap <- uwot::umap(dat_pca$ind$coord, 
                       n_neighbors = 20, 
                       n_components = 2, 
                       pca = NULL, 
                       verbose = FALSE, 
                       init = "normlaplacian")

res_umap <- data.frame(Dim.1 = res_umap[,1], Dim.2 = res_umap[,2])


class_name <- "SequencerSerial"
res_umap$category <- class[[class_name]]
ggplot(res_umap, aes(Dim.1, Dim.2, color = category)) + geom_point(shape = "+", size = 3) + 
  theme_bw() + scale_color_brewer(palette = "Dark2") + 
  labs(x = "Z1", y = "Z2", color = class_name)
ggplot(res_umap, aes(Dim.1, Dim.2, color = category == "SN587")) + geom_point(shape = "+", size = 3) + 
  theme_bw() + scale_color_brewer(palette = "Dark2") + 
  labs(x = "Z1", y = "Z2", color = class_name)

res_umap$size <- scanb_y[ss_filter, "Size.mm"][[1]]
ggplot(res_umap, aes(Dim.1, Dim.2, color = log10(size))) + geom_point(shape = "+", size = 3) + 
  theme_bw() + scale_color_distiller(palette = "RdBu") + 
  labs(x = "Z1", y = "Z2", color = "Size log10(mm)")


res_umap$subtype <- factor(scanb_y[ss_filter, "NCN.PAM50"][[1]])
ggplot(res_umap, aes(Dim.1, Dim.2, color = subtype)) + geom_point(shape = "+", size = 3) + 
  theme_bw() + scale_color_brewer(palette = "Dark2") + 
  labs(x = "Z1", y = "Z2", color = "NCN PAM50")

res_umap$pool <- factor(scanb_y[ss_filter, "PoolName"][[1]])
ggplot(res_umap, aes(Dim.1, Dim.2, color = pool)) + geom_point(shape = "+", size = 3) + 
  theme_bw() + #scale_color_brewer(palette = "Dark2") + 
  labs(x = "Z1", y = "Z2") + theme(legend.position="none")

res_umap$libbatch <- factor(scanb_y[ss_filter, "LibraryBatchNo"][[1]])
ggplot(res_umap, aes(Dim.1, Dim.2, color = libbatch)) + geom_point(shape = "+", size = 3) + 
  theme_bw() + #scale_color_brewer(palette = "Dark2") + 
  labs(x = "Z1", y = "Z2") + theme(legend.position="none")

table(table(scanb_y$LibraryBatchNo, scanb_y$PoolName))

# Remove samples with duplicated ids
duplicated_ids <- unique(scanb_y[["Patient"]][duplicated(scanb_y[["Patient"]])])
ss_filter <- scanb_y[["SequencerSerial"]] %in% c("NB501986", "NS500281", "SN587")
di_ss_filter <- ss_filter & (!scanb_y[["Patient"]] %in% duplicated_ids)


if (FALSE) {
  scanb_corrected1 <- sva::ComBat_seq(fixed_matrix[,di_ss_filter], 
                                      batch = factor(scanb_y[di_ss_filter, "SequencerSerial"][[1]]))
  scanb_corrected2 <- sva::ComBat_seq(scanb_corrected1, 
                                      batch = factor(scanb_y[di_ss_filter, "LibraryProtocol"][[1]]))
  scanb_corrected <- scanb_corrected2
} else {
  scanb_batch <- interaction(scanb_y[di_ss_filter,"LibraryProtocol"][[1]], 
                             scanb_y[di_ss_filter,"SequencerSerial"][[1]])
  scanb_corrected <- sva::ComBat_seq(fixed_matrix[,di_ss_filter], batch = droplevels(scanb_batch))
}

dge <- edgeR::DGEList(counts=scanb_corrected)
dge <- edgeR::calcNormFactors(dge, method = "upperquartile")
scanb_normalized <- edgeR::cpm(dge, normalized.lib.sizes = TRUE)

# run batch-association analysis again
dat <- log2(scanb_normalized+1)
class <- scanb_y[di_ss_filter,n_unique <= 600 & n_unique > 1]
n_pc_max <- 20

if (is.null(dim(class))) class <- data.frame(class = as.character(class))
class <- lapply(class, as.factor)
class <- class[sapply(class, function(x) length(levels(x)))>1]

pca_silh <- c()
pca_reg <- list()
DSC_res <- c()

dat_pca <- FactoMineR::PCA(t(dat),
                           scale.unit = FALSE,
                           ncp = min(c(n_pc_max, dim(dat))),
                           graph = FALSE)

for (i in 1:length(class)) {
  # PCA based
  pca_silh[i] <- kBET::batch_sil(list(x = dat_pca$ind$coord[!is.na(class[[i]]),]),
                                 class[[i]][!is.na(class[[i]])],
                                 nPCs = min(c(n_pc_max, dim(dat))))
  pca_reg[[i]] <- try(kBET::pcRegression(list(x = dat_pca$ind$coord[!is.na(class[[i]]),], 
                                              sdev = sqrt(dat_pca$eig[1:n_pc_max,"eigenvalue"])),
                                         class[[i]][!is.na(class[[i]])],
                                         n_top = min(c(n_pc_max, dim(dat)))))
  
  # Other
  DSC_res[i] <- COPS::DSC(dat[,!is.na(class[[i]])], class[[i]][!is.na(class[[i]])])
}

PC_R2_max = sapply(pca_reg, function(a) ifelse(class(a) != "try-error", a[["maxR2"]], NA))

data.frame(var = names(class)[order(DSC_res)], DSC = sort(DSC_res))
data.frame(var = names(class)[order(pca_silh)], pca_silh = sort(pca_silh))
data.frame(var = names(class)[order(PC_R2_max)], PC_R2_max = sort(PC_R2_max, na.last = FALSE))

res_umap <- uwot::umap(dat_pca$ind$coord, 
                       n_neighbors = 20, 
                       n_components = 2, 
                       pca = NULL, 
                       verbose = FALSE, 
                       init = "normlaplacian")

res_umap <- data.frame(Dim.1 = res_umap[,1], Dim.2 = res_umap[,2])

res_umap <- data.frame(res_umap, class)
require(ggplot2)
ggplot(res_umap, aes(Dim.1, Dim.2, color = SequencerSerial)) + geom_point(shape = "+", size = 3) + 
  theme_bw() + scale_color_brewer(palette = "Dark2") + 
  labs(x = "Z1", y = "Z2", color = "SequencerSerial")

ggplot(res_umap, aes(Dim.1, Dim.2, color = PoolName)) + geom_point(shape = "+", size = 3) + 
  theme_bw() + #scale_color_brewer(palette = "Dark2") + 
  labs(x = "Z1", y = "Z2", color = "PoolName") + theme(legend.position="none")

ggplot(res_umap, aes(Dim.1, Dim.2, color = LibraryBatchNo)) + geom_point(shape = "+", size = 3) + 
  theme_bw() + #scale_color_brewer(palette = "Dark2") + 
  labs(x = "Z1", y = "Z2", color = "LibraryBatchNo") + theme(legend.position="none")

ggplot(res_umap, aes(Dim.1, Dim.2, color = LibraryBatchNo)) + geom_point(shape = "+", size = 3) + 
  theme_bw() + #scale_color_brewer(palette = "Dark2") + 
  labs(x = "Z1", y = "Z2", color = "LibraryBatchNo") + theme(legend.position="none")

###########################################
# Re-apply batch correction
###########################################
# Subtype distribution in batches
subtype_dist <- sweep(table(res_umap$PoolName, res_umap$SSP.Subtype), 1, 
                      apply(table(res_umap$PoolName, res_umap$SSP.Subtype), 1, sum), FUN = "/")
subtype_dist <- sweep(table(res_umap$LibraryBatchNo, res_umap$SSP.Subtype), 1, 
                      apply(table(res_umap$LibraryBatchNo, res_umap$SSP.Subtype), 1, sum), FUN = "/")

hist(subtype_dist[,1])

pool_name_table <- table(scanb_y[di_ss_filter,"PoolName"][[1]])
pool_filter <- scanb_y[di_ss_filter, "PoolName"][[1]] %in% names(pool_name_table)[pool_name_table > 2]

batch_table <- table(scanb_y[di_ss_filter,"LibraryBatchNo"][[1]])
batch_filter <- scanb_y[di_ss_filter, "LibraryBatchNo"][[1]] %in% names(batch_table)[batch_table > 2]

second_filter <- pool_filter & batch_filter
scanb_corrected2 <- sva::ComBat_seq(scanb_corrected[,second_filter], 
                                    batch = factor(scanb_y[di_ss_filter,"PoolName"][[1]][second_filter]),
                                    group  = factor(scanb_y[di_ss_filter,"SSP.Subtype"][[1]][second_filter]))
scanb_corrected3 <- sva::ComBat_seq(scanb_corrected[,second_filter], #scanb_corrected2, 
                                    batch = factor(scanb_y[di_ss_filter,"LibraryBatchNo"][[1]][second_filter]),
                                    group  = factor(scanb_y[di_ss_filter,"SSP.Subtype"][[1]][second_filter]),
                                    shrink = TRUE, shrink.disp = TRUE)

dge <- edgeR::DGEList(counts=scanb_corrected3)
dge <- edgeR::calcNormFactors(dge, method = "upperquartile")
scanb_normalized <- edgeR::cpm(dge, normalized.lib.sizes = TRUE)

# run batch-association analysis again
dat <- log2(scanb_normalized+1)
class <- scanb_y[di_ss_filter,][second_filter,][,n_unique <= 600 & n_unique > 1]
n_pc_max <- 20

if (is.null(dim(class))) class <- data.frame(class = as.character(class))
class <- lapply(class, as.factor)
class <- class[sapply(class, function(x) length(levels(x)))>1]

pca_silh <- c()
pca_reg <- list()
DSC_res <- c()

dat_pca <- FactoMineR::PCA(t(dat),
                           scale.unit = FALSE,
                           ncp = min(c(n_pc_max, dim(dat))),
                           graph = FALSE)

for (i in 1:length(class)) {
  # PCA based
  pca_silh[i] <- kBET::batch_sil(list(x = dat_pca$ind$coord[!is.na(class[[i]]),]),
                                 class[[i]][!is.na(class[[i]])],
                                 nPCs = min(c(n_pc_max, dim(dat))))
  pca_reg[[i]] <- try(kBET::pcRegression(list(x = dat_pca$ind$coord[!is.na(class[[i]]),], 
                                              sdev = sqrt(dat_pca$eig[1:n_pc_max,"eigenvalue"])),
                                         class[[i]][!is.na(class[[i]])],
                                         n_top = min(c(n_pc_max, dim(dat)))))
  
  # Other
  DSC_res[i] <- COPS::DSC(dat[,!is.na(class[[i]])], class[[i]][!is.na(class[[i]])])
}

PC_R2_max = sapply(pca_reg, function(a) ifelse(class(a) != "try-error", a[["maxR2"]], NA))

data.frame(var = names(class)[order(DSC_res)], DSC = sort(DSC_res))
data.frame(var = names(class)[order(pca_silh)], pca_silh = sort(pca_silh))
data.frame(var = names(class)[order(PC_R2_max)], PC_R2_max = sort(PC_R2_max, na.last = FALSE))

res_umap <- uwot::umap(dat_pca$ind$coord, 
                       n_neighbors = 20, 
                       n_components = 2, 
                       pca = NULL, 
                       verbose = FALSE, 
                       init = "normlaplacian")

res_umap <- data.frame(Dim.1 = res_umap[,1], Dim.2 = res_umap[,2])

res_umap <- data.frame(res_umap, class)
require(ggplot2)
ggplot(res_umap, aes(Dim.1, Dim.2, color = SequencerSerial)) + geom_point(shape = "+", size = 3) + 
  theme_bw() + scale_color_brewer(palette = "Dark2") + 
  labs(x = "Z1", y = "Z2")

ggplot(res_umap, aes(Dim.1, Dim.2, color = SSP.Subtype)) + geom_point(shape = "+", size = 3) + 
  theme_bw() + scale_color_brewer(palette = "Dark2") + 
  labs(x = "Z1", y = "Z2")

ggplot(res_umap, aes(Dim.1, Dim.2, color = PoolName)) + geom_point(shape = "+", size = 3) + 
  theme_bw() + #scale_color_brewer(palette = "Dark2") + 
  labs(x = "Z1", y = "Z2") + theme(legend.position="none")

ggplot(res_umap, aes(Dim.1, Dim.2, color = LibraryBatchNo)) + geom_point(shape = "+", size = 3) + 
  theme_bw() + #scale_color_brewer(palette = "Dark2") + 
  labs(x = "Z1", y = "Z2") + theme(legend.position="none")
