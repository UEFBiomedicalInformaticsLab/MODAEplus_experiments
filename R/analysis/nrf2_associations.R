script_root <- Sys.getenv("MODAE_SCRIPT_PATH")
if (script_root == "") {
  stop("Please define the MODAE_SCRIPT_PATH environment variable.")
}
source(paste0(script_root, "R/setup.R"))
#Sys.setenv(LD_LIBRARY_PATH = paste0("~/JAGS/lib:", Sys.getenv("LD_LIBRARY_PATH")))
#Sys.setenv(PKG_CONFIG_PATH = paste0("~/JAGS/lib/pkgconfig:", Sys.getenv("PKG_CONFIG_PATH")))
#Sys.setenv(JAGS_PREFIX = "~/JAGS")
#install.packages("rjags", configure.args="--with-jags-prefix")

for (save_dir in save_dirs) {
  if(exists("var_list")) try(detach(var_list))
  var_list <- get_paths_and_parameters(save_dir, remove_incomplete = FALSE)
  attach(var_list)
  best_task <- readLines(paste0(save_path, "best_task.txt"))
  param_best_task_ind <- match(best_task, parameters[["task"]])
  shared_embedding_names <- paste0("z", 1:parameters[param_best_task_ind, "bottle_neck"])
  
  # Final predictions
  prediction_paths <- paste0(save_path, "external_evaluation/")
  if (dir.exists(paste0(save_path, "reverse_external_evaluation"))) {
    prediction_paths[2] <- paste0(save_path, "reverse_external_evaluation/")
  }
  for (prediction_path in prediction_paths) {
    fn <- paste0(prediction_path, "internal/unified_patient_output_table.csv.gz")
    tcga_final_output_table <- readr::read_csv(fn)
    fn <- paste0(plot_path, "best_model_drugwise_performance.csv")
    dr_perf_table <- readr::read_csv(fn)
    
    sensitivity_cols <- grep("_sensitivity$", colnames(tcga_final_output_table))
    sensitivity_cols <- colnames(tcga_final_output_table)[sensitivity_cols]
    # should be identical order to names saved in file
    sensitivity_drug_names <- get_ctrp_drugs(n_drugs = length(sensitivity_cols))
    sensitivity_drug_names <- process_ctrp_drug_names(sensitivity_drug_names)
    #data.frame(gsub("_sensitivity$", "", sensitivity_cols), sensitivity_drug_names)
    sens_fractions <- lapply(
      sensitivity_cols, 
      function(x) tapply(
        tcga_final_output_table[[x]] > 0,
        tcga_final_output_table[["type"]], 
        mean
      )
    )
    sens_counts <- lapply(
      sensitivity_cols, 
      function(x) tapply(
        tcga_final_output_table[[x]] > 0,
        tcga_final_output_table[["type"]], 
        sum
      )
    )
    names(sens_fractions) <- sensitivity_cols
    names(sens_counts) <- sensitivity_cols
    if (FALSE) {
      fn <- paste0(prediction_path, "internal/NRF2_global_sens_cor.csv.gz")
      global_correlation_df <- readr::read_csv(fn)
      fn <- paste0(prediction_path, "internal/NRF2_cancer_sens_cor.csv.gz")
      cancer_correlation_df <- readr::read_csv(fn)
      fn <- paste0(prediction_path, "internal/NRF2_cancer_sens_betareg.csv.gz")
      cancer_betareg_df <- readr::read_csv(fn)
    }
    global_correlation_df <- data.frame()
    cancer_correlation_df <- data.frame()
    cancer_betareg_df <- data.frame()
    cancer_t_test_df <- data.frame()
    
    nrf2_cor_path <- paste0(plot_path, "sensitivity_nrf2_correlation/")
    for (i in pan_cancer_types) {
      dir.create(paste0(nrf2_cor_path, i), recursive = TRUE)
    }
    
    for (i in sensitivity_cols) { #[which(sensitivity_cols=="nilotinib_sensitivity"):length(sensitivity_cols)]) {
      #i <- "nilotinib_sensitivity"
      y <- tcga_final_output_table[[i]]
      if (sum(y>0) < 2) next
      temp <- data.frame(
        x = tcga_final_output_table[["NRF2_score"]], 
        y = y, 
        z = factor(tcga_final_output_table[["type"]])
      )
      temp <- temp[!is.na(temp[["x"]]),]
      temp <- temp[!is.na(temp[["y"]]),]
      if (FALSE) {
        zoib_mod <- zoib::zoib(
          y ~ x * z | 1 | x * z | 1, 
          data = temp, 
          zero.inflation = TRUE, 
          one.inflation = FALSE, 
          random = 0, 
          n.chain = 2, 
          n.iter = 1000, 
          n.burn = 100, 
          n.thin = 5
        )
        summary(zoib_mod$coeff)
      }
      if (FALSE) {
        betareg_mod <- betareg::betareg(y ~ x*z, temp)
        summary(betareg_mod)#$coefficients$mean[2,"Pr(>|z|)"]
      }
      if (FALSE) {
        test_to_df <- function(x) data.frame(rho = x$estimate, p.value = x$p.value)
        cor_df <- data.frame()
        for (cor_method in c("spearman", "pearson", "kendall")) {
          cor_df <- rbind(
            cor_df, 
            data.frame(
              treatment = i, 
              method = cor_method, 
              test_to_df(cor.test(temp[["x"]], temp[["y"]], method = cor_method))
            )
          )
        }
        global_correlation_df <- rbind(global_correlation_df, cor_df)
      }
      if (TRUE) { # cancer specific
        for (canceri in names(sens_counts[[i]])[sens_counts[[i]] > 1]) {
          temp_canceri <- temp[temp[["z"]] == canceri,]
          if (FALSE) {
            # correlation
            test_to_df <- function(x) data.frame(rho = x$estimate, p.value = x$p.value)
            cor_df <- data.frame()
            for (cor_method in c("spearman", "pearson", "kendall")) {
              cor_df <- rbind(
                cor_df, 
                data.frame(
                  cancer = canceri, 
                  treatment = i, 
                  method = cor_method, 
                  test_to_df(cor.test(temp_canceri[["x"]], temp_canceri[["y"]], method = cor_method))
                )
              )
            }
            cancer_correlation_df <- rbind(cancer_correlation_df, cor_df)
          }
          if (FALSE) {
            # beta regression
            try({
              betareg_mod <- betareg::betareg(y ~ x, temp_canceri)
              betareg_sum <- summary(betareg_mod)#$coefficients$mean[2,"Pr(>|z|)"]
              cancer_betareg_df <- rbind(
                cancer_betareg_df, 
                data.frame(
                  cancer = canceri, 
                  treatment = i, 
                  effect = betareg_sum$coefficients$mu["x", "Estimate"], 
                  p.value = betareg_sum$coefficients$mu["x", "Pr(>|z|)"]
                )
              )
            })
          }
          if (TRUE) {
            # group and t-test
            temp_canceri[["g"]] <- ifelse(
              temp_canceri[["y"]] > mean(temp_canceri[["y"]]), 
              "sensitive", 
              "resistant"
            )
            if (all(table(temp_canceri[["g"]])>1)) {
              res <- t.test(x ~ g, data = temp_canceri)
              
              p1 <- ggplot(temp_canceri, aes(x = x, y = y, color = g)) + 
                geom_point(shape = "+", size = 3) + 
                theme_bw() + 
                scale_color_brewer(palette = "Set1") +
                labs(x = "NRF2_score", y = i, color = "group")
              save_figure_safe(
                ggExtra::ggMarginal(p1, type = "boxplot", margins = "both", groupFill = TRUE), 
                png, 
                paste0(nrf2_cor_path, canceri, "/", i, "_nrf2_cor.png"), 
                width = plot_width * 1, 
                height = plot_width * 1, 
                res = plot_res, 
                units = plot_units
              )
              
              cancer_t_test_df <- rbind(
                cancer_t_test_df, 
                data.frame(
                  cancer = canceri, 
                  treatment = i, 
                  resistant_fraction = mean(temp_canceri[["g"]] == "resistant"), 
                  mean_sensitivity_in_resistant = with(temp_canceri, mean(y[g == "resistant"])), 
                  mean_sensitivity_in_sensitive = with(temp_canceri, mean(y[g == "sensitive"])), 
                  t(res$estimate), 
                  t.test.p.val = res$p.value
                )
              )
            }
          }
        }
      }
    }
    cancer_t_test_df[["effect"]] <- with(cancer_t_test_df, mean.in.group.sensitive - mean.in.group.resistant)
    hist(cancer_t_test_df[["resistant_fraction"]])
    ggplot(cancer_t_test_df, aes(x = effect, y = -log10(t.test.p.val))) + 
      geom_point(shape = "+", size = 3) + 
      theme_bw() + 
      facet_wrap(cancer ~., scales = "free") 
    
    fnw <- paste0(prediction_path, "internal/NRF2_global_sens_cor.csv.gz")
    readr::write_csv(global_correlation_df, file = gzfile(fnw))
    fnw <- paste0(prediction_path, "internal/NRF2_cancer_sens_cor.csv.gz")
    readr::write_csv(cancer_correlation_df, file = gzfile(fnw))
    fnw <- paste0(prediction_path, "internal/NRF2_cancer_sens_betareg.csv.gz")
    readr::write_csv(cancer_betareg_df, file = gzfile(fnw))
    fnw <- paste0(prediction_path, "internal/NRF2_cancer_sens_group_ttest.csv.gz")
    readr::write_csv(cancer_t_test_df, file = gzfile(fnw))
    ggplot(global_correlation_df, aes(x = rho, color = method)) + geom_density()
    ggplot(cancer_correlation_df, aes(x = rho, color = method)) + geom_density()
    ggplot(cancer_betareg_df, aes(x = effect)) + geom_density()
    df_filter <- with(cancer_betareg_df, cancer %in% c("LUAD", "LUSC"))
    df_filter <- df_filter & with(cancer_betareg_df, treatment %in% unique(treatment[df_filter & abs(effect) > 0.05]))
    ggplot(cancer_betareg_df[df_filter,], aes(x = treatment, y = effect, fill = cancer)) + 
      geom_bar(stat = "identity") + 
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90))
    
    
  }
}