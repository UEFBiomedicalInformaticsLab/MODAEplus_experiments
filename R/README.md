## R-based analyses

These have been separated into different folders and scripts. These scripts 
have been created to handle pre-processing (and retrieval) of the data as 
well as the post-hoc analysis of the results after training and main 
evaluation python scripts have run. 

* Hyperparameter tuning summary and parameter selection. 
* Performance summaries comparing to baselines. 
* External validation summaries.
* Embedding accuracy analysis.
* Survival based drug-response analysis.
* Differential drug-target expression analysis. 

### setup.R

This file defines many helpeor functions and sets many of the paths based 
on the environment variables MODAE_DATA_PATH and MODAE_OUTPUT_PATH. This 
works as long as the data is organized in the same way as in Zenodo 
(or group folder). 

## analysis

Various scripts for analysing the results. All plotting has been wrapped 
with the function save_figure_safe, so all plots can be found with:
```
grep -r . -ne save_figure_safe
```
Plots are saved to MODAE_OUTPUT_PATH/yyyymmdd_random_search/plots. 
An overview of scripts and their output:

### Important analysis scripts
* analyse_results.R combines all the hyper-parameter tuning results and
  determines the best search iteration (a.k.a. "task").
  * Pareto plots named "objectives_plot*.png" with a few variations, main
    one being objectives_plot_colored_grad.
  * Plot of evaluation metrics named "metrics_plot.png" at training,
    validation, retraining, and testing.
  * Plot of drug prediction correlations named
    "drug_prediction_correlations_example.pdf". 
  * Plot of training losses over time named "diagnostic_plot_*.pdf".
    If the --no_diagnostics flag has been set these will not be
    generated. In the present script it has been set for
    hyperparameter tuning, but not for final training (plotted later). 
  * Legacy plots of embedding UMAPs, PCs, and t-SNEs using hyperparameter
    tuning results by combining training and test embeddings (i.e.,
    not from final training). Named "embedding_*_brca_classes.png", 
    "embedding_*_cancer_types.png", and "embedding_*_cancer_types.png".
  * Plots of pairwise embeddings are currently disabled as they were
    not particularly informative and take a long time to run.
    Saved in "pairwise_embeddings/" if run. 
* compare_sensitivity_predictions.R can be used to compare the drug
  sensitivity predictions of two results.
  * Plot of correlations between predictions for selected drugs. Saved
    in "patient_sensitivity_cor/".
  * "sensitivity_difference_gene_correlation_histogram.png"
  * "patient_pred_cor.csv"
* diagnostics.R plots diagnostics from final training run as
  "external_diagnostic_plot_*.pdf".
* drug_responders_cancer_specific.R visualizes CTRDB TCGA patient drug
  response vs predicted drug sensitivity.
  * Saved in "treatment_response/".
  * Box- and scatter-plots.
  * Pan-cancer version in "drug_responders.R"
* drug_responders_nontcga.R visualizes CTRDB non-TCGA patient drug
  responses vs predicted drug sensitivity.
  * Caveat: many of these show signs of bad applicability domain.
  * Embedding UMAPs combined with original data saved in
    "treatment_response/" as "*_embeddings.png" and 
    "*_scatter.png".
* drug_targets.R performs t-tests and regression analysis between
  drug-target gene expression and predicted drug-sensitivity.
  * Saved in "drug_target_expression/".
  * Beta-regression for data bounded in (0,1). Includes code for
    zero-inflated models that should work better with [0,1],
    but running these models is too slow.
  * Resistant vs. sensitive boxplots of each drug's target
    genes expression. Pan-cancer and cancer-specific plots. 
  * Scatter plots of drug-target expression and predicted drug-sensitivity.
  * Tables in folders named by target source:
    "t_tested_differences_by_cancer.csv" and "t_tested_differences.csv".
* drug_target_summary.R summarizes drug-target expression t-test.
  * Saved in "drug_target_expression/". 
  * Summary tables: "t_test_significance_summary*.csv".
  * Bar-plots: "*_drug_target_expression_*_t_test*.png.
* embedding_plots.R creates visualizations of the embeddings as well as
  a comparison plot comparing final model performance vs baselines.
  * Plots: "*final_embedding_OT_level1.png",
    "performance_composite_plot.png".
  * CV result from hyper-parameter tuning (using final parameters).
* embedding_sensitivity_plots.R creates visualizations of the
  embeddings coloured by predicted drug-sensitivity.
  * Saved in "sensitivity_umaps/".
  * Pan-cancer and cancer-specific plots. 
* drug_treatment_survival.R visualizes the survival comparisons between 
  predicted resistant and sensitive patient groups. 
  * Saved in "treatment_survival_km".
  * Survival curves and Cox PH hazard ratios.
  * Generate equivalent plots for CODE-AE prediction for comparison.
  * Run survival::cox.zph to test the proportionality assumption 
    and plot them in "tcga_survival_assumptions*.pdf". 
  * Format results into summary tables. TODO: save tables.
* result_exploration.R plot drug sensitivity heatmaps with clnical
  information.
  * Saved in "sensitivity_heatmaps/".
* tissue_classifier_comparison.R combines tissue classification metrics
  from kNN classifiers for different embeddings including CODE-AE,
  Celligner, ComBat, and ablated MODAE models.
  * Tables: "combined_tissue_metrics.csv", "tissue_metrics.csv".
* upload_ready_outputs.R combines embeddings, predictions, and
  several clinical information sources to create a single output
  table for all TCGA patients.
  * Saved in "../pancan_test/external_evaluation/internal/".
  * Table: "unified_patient_output_table.csv.gz"
  * Also formats best model drug performance into a table:
    "best_model_drugwise_performance.csv" (saved in "plots/").
* validation_analysis.R summarizes external validation performance.
  * Table: "external_drug_response_performance.csv".
  * Plot: "external_validation_drug_performance_text.png". 

### Less important analysis scripts (prototypes to be deleted or completed)
* lung_cancer_summary.R visualizes LUAD and LUSC patient
  drug-sensitivity and gene-expression with heatmaps.
* nrf2_associations.R regression modeling and sensitivity group
  t-tests for NRF2-score (Härkönen et al. 2023) associations.
* prediction_clustering.R runs random-forest based clustering on
  drug-sensitivity predictions and clinical annotations.
  * Incomplete.

## baseline_models
Code for batch correction baselines. 
* baseline_batch_correction.R

## data_exploration

## data_preprocessing

## plots
This folder is currently underutilized and some of the scripts in 
analysis should be moved here. Currently it contains:
* combined_embedding_plot.R which combines the UMAPs of 
  patient and cell-line tissue overlap.



