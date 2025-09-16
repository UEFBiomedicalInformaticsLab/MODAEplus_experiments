#!/usr/bin/env nextflow
process multi_objective_selection {
    
    publishDir 'results', mode: 'copy'
    container '../../container/R.sif'

    input:
    path metric_files
    path batch_effect_files
    path drug_sensitivity_files
    path tissue_classifier_files

    output:
    path 'best_task.txt'
    path '*.png'

    script:
    """
    Rscript --no-save --no-restore /scripts/R/analysis/analyse_results.R
    """
}