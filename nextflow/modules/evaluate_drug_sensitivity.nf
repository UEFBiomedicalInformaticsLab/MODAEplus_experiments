#!/usr/bin/env nextflow
process evaluate_drug_sensitivity {
    
    publishDir 'results', mode: 'copy'
    container '../../container/modae.sif'

    input:
    path embedding_files
    path data_spec_file

    output:
    path '*drug_response_mse.csv.gz'
    path '*drug_response_r2.csv.gz'
    path '*drug_response_cor.csv.gz'

    script:
    """
    python \
        /scripts/python/hyperparameter_evaluation/drug_sensitivity_performance.py \
        2> >(grep -v " I " >&2)
    """
}