#!/usr/bin/env nextflow
process evaluate_tissue_classification {
    
    publishDir 'results', mode: 'copy'
    container '../../container/modae.sif'

    input:
    path prediction_files
    path data_spec_file

    output:
    path '*tissue_classifier_performance.csv.gz'

    script:
    """
    python \
        /scripts/python/hyperparameter_evaluation/tissue_prediction_performance.py \
        2> >(grep -v " I " >&2)
    """
}