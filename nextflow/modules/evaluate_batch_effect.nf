#!/usr/bin/env nextflow
process evaluate_batch_effect {
    
    publishDir 'results', mode: 'copy'
    container '../../container/modae.sif'

    input:
    path prediction_files
    path data_spec_file

    output:
    path '*batch_predictions.csv'

    script:
    """
    python \
        /scripts/python/hyperparameter_evaluation/batch_effect_analysis.py \
        2> >(grep -v " I " >&2)
    """
}