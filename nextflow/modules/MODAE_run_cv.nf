#!/usr/bin/env nextflow
process MODAE_run_cv {
    
    publishDir 'results', mode: 'copy'
    container '../../container/modae.sif'

    input:
    path model_config_file
    path data_spec_file
    val task_id

    output:
    path '*embeddings*.csv.gz'
    path '*predictions*.csv.gz'
    path '*metrics*.csv'

    script:
    """
    python -m modae.training_program \
        --config='$model_config_file' \
        --data-spec='$data_spec_file' \
        --task_id='$task_id'
    """
}