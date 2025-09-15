#!/usr/bin/env nextflow
process MODAE_setup {
    
    publishDir 'results', mode: 'copy'
    container '../../container/modae.sif'

    input:
    path data_config_file

    output:
    path '*.tfrecords'
    path '*serialized_data_spec.json', emit: data_spec_file

    script:
    """
    python -m modae --config='$data_config_file'
    """
}