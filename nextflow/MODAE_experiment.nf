#!/usr/bin/env nextflow
/*
 * Pipeline parameters
 */
 
 params.data_config_file = 'data_config.json'
 params.model_config_file = 'model_config.json'
 params.search_iterations = 100

 // Include modules
include { MODAE_setup } from './modules/MODAE_setup.nf'
include { MODAE_run_cv } from './modules/MODAE_run_cv.nf'

workflow {

    // Setup data files
    MODAE_setup( params.data_config_file )

    // Setup search iteration channel
    search_channel = Channel.of( 1..params.search_iterations )

    // Run iterations
    MODAE_run_cv(
        params.model_config_file, 
        MODAE_setup.out.data_spec_file, 
        search_channel
    )
}