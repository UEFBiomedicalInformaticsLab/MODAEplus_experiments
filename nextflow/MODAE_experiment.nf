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
include { evaluate_batch_effect } from './modules/evaluate_batch_effect.nf'
include { evaluate_drug_sensitivity } from './modules/evaluate_drug_sensitivity.nf'
include { evaluate_tissue_classification } from './modules/evaluate_tissue_classification.nf'
include { multi_objective_selection } from './modules/multi_objective_selection.nf'

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

    // Run additional evaluation
    evaluate_batch_effect(
        MODAE_run_cv.out.embedding_files.collect(),
        MODAE_setup.out.data_spec_file
    )
    evaluate_drug_sensitivity(
        MODAE_run_cv.out.prediction_files.collect(),
        MODAE_setup.out.data_spec_file
    )
    evaluate_tissue_classification(
        MODAE_run_cv.out.prediction_files.collect(),
        MODAE_setup.out.data_spec_file
    )

    // Compute multi-objective performance and identify best setting
    multi_objective_selection(
        MODAE_run_cv.out.metric_files.collect(),
        evaluate_batch_effect.out.collect(), 
        evaluate_drug_sensitivity.out.collect(), 
        evaluate_tissue_classification.out.collect()
    )
    
    // Final training 
    
}