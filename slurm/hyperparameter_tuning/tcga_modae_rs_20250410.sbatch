#!/bin/bash
#SBATCH --job-name=MODAE_HPO
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=5
##SBATCH --partition=gpu
#SBATCH --partition=small
#SBATCH --array=1-100%50
#SBATCH --mem=25G
##SBATCH --gres=gpu:1
#SBATCH -t 3-00:00:00
#SBATCH --output=log/MODAE_HPO.%A_%a.out
#SBATCH --error=log/MODAE_HPO.%A_%a.err
##SBATCH --mail-user=
#SBATCH --mail-type=NONE

# Host paths
source path.conf

# Container internal paths
export SAVE_PATH=/output/20250410_random_search/

#export SAVE_PATH=${SAVE_PATH}lung_test/
#export NAME=lung_test_

#export SAVE_PATH=${SAVE_PATH}brca_test/
#export NAME=brca_test_

#export SAVE_PATH=${SAVE_PATH}brca_test_noclfilter_pretrain/
#export NAME=brca_test_noclfilter_pretrain_

#export SAVE_PATH=${SAVE_PATH}scanb_test_noclfilter/
#export NAME=scanb_test_noclfilter_

export SAVE_PATH=${SAVE_PATH}pancan_test/
export NAME=pancan_test_

apptainer run \
  --bind $MODAE_OUTPUT_PATH:/output \
  ../container/modae.sif \
    -l mrna \
    -o $NAME \
    -p $SAVE_PATH \
    --max_cv_repeats=1 \
    --epochs=150 \
    --no_pre_train_ae \
    --no_pre_train_cl \
    --no_pre_train_sr \
    --no_pre_train_bd \
    --no_pre_train_bc \
    --no_pre_train_dr \
    --fast_adaptive_multitask_optimization \
    --famo_weight_decay=0.00005-0.005 \
    --decoder_split_loss \
    --classifier_split_loss \
    --hidden_activation=selu \
    --minibatchsize=32 \
    --mse \
    --survival_model \
    --drug_response_model \
    --batch_correction \
    --classifier \
    --batch_loss=wasserstein \
    --encoder_layers=512-128,128-32,64-8 \
    --decoder_layers=128-32,128-512 \
    --survival_model_layers=8-16 \
    --drug_response_model_layers=32-128 \
    --drug_response_model_drugwise_layers=2-8 \
    --batch_detector_layers=32-128,16-64 \
    --classifier_layers=32-128,16-64 \
    --deconfounder_layers_per_batch=2-8 \
    --batch_adversarial_gradient_penalty=1. \
    --deconfounder_norm_penalty=1. \
    --drug_response_model_output_activation=relu \
    --learning_rate=0.000001-0.0001 \
    --adversarial_learning_rate_multi=10. \
    --adversarial_pre_training_optimizer=adam \
    --adam_beta1=0. \
    --adam_beta2=0.9 \
    --reg_weights_type= \
    --reg_weights=0. \
    --reg_noise_sd=0. \
    --reg_dropout_rate_input=0.-0.5 \
    --reg_dropout_rate_autoencoder=0. \
    --reg_dropout_rate_survival=0.-0.5 \
    --reg_dropout_rate_classifier=0.-0.25 \
    --reg_dropout_rate_batch=0.-0.25 \
    --reg_dropout_rate_drug_response=0.-0.5 \
    --data_var_cutoff_rank=-1 \
    --serialize_data \
    --parallel \
    --nthreads=1 \
    --nthreads_interop=1 \
    --nprocess=5 \
    --ptimeout=0 \
    --cache_data \
    --no_diagnostics \
    2> >(grep -v " I " >&2)

