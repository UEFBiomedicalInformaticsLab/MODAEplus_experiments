#!/bin/bash
#SBATCH --job-name=TCGA_SUPERAE
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=5
##SBATCH --partition=gpu
#SBATCH --partition=small
##SBATCH --array=1
#SBATCH --mem=25G
##SBATCH --gres=gpu:v100:1
##SBATCH --gres=gpu:1
#SBATCH -t 3-00:00:00
#SBATCH --output=/scratch/users/rintala/log/TCGA_SUPERAE.%A_%a.out
#SBATCH --error=/scratch/users/rintala/log/TCGA_SUPERAE.%A_%a.err
#SBATCH --mail-user=teemu.rintala@uef.fi
#SBATCH --mail-type=NONE
##SBATCH --account=project_2005873
#SBATCH --exclude=kudos[12]
#sleep 5

klist -l
#export PROJ_PATH=/scratch/project_2005873
export PROJ_PATH=/research/work/rintala
export SAVE_PATH=$PROJ_PATH/superAE_HPO/20250410_random_search/
touch ${SAVE_PATH}touch_file

export SAVE_PATH=${SAVE_PATH}pancan_ablation_test/
export NAME=pancan_test_

#module load tensorflow/2.13
#source $HOME/tensorflow_23/bin/activate
source /userapps/rintala/tensorflow/bin/activate
export TF_CPP_MIN_LOG_LEVEL=2
#export PROG_PATH=/projappl/project_2005873/superAE/src/sae/training_program.py
#export PROG_PATH=$HOME/superAE/src/sae/training_program.py
if false; then
python -m sae.training_program \
    -l mrna \
    -o $NAME \
    -p $SAVE_PATH \
    --max_cv_repeats=1 \
    --epochs=150 \
    --epochs_pre_ae=1 \
    --no_pre_train_ae \
    --epochs_pre_cl=1 \
    --no_pre_train_cl \
    --epochs_pre_sr=1 \
    --no_pre_train_sr \
    --epochs_pre_bd=1 \
    --no_pre_train_bd \
    --epochs_pre_bc=1 \
    --no_pre_train_bc \
    --epochs_pre_dr=1 \
    --no_pre_train_dr \
    --fast_adaptive_multitask_optimization \
    --famo_weight_decay=0.0003551639396349761 \
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
    --encoder_layers=226,106,16 \
    --decoder_layers=95,445 \
    --survival_model_layers=16 \
    --drug_response_model_layers=67 \
    --drug_response_model_drugwise_layers=6 \
    --batch_detector_layers=36,64 \
    --classifier_layers=105,31 \
    --deconfounder_layers_per_batch=5 \
    --batch_adversarial_gradient_penalty=1. \
    --deconfounder_norm_penalty=1. \
    --drug_response_model_output_activation=relu \
    --pre_learning_rate_ae=0.00008590623250698555 \
    --pre_learning_rate_cl=0.00008590623250698555 \
    --pre_learning_rate_sr=0.00008590623250698555 \
    --pre_learning_rate_bd=0.00008590623250698555 \
    --pre_learning_rate_bc=0.00008590623250698555 \
    --pre_learning_rate_dr=0.00008590623250698555 \
    --learning_rate=0.00008590623250698555 \
    --adversarial_learning_rate_multi=10. \
    --adversarial_pre_training_optimizer=adam \
    --batch_correction_pre_training_optimizer=adam \
    --adam_beta1=0. \
    --adam_beta2=0.9 \
    --reg_weights_type= \
    --reg_weights=0. \
    --reg_noise_sd=0. \
    --reg_dropout_rate_input=0.29179563810929937 \
    --reg_dropout_rate_autoencoder=0. \
    --reg_dropout_rate_survival=0.03478547730630027 \
    --reg_dropout_rate_classifier=0.21685112099827208 \
    --reg_dropout_rate_batch=0.033310129812936934 \
    --reg_dropout_rate_drug_response=0.08906233077974918 \
    --data_var_cutoff_rank=-1 \
    --serialize_data \
    --parallel \
    --nthreads=1 \
    --nthreads_interop=1 \
    --nprocess=5 \
    --ptimeout=0 \
    --gpu_memory=4800 \
    --cache_data \
    --no_diagnostics \
    2> >(grep -v " I " >&2)
fi

if true; then
export NAME=pancan_test_no_private_without_deconfounding_
python -m sae.training_program \
    -l mrna \
    -o $NAME \
    -p $SAVE_PATH \
    --max_cv_repeats=1 \
    --epochs=150 \
    --epochs_pre_ae=1 \
    --no_pre_train_ae \
    --epochs_pre_cl=1 \
    --no_pre_train_cl \
    --epochs_pre_sr=1 \
    --no_pre_train_sr \
    --epochs_pre_bd=1 \
    --no_pre_train_bd \
    --epochs_pre_bc=1 \
    --no_pre_train_bc \
    --epochs_pre_dr=1 \
    --no_pre_train_dr \
    --fast_adaptive_multitask_optimization \
    --famo_weight_decay=0.0003551639396349761 \
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
    --encoder_layers=226,106,16 \
    --decoder_layers=95,445 \
    --survival_model_layers=16 \
    --drug_response_model_layers=67 \
    --drug_response_model_drugwise_layers=6 \
    --batch_detector_layers=36,64 \
    --classifier_layers=105,31 \
    --deconfounder_layers_per_batch=5 \
    --batch_adversarial_gradient_penalty=1. \
    --deconfounder_norm_penalty=0. \
    --drug_response_model_output_activation=relu \
    --pre_learning_rate_ae=0.00008590623250698555 \
    --pre_learning_rate_cl=0.00008590623250698555 \
    --pre_learning_rate_sr=0.00008590623250698555 \
    --pre_learning_rate_bd=0.00008590623250698555 \
    --pre_learning_rate_bc=0.00008590623250698555 \
    --pre_learning_rate_dr=0.00008590623250698555 \
    --learning_rate=0.00008590623250698555 \
    --adversarial_learning_rate_multi=10. \
    --adversarial_pre_training_optimizer=adam \
    --batch_correction_pre_training_optimizer=adam \
    --adam_beta1=0. \
    --adam_beta2=0.9 \
    --reg_weights_type= \
    --reg_weights=0. \
    --reg_noise_sd=0. \
    --reg_dropout_rate_input=0.29179563810929937 \
    --reg_dropout_rate_autoencoder=0. \
    --reg_dropout_rate_survival=0.03478547730630027 \
    --reg_dropout_rate_classifier=0.21685112099827208 \
    --reg_dropout_rate_batch=0.033310129812936934 \
    --reg_dropout_rate_drug_response=0.08906233077974918 \
    --data_var_cutoff_rank=-1 \
    --serialize_data \
    --parallel \
    --nthreads=1 \
    --nthreads_interop=1 \
    --nprocess=5 \
    --ptimeout=0 \
    --gpu_memory=4800 \
    --cache_data \
    --no_diagnostics \
    2> >(grep -v " I " >&2)
fi

if false; then
export NAME=pancan_test_no_deconfounding_
python -m sae.training_program \
    -l mrna \
    -o $NAME \
    -p $SAVE_PATH \
    --max_cv_repeats=1 \
    --epochs=150 \
    --epochs_pre_ae=1 \
    --no_pre_train_ae \
    --epochs_pre_cl=1 \
    --no_pre_train_cl \
    --epochs_pre_sr=1 \
    --no_pre_train_sr \
    --epochs_pre_bd=1 \
    --no_pre_train_bd \
    --epochs_pre_bc=1 \
    --no_pre_train_bc \
    --epochs_pre_dr=1 \
    --no_pre_train_dr \
    --fast_adaptive_multitask_optimization \
    --famo_weight_decay=0.0003551639396349761 \
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
    --encoder_layers=226,106,21 \
    --decoder_layers=95,445 \
    --survival_model_layers=16 \
    --drug_response_model_layers=67 \
    --drug_response_model_drugwise_layers=6 \
    --batch_detector_layers=36,64 \
    --classifier_layers=105,31 \
    --deconfounder_layers_per_batch=0 \
    --batch_adversarial_gradient_penalty=1. \
    --deconfounder_norm_penalty=0. \
    --drug_response_model_output_activation=relu \
    --pre_learning_rate_ae=0.00008590623250698555 \
    --pre_learning_rate_cl=0.00008590623250698555 \
    --pre_learning_rate_sr=0.00008590623250698555 \
    --pre_learning_rate_bd=0.00008590623250698555 \
    --pre_learning_rate_bc=0.00008590623250698555 \
    --pre_learning_rate_dr=0.00008590623250698555 \
    --learning_rate=0.00008590623250698555 \
    --adversarial_learning_rate_multi=10. \
    --adversarial_pre_training_optimizer=adam \
    --batch_correction_pre_training_optimizer=adam \
    --adam_beta1=0. \
    --adam_beta2=0.9 \
    --reg_weights_type= \
    --reg_weights=0. \
    --reg_noise_sd=0. \
    --reg_dropout_rate_input=0.29179563810929937 \
    --reg_dropout_rate_autoencoder=0. \
    --reg_dropout_rate_survival=0.03478547730630027 \
    --reg_dropout_rate_classifier=0.21685112099827208 \
    --reg_dropout_rate_batch=0.033310129812936934 \
    --reg_dropout_rate_drug_response=0.08906233077974918 \
    --data_var_cutoff_rank=-1 \
    --serialize_data \
    --parallel \
    --nthreads=1 \
    --nthreads_interop=1 \
    --nprocess=5 \
    --ptimeout=0 \
    --gpu_memory=4800 \
    --cache_data \
    --no_diagnostics \
    2> >(grep -v " I " >&2)
fi

if false; then
export NAME=pancan_test_no_famo_
python -m sae.training_program \
    -l mrna \
    -o $NAME \
    -p $SAVE_PATH \
    --max_cv_repeats=1 \
    --epochs=150 \
    --epochs_pre_ae=1 \
    --no_pre_train_ae \
    --epochs_pre_cl=1 \
    --no_pre_train_cl \
    --epochs_pre_sr=1 \
    --no_pre_train_sr \
    --epochs_pre_bd=1 \
    --no_pre_train_bd \
    --epochs_pre_bc=1 \
    --no_pre_train_bc \
    --epochs_pre_dr=1 \
    --no_pre_train_dr \
    --decoder_split_loss \
    --classifier_split_loss \
    --reconstruction_weight=0.34593765 \
    --reconstruction_cl_weight=0.45876225 \
    --classifier_weight=0.01483179 \
    --classifier_cl_weight=0.01327784 \
    --survival_risk_weight=0.01474393 \
    --domain_adaptation_weight=0.01309518 \
    --drug_regression_weight=0.10831788 \
    --deconfounder_norm_penalty=0.03103348 \
    --objective_weight_norm_order=-1. \
    --hidden_activation=selu \
    --minibatchsize=32 \
    --mse \
    --survival_model \
    --drug_response_model \
    --batch_correction \
    --classifier \
    --batch_loss=wasserstein \
    --encoder_layers=226,106,16 \
    --decoder_layers=95,445 \
    --survival_model_layers=16 \
    --drug_response_model_layers=67 \
    --drug_response_model_drugwise_layers=6 \
    --batch_detector_layers=36,64 \
    --classifier_layers=105,31 \
    --deconfounder_layers_per_batch=5 \
    --batch_adversarial_gradient_penalty=1. \
    --drug_response_model_output_activation=relu \
    --pre_learning_rate_ae=0.00008590623250698555 \
    --pre_learning_rate_cl=0.00008590623250698555 \
    --pre_learning_rate_sr=0.00008590623250698555 \
    --pre_learning_rate_bd=0.00008590623250698555 \
    --pre_learning_rate_bc=0.00008590623250698555 \
    --pre_learning_rate_dr=0.00008590623250698555 \
    --learning_rate=0.00008590623250698555 \
    --adversarial_learning_rate_multi=10. \
    --adversarial_pre_training_optimizer=adam \
    --batch_correction_pre_training_optimizer=adam \
    --adam_beta1=0. \
    --adam_beta2=0.9 \
    --reg_weights_type= \
    --reg_weights=0. \
    --reg_noise_sd=0. \
    --reg_dropout_rate_input=0.29179563810929937 \
    --reg_dropout_rate_autoencoder=0. \
    --reg_dropout_rate_survival=0.03478547730630027 \
    --reg_dropout_rate_classifier=0.21685112099827208 \
    --reg_dropout_rate_batch=0.033310129812936934 \
    --reg_dropout_rate_drug_response=0.08906233077974918 \
    --data_var_cutoff_rank=-1 \
    --serialize_data \
    --parallel \
    --nthreads=1 \
    --nthreads_interop=1 \
    --nprocess=5 \
    --ptimeout=0 \
    --gpu_memory=4800 \
    --cache_data \
    --no_diagnostics \
    2> >(grep -v " I " >&2)
fi

if false; then
export NAME=pancan_test_no_classifier_
python -m sae.training_program \
    -l mrna \
    -o $NAME \
    -p $SAVE_PATH \
    --max_cv_repeats=1 \
    --epochs=150 \
    --epochs_pre_ae=1 \
    --no_pre_train_ae \
    --epochs_pre_cl=1 \
    --no_pre_train_cl \
    --epochs_pre_sr=1 \
    --no_pre_train_sr \
    --epochs_pre_bd=1 \
    --no_pre_train_bd \
    --epochs_pre_bc=1 \
    --no_pre_train_bc \
    --epochs_pre_dr=1 \
    --no_pre_train_dr \
    --fast_adaptive_multitask_optimization \
    --famo_weight_decay=0.0003551639396349761 \
    --decoder_split_loss \
    --classifier_split_loss \
    --hidden_activation=selu \
    --minibatchsize=32 \
    --mse \
    --survival_model \
    --drug_response_model \
    --batch_correction \
    --batch_loss=wasserstein \
    --encoder_layers=226,106,16 \
    --decoder_layers=95,445 \
    --survival_model_layers=16 \
    --drug_response_model_layers=67 \
    --drug_response_model_drugwise_layers=6 \
    --batch_detector_layers=36,64 \
    --classifier_layers=105,31 \
    --deconfounder_layers_per_batch=5 \
    --batch_adversarial_gradient_penalty=1. \
    --deconfounder_norm_penalty=1. \
    --drug_response_model_output_activation=relu \
    --pre_learning_rate_ae=0.00008590623250698555 \
    --pre_learning_rate_cl=0.00008590623250698555 \
    --pre_learning_rate_sr=0.00008590623250698555 \
    --pre_learning_rate_bd=0.00008590623250698555 \
    --pre_learning_rate_bc=0.00008590623250698555 \
    --pre_learning_rate_dr=0.00008590623250698555 \
    --learning_rate=0.00008590623250698555 \
    --adversarial_learning_rate_multi=10. \
    --adversarial_pre_training_optimizer=adam \
    --batch_correction_pre_training_optimizer=adam \
    --adam_beta1=0. \
    --adam_beta2=0.9 \
    --reg_weights_type= \
    --reg_weights=0. \
    --reg_noise_sd=0. \
    --reg_dropout_rate_input=0.29179563810929937 \
    --reg_dropout_rate_autoencoder=0. \
    --reg_dropout_rate_survival=0.03478547730630027 \
    --reg_dropout_rate_classifier=0.21685112099827208 \
    --reg_dropout_rate_batch=0.033310129812936934 \
    --reg_dropout_rate_drug_response=0.08906233077974918 \
    --data_var_cutoff_rank=-1 \
    --serialize_data \
    --parallel \
    --nthreads=1 \
    --nthreads_interop=1 \
    --nprocess=5 \
    --ptimeout=0 \
    --gpu_memory=4800 \
    --cache_data \
    --no_diagnostics \
    2> >(grep -v " I " >&2)
fi





export NAME=pancan_test_





