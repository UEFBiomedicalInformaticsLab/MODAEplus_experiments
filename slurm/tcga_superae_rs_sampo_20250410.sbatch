#!/bin/bash
#SBATCH --job-name=TCGA_SUPERAE
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=5
##SBATCH --partition=gpu
#SBATCH --partition=small
#SBATCH --array=101-400%50
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

#module load tensorflow/2.13
#source $HOME/tensorflow_23/bin/activate
source /userapps/rintala/tensorflow/bin/activate
export TF_CPP_MIN_LOG_LEVEL=2
#export PROG_PATH=/projappl/project_2005873/superAE/src/sae/training_program.py
#export PROG_PATH=$HOME/superAE/src/sae/training_program.py
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
    --pre_learning_rate_ae=0.00001-0.0001 \
    --pre_learning_rate_cl=0.00001-0.0001 \
    --pre_learning_rate_sr=0.00001-0.0001 \
    --pre_learning_rate_bd=0.00001-0.001 \
    --pre_learning_rate_bc=0.00001-0.0001 \
    --pre_learning_rate_dr=0.0001-0.001 \
    --learning_rate=0.000001-0.0001 \
    --adversarial_learning_rate_multi=10. \
    --adversarial_pre_training_optimizer=adam \
    --batch_correction_pre_training_optimizer=adam \
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
    --gpu_memory=4800 \
    --cache_data \
    --no_diagnostics \
    2> >(grep -v " I " >&2)

