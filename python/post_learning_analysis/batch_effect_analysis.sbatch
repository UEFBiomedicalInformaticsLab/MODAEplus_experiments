#!/bin/bash
#SBATCH --job-name=TCGA_SUPERAE_EVAL
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=96
#SBATCH --partition=small
##SBATCH --partition=gpu
##SBATCH --gres=gpu:1
#SBATCH --mem=100G
#SBATCH -t 3-00:00:00
#SBATCH --output=/scratch/users/rintala/log/TCGA_SUPERAE_EVAL.%A.out
#SBATCH --error=/scratch/users/rintala/log/TCGA_SUPERAE_EVAL.%A.err
#SBATCH --mail-user=teemu.rintala@uef.fi
#SBATCH --mail-type=NONE
#SBATCH --exclude=kudos12

#source $HOME/tensorflow_23/bin/activate
source /userapps/rintala/tensorflow/bin/activate
export TF_CPP_MIN_LOG_LEVEL=2

python post_learning_analysis/batch_effect_analysis.py 2> >(grep -v " I " >&2)
