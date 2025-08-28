#!/bin/bash
#SBATCH --job-name=TCGA_SUPERAE_EVAL
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=96
#SBATCH --partition=small
#SBATCH --mem=400G
#SBATCH -t 3-00:00:00
#SBATCH --output=/scratch/users/rintala/log/TCGA_SUPERAE_EVAL.%A.out
#SBATCH --error=/scratch/users/rintala/log/TCGA_SUPERAE_EVAL.%A.err
#SBATCH --mail-user=teemu.rintala@uef.fi
#SBATCH --mail-type=NONE
#SBATCH --exclude=kudos12

#source $HOME/tensorflow_23/bin/activate
source /userapps/rintala/tensorflow/bin/activate
export TF_CPP_MIN_LOG_LEVEL=2

python post_learning_analysis/cancer_specific_drug_sensitivity_performance.py 2> >(grep -v " I " >&2)
