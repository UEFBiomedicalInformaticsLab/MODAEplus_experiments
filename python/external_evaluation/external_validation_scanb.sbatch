#!/bin/bash
#SBATCH --job-name=TCGA_SUPERAE_EVAL
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=10G
#SBATCH -t 0-02:00:00
#SBATCH --output=log/TCGA_SUPERAE_EVAL.%A.out
#SBATCH --error=log/TCGA_SUPERAE_EVAL.%A.err
#SBATCH --mail-user=teemu.rintala@uef.fi
#SBATCH --mail-type=NONE

source $HOME/tensorflow_23/bin/activate
export TF_CPP_MIN_LOG_LEVEL=2

python external_evaluation/external_validation_scanb.py 2> >(grep -v " I " >&2)
