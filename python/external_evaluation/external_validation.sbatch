#!/bin/bash
#SBATCH --job-name=TCGA_MODAE_EVAL
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=2
#SBATCH --partition=small
##SBATCH --gres=gpu:1
#SBATCH --mem=40G
#SBATCH -t 3-00:00:00
#SBATCH --output=log/TCGA_MODAE_EVAL.%A.out
#SBATCH --error=log/TCGA_MODAE_EVAL.%A.err
#SBATCH --mail-user=teemu.rintala@uef.fi
#SBATCH --mail-type=NONE

export TF_CPP_MIN_LOG_LEVEL=2
apptainer execute modae.sif python \
    external_evaluation/external_validation.py \
    2> >(grep -v " I " >&2)
