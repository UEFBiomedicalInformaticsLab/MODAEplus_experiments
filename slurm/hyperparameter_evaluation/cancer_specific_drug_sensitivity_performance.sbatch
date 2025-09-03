#!/bin/bash
#SBATCH --job-name=TCGA_SUPERAE_EVAL
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=96
#SBATCH --partition=small
#SBATCH --mem=400G
#SBATCH -t 3-00:00:00
#SBATCH --output=log/TCGA_SUPERAE_EVAL.%A.out
#SBATCH --error=log/TCGA_SUPERAE_EVAL.%A.err
##SBATCH --mail-user=teemu.rintala@uef.fi
##SBATCH --mail-type=NONE
#SBATCH --exclude=kudos12

apptainer exec \
  --bind /research/groups/fortino/teemu/MODAE_data:/data \
  --bind $HOME/MODAEplus_experiments:/scripts \
  modae.sif python \
  /scripts/python/hyperparameter_evaluation/cancer_specific_drug_sensitivity_performance.py \
  2> >(grep -v " I " >&2)
