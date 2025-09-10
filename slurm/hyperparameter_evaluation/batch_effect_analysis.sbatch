#!/bin/bash
#SBATCH --job-name=MODAE_EVALUATION
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=96
#SBATCH --partition=small
#SBATCH --mem=100G
#SBATCH -t 3-00:00:00
#SBATCH --output=log/MODAE_EVALUATION.%A.out
#SBATCH --error=log/MODAE_EVALUATION.%A.err
##SBATCH --mail-user=
#SBATCH --mail-type=NONE

# Host paths
source path.conf

apptainer exec \
  --bind $MODAE_DATA_PATH:/data \
  --bind $MODAE_SCRIPT_PATH:/scripts \
  --bind $MODAE_OUTPUT_PATH:/output \
  ../container/modae.sif python \
  /scripts/python/hyperparameter_evaluation/batch_effect_analysis.py \
  2> >(grep -v " I " >&2)
