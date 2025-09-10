#!/bin/bash
#SBATCH --job-name=MODAE_FINAL
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=2
#SBATCH --partition=small
##SBATCH --gres=gpu:1
#SBATCH --mem=40G
#SBATCH -t 3-00:00:00
#SBATCH --output=log/MODAE_FINAL.%A.out
#SBATCH --error=log/MODAE_FINAL.%A.err
##SBATCH --mail-user=
#SBATCH --mail-type=NONE

# Host paths
source path.conf

apptainer exec \
  --bind $MODAE_DATA_PATH:/data \
  --bind $MODAE_SCRIPT_PATH:/scripts \
  --bind $MODAE_OUTPUT_PATH:/output \
  ../container/modae.sif python \
    /scripts/python/final_training/final_training.py \
    2> >(grep -v " I " >&2)
