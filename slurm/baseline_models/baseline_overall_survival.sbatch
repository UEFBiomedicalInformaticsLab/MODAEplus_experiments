#!/bin/bash
#SBATCH --job-name=TCGA_baseline
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=40
#SBATCH --partition=small
#SBATCH --mem=100G
#SBATCH -t 3-00:00:00
#SBATCH --output=log/TCGA_baseline.%A.out
#SBATCH --error=log/TCGA_baseline.%A.err
##SBATCH --mail-user=
#SBATCH --mail-type=NONE

# Host paths
source path.conf

apptainer exec \
  --bind $MODAE_DATA_PATH:/data \
  --bind $MODAE_SCRIPT_PATH:/scripts \
  --bind $MODAE_OUTPUT_PATH:/output \
  ../container/modae.sif python \
  /scripts/python/baseline_models/baseline_overall_survival.py
