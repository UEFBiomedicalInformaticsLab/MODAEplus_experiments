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
#SBATCH --mail-user=teemu.rintala@uef.fi
#SBATCH --mail-type=NONE
#SBATCH --exclude=kudos[4,12]

apptainer exec \
  --bind /research/groups/fortino/teemu/MODAE_data:/data \
  --bind $HOME/MODAEplus_experiments:/scripts \
  ../container/modae.sif python \
  /scripts/python/baseline_models//baseline_overall_survival.py
