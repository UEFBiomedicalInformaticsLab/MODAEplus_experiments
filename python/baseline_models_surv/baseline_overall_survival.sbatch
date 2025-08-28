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

#source $HOME/tensorflow_23/bin/activate 
source /userapps/rintala/tensorflow/bin/activate
python baseline_models_surv/baseline_overall_survival.py
