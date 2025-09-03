#!/bin/bash
#SBATCH --job-name=CTRP_feature_selection
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=20
#SBATCH --partition=serial
#SBATCH --mem=100G
#SBATCH -t 3-00:00:00
#SBATCH --output=log/CTRP_feature_selection.%A.out
#SBATCH --error=log/CTRP_feature_selection.%A.err
##SBATCH --mail-user=teemu.rintala@uef.fi
##SBATCH --mail-type=NONE

apptainer exec \
  --bind /research/groups/fortino/teemu/MODAE_data:/data \
  --bind $HOME/MODAEplus_experiments:/scripts \
  ../container/modae.sif python \
  /scripts/python/baseline_models/baseline_drug_sensitivity_main.py \
    --root_dir=/research/work/rintala/ctrp/ \
    --result_file=forward_selection_all_drugs.csv \
    --identical_data \
    --threads=20 \
    --feature_selection_only \
    --multivariate_fs
