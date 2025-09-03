#!/bin/bash
#SBATCH --job-name=CTRP_baseline
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=10
#SBATCH --partition=serial
#SBATCH --mem=30G
#SBATCH -t 3-00:00:00
#SBATCH --output=log/CTRP_baseline.%A.out
#SBATCH --error=log/CTRP_baseline.%A.err
##SBATCH --mail-user=teemu.rintala@uef.fi
##SBATCH --mail-type=NONE

apptainer exec \
  --bind /research/groups/fortino/teemu/MODAE_data:/data \
  --bind $HOME/MODAEplus_experiments:/scripts \
  ../container/modae.sif python \
  /scripts/python/baseline_models/baseline_drug_sensitivity_main.py \
    --root_dir=/data/drug_sensitivity_xia/ \
    --result_file=pca10_elasticnet_aac_identical.csv \
    --identical_data \
    --standardize \
    --use_pca \
    --pca_npc=10 \
    --gene_preselection \
    --model=elasticnet \
    --threads=10
