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
#SBATCH --mail-user=teemu.rintala@uef.fi
#SBATCH --mail-type=NONE

source $HOME/tensorflow_23/bin/activate 
python baseline_models_dr/baseline_drug_sensitivity_main.py \
  --root_dir=/research/work/rintala/drug_response_dataset/ \
  --result_file=pca100_elasticnet_aac_identical.csv \
  --identical_data \
  --gene_preselection \
  --cl_exp_root=/research/users/rintala/ccle/ \
  --cl_exp_file=CCLE_expression.csv \
  --cl_dr_root=/research/users/rintala/ctrp/ \
  --cl_dr_file=screen_aac.csv.gz \
  --cl_dr_row_info=screen_rowinfo.csv.gz \
  --cl_dr_row_map=ctrp_ccle_clid_map.csv.gz \
  --pca_npc=100 \
  --model=elasticnet \
  --threads=10 \
  --use_pca #\
#  --drug_features=dragon_pfp 
