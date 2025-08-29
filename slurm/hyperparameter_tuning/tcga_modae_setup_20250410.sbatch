#!/bin/bash
#SBATCH --job-name=TCGA_SUPERAE_SETUP
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --partition=small
##SBATCH --partition=serial
##SBATCH --partition=gpu
#SBATCH --mem=10G
#SBATCH -t 0-00:15:00
#SBATCH --output=/scratch/users/rintala/log/TCGA_SUPERAE_SETUP.%A.out
#SBATCH --error=/scratch/users/rintala/log/TCGA_SUPERAE_SETUP.%A.err
#SBATCH --mail-user=teemu.rintala@uef.fi
#SBATCH --mail-type=NONE
##SBATCH --account=project_2005873
#SBATCH --exclude=kudos[12]

export PROJ_PATH=/research/work/rintala
export SAVE_PATH=$PROJ_PATH/superAE_HPO/20250410_random_search/

export SOLID_CANCER_LIST=ACC,BLCA,BRCA,CESC,CHOL,COAD,DLBC,ESCA,GBM,HNSC,KICH,KIRC,KIRP,LGG,LIHC,LUAD,LUSC,MESO,OV,PAAD,PCPG,PRAD,READ,SARC,SKCM,STAD,TGCT,THCA,THYM,UCEC,UCS,UVM

#export SAVE_PATH=${SAVE_PATH}lung_test/
#export NAME=lung_test_
#export CANCER_LIST=LUAD,LUSC
##export ONCOTREE_LEVEL=level_3
##export ONCOTREE_LIST=LUAD,LUSC
#export FILTER_FILE=Model_augmented.csv
#export FILTER_COLUMN=level_3
#export FILTER_INCLUDE=LUAD,LUSC

#export SAVE_PATH=${SAVE_PATH}brca_test/
#export NAME=brca_test_
#export CANCER_LIST=BRCA
##export ONCOTREE_LEVEL=level_2
##export ONCOTREE_LIST=BRCA
#export FILTER_FILE=Model_augmented.csv
#export FILTER_COLUMN=level_2
#export FILTER_INCLUDE=BRCA

#export SAVE_PATH=${SAVE_PATH}brca_test_noclfilter_pretrain/
#export NAME=brca_test_noclfilter_pretrain_
#export CANCER_LIST=BRCA
##export ONCOTREE_LEVEL=
##export ONCOTREE_LIST=
#export FILTER_FILE=Model_augmented.csv
#export FILTER_COLUMN=solid
#export FILTER_INCLUDE=solid

export SAVE_PATH=${SAVE_PATH}pancan_test/
export NAME=pancan_test_
export CANCER_LIST=$SOLID_CANCER_LIST
export SOLID_CANCER_LIST=
##export ONCOTREE_LEVEL=level_1
##export ONCOTREE_LIST=
export FILTER_FILE=Model_augmented.csv
export FILTER_COLUMN=solid
export FILTER_INCLUDE=solid

#export SAVE_PATH=${SAVE_PATH}scanb_test_noclfilter/
#export NAME=scanb_test_noclfilter_
#export CANCER_LIST=.
##export ONCOTREE_LEVEL=
##export ONCOTREE_LIST=
#export FILTER_FILE=Model_augmented.csv
#export FILTER_COLUMN=solid
#export FILTER_INCLUDE=solid

#module load tensorflow/2.13
#source $HOME/tensorflow_24/bin/activate
source /userapps/rintala/tensorflow/bin/activate
export TF_CPP_MIN_LOG_LEVEL=2
#export PROG_PATH=/projappl/project_2005873/superAE/src/sae/training_program.py
#export PROG_PATH=$HOME/superAE/src/sae/training_program.py
python -m sae.training_program \
    -l mrna \
    -o $NAME \
    -p $SAVE_PATH \
    --survival_model \
    --drug_response_model \
    --classifier \
    --data_var_cutoff_rank=-1 \
    --serialize_data \
    --patient_expression_root_dir=$PROJ_PATH/tcga/pan_cancer/data_full/ \
    --patient_expression_cancer_list=$CANCER_LIST \
    --patient_expression_file_name=mrna.csv.gz \
    --patient_expression_sample_file=sample_info.csv.gz \
    --patient_expression_sample_id_col=colname \
    --patient_expression_sample_type_col=type \
    --patient_expression_sample_type_include=1 \
    --patient_gene_mapping_file=TCGA_gene_mapping.csv \
    --patient_expression_mean_cut=1. \
    --patient_class_file=oncotree_level1.csv \
    --patient_class_col=level_1 \
    --patient_id_col=bcr_patient_barcode \
    --patient_redaction_col=Redaction \
    --patient_survival_file=survival.csv.gz \
    --patient_survival_event_col=OS \
    --patient_survival_time_col=OS.time \
    --patient_survival_covar_cols=age_at_initial_pathologic_diagnosis \
    --patient_survival_covar_onehot=0 \
    --cell_line_expression_root_dir=$PROJ_PATH/ccle/ \
    --cell_line_expression_file=CCLE_expression.csv \
    --cell_line_gene_mapping_file=CCLE_gene_mapping.csv \
    --cell_line_expression_mean_cut=1. \
    --cell_line_class_file=oncotree_level1.csv \
    --cell_line_class_col=level_1 \
    --cell_line_filter_mapping_file=$FILTER_FILE \
    --cell_line_filter_mapping_column=$FILTER_COLUMN \
    --cell_line_filter_inclusion_list=$FILTER_INCLUDE \
    --cell_line_drug_response_root_dir=$PROJ_PATH/drug_response_dataset/ \
    --cell_line_drug_response_file=xia_ctrp_screen_dss.csv.gz \
    --cell_line_drug_response_row_info_file=xia_ctrp_screen_info.csv.gz \
    --cell_line_drug_response_row_info_sample_col=CELL \
    --cell_line_drug_response_row_info_treat_col=DRUG \
    --cell_line_drug_response_row_map_file=xia_id_ccle_clid_map.csv.gz \
    --cell_line_drug_response_row_map_column=depmap_id \
    --cell_line_drug_response_target_column=DSS1 \
    --complete_data_shuffle_seed=0 \
    --data_standardize \
    --data_setup_only \
    --debug_data \
    2> >(grep -v " I " >&2)
