# MODAE experiments

## Step 1 setup files (TODO: update URLs when published)

```
git clone https://github.com/trintala/MODAEplus.git
git clone https://github.com/trintala/MODAEplus_experiments.git
```

Download data from Zenodo (not uploaded at the moment, found in group folder).

## Step 2 build Apptainer container (optional)

```
cd MODAEplus_experiments/container
sudo apptainer build modae.sif modae.def
```
The image can also be found on Zenodo (group folder).

## Step 3 SLURM config

Check that files in slurm are configured appropriately, you may need to adjust 
image, data, script, and output paths. 

## Step 4 run hyper-parameter search

A setup script must be run to process the data: 
[setup.sbatch](slurm/hyperparameter_tuning/tcga_modae_setup_20250410.sbatch)

A SLURM array job can be used to run multiple iterations of random search: 
[rs.sbatch](slurm/hyperparameter_tuning/tcga_modae_rs_20250410.sbatch)

## Step 5 run additional hyper-parameter search evaluation

MODAE saves some metrics while running, but there are several additional metrics 
that haven't been integrated into the main code. 

See scripts in [slurm/hyperparameter_evaluation](slurm/hyperparameter_evaluation). 

## Step 6 identify the best setting

This experimental pipeline uses R to process and visualize MODAE results. 
[R/analysis/analyse_results.R](R/analysis/analyse_results.R) calculates summaries 
of evaluation metrics for all hyper-parameter settings and saves several figures 
to the output plot path. Adjusting [R/setup.R](R/setup.R) is necessary when 
working on multiple experiments. 

## Step 7 run final training and external validation

Adjust [python/final_training/final_training.py](python/final_training/final_training.py) 
and run [slurm/final_training/final_training.sbatch](slurm/final_training/final_training.sbatch) 
to train a final model and external validation. Some post-processing is necessary 
and can be done with [python/final_training/external_validation_post_processing.py](python/final_training/external_validation_post_processing.py). 

## Step 8 analyse final results

See various scripts in [R/analysis/](R/analysis/).
