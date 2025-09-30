# MODAE experiments

## Step 1 setup files and environment (TODO: update URLs when published)

```
git clone https://github.com/trintala/MODAEplus.git
git clone https://github.com/trintala/MODAEplus_experiments.git
```

Download data from Zenodo (not uploaded at the moment, found in group folder).

Note that the scripts require environmental variables that define paths to 
data, scripts, and output. To keep the code portable and updateable, it is 
recommended to use these variables instead of editing the scripts directly. 
They can be setup like this:

```
export MODAE_DATA_PATH=$HOME/MODAE_data/
export MODAE_SCRIPT_PATH=$HOME/MODAEplus_experiments/
export MODAE_OUTPUT_PATH=$HOME/MODAE_output/
```

## Step 2 build Apptainer container (optional)

```
cd MODAEplus_experiments/container
sudo apptainer build modae.sif modae.def
```
The image can also be found on Zenodo (group folder).

## Step 3 SLURM config

Check that files in slurm are configured appropriately, adjust 
image, data, script, and output paths as necessary. The paths 
can be configured in [slurm/path.conf](slurm/path.conf) which defaults 
to subdirectories in $HOME. 

To execute the analyses in the following steps run the sbatch command in 
the slurm directory. 

## Step 4 run hyper-parameter search

A setup script must be run to process the data: 
[setup.sbatch](slurm/hyperparameter_tuning/tcga_modae_setup_20250410.sbatch)
* It processes the selected datasets.
* It splits the data into cross-validation folds for both parameter search and testing.
* It serializes the data to files to minimize the memory footprint.

A SLURM array job can be used to run multiple iterations of random search: 
[rs.sbatch](slurm/hyperparameter_tuning/tcga_modae_rs_20250410.sbatch)
* For each run it generates new random hyper-parameter settings.
* It runs CV once with 3-to-1 train-test for the parameter search and once with
  (3+1)-to-1 train-test for the testing. This corresponds to CV with a nested train-test
  split for hyper-parameter selection. 
* The CV is run using multiprocessing and can share one GPU between multiple processes.
* The model is relatively small and often CPU-based training is more cost-effective. 

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

See various scripts in [R/analysis/](R/analysis/) and [R/plotting/](R/plotting/). 
Futher information about the produced plots is listed in [R/README.md](R/README.md). 
