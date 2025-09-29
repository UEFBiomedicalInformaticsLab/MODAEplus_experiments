## R-based analyses

These have been separated into different folders and scripts. 

### setup.R

This file defines many helpeor functions and sets many of the paths based 
on the environment variables MODAE_DATA_PATH and MODAE_OUTPUT_PATH. This 
works as long as the data is organized in the same way as in the group 
folder (Zenodo). 

## analysis

Various scripts for analysing the results. 

* analyse_results.R combines all the hyper-parameter tuning results and
  determines the best search iteration (a.k.a. "task"). 
* embedding_plots.R creates visualizations of the embeddings
  * See 
