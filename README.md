# AutoencoderABM
Repository to reproduce results and figures for [proceedings submission](https://doi.org/10.64898/2026.03.13.711699). This repository consists of four parts: mapping spatial data to ABM grid, optimization of autoencoder, optimization of ABM and reproducing proceedings figures. 

## Installation 
The library [mesa](https://mesa.readthedocs.io/latest/) is used to implement the ABM in Python. The following is also provided by them to install the latest stable release:

```ruby
pip install -U mesa
```

And to install the recommend dependencies: 

```ruby
pip install -U mesa[rec]
```

If you experience errors running the code, it could be that an older version is required. The following version has been tested to work for certain with the code: 

```ruby
pip install mesa==3.0.0a5 
```

Particle swarm optimization (PSO) is used for ABM parameter optimization, using the library [PySwarms](https://pyswarms.readthedocs.io/en/latest/installation.html). The library can be installed as followed:

```ruby
pip install pyswarms
```

For the TCGA specifically, patchify was used to divide the images into smaller sections:
```ruby
pip install patchify
```

## Data availability 
Note that some example files are present in the repository, although the data used for the proceedings need to be downloaded from [Zenodo](https://zenodo.org/records/19022344). Zenodo includes the following files:
* original tumoroid images
* mapped tumoroid and tcga images (both CSV file for ABM input and image), and generated synthetic images
* optimized autoencoders and encoders for three data sets: synthetic, tumoroid and TCGA
* simulation images with optimized ABM parameters

For each of the sections below, there will be an indication if additional data is needed from Zenodo to reproduce the proceedings results. 

## 00_prep
*Download Zenodo folders 00_original/M07 and 00_original/M10 for the original tumoroid images* 

This folder contains the scripts to convert TCGA slides, using [SpoTLighT](https://pubmed.ncbi.nlm.nih.gov/39511301/) outputs, to ABM grid, or tumoroid microscopy images to ABM grid. There are two main scripts (`tcga_patches.py` and `tumoroid_conversion_pipeline.py`) with each a separate utils script with additional functions. 

The tumoroid microscopy images are from [Liao et al.](https://pmc.ncbi.nlm.nih.gov/articles/instance/11322607/bin/42003_2024_6682_MOESM4_ESM.avi). Images used specifically for the proceedings are also added to Zenodo. The TCGA images are from the [TCGA SKCM](https://www.cancer.gov/tcga) cohort, with corresponding [gene expression data](https://gdac.broadinstitute.org). 

## 01_autoencoder_optimization
*See Zenodo folders "01_mapped/TCGA" for mapped TCGA patches. Mapped tumoroid of two experiments can be found in "01_mapped/data_M07" and "01_mapped/data_M10. The synthetic data are available in 01_mapped/synthetic_data*

A general script (`autoencoder.py`) to optimize loaded images. For each data set, a separate loading function and PCA visualization function are defined (`utils.py`).

## 02_ABM_optimization
*See Zenodo folder "02_autoencoder" for the autoencoder and encoder models for all three data sets (synthetic, tumoroid and TCGA)* 
The agent-based model was first adapted from [Kather et al.](https://pubmed.ncbi.nlm.nih.gov/28923860/) to Python, in which we also added additional processes focusing on the tumor-lymphocyte interactions. 

The following main scripts are present in this folder for the agent-based model:
- `Parameters.py` : contains model parameters and is the script in which you can change these values
- `Logging.py` : contains a function to write the Parameters.py into a textfile
- `CancerModel.py` : main class for the ABM, the script describes what needs to be done in each step of the simulation. For each cell in the simulation, a new class instance will be created (TumorCell.py or Lymphocyte.py)
- `TumorCell.py` : class for tumor cells
- `Lymphocyte.py` : class for lymphocytes
- `SimResults.py` : a class to save statistics at each step of the simulation (e.g. tumor cell count) and to visualize or calculate additional statistics 
- `Visualization.py` : contains a function to create a visualization of the cells within the simulation
- `ModelOptimization.py` : script with functions used for optimization, such as the user-defined loss function and a special simulation function to inspect final model
   
There are several main scripts that can be run, either for simulation only (e.g. with optimized parameters, `00_simulation.py`) or for optimization (one script for each data set). In the folder `parameters`, the contents of the script for each data set can be copied into the file with model parameters `Parameters.py` to reproduce the results from the proceedings. 

## 03_analysis
*See Zenodo folder "03_ABMoptimization" for optimization results of the ABM which can be further analyzed. Simulations using optimized parameters are also included.*

The data of each figure in the proceedings can be replicated using the R scripts and data provided in this folder. 



