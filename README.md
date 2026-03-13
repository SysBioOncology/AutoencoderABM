# AutoencoderABM
Repository to reproduce results and figures for proceedings submission DOI TO BE ADDED

This repository consists of four parts: mapping spatial data to ABM grid, optimization of autoencoder, optimization of ABM and reproducing proceedings figures. 

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

## Data availability 
Note that some example files are present in the repository, although the data used for the proceedings need to be downloaded from Zenodo [ADD LINK]. Zenodo includes the following files:
* original tumoroid and tcga images
* mapped tumoroid and tcga images (both CSV file for ABM input and image)
* optimized autoencoders and encoders for three data sets: synthetic, tumoroid and TCGA
* simulation images with optimized ABM parameters

For each of the sections below, there will be an indication if additional data is needed from Zenodo to reproduce the proceedings results. 

## 00_prep
*Download Zenodo folder "00_original/TCGA". These are the original TCGA images* 

This folder contains the scripts to convert TCGA slides, using [SpoTLighT](https://pubmed.ncbi.nlm.nih.gov/39511301/) outputs, to ABM grid, or tumoroid microscopy images to ABM grid. There are two main scripts (`tcga_patches.py` and `tumoroid_conversion_pipeline.py`) with each a separate utils script with additional functions. 

## 01_autoencoder_optimization
*See Zenodo folders "01_mapped/TCGA" and "01_mapped/tumoroid" for all the mapped images. See Zenodo folder "02_autoencoder" for the autoencoder and encoder models for all three data sets (synthetic, tumoroid and TCGA)*

A general script (`autoencoder.py`)

