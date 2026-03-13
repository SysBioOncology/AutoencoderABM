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

## 00_prep
First

