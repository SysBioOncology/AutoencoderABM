###########################%% TUMOROID %%###################################
#%% Get all libraries we need 

# agent based modeling 
import mesa 

# misc. libraries 
import numpy as np 
from collections import deque
import random
import math 

# visualization 
import seaborn as sns 
import matplotlib.pyplot as plt

#%% General parameters 
oneStepDuration = 1 # hours
height = 100
width = 100
fontsize = 6 # fontsize of annotations on simulation plots 
N_dict = {'tumor':10, 'lymphocytes':2} # starting number of cells if random initialization

nSteps = 72

#%% Optimization parameters
optim_repeats = 1 # number of times to do optimization 
sim_iter = 1 # number of times to perform simulation with same parameter set (mean is used for loss function)

# parameters for PSO 
optim_steps = 3 #30 # number of PSO iterations in one optimization  
options = {'c1': 2, 'c2': 1, 'w': 0.5} 
n_particles = 2

# values indicating different cell types, used during autoencoder training
tumor_value = 106 / 255
lymph_value = 156 / 255
unoccupied_value = 242 / 255

# size of encoder input, used to resize image if grid height / width doesn't match
encoder_input = 100

# thresholds to remove blurriness if image was upscaled 
thresholdA = 0.55
thresholdB = 0.65

#%% Simulation of therapy 
# each entry in the list is a therapy 
therapy_administration = None #[10, 10] # simulation step at which therapy is administered  
therapy_effect = ['lymphocyte', 'system'] # indication of which elements of the model the therapy affects (lymphocyte, system, tumor)
therapy_par = [{'IMpkill':20}, {'IMinfluxProb':10, 'IMinfluxRate':3}] # parameters that are affected with corresponding factor to change value (par_name : factor)


#%% Tumor cell parameters 
TUpprol = 0.15 #0.5055 # leiden = 0.21
TUpmig = 0.0005 # default 0.35, leiden = 0.05 
TUpdeath = 0.05 # 0.1216 for Leiden, 0.3 default
TUpmax = 10 
TUdamageThresh = 2 
TUps = 0
TUstem = 0 
tumor_pars = ['TUpprol', 'TUpmig', 'TUpdeath', 'TUpmax', 'TUdamageThresh', 'TUps', 'TUpmut', 'TUstem']


#%% Lymphocytes parameters 
IMkmax = 15
IMpmax = 5
IMpmig = 0.9
IMpkill = 0.9
IMpprol = 4.6289e-04
IMpdeath = 0.2 #1.5173e-04
IMinfluxProb = 0.2 #0.1 for directed migration, 0.33 default 
IMinfluxRate = 1
IMrateDynamic = 0.0075
IMrwalk = 0
engagementDuration = 48
IMinfluxEdge = False # influx only at edge if True, else random through entire grid
IMdirectedWidth = 0

# decay function for directed migration 
IMdecayWeight = 0
IMdecay = np.array([math.exp(-IMdecayWeight * xi) for xi in list(range(0, height))])

lymph_pars = ['IMkmax', 'IMpmax', 'IMpmig', 'IMpkill', 'IMrwalk', 'IMspeed', 'IMpprol', 'IMpdeath', 'IMrwalk',
              'IMinfluxEdge', 'IMdirectedWidth', 'IMdecayWeight', 'IMdecay']
system_pars = ['IMinfluxProb', 'IMinfluxRate', 'IMrateDynamic', 'engagementDuration']


#%% Set initial parameter values 
init_pars = {'TUpprol':TUpprol, 'TUpmig':TUpmig, 
              'TUpdeath':TUpdeath, 'TUpmax':TUpmax,
              'TUps':TUps, 'IMkmax':IMkmax, 'IMpmax':IMpmax,
              'IMrwalk':IMrwalk,
              'IMpmig':IMpmig, 'IMpkill':IMpkill,
              'IMpprol':IMpprol, 'IMpdeath':IMpdeath,
              'IMinfluxRate':IMinfluxRate, 
              'IMrateDynamic':IMrateDynamic,
              'IMinfluxProb':IMinfluxProb,
              'IMinfluxEdge':IMinfluxEdge,
              'IMdirectedWidth':IMdirectedWidth,
              'IMdecayWeight':IMdecayWeight,
              'IMdecay':IMdecay}
