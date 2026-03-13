# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 09:48:33 2024

@author: 20192020
"""
#%% Load libraries 
from ModelOptimization import loss_function_autoencoder, simulation, load_img
from CancerModel import CancerModel
import Parameters as P
import numpy as np 
import pandas as pd 
import pickle 
from Logging import write_python_file

# Import PySwarms
from pyswarms.single.global_best import GlobalBestPSO
from pyswarms.single import LocalBestPSO

import tensorflow as tf
 
#%% Prepare optimization 
init_pars = P.init_pars
to_optimize = ['TUpprol', 'IMpkill', 'IMinfluxProb'] # parameters to optimize with PSO
model_outputs = ['final_T', 'final_L'] # which outputs to keep track off 

# load encoder and define whether intermediate images are saved 
folder = 'tcga'
all_patchnames = pd.read_csv('../data/' + folder + '/00_filenames_annots.csv').head(5) # only select 5 example images
encoder = tf.keras.models.load_model('../data/autoencoder/encoder_synthetic.keras')
save_img = False # save images of simulations done during optimization 

# set parameter boundaries, order corresponds with to_optimize 
x_min = np.array([0.1, 0.1, 0.1])
x_max = np.array([1, 1, 1]) 
bounds = (x_min, x_max)
    
# load init configuration if applicable 
filename = '../data/init_config/25x25_proceedings.csv'
init_config = pd.read_csv(filename)
init_config['x'] += 3 
init_config['y'] += 3 

outputFolder = 'tcga'
startSample = 0
nSamples = all_patchnames.shape[0]

write_python_file("Parameters.py", outputFolder)


#%% Run optimization
for sample in range(startSample, nSamples):
    # define y train 
    X_train = np.zeros((1, 100, 100))
    one_patch = all_patchnames.iloc[sample,][['file']].values[0]
    
    print(f'START!! Sample {sample} with name {one_patch} has taken off ========')
    X_train[0,] = load_img('../data/'+folder+'/'+one_patch+'.png',
                           x1 = 946, x2 = 3902, y1 = 610, y2 = 3566)
    y_train = encoder.predict(X_train).flatten()

    # collect everything for the optimization 
    kwargs={"y_train":y_train, 
            "init_pars":init_pars,
            "to_optimize":to_optimize,
            "model_outputs":model_outputs,
            "init_config":init_config,
            "encoder":encoder, 
            "save_img":save_img,
            "title":None,
            "outputFolder":outputFolder,
            "sample":sample}

    print(f"======== START OPTIMIZATION: {P.optim_repeats} repeats, {P.sim_iter} iterations, {P.nSteps} steps per simulation ============")
    # optimization 
    for i in range(P.optim_repeats):   
        kwargs['title'] = outputFolder + '/Step1'
        optimizer = GlobalBestPSO(n_particles=P.n_particles, dimensions=len(to_optimize), 
                                  options=P.options, bounds=bounds)
        cost, pos = optimizer.optimize(loss_function_autoencoder, P.optim_steps, verbose=True, **kwargs)
        
        
    










