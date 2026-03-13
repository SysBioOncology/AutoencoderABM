#%% Load libraries 
from ModelOptimization import loss_function_autoencoder, simulation, load_img
from CancerModel import CancerModel
import Parameters as P
from Logging import write_python_file
import numpy as np 
import pandas as pd 
import pickle 

# Import PySwarms
from pyswarms.single.global_best import GlobalBestPSO

import tensorflow as tf
 
#%% Prepare optimization 
init_pars = P.init_pars
to_optimize = ['IMpkill']
model_outputs = ['final_T', 'final_L'] # which outputs to keep track off 

# load encoder and define whether intermediate images are saved 
encoder = tf.keras.models.load_model('../data/autoencoder/encoder_tumoroid.keras')
save_img = False
n_samples = 5
all_samples = range(1,  n_samples+1) 

# set parameter boundaries, order corresponds with to_optimize 
x_min = np.array([0.001])
x_max = np.array([0.8]) 
bounds = (x_min, x_max)
    
experiment = 'M07'
inputFolder = 'tumoroid/' + experiment
outputFolder = 'tumoroid'
write_python_file("Parameters.py", outputFolder)


#%% Run optimization
for sample in all_samples:
    # define y train 
    X_train = np.zeros((1, 100, 100))
    X_train[0,] = load_img(f"../data/{inputFolder}/sample{sample}_train.png")
    y_train = encoder.predict(X_train).flatten()

    # load init configuration 
    filename = f"../data/{inputFolder}/sample{sample}_init.csv"
    init_config = pd.read_csv(filename)

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
    
    print(f"======================= START OPTIMIZATION: sample {sample} ==========================================")
    for i in range(P.optim_repeats):   
        kwargs['title'] = outputFolder + '/Sample' + str(sample) + '_Step1'
        
        optimizer = GlobalBestPSO(n_particles=P.n_particles, dimensions=len(to_optimize), 
                                  options=P.options, bounds=bounds)
        cost, pos = optimizer.optimize(loss_function_autoencoder, P.optim_steps, verbose=True, **kwargs)









