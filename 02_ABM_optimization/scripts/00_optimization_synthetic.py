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
to_optimize = ['IMpkill', 'TUpprol', 'IMrwalk'] # parameters to optimize with PSO
model_outputs = ['final_T', 'final_L'] # which outputs to keep track off 

# load encoder and define whether intermediate images are saved 
folder = 'synthetic'
encoder = tf.keras.models.load_model('../data/autoencoder/encoder_synthetic.keras')
save_img = False

# set parameter boundaries, order corresponds with to_optimize 
x_min = np.array([0.1, 0.1, 0])
x_max = np.array([0.8, 0.5, 1]) 
bounds = (x_min, x_max)
    
# load init configuration if applicable 
init_config_sample = 1 
filename = '../data/init_config/25x25_proceedings.csv'
init_config = pd.read_csv(filename)
init_config['x'] += 40
init_config['y'] += 43 
print(f"Initial #tumor: {init_config['cell'].value_counts().get('tumor', 0)}, initial lymphocyte: {init_config['cell'].value_counts().get('lymphocyte', 0)}")

outputFolder = 'synthetic'
startSample = 0
nSamples = 1

# if only optimize selected samples (e.g. specific TUpprol values)
#samples_df = pd.read_csv('../data/synthetic_proceedings_specific_range.csv')
#TUpprol = 0.5
#init_pars['TUpprol'] = TUpprol
#selected_samples = list(samples_df[samples_df['TUpprol']==TUpprol]['IDX'])
selected_samples = []

write_python_file("Parameters.py", outputFolder)
temp_par = pd.read_csv('../data/synthetic/synthetic_proceedings.csv')


#%% Run optimization
for sample in range(startSample, nSamples):
    print(f'START!! Sample {sample} has taken off ==============================')
    
    if len(selected_samples) > 0:
        sample = selected_samples[sample]
        print(f">> Correction, sample {sample} :)")
        print(f">> TUpprol value of {temp_par.iloc[sample,][['TUpprol']].values[0]}")
        print(f">> IMpkill value of {temp_par.iloc[sample,][['IMpkill']].values[0]}")
        print(f">> Check: TUpprol = {init_pars['TUpprol']}")
    
    # define y train 
    X_train = np.zeros((1, 100, 100))
    X_train[0,] = load_img('../data/'+folder+'/'+str(sample)+'_1_TumorCells_ImmuneCells.png')
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
    
    all_par_values = []
    all_cost = []
    
    print(f"======== START OPTIMIZATION: {P.optim_repeats} repeats, {P.sim_iter} iterations, {P.nSteps} steps per simulation ============")
    # optimize ABM parameters
    for i in range(P.optim_repeats):   
        kwargs['title'] = outputFolder + '/Step1'
        
        # use PSO for parameter estimation 
        optimizer = GlobalBestPSO(n_particles=P.n_particles, dimensions=len(to_optimize), 
                                  options=P.options, bounds=bounds)
        cost, pos = optimizer.optimize(loss_function_autoencoder, P.optim_steps, verbose=True, **kwargs)

