# Load modules and libraries 
from CancerModel import CancerModel
from Visualization import get_current_occupation
from Logging import write_python_file, calculate_ci

from PIL import Image
import pandas as pd 
import pickle
import os 
import numpy as np
import sys
#from line_profiler import LineProfiler

# load script with parameter values
from Parameters import plt
import Parameters as P

#%% User-defined variables related to file names, logging etc.
# set output folder to save results 
outputFolder = 'tcga_simulation'

# indicate whether intermediate frames and / or final GIF need to be saved 
save_frames = True
save_gif = False
verbose = 0
    # 0 = only print sample + iteration + saving / deleting visuals
    # 1 = print intermediate step indicators
    # 2 = print information during simulation

# initial parameter values 
init_pars = P.init_pars
#par_values = pd.DataFrame([list(init_pars.values())], columns=list(init_pars.keys()))

# if specific parameter values need to be simulated 
par_values = pd.read_csv('../data/optimized/00_estimated_tcga.csv')
parameter_names = list(par_values.columns)

# save model parameter values 
write_python_file("Parameters.py", outputFolder)
 
# load init configuration if applicable 
filename = '../data/init_config/25x25_proceedings.csv'
init_config = pd.read_csv(filename)
#init_config = None

tumoroid = False
# move configuration to middle (TCGA)
init_config['x'] += 3 
init_config['y'] += 3  

# synthetic
#init_config['x'] += 40
#init_config['y'] += 43 

# tumoroid
experiment = 'M07'
#tumoroid = True

# saving visualizations 
if save_gif:
    filenames = {}
    filenames_Ncount = {}
    
# saving model outputs 
outputs = ['exhaust_L', 'killed_T', 'total_T', 'total_L', 'final_T', 'final_L', 'ratio_LT']
all_output_log = {x:[] for x in parameter_names + outputs}
# note that these outputs are calculated using specific functions in SimResults class
# user can adjust this in the next section 

# set number of iterations for this simulation (how often simulation with same 
    # parameter values is run)
nIter_start = 0
nIter = 5
row_start = 0
row_end = 1 #len(par_values.index)


#%% Simulation with predefined model parameter values 
all_models = []
idx = 0 
for row in range(row_start, row_end):
    print('')
    print(f'============ START : SAMPLE {row} ============')
    
    if tumoroid:
        filename = '../data/leiden/' + experiment + '/sample'+str(row+1)+'_init.csv' 
        init_config = pd.read_csv(filename)
    
    # set parameter values 
    for par in parameter_names:
        init_pars[par] = par_values.iloc[row,][[par]].values[0]
        print(f"{par}={par_values.iloc[row,][[par]].values[0]}")
    
    all_models.append([])
    for test_id in range(nIter_start, nIter):
        
        if verbose > 0:
            print(f'------------ : ITERATION {test_id} ------------')
        
        # create ABM 
        model = CancerModel(P.N_dict, P.width, P.height, init_pars, verbose=max(0, verbose-1),
                            init_config = init_config) 
        
        # dictionary of functions related to different types of outputs 
        output_functions = {'exhaust_L':model.sim_results.calculate_exhaust_lymp_ratio,
                            'killed_T':model.sim_results.calculate_killed_tumor_ratio,
                            'total_T':model.sim_results.get_total_tumor,
                            'total_L':model.sim_results.get_total_lymp,
                            'final_T':model.sim_results.get_final_tumor,
                            'final_L':model.sim_results.get_final_lymp,
                            'ratio_LT':model.sim_results.calculate_lymp_tumor_ratio,
                           'ratio_stemT':model.sim_results.calculate_stem_tumor_ratio
                            }

        # perform multiple simulations steps with teh model 
        for i in range(1, (P.nSteps+1)):
            # only continue if there are still tumor cells present 
            if model.N_tumor > 0:
                model.step()
                
                if verbose > 0:
                    print('--'*25)
                    print(f'------------ Step: {i} ------------ ')
                
            else:
                print('Simulation stopped, all tumor cells are gone!')
                break
            
        if save_frames or save_gif:
            # plot initial configuration 
            plot_title = f'{row}_{test_id}_TumorCells_ImmuneCells'
            agent_counts, labels = get_current_occupation(model, save = True, annot_id = False,
                                   plot_title = plot_title, outputFolder = outputFolder,
                                   plot_label = False)
            plt.close('all')
            
            if save_gif: filenames_Ncount[test_id].append(plot_title)
            
        # record model outputs (all)
        for par in parameter_names:
            all_output_log[par].append(init_pars[par])
            
        for output_idx in range(len(outputs)):
            temp = output_functions[outputs[output_idx]]()
            all_output_log[outputs[output_idx]].append(temp)
        
    output_csv = pd.DataFrame(all_output_log)
    output_csv.to_csv(f'../output/{outputFolder}/00_all_outputs_{row_start}_{row}.csv', index=False)
    
    # remove previously generated model outputs CSV files and saved models 
    if row > row_start:
        os.remove(f'../output/{outputFolder}/00_all_outputs_{row_start}_{row-1}.csv')
            
    idx += 1 
        
    
plt.close()

print('=============== DONE! ===============')






