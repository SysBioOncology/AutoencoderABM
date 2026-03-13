import Parameters as P
from CancerModel import CancerModel
from Visualization import get_current_occupation
import numpy as np
import os 
import pandas as pd
import matplotlib.pyplot as plt 
import cv2
       

def loss_function(x, y_train, y_train_sd, init_pars, to_optimize, 
                  model_outputs, init_config): 
    """
    Calculate loss function 

    Parameters
    ----------
    x : numpy array 
        Value of parameters to optimize 
    y_train, y_train_sd : numpy array 
        Training values for model outputs 
    init_pars : dictionary 
        Dictionary to keep track of every model parameter and current value 
    to_optimize : list
        List of parameters to optimize
    model_outputs : list
        List of model outputs to use for optimization and given in y_train
    init_config : dataframe
        Dataframe for set initial configuration of model, can also be randomized

    Returns
    -------
    Loss function value for PSO 

    """
    f_tot = 0 
    
    # go over each particle
    for i in range(np.shape(x)[0]):
        
        # adjust parameters before feeding it to the model 
        pars = init_pars.copy()
        for idx in range(len(to_optimize)):
            # get parameter value 
            pars[to_optimize[idx]] = x[i, idx]
        
        all_properties = {x:[] for x in model_outputs}
        
        for j in range(P.sim_iter):
            model = CancerModel(P.N_dict, P.width, P.height, pars, verbose=0,
                                init_config = init_config)
            
            # change this dictionary if other functions are needed for output calculation
            output_functions = {'exhaust_L':model.sim_results.calculate_exhaust_lymp_ratio,
                                'killed_T':model.sim_results.calculate_killed_tumor_ratio,
                                'total_T':model.sim_results.get_total_tumor,
                                'total_L':model.sim_results.get_total_lymp,
                                'final_T':model.sim_results.get_final_tumor,
                                'final_L':model.sim_results.get_final_lymp,
                                'ratio_LT':model.sim_results.calculate_lymp_tumor_ratio,
                                'ratio_stemT':model.sim_results.calculate_stem_tumor_ratio
                                }
            
            # run simulation for certain number of steps 
            for step in range(P.nSteps):            
                # check if there are still tumor cells
                if model.N_tumor > 0:
                    model.step()
                    
                else:
                    print('Simulation stopped, all tumor cells are gone!')
                    break
                
            # collect model outputs 
            for key in all_properties.keys():
                all_properties[key].append(output_functions[key]())
        
        # save mean output value to use for optimization 
        for key in all_properties.keys():
            all_properties[key] = np.mean(all_properties[key])
        
        # predicted output value 
        temp_outputs = []
        for output in model_outputs:
            temp_outputs.append(all_properties[output])
        y_pred = np.array(temp_outputs)
           
        # calculate loss function
        f_tot += np.sum(((y_pred - y_train) / y_train_sd)**2)**0.5
        
    return(f_tot / np.shape(x)[0])



def simulation(model, filenames, filenames_Ncount, outputFolder, save = False,
               verbose = 0, nSteps = 20):
    """Run model simulation using optimized parameters"""
    filenames = []
    filenames_Ncount = []
    
    # run simulation for certain number of steps 
    for i in range(nSteps):
        # check if there are still tumor cells
        if model.N_tumor > 0:
            model.step()
            
            if verbose:
                print('--'*50)
                print(f'Step: {i}')
                
            plot_title = 'TumorCells_ImmuneCells_0000'
            filenames.append('../output/' + outputFolder + '/' + plot_title + '.png')
            
            plot_title_Ncounts = f'../output/{outputFolder}/NCounts_Figure_0000.png'
            filenames_Ncount.append(plot_title_Ncounts)
            
            if save:
                # plot initial configuration 
                agent_counts, labels = get_current_occupation(model, save = True, annot_id = False,
                                       plot_title = plot_title, outputFolder = outputFolder)
                
                # plot initial Ncounts
                model.sim_results.visualize_cell_counts(plot_title=plot_title_Ncounts)
                plt.close('all')
            
        else:
            print('Simulation stopped, all tumor cells are gone!')
            break
        
    return(model, filenames, filenames_Ncount)

def get_spatial_output(model):
    """Get array containing spatial configuration of current model state"""
    agent_map =  np.zeros((model.grid.width, model.grid.height))
    for cell_content, (x, y) in model.grid.coord_iter():
        if cell_content != None:
            # if tumor cell 
            if cell_content.type == 'tumor':
                agent_map[x][y] = P.tumor_value 
            # if T cell
            elif cell_content.type == 'lymphocyte':
                agent_map[x][y] = P.lymph_value 
        # if unoccupied
        else:
            agent_map[x][y] = P.unoccupied_value 
                
    return agent_map

def load_img (filename, x1 = 936, x2 = 3912, y1=600, y2=3576):    
    """Specific function written to load ABM images into array for encoder"""
    imarray = cv2.imread(filename, cv2.IMREAD_GRAYSCALE)
    imarray = imarray[x1:x2, y1:y2]
    imarray = cv2.resize(imarray, dsize=(P.encoder_input, P.encoder_input))
    imarray = imarray.astype(float) / 255
    
    return imarray

def loss_function_autoencoder(x, y_train, init_pars, to_optimize, 
                  model_outputs, init_config, encoder, save_img, title,
                  outputFolder, sample): 
    """
    Calculate loss function 

    Parameters
    ----------
    x : numpy array 
        Value of parameters to optimize 
    y_train : numpy array 
        Array with latent feature values for training image
    init_pars : dictionary 
        Dictionary to keep track of every model parameter and current value 
    to_optimize : list
        List of parameters to optimize
    model_outputs : list
        List of model outputs to use for optimization and given in y_train
    init_config : dataframe
        Dataframe for set initial configuration of model, can also be randomized
    encoder : Keras model
        Encoder layers from a trained autoencoder 
    save_img : Boolean
        Boolean to indicate whether to save images of intermediate optimization steps
    outputFolder : string
        If image is saved, this folder will be used for the file location
         

    Returns
    -------
    Loss function value for PSO with autoencoder integrated 

    """
    f_tot = 0 
    output_error = {'particle':[], 'sim_iter':[], 'error':[]}
    
    # add parameter values 
    for par in to_optimize:
        output_error[par] = []
    
    # go over each particle
    for i in range(np.shape(x)[0]):
        
        # adjust parameters before feeding it to the model 
        pars = init_pars.copy()
        for idx in range(len(to_optimize)):
            # get parameter value 
            pars[to_optimize[idx]] = x[i, idx]
        
        for j in range(P.sim_iter):
            model = CancerModel(P.N_dict, P.width, P.height, pars, verbose=0,
                                init_config = init_config)
            
            # run simulation for certain number of steps 
            for step in range(P.nSteps):            
                # check if there are still tumor cells
                if model.N_tumor > 0:
                    model.step()
                    
                else:
                    print('Simulation stopped, all tumor cells are gone!')
                    break
                
            # get encoded features 
            X_pred = np.zeros((1, P.width, P.height))
            X_pred[0,] = get_spatial_output(model)
            
            # resize image if ABM grid doesn't match with encoder input
            if X_pred.shape[0] != P.encoder_input:
                X_pred = cv2.resize(X_pred[0], dsize=(P.encoder_input, P.encoder_input), interpolation=cv2.INTER_NEAREST)
                
                # process image a bit to remove blurriness caused by upscaling resolution
                X_pred[X_pred >= P.thresholdB] = P.unoccupied_value * 255
                X_pred[X_pred < P.thresholdA] = P.tumor_value * 255
                X_pred[(X_pred < P.thresholdB) & (X_pred >= P.thresholdA)] = P.lymph_value * 255
                X_pred = X_pred.astype(float) / 255
                
                X_pred = X_pred.reshape((1, P.encoder_input, P.encoder_input))
                
            y_pred = encoder.predict(X_pred).flatten()
            
            current_error = np.sum((y_pred - y_train)**2)**0.5
            f_tot += current_error
            
            # save to dictionary which will be saved 
            output_error['particle'].append(i)
            output_error['sim_iter'].append(j)
            output_error['error'].append(current_error)
            
            for idx in range(len(to_optimize)):
                par = to_optimize[idx]
                output_error[par].append(x[i, idx])
            
    # save one image from this optimization step 
    # get suitable plot title to prevent overiding previous image 
    if save_img:
        current_step = 1
        while os.path.exists(f"../output/{outputFolder}/Sample{sample}_step{current_step}.png"):
            current_step += 1 
        img_title = f"Sample{sample}_step{current_step}"
        
        get_current_occupation(model, save = True, annot_id = False, 
                                   plot_title = img_title, outputFolder = outputFolder,
                                   plot_label = False)
        plt.close('all')
    
    # save error value from this step
    current_step = 1
    

    while os.path.exists(f"../output/{outputFolder}/Sample{sample}_step{current_step}_error.csv"):
        current_step += 1 
    
    pd.DataFrame.from_dict(output_error).to_csv(f"../output/{outputFolder}/Sample{sample}_step{current_step}_error.csv",
                                                index = False)


    return(f_tot / (np.shape(x)[0] * P.sim_iter))







