import pandas as pd
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt 

from utils import convert_spotlight2abm, get_patches

#%% Get SpoTlighT output 
val_pred = pd.read_csv('TCGA/tcga_validation_tile_predictions_proba.csv',
                       delimiter = '\t') 
patient_id = pd.read_csv('subset_tcga_skcm.csv') # only get desert or immune-enriched
val_pred = patient_id.merge(val_pred, how='left', left_on='patient', right_on='TCGA_patient_ID')

samples = val_pred['slide_submitter_id'].unique()
x_coords = val_pred['Coord_X']
y_coords = val_pred['Coord_Y']

#%% Prepare some settings for deriving initial configuration 
x_limits = (20000, 150000) 
y_limits = (2000, 50000) 
size_new_coordinates = 5000 # multiple of 10 to transform sparse matrix into a grid for ABM
threshold_prob = 0.5 # minimum probability for a tile to get the label of a specific cell type 
key={'tumor':1, 'lymphocyte':2}
outputFolder = "TCGA/example_data"
limits = True
take_half = True

#%% Convert spotlight to ABM
tile_size_x, tile_size_y = convert_spotlight2abm(val_pred, samples, x_limits, y_limits, 
                          size_new_coordinates, threshold_prob, key, limits,
                          outputFolder, take_half)

#%% Create patches
folder = 'TCGA/example_data' # folder with converted images
outputFolder = folder + '/patches' # folder to save patches 
patch_size = (1200, 1200)
convert_dict = {'tumor':106, 'lymphocyte':156, 'unoccupied':242}
key = {'tumor':1, 'lymphocyte':2}
filled_ratio = 0.8

chosen_filenames = get_patches(folder, samples, patch_size, convert_dict, key, outputFolder, filled_ratio, 
                  (100, 100), patient_id)





