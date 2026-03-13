#%% Load libraries
import pandas as pd
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt 
import cv2
import random
from PIL import Image
from tumoroid_utils import read_plot_image, resize_image, map2inputdf

plt.close('all')

#%% Run for multiple files 
# indicate which images to map 
folder = "tumoroid/example_data"
idx = '07'
all_zstacks = [4, 7]
n_time_points = 5

# options for mapping
new_size = (100, 100)
interp = cv2.INTER_CUBIC
noise_boundary = 50
fill_middle = False

# variables to save conflict information 
all_conflict = []
all_no_conflict = []
conflict_ratio = np.zeros((len(all_zstacks), n_time_points))

for z in all_zstacks:
    z_idx = 0
    
    for t in range(1, n_time_points+1):
        # open files 
        im_tumor, imarray_tumor = read_plot_image(f"{folder}/20200925_spheroid_72h time lapsexy{idx}c1z{str(z).zfill(2)}t{str(t).zfill(2)}.tif",
                                                  'Tumor cell')
        im_Tcell, imarray_Tcell = read_plot_image(f"{folder}/20200925_spheroid_72h time lapsexy{idx}c2z{str(z).zfill(2)}t{str(t).zfill(2)}.tif",
                                                  'T cell')

        # resizing files 
        res_tumor = resize_image(imarray_tumor, new_size, interp, 
                                 "Tumor cell (resized from 512x512 to 100x100)")
        
        res_Tcell = resize_image(imarray_Tcell, new_size, interp, 
                                 "T cell (resized from 512x512 to 100x100)")
        plt.close("all")
        
        # convert map to input dataframe 
        output_filename = f"{folder}_mapped/z{str(z).zfill(2)}t{str(t).zfill(2)}.csv"
        conflict, no_conflict = map2inputdf(res_tumor, res_Tcell, noise_boundary, output_filename,
                                            fill_middle = fill_middle)
        
        # save the number of conflicts that occured 
        all_conflict.append(conflict)
        all_no_conflict.append(no_conflict)
        
        conflict_ratio[z_idx, t-1] = conflict / (no_conflict + conflict)
        
    z_idx += 1
        
        
#%% Plot conflict ratio 
plt.figure()
g = sns.heatmap(conflict_ratio, cbar=True, vmin=0)
g.set(title='Conflict ratio',
      xlabel='time points', ylabel='z-stacks');
g.set_yticklabels(all_zstacks);
        















