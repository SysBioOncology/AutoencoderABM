import pandas as pd
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt 
import cv2
import random
from PIL import Image

def read_plot_image(filepath, plot_title, plot = False):
    """Open and plot image, returns image"""
    
    im = Image.open(filepath)
    imarray = np.array(im)
    
    if plot:
        plt.figure()
        g = sns.heatmap(imarray)
        g.set(xticks=[], yticks=[])
        g.set(title=plot_title);
    
    return im, imarray

def resize_image(imarray, new_size, interp, plot_title, plot = False):
    """Resize image using specific interpolation approach"""
    res = cv2.resize(imarray, dsize=new_size, interpolation=interp)
    
    if plot:
        plt.figure()
        g = sns.heatmap(res)
        g.set(xticks=[], yticks=[])
        g.set(title=plot_title);
    
    return res

def map2inputdf(res_tumor, res_Tcell, noise_boundary, 
                output_filename, fill_middle = False):
    """
    Map tumoroid fluorescent images onto ABM grid 

    Parameters
    ----------
    res_tumor : array
        Image array of tumor cell signals
    res_Tcell : array
        Image array of T cell signals 
    noise_boundary : integer
        Threshold to remove background signal, only signals above this 
        value will be considered 
    output_filename : string
        Filepath to location where CSV with mapped image can be saved
    fill_middle : Boolean, optional
        Whether to use edge detection + morphology to fill in middle, 
        if there is an assumption that the tumor mass is solid while 
        the microscopy images don't show this due to technical limtiations. 
        The default is False.

    Returns
    -------
    conflict, no_conflict : integer
        Value representing how many times both tumor and T cell were mapped 
        to same grid cell (or when they didn't have this issue)
    """
    conflict = 0
    no_conflict = 0
    converted_array = np.zeros(res_tumor.shape, dtype="uint8")
    
    # go over every value within the array to assign tumor and T cells to 
    # correct grid cell 
    output_df = {'cell':[], 'x':[], 'y':[], 'stem':[]}
    for row in range(res_tumor.shape[0]):
        for col in range(res_tumor.shape[1]):
            # get signal value 
            value_tumor = res_tumor[row, col]
            value_Tcell = res_Tcell[row, col]
            
            # check if signal is above noise boundary 
            if value_tumor > noise_boundary or value_Tcell > noise_boundary:
                if value_tumor > value_Tcell:
                    cell = 'tumor'
                    converted_array[row, col] = 106
                elif value_tumor < value_Tcell:
                    cell = 'lymphocyte'
                    converted_array[row, col] = 156
                else: # if tumor and lymphocyte are equal 
                    cell = random.choice(['tumor', 'lymphocyte'])
                    if cell == 'tumor': 
                        converted_array[row, col] = 106
                    else: 
                        converted_array[row, col] = 156
                    
                # get an idea of how many times there's a 'conflict' 
                if value_tumor > noise_boundary and value_Tcell > noise_boundary:
                    conflict += 1 
                else:
                    no_conflict += 1 
        
                add = True
            else: # if no signal is found above noise boundary (= unoccupied grid cell)
                add = False 
                converted_array[row, col] = 242
                
            
            if add:
                stem = 0
                output_df['cell'].append(cell)
                output_df['x'].append(row)
                output_df['y'].append(col)
                output_df['stem'].append(stem)
    
    if fill_middle:
        # create mask of filled in area 
        mask = get_mask(converted_array)
        
        # get indices of mask 
        coordinates = np.argwhere(mask == 255)
        for row in coordinates:
            # check if empty 
            if converted_array[row[0], row[1]] == 242:
                # seed tumor cell 
                converted_array[row[0], row[1]] = 106
                
                output_df['cell'].append('tumor')
                output_df['x'].append(row[0])
                output_df['y'].append(row[1])
                output_df['stem'].append(0)
    
    # save mapped image as CSV 
    output_df = pd.DataFrame.from_dict(output_df)
    output_df.to_csv(output_filename, index=False)
    
    return conflict, no_conflict


def get_mask(converted_array):
    """Fills middle of a structure"""
    
    mask = np.zeros(converted_array.shape, dtype="uint8")
    mask[converted_array == 106] = 255
    edged = cv2.Canny(mask, 30, 200)
    
    # use morphology to close figure
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (35,35))
    morph = cv2.morphologyEx(edged, cv2.MORPH_CLOSE, kernel, )
    
    # find contours and bounding boxes
    mask = np.zeros_like(edged)
    contours = cv2.findContours(morph, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    contours = contours[0] if len(contours) == 2 else contours[1]
    for cntr in contours:
        cv2.drawContours(mask, [cntr], 0, 255, -1)
        
    return mask 