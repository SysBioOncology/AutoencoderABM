import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import seaborn as sns
from patchify import patchify
import cv2

def convert_spotlight2abm(val_pred, samples, x_limits, y_limits, 
                          size_new_coordinates, threshold_prob, key, limits, 
                          outputFolder, take_half):
    """
    Convert spotlight outputs to ABM initial configuration coordinates

    Parameters
    ----------
    val_pred : TYPE
        DESCRIPTION.
    samples : TYPE
        DESCRIPTION.
    x_limits : TYPE
        DESCRIPTION.
    y_limits : TYPE
        DESCRIPTION.
    size_new_coordinates : TYPE
        DESCRIPTION.
    threshold_prob : TYPE
        DESCRIPTION.
    key : TYPE
        DESCRIPTION.
    limits : TYPE
        DESCRIPTION.
    outputFolder : TYPE
        DESCRIPTION.
    take_half : TYPE
        DESCRIPTION.

    Returns
    -------
    tile_size_x : TYPE
        DESCRIPTION.
    tile_size_y : TYPE
        DESCRIPTION.

    """
    tile_size_x = dict()
    tile_size_y = dict()
    
    # go over all images 
    for i in range(0, len(samples)):
        print(f">> Start sample {i} =====================================")
        
        # 1) get small section of spotlight output 
        single_id = samples[i]
        plot_test_data = val_pred[val_pred['slide_submitter_id']==single_id]
        
        # crop image if required 
        if limits:
            plot_test_data = plot_test_data[plot_test_data['Coord_X'] < x_limits[1]]
            plot_test_data = plot_test_data[plot_test_data['Coord_X'] > x_limits[0]]
            plot_test_data = plot_test_data[plot_test_data['Coord_Y'] < y_limits[1]]
            plot_test_data = plot_test_data[plot_test_data['Coord_Y'] > y_limits[0]]
        elif take_half:
            xlimit = round(max(plot_test_data['Coord_X']) / 2)
            plot_test_data = plot_test_data[plot_test_data['Coord_X'] < xlimit]
            
        # 2) normalize tile coordinates and map them to new grid coordinates
        if plot_test_data.index.size > 0:            
            plot_test_data['Coord_X'] = round((plot_test_data['Coord_X'] - min(plot_test_data['Coord_X'])) / (
                max(plot_test_data['Coord_X']) - min(plot_test_data['Coord_X'])) * size_new_coordinates)
            plot_test_data['Coord_Y'] = round((plot_test_data['Coord_Y'] - min(plot_test_data['Coord_Y'])) / (
                max(plot_test_data['Coord_Y']) - min(plot_test_data['Coord_Y'])) * size_new_coordinates)
            
            # remove rows if there are NaN values 
            plot_test_data = plot_test_data.dropna() 
            
            # get spacing between tiles 
            all_x = plot_test_data['Coord_X'].drop_duplicates().sort_values()
            all_y = plot_test_data['Coord_Y'].drop_duplicates().sort_values()
            tile_size_x[single_id] = all_x.diff()
            tile_size_y[single_id] = all_y.diff()
            
            # 3) retrieve cell type information per tile 
            plot_test_array =  np.zeros((size_new_coordinates+1, size_new_coordinates+1)) # values for visualization
            probabilities = np.zeros((size_new_coordinates+1, size_new_coordinates+1)) # actual probability values
            df_x, df_y, df_cell, df_stem = [], [], [], [] # lists that will be columns in output for ABM
            
            for row in range(0, plot_test_data.index.size): 
                # get (transformed) x and y coordinates
                one_x = int(plot_test_data.iloc[row,]['Coord_X'])
                one_y = int(plot_test_data.iloc[row,]['Coord_Y'])
                
                # determine what type of cell we're dealing with 
                tumor = plot_test_data.iloc[row,]['tumor_purity']
                lymphocytes = plot_test_data.iloc[row,]['T_cells']
                if tumor > threshold_prob or lymphocytes > threshold_prob:
                    df_x.append(one_y)
                    df_y.append(one_x) # x and y need to be switched when giving to ABM
                    df_stem.append(0) # no data whether something is a stem cell so a standard value of 0 
                    
                    if tumor >= lymphocytes:
                        label = 'tumor' 
                        value = key[label]
                        prob = tumor
                        df_cell.append(label)
                        
                    elif lymphocytes > tumor and lymphocytes > threshold_prob:
                        label = 'lymphocyte' 
                        value = key[label]
                        prob = lymphocytes
                        df_cell.append(label)
                        
                else: # probability not high enough to assign cell type to this cell 
                    value = 0
                    prob = 0
                
                # determine area that will get same label and probability value 
                fill_area_x = round(tile_size_x[single_id].mean(numeric_only=True) / 2)
                fill_area_y = round(tile_size_y[single_id].mean(numeric_only=True) / 2)
                
                fill_x_min = max(one_x - fill_area_x, 0)
                fill_x_max = min(size_new_coordinates+1, one_x + fill_area_x) + 1
                fill_y_min = max(one_y - fill_area_y, 0)
                fill_y_max = min(size_new_coordinates+1, one_y + fill_area_y) + 1
                
                # update variables 
                plot_test_array[fill_x_min:fill_x_max,fill_y_min:fill_y_max] = value
                probabilities[fill_x_min:fill_x_max,fill_y_min:fill_y_max] = prob
            
            # save CSV that can be used to create initial configuration in the ABM
            abm_output = pd.DataFrame({'cell':df_cell, 'x':df_x, 'y':df_y, 'stem':df_stem})
            abm_output.to_csv(f'{outputFolder}/{single_id}.csv', index=False)
            
            # 4) plot visualization of ABM initial configuration 
            plt.figure()
            n = len(key) + 1
    
            cmap_dict = {0: '#F2F2F2', 1: '#A75443', 2: '#54BAC2'}
            cmap = ListedColormap([cmap_dict[i] for i in range(len(cmap_dict))])
    
            g = sns.heatmap(plot_test_array, cmap=cmap, cbar=True,
                        square=True, fmt='', cbar_kws={"shrink": 0.5},
                        vmin=0, vmax=len(key) + 1)
    
            colorbar = g.collections[0].colorbar
            r = colorbar.vmax - colorbar.vmin
            colorbar.set_ticks([colorbar.vmin + 0.5 * r / (n) + r * i / (n) for i in range(n)])
            colorbar.set_ticklabels(['Unoccupied', 'Tumor', 'Lymphocyte'])
    
            # configure panel border and background 
            g.axhline(y=0, color='k',linewidth=2)
            g.axhline(y=plot_test_array.shape[1], color='k',linewidth=2)
            g.axvline(x=0, color='k',linewidth=2)
            g.axvline(x=plot_test_array.shape[0], color='k',linewidth=2)
    
            g.set(title=f"id = {single_id}");
    
            g.figure.set_size_inches(8, 8)
            g.set(xticks=[], yticks=[])
            
            g.figure.savefig(f'{outputFolder}/{single_id}.png',dpi=600)
            plt.close('all')
            
        else:
            print(f'Given x and y limits resulted in an empty array for {single_id}!')
            
    return tile_size_x, tile_size_y

def get_patches(folder, filenames, patch_size, convert_dict, key, outputFolder, 
                threshold_ratio, new_size, patient_id):
    """
    Divide image into patches 

    Parameters
    ----------
    folder : TYPE
        DESCRIPTION.
    filenames : TYPE
        DESCRIPTION.
    patch_size : TYPE
        DESCRIPTION.
    convert_dict : TYPE
        DESCRIPTION.
    key : TYPE
        DESCRIPTION.
    outputFolder : TYPE
        DESCRIPTION.
    threshold_ratio : TYPE
        DESCRIPTION.
    new_size : TYPE
        DESCRIPTION.
    patient_id : TYPE
        DESCRIPTION.

    Returns
    -------
    chosen_filenames : TYPE
        DESCRIPTION.

    """
    chosen_filenames = {'file':[], 'filled_ratio':[], 'n_tumor':[], 'n_lymph':[]}
    
    for file in filenames:
        # read image and crop accordingly 
        imarray = cv2.imread(folder+file+'.png', cv2.IMREAD_GRAYSCALE)
        imarray = imarray[946:3902, 610:3566]
        patches = patchify(imarray, patch_size, patch_size[0])
        patch_id = 1
        
        # go over all patches created 
        for i in range(patches.shape[0]): 
            for j in range(patches.shape[1]): 
                patch = patches[i, j]
                
                # get some properties of this patch 
                unique, counts = np.unique(patch, return_counts=True)
                filled_ratio = (counts[unique==convert_dict['tumor']] + counts[unique==convert_dict['lymphocyte']]) / (sum(counts))
                
                # check if patch meets requirement 
                if filled_ratio > threshold_ratio and counts[unique==convert_dict['tumor']] > 50:
                    plot_data = patch.copy()
                    plot_data[plot_data==convert_dict['unoccupied']] = 0
                    plot_data[plot_data==convert_dict['tumor']] = key['tumor']
                    plot_data[plot_data==convert_dict['lymphocyte']] = key['lymphocyte']
                    
                    # plot visualization of ABM initial configuration 
                    plt.figure()
                    n = len(key) + 1
            
                    cmap_dict = {0: '#F2F2F2', 1: '#A75443', 2: '#54BAC2'}
                    cmap = ListedColormap([cmap_dict[i] for i in range(len(cmap_dict))])
            
                    g = sns.heatmap(plot_data, cmap=cmap, cbar=True,
                                square=True, fmt='', cbar_kws={"shrink": 0.5},
                                vmin=0, vmax=len(key) + 1)
            
                    colorbar = g.collections[0].colorbar
                    r = colorbar.vmax - colorbar.vmin
                    colorbar.set_ticks([colorbar.vmin + 0.5 * r / (n) + r * i / (n) for i in range(n)])
                    colorbar.set_ticklabels(['Unoccupied', 'Tumor', 'Lymphocyte'])
            
                    # configure panel border and background 
                    g.axhline(y=0, color='k',linewidth=2)
                    g.axhline(y=plot_data.shape[1], color='k',linewidth=2)
                    g.axvline(x=0, color='k',linewidth=2)
                    g.axvline(x=plot_data.shape[0], color='k',linewidth=2)
            
                    mfp = patient_id.loc[patient_id['patient']==file[:-11], ['MFP']].values[0]
                    g.set(title=f"id = {file}; MFP = {mfp}; filled ratio = {round(float(filled_ratio), 2)}");
            
                    g.figure.set_size_inches(8, 8)
                    g.set(xticks=[], yticks=[])
                    
                    g.figure.savefig(f'{outputFolder}/{file}_patch={patch_id}.png',dpi=600)
                    
                    # convert patch to dataframe for an output CSV
                    patch = cv2.resize(patch, dsize=new_size)
                    convert_to_df(patch, 
                                  convert_dict, f'{outputFolder}/{file}_patch={patch_id}.csv')
                    
                    # update output dictionary 
                    unique, counts = np.unique(patch, return_counts=True)
                    chosen_filenames['file'].append(f"{outputFolder}/{file}_patch={patch_id}")
                    chosen_filenames['filled_ratio'].append(filled_ratio)
                    chosen_filenames['n_tumor'].append(counts[unique==convert_dict['tumor']])
                    chosen_filenames['n_lymph'].append(counts[unique==convert_dict['lymphocyte']])
                    
                    patch_id += 1
                    
                    print(f">> Added new sample! We currently have {len(chosen_filenames['file'])} images!")
                    
                plt.close('all')
         
    pd.DataFrame(chosen_filenames).to_csv(f"{outputFolder}/00_filenames.csv", index=False)
    return chosen_filenames 
            
    
            
def convert_to_df(array, conv_dict, outputfile):      
    """
    Conversion image array to dataframe which can be used as initial 
    configuration of ABM

    Parameters
    ----------
    array : TYPE
        DESCRIPTION.
    conv_dict : TYPE
        DESCRIPTION.
    outputfile : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    output_dict = {'cell':[], 'x':[], 'y':[], 'stem':[]}
    for row in range(array.shape[0]):
        for col in range(array.shape[1]):
            value = array[row, col]
            
            if value == conv_dict['tumor']:
                output_dict['cell'].append('tumor')
                output_dict['stem'].append(0)
                output_dict['x'].append(row)
                output_dict['y'].append(col)
                
            elif value == conv_dict['lymphocyte']:
                output_dict['cell'].append('lymphocyte')
                output_dict['stem'].append(0)
                output_dict['x'].append(row)
                output_dict['y'].append(col)
                
    pd.DataFrame(output_dict).to_csv(outputfile, index=False)




















            
