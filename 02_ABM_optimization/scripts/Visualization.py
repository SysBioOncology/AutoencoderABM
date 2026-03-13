# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 10:13:43 2024

@author: 20192020
"""
import Parameters as P
from Parameters import np, sns, plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import ListedColormap

def get_current_occupation(model, save = False, annot_id = True, 
                           key={'tumor':1, 'lymphocyte':2}, plot_title='', 
                           outputFolder = 'test_runs', plot_label = True):
    """Visualizes grid of model and shows which coordinates are occupied"""
    # create arrays used for plotting of the cells 
    labels = np.full((model.grid.width, model.grid.height), '', dtype='U50')
    agent_counts =  np.zeros((model.grid.width, model.grid.height))
    for cell_content, (x, y) in model.grid.coord_iter():
        if cell_content != None:
            
            if annot_id:
                agent_count = cell_content.unique_id # get unique ID
                labels[x][y] = agent_count
            else:
                # if tumor cell is stem cell 
                if cell_content.type == 'tumor' and cell_content.is_stem:
                    labels[x][y] = '+'
                # if T cell is exhausted
                elif cell_content.type == 'lymphocyte' and (
                        cell_content.n_killed >= cell_content.IMkmax):
                    labels[x][y] ='/'
                # if cell is engaged
                elif cell_content.engaged:
                    labels[x][y] = 'o'
                else:
                    labels[x][y] = ''
                
            agent_counts[x][y] = key[cell_content.type]
            
            
    # Plot using seaborn, with a size of 5x5
    plt.figure()
    n = len(key) + 1
    
    cmap_dict = {0: '#F2F2F2', 1: '#A75443', 2: '#54BAC2'}
    cmap = ListedColormap([cmap_dict[i] for i in range(len(cmap_dict))])

    if plot_label:
        g = sns.heatmap(agent_counts, cmap=cmap, annot=labels, cbar=True,
                    square=True, annot_kws={"fontsize":P.fontsize}, fmt='', cbar_kws={"shrink": 0.5},
                    vmin=0, vmax=len(key) + 1)
    else:
        g = sns.heatmap(agent_counts, cmap=cmap, cbar=True,
                    square=True, annot_kws={"fontsize":P.fontsize}, fmt='', cbar_kws={"shrink": 0.5},
                    vmin=0, vmax=len(key) + 1)
    
    colorbar = g.collections[0].colorbar
    r = colorbar.vmax - colorbar.vmin
    colorbar.set_ticks([colorbar.vmin + 0.5 * r / (n) + r * i / (n) for i in range(n)])
    colorbar.set_ticklabels(['Unoccupied', 'Tumor', 'Lymphocyte'])
    
    # configure panel border and background 
    g.axhline(y=0, color='k',linewidth=2)
    g.axhline(y=agent_counts.shape[1], color='k',linewidth=2)
    g.axvline(x=0, color='k',linewidth=2)
    g.axvline(x=agent_counts.shape[0], color='k',linewidth=2)
    
    g.figure.set_size_inches(8, 8)
    g.set(xticks=[], yticks=[])
    
    if plot_label:
        g.set(title=f"Agents on grid ({model.grid.height}x{model.grid.width}) after step {model.nSteps} \n o = engaged cell, + = tumor stem cell, / = exhausted lymphocyte");
    else:
        g.set(title=f"Agents on grid ({model.grid.height}x{model.grid.width}) after step {model.nSteps}");
    
    if save:
        g.figure.savefig(f'../output/{outputFolder}/{plot_title}.png',dpi=600)
        
    return agent_counts, labels


