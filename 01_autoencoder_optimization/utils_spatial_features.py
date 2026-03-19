import cv2 
import numpy as np
import math 
import networkx as nx
import pandas as pd
import copy
from sklearn.feature_extraction.image import grid_to_graph

def calc_complexity(imarray, value2cell, int2cell):
    """
    Calculate complexity score of spherical structure 

    Parameters
    ----------
    imarray : array
        Image to be converted
    value2cell : dictionary
        Converting float values to string (which value corresponds with 
                                               which cell type)
    int2cell : dictionary
        Converting string value to integers (which cell type corresponds with 
                                             which integer)

    Returns
    -------
    area, perimeter, complexity : float
        Return area, perimeter and resulting complexity score for image
    """
    # prepare array with correct values 
    imarray = imarray.copy()
    for key in value2cell:
        imarray[imarray==float(key)] = int2cell[value2cell[key]]
    
    imarray = cv2.GaussianBlur(imarray,(3, 3),0)
    imarray[imarray>115] = 255
    
    # create mask for tumor cells
    mask = np.zeros(imarray.shape[:2], dtype="uint8")
    mask[imarray == 106] = 255

    # edge detection
    edged = cv2.Canny(mask, 30, 200)
    hh, ww = edged.shape[:2]
    
    # use morphology to close figure
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (35,35))
    morph = cv2.morphologyEx(edged, cv2.MORPH_CLOSE, kernel, )

    # find contours and bounding boxes
    mask = np.zeros_like(edged)
    contours = cv2.findContours(morph, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    contours = contours[0] if len(contours) == 2 else contours[1]

    # get area, perimeter and complexity (if possible from image)
    if len(contours) > 0:
        area = cv2.contourArea(contours[0])
        perimeter = cv2.arcLength(contours[0],True)
        
        if area > 0 and perimeter > 0:
            complexity = perimeter**2 / (4*math.pi*area)
        else:
            complexity = float("nan")
    else:
        area = float("nan")
        perimeter = float("nan")
        complexity = float("nan")
        
    return area, perimeter, complexity

def get_edges_df(imarray, value2cell):
    """
    Convert array of image into network, assigning edges if two cells are neighbors

    Parameters
    ----------
    imarray : array
        Image to be converted, consisting of only three values representing 
        cell types as provided in value2cell
    value2cell : dictionary
        Conversion of numerical values to cell type 

    Returns
    -------
    Dataframe with edges of the network as derived from the image and a dictionary
    pos with original positions in the ABM grid
    """
    output_dict = {0:[], 1:[]}
    pos = {}
    idx = 0
    
    #1) go over every cell 
    for x in range(imarray.shape[0]):
        neighbors_x = [x-1, x-1, x-1, x, x, x+1, x+1, x+1]
        
        for y in range(imarray.shape[1]):  
            if value2cell[str(imarray[x, y])] != 'unoccupied':
                #2) get neighbors 
                neighbors_y = [y-1, y, y+1, y-1, y+1, y-1, y, y+1]
                neighbors_idx = [idx-imarray.shape[1]-1, idx-imarray.shape[1], 
                                 idx-imarray.shape[1]+1, idx-1, idx+1, 
                                 idx+imarray.shape[1]-1, idx+imarray.shape[1],
                                 idx+imarray.shape[1]+1]
                
                #3) go over all neighbors and check if there is a cell present 
                for (xi, yi, idx_i) in zip(neighbors_x, neighbors_y, neighbors_idx):
                    
                    # use try-catch to make sure we're still within the array 
                    try:
                        if xi>0 and yi>0 and value2cell[str(imarray[xi, yi])] != 'unoccupied':
                            output_dict[0].append(idx)
                            output_dict[1].append(idx_i)
                            
                            # make sure we don't add pos of isolated node
                            # also mirror y as in the image, y=0 is at the top
                            pos[idx]=(y, abs(x-imarray.shape[0]+1))
                            pos[idx_i]=(yi, abs(xi-imarray.shape[0]+1))
                    except IndexError:
                        pass 
                        
            idx += 1 
                    
    return pd.DataFrame.from_dict(output_dict), pos

def calculate_direct_neighbors(G, cells):
    """
    Derive neighbors of each cell within the network 

    Parameters
    ----------
    G : networkx graph  
        Networkx graph representing spatial configuration of ABM grid
    cells : list
        List of cell types that are present in the network 

    Returns
    -------
    neighbors_dict : dictionary
        Dictionary that keeps track of all the neighbor counts of every agent 
        in the simulation. For this nested dictionary, the first keys are the 
        source node, the second the target node. So the list of values of e.g. 
        neighbors_dict['tumor']['lymphocyte'] represents how many lymphocytes 
        are next to tumor cells. 
    neighbors_dict_mean : dictionary
        Mean over all the values for each cell-cell interaction of the neighbors_dict

    """
    
    # get all completely connected components 
    components = nx.connected_components(G)

    # initialize output dictionary 
    neighbors_dict = {x:{y:[] for y in cells} for x in cells}

    # go over every component to get neighbors 
    for s in components:
        # go over each node 
        for n in s:
            source_type = G.nodes[n]["cell type"]
            neighbor_type = [G.nodes[x]["cell type"] for x in G.neighbors(n)]
            
            for cell_type in cells:
                neighbors_dict[source_type][cell_type].append(neighbor_type.count(cell_type))
            
    # get average
    neighbors_dict_mean = copy.deepcopy(neighbors_dict)

    for source_type in cells:
        for neighbor_type in cells:
            neighbors_dict_mean[source_type][neighbor_type] = np.mean(
                neighbors_dict_mean[source_type][neighbor_type])
    
    return neighbors_dict, neighbors_dict_mean


def get_neighbors(imarray, cells, value2cell):
    """
    Transform ABM grid into graph and calculate number of cell type-specific 
    neighbors of each agent type 

    Parameters
    ----------
    imarray : array
        Image for which neighbors are going to be calculated
    cells : list
        List of cell types present in the ABM 
    value2cell : dictionary
        Dictionary to map image value to cell type 

    Returns
    -------
    neighbors_dict : dictionary
        Dictionary that keeps track of all the neighbor counts of every agent 
        in the simulation. For this nested dictionary, the first keys are the 
        source node, the second the target node. So the list of values of e.g. 
        neighbors_dict['tumor']['lymphocyte'] represents how many lymphocytes 
        are next to tumor cells. 
    neighbors_dict_mean : dictionary
        Mean over all the values for each cell-cell interaction of the neighbors_dict

    """
    X_tumor = np.reshape(imarray, (-1, 1))
    connectivity = grid_to_graph(*imarray.shape)

    # remove unoccupied from connectivity 
    edges_df = pd.DataFrame(np.concatenate((connectivity.row, connectivity.col, connectivity.data)).reshape(3, -1).T)
    edges_df = edges_df[X_tumor[edges_df[0]] != list(value2cell.values())[0]] # remove links with unoccupied cells
    edges_df = edges_df[X_tumor[edges_df[1]] != list(value2cell.values())[0]]
    edges_df = edges_df[edges_df[0] != edges_df[1]] # remove self-loops

    # add idx of cells so can use for mapping to coordinates 
    edges_df = edges_df.assign(cell_id = edges_df[[0, 1]].agg(tuple, axis=1))
    
    # get edges
    edges_df, all_pos = get_edges_df(imarray, value2cell)
    
    # create graph 
    G = nx.from_pandas_edgelist(edges_df, 0, 1)
    
    # set node attributes 
    celltype_dict = {x:{'cell type':value2cell[str(val)]} for (x, val) in enumerate(X_tumor[:,0])}
    celltype_dict = {key:val for key, val in celltype_dict.items() if val['cell type'] != 'unoccupied'}
    nx.set_node_attributes(G, celltype_dict)
    
    neighbors_dict, neighbors_dict_mean = calculate_direct_neighbors(G, cells)    
            
    return neighbors_dict_mean, neighbors_dict