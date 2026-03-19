import cv2
import numpy as np 
from utils_spatial_features import get_neighbors, calc_complexity 

n_images = 1
new_size = (100, 100) # ABM grid size, height x width 
X = np.zeros((n_images, new_size[0], new_size[1])) # prepare array with all images

# lines below can be put into for-loop to go over multiple images 
example_img = 'data/tumoroid/data_M07/z04t22_TumorCells_ImmuneCells_0000.png'
imarray = cv2.imread(example_img, cv2.IMREAD_GRAYSCALE)
imarray = imarray[936:3912, 600:3576]
imarray = cv2.resize(imarray, dsize=new_size)
X[0,] = imarray.astype(float) / 255 

# variables for complexity and neighbors calculation 
cells = ['tumor', 'lymphocyte']
int2cell = {'unoccupied':242, 'tumor':106, 'lymphocyte':156}
value2cell = {str(round(242/255,2)):'unoccupied', str(round(106/255, 2)):'tumor', str(round(156/255, 2)):'lymphocyte'}

# code below can (same as above) be put into for-loop to go over multiple images
# neighbors
train_img = X[0,].round(2)
avg_neighbors, all_neighbors = get_neighbors(train_img, cells, value2cell)
area, perimeter, complexity = calc_complexity(train_img, value2cell, int2cell)
