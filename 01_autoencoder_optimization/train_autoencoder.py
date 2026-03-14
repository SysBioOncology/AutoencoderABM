#%% Import necessary libraries for loading and processing data
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA

import matplotlib.pyplot as plt
import plotly.express as px

# Import Keras for implementing autoencoders
import keras
from keras.models import Sequential
from keras.layers import Input, Dense, Conv2D, MaxPooling2D, Conv2DTranspose, Flatten, Reshape
from PIL import Image

import tensorflow as tf
import cv2

# load some other functions 
from utils import load_data_synthetic, load_data_tcga, load_data_tumoroid
from utils import vis_synthetic, vis_tcga, vis_tumoroid

#%% Set some variables 
# folder to load data / save data to 
input_folder = 'data/synthetic'

# filenames to save autoencoder and encoder 
model_filename = 'encoder_architecture.keras'
autoencoder_filename = 'autoencoder_architecture.keras'

# data set 
data_set = 'synthetic'

# variables for loading data 
# synthetic 
IMrwalk = [0, 0.5, 1] 
duplicates = 10 
samples = 5
new_size = (100, 100)

# tcga 
#samples = 15 

# training autoencoder
epochs = 200
batch = 128

# PCA
n_components = 25


#%% Load data
if data_set == 'synthetic':
    X_train, y_train = load_data_synthetic(input_folder, samples, duplicates, 
                                           IMrwalk, new_size)
elif data_set == 'tumoroid':
    # get data from first experiment 
    one_folder = input_folder + '/data_M07'
    n_zstacks = 1
    start = 4
    n_time = 72
    n_zstacks_m10 = 1
    start_m10 = 3
    a = n_zstacks * n_time + n_zstacks_m10 * n_time # total number of images
    
    X_train, all_z, all_t, idx = load_data_tumoroid(one_folder, n_zstacks, a,
                           start, n_time, X_train = None,
                           all_z = None, all_t = None, idx = 0,
                           new_size = (100,100))
    
    # get data from second experiment 
    one_folder = input_folder + '/data_M10'
    X_train, all_z, all_t, _ = load_data_tumoroid(one_folder, n_zstacks_m10, a,
                           start_m10, n_time, X_train = X_train,
                           all_z = all_z, all_t = all_t, idx=idx,
                           new_size = (100,100))
    
    y_train = np.concatenate((np.array(all_z).reshape(a,1),np.array(all_t).reshape(a,1)), axis=1)
    
elif data_set == 'tcga':
    X_train, y_train = load_data_tcga(input_folder, new_size, samples)


# split into train and validation/test sets 
X_train, X_test, y_train, y_test = train_test_split(X_train, y_train, test_size=0.1, random_state=42)
        

#%% Setup autoencoder 
# Define the autoencoder architecture
input_dim = X_train.shape[1]
encoding_dim = 2

input_layer = keras.layers.Input(shape=(input_dim, input_dim, 1))

# Encoder
x = Conv2D(1, (9, 9), activation="relu", padding="same")(input_layer)
x = MaxPooling2D((2, 2), padding="same")(x)
x = Conv2D(1, (9, 9), activation="relu", padding="same")(x)
x = MaxPooling2D((2, 2), padding="same")(x)
x = Flatten()(x)

# Decoder
x = Reshape((25, 25, 1))(x)
x = Conv2DTranspose(1, (9, 9), strides=2, activation="relu", padding="same")(x)
x = Conv2DTranspose(1, (9, 9), strides=2, activation="relu", padding="same")(x)
x = Conv2D(1, (3, 3), activation="sigmoid", padding="same")(x)

# Autoencoder
autoencoder = keras.Model(input_layer, x)
autoencoder.compile(optimizer="adam", loss="mse")
autoencoder.summary()


#%% Train autoencoder
fitting_results = autoencoder.fit(X_train, X_train, epochs=epochs, batch_size = batch, 
                                  shuffle = True)
y_pred = autoencoder.predict(X_train)


#%% Save autoencoder and get latent features 
encoder = keras.Model(inputs=autoencoder.input, outputs=autoencoder.layers[5].output)
encoded_features_train = encoder.predict(X_train)

encoder.save(f"{input_folder}/{model_filename}")
autoencoder.save(f"{input_folder}/{autoencoder_filename}")


#%% Visualize latent space 
# perform PCA
pca = PCA(n_components=n_components)
pca_results = pca.fit_transform(encoded_features_train)
encoded_features_test = encoder.predict(X_test)
pca_results_test = pca.transform(encoded_features_test)
explained_variance = pca.explained_variance_ratio_
print(f"Train; explained variance: {explained_variance[0] + explained_variance[1]}")

# visualize PCA
if data_set == 'synthetic':
    vis_synthetic(pca_results, pca_results_test, y_train, 
                      X_test, y_test,
                      explained_variance, samples, epochs, batch, 
                      input_folder)
elif data_set == 'tumoroid':
    vis_tumoroid(pca_results, y_train, explained_variance,
                     input_folder, samples, epochs, batch, 
                     pca_results_test, y_test)
elif data_set == 'tcga':
    vis_tcga(y_train, y_test, pca_results, pca_results_test, 
             explained_variance, input_folder, 
             samples, epochs, batch)






