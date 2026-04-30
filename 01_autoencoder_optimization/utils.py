import pandas as pd 
import numpy as np 
import cv2

import matplotlib.pyplot as plt 


def load_data_synthetic(input_folder, samples, duplicates, IMrwalk,
                        new_size):
    """Prepare X_train and y_train of synthetic data"""
    
    # prepare labels (for later plotting with latent space)
    one_y = pd.read_csv(f"{input_folder}/synthetic_proceedings.csv").head(samples)
    one_y = pd.DataFrame(np.repeat(one_y.values, repeats=duplicates, axis=0),
                         columns=one_y.columns)
    y_train = one_y.copy()
    y_train['IMrwalk'] = IMrwalk[0]

    # get y_train of every IMrwalk folder
    for r in IMrwalk[1:]:
        temp = one_y.copy()
        temp['IMrwalk'] = r
        y_train = pd.concat([y_train, temp])
    y_train = np.array(y_train)
    
    # prepare images 
    X_train = np.zeros((int(len(IMrwalk)*samples*duplicates), new_size[0], new_size[1]))
    idx = 0 
    idx_one_sample = 0
    for rwalk in IMrwalk:    
        idx += samples
        for one_sample in range(samples):
            print(f"Loading IMrwalk={str(rwalk)}, sample {one_sample}")
            for one_dup in range(0, duplicates):
                filename = f"{input_folder}/IMrwalk={str(rwalk)}/{one_sample}_{one_dup}_TumorCells_ImmuneCells.png"
                imarray = cv2.imread(filename, cv2.IMREAD_GRAYSCALE) # read image
                imarray = imarray[936:3912, 600:3576] # crop 
                imarray = cv2.resize(imarray, dsize=new_size) # resize to correct height and width
                X_train[idx_one_sample,] = imarray.astype(float) / 255 # convert to float for autoencoder
                idx_one_sample += 1 
                
    return X_train, y_train


def load_data_tumoroid(input_folder, n_zstacks, a,
                       start, n_time, X_train = None, 
                       all_z = None, all_t = None, idx = 0, 
                       new_size = (100,100)):
    """Prepare X_train and y_train of tumoroid data"""
    # intialize some variables if needed 
    if X_train is None:
        X_train = np.zeros((a, new_size[0], new_size[1]))

        all_z = []
        all_t = []

    # load all input arrays 
    for z in range(start, start+n_zstacks):
        for t in range(1, n_time+1):
            all_z.append(z)
            all_t.append(t)
            
            filename = f"{input_folder}/z{str(z).zfill(2)}t{str(t).zfill(2)}_TumorCells_ImmuneCells_0000.png"
            imarray = cv2.imread(filename, cv2.IMREAD_GRAYSCALE)
            imarray = imarray[936:3912, 600:3576]
            imarray = cv2.resize(imarray, dsize=new_size)
            X_train[idx,] = imarray.astype(float) / 255
            idx += 1 
                    
    return X_train, all_z, all_t, idx


def load_data_tcga(input_folder, new_size, nSamples):
    """Prepare X_train and y_train of TCGA patches"""
    # load patient filenames
    filenames = pd.read_csv(f"{input_folder}/00_filenames_annots.csv")['file'].to_list()

    # load all input arrays 
    X_train = np.zeros((len(filenames), new_size[0], new_size[1]))
    idx = 0
    for file in filenames[:nSamples]:
        print(f">> Loading sample {file}")
        filename = f"{input_folder}/{file}.png"
        imarray = cv2.imread(filename, cv2.IMREAD_GRAYSCALE)
        imarray = imarray[946:3902, 610:3566] 
        imarray = cv2.resize(imarray, dsize=new_size)
        X_train[idx,] = imarray.astype(float) / 255
        idx += 1 
    
    # prepare y_train based on immune subtype 
    y_train = np.array(pd.read_csv(f"{input_folder}/00_filenames_annots.csv")['immune subtype']).flatten()
    
    return X_train, y_train

def vis_synthetic(pca_results, pca_results_test, y_train, 
                  X_test, y_test, explained_variance, samples, 
                  epochs, batch, output_folder
                  ):
    """Visualize PCA of synthetic data, by coloring the principal component 
    plots based on TUpprol, IMpkill or IMrwalk values"""
    
    # TUpprol of training set
    fig, ax = plt.subplots(figsize=(6, 5))
    pca_fig = ax.scatter(pca_results[:,0], pca_results[:,1], c=y_train[:,0])
    cb = fig.colorbar(pca_fig, ax=ax)
    cb.set_label('TUpprol')
    ax.set_xlabel('PC1');
    ax.set_ylabel('PC2');
    title_variance = explained_variance[0] + explained_variance[1]
    title_variance = round(title_variance, 4)
    ax.set_title(f"Explained variance: {title_variance:.4f}")
    fig.savefig(f"{output_folder}/train_encoded_TUpprol_nSamples={samples}_epochs={epochs}_batch={batch}.png")

    # TUpprol of validation set
    fig, ax = plt.subplots(figsize=(6, 5))
    pca_fig = ax.scatter(pca_results_test[:,0], pca_results_test[:,1], c=y_test[:,0])
    cb = fig.colorbar(pca_fig, ax=ax)
    cb.set_label('TUpprol')
    ax.set_xlabel('PC1');
    ax.set_ylabel('PC2');
    title_variance = explained_variance[0] + explained_variance[1]
    title_variance = round(title_variance, 4)
    ax.set_title(f"Validation; explained variance: {title_variance:.4f}")
    fig.savefig(f"{output_folder}/val_encoded_TUpprol_nSamples={samples}_epochs={epochs}_batch={batch}.png")

    # IMpkill of training set
    fig, ax = plt.subplots(figsize=(6, 5))
    pca_fig = ax.scatter(pca_results[:,0], pca_results[:,1], c=y_train[:,1])
    cb = fig.colorbar(pca_fig, ax=ax)
    cb.set_label('IMpkill')
    ax.set_xlabel('PC1');
    ax.set_ylabel('PC2');
    title_variance = explained_variance[0] + explained_variance[1]
    title_variance = round(title_variance, 4)
    ax.set_title(f"Explained variance: {title_variance:.4f}")
    fig.savefig(f"{output_folder}/train_encoded_IMpkill_nSamples={samples}_epochs={epochs}_batch={batch}.png")

    # IMpkill of validation set
    fig, ax = plt.subplots(figsize=(6, 5))
    pca_fig = ax.scatter(pca_results_test[:,0], pca_results_test[:,1], c=y_test[:,1])
    cb = fig.colorbar(pca_fig, ax=ax)
    cb.set_label('IMpkill')
    ax.set_xlabel('PC1');
    ax.set_ylabel('PC2');
    title_variance = explained_variance[0] + explained_variance[1]
    title_variance = round(title_variance, 4)
    ax.set_title(f"Validation; explained variance: {title_variance:.4f}")
    fig.savefig(f"{output_folder}/val_encoded_IMpkill_nSamples={samples}_epochs={epochs}_batch={batch}.png")

    # IMrwalk of training set 
    fig, ax = plt.subplots(figsize=(6, 5))
    pca_fig = ax.scatter(pca_results[:,0], pca_results[:,1], c=y_train[:,2])
    cb = fig.colorbar(pca_fig, ax=ax)
    cb.set_label('IMrwalk')
    ax.set_xlabel('PC1');
    ax.set_ylabel('PC2');
    title_variance = explained_variance[0] + explained_variance[1]
    title_variance = round(title_variance, 4)
    ax.set_title(f"Explained variance: {title_variance:.4f}")
    fig.savefig(f"{output_folder}/train_encoded_IMrwalk_nSamples={samples}_epochs={epochs}_batch={batch}.png")

    # IMrwalk of validation set 
    fig, ax = plt.subplots(figsize=(6, 5))
    pca_fig = ax.scatter(pca_results_test[:,0], pca_results_test[:,1], c=y_test[:,2])
    cb = fig.colorbar(pca_fig, ax=ax)
    cb.set_label('IMrwalk')
    ax.set_xlabel('PC1');
    ax.set_ylabel('PC2');
    title_variance = explained_variance[0] + explained_variance[1]
    title_variance = round(title_variance, 4)
    ax.set_title(f"Validation; explained variance: {title_variance:.4f}")
    fig.savefig(f"{output_folder}/val_encoded_IMrwalk_nSamples={samples}_epochs={epochs}_batch={batch}.png")


def vis_tumoroid(pca_results, y_train, explained_variance,
                 output_folder, samples, epochs, batch, 
                 pca_results_test, y_test):
    """Visualize PCA of tumoroid data, by coloring the principal component 
    plots based on time or z stack"""
    
    # z stacks, training
    fig, ax = plt.subplots(figsize=(6, 5))
    pca_fig = ax.scatter(pca_results[:,0], pca_results[:,1], c=y_train[:,0])
    cb = fig.colorbar(pca_fig, ax=ax)
    cb.set_label('z stacks')
    ax.set_xlabel('PC1');
    ax.set_ylabel('PC2');
    title_variance = explained_variance[0] + explained_variance[1]
    title_variance = round(title_variance, 4)
    ax.set_xlim(14.5, 16.5)
    ax.set_title(f"Explained variance: {title_variance:.4f}")
    fig.savefig(f"{output_folder}/train_z_encoded_nSamples={samples}_epochs={epochs}_batch={batch}.png")

    # time points, training
    fig, ax = plt.subplots(figsize=(6, 5))
    pca_fig = ax.scatter(pca_results[:,0], pca_results[:,1], c=y_train[:,1])
    cb = fig.colorbar(pca_fig, ax=ax)
    cb.set_label('time points')
    ax.set_xlabel('PC1');
    ax.set_ylabel('PC2');
    title_variance = explained_variance[0] + explained_variance[1]
    title_variance = round(title_variance, 4)
    ax.set_xlim(14.5, 16.5)
    ax.set_title(f"Explained variance: {title_variance:.4f}")
    fig.savefig(f"{output_folder}/train_t_encoded_nSamples={samples}_epochs={epochs}_batch={batch}.png")

    # z stacks, validation
    fig, ax = plt.subplots(figsize=(6, 5))
    pca_fig = ax.scatter(pca_results_test[:,0], pca_results_test[:,1], c=y_test[:,0])
    cb = fig.colorbar(pca_fig, ax=ax)
    cb.set_label('z stacks')
    ax.set_xlabel('PC1');
    ax.set_ylabel('PC2');
    title_variance = explained_variance[0] + explained_variance[1]
    title_variance = round(title_variance, 4)
    ax.set_xlim(14.5, 16.5)
    ax.set_title(f"Explained variance: {title_variance:.4f}")
    fig.savefig(f"{output_folder}/val_z_encoded_nSamples={samples}_epochs={epochs}_batch={batch}.png")

    # time points, validation
    fig, ax = plt.subplots(figsize=(6, 5))
    pca_fig = ax.scatter(pca_results_test[:,0], pca_results_test[:,1], c=y_test[:,1])
    cb = fig.colorbar(pca_fig, ax=ax)
    cb.set_label('time points')
    ax.set_xlabel('PC1');
    ax.set_ylabel('PC2');
    title_variance = explained_variance[0] + explained_variance[1]
    title_variance = round(title_variance, 4)
    ax.set_xlim(14.5, 16.5)
    ax.set_title(f"Explained variance: {title_variance:.4f}")
    fig.savefig(f"{output_folder}/val_t_encoded_nSamples={samples}_epochs={epochs}_batch={batch}.png")


def vis_tcga(y_train, y_test, pca_results, pca_results_test, 
             explained_variance, output_folder, 
             samples, epochs, batch):
    """Visualize PCA of TCGA patches, by coloring the principal component 
    plots based on immune subtype"""
    
    # convert immune subtypes to integer
    y_train[y_train=='D']=0
    y_train[y_train=='IE']=1
    y_test[y_test=='D']=0
    y_test[y_test=='IE']=1

    # training 
    samples = y_train.shape[0]
    fig, ax = plt.subplots(figsize=(6, 5))
    pca_fig = ax.scatter(pca_results[:,0], pca_results[:,1], c=y_train)
    cb = fig.colorbar(pca_fig, ax=ax)
    cb.set_label('immune subtype')
    ax.set_xlabel('PC1');
    ax.set_ylabel('PC2');
    title_variance = explained_variance[0] + explained_variance[1]
    title_variance = round(title_variance, 4)
    ax.set_title(f"Explained variance: {title_variance:.4f}\ndesert=0, immune-enriched=1")
    fig.savefig(f"{output_folder}/train_encoded_nSamples={samples}_epochs={epochs}_batch={batch}.png")

    # validation 
    fig, ax = plt.subplots(figsize=(6, 5))
    pca_fig = ax.scatter(pca_results_test[:,0], pca_results_test[:,1], c=y_test)
    cb = fig.colorbar(pca_fig, ax=ax)
    cb.set_label('immune subtype')
    ax.set_xlabel('PC1');
    ax.set_ylabel('PC2');
    title_variance = explained_variance[0] + explained_variance[1]
    title_variance = round(title_variance, 4)
    ax.set_title(f"Validation; explained variance: {title_variance:.4f}\ndesert=0, immune-enriched=1")
    fig.savefig(f"{output_folder}/val_encoded_nSamples={samples}_epochs={epochs}_batch={batch}.png")







        