# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 08:52:10 2024

@author: 20192020
"""
import Parameters as P
import numpy as np 
import scipy.stats as stats

def write_python_file(filename, outputFolder):
    """ Save the parameters script as a textfile """
    
    # retrieve file with parameter values 
    with open(filename) as f:
        data = f.read()
        f.close()

    # save Python script as textfile 
    with open(f"../output/{outputFolder}/00_default_parameters.txt", mode="w") as f:
        f.write(data)
        f.close()

def calculate_ci(df):        
    m, s, n = df['data'].mean(), df['data'].std(ddof=1), len(df)
    t = stats.t.ppf(0.975, df=n-1)  
    e = t * (s / np.sqrt(n))  
    
    return (m-e, m+e)