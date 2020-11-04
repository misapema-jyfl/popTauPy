#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 14:16:24 2020

@author: miha
"""
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# Set font for plots
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 12}

matplotlib.rc('font', **font)


'''
Give the data file file paths. Must be in numerical order!
'''
dataFilePaths = [
     "/home/miha/uni/research/population_confinement_time/data/2020 06 29 measurement wo noise 2/K1 wo noise.csv",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 measurement wo noise 2/K2 wo noise.csv",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 measurement wo noise 2/K3 wo noise.csv",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 measurement wo noise 2/K4 wo noise.csv",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 measurement wo noise 2/K5 wo noise.csv",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 measurement wo noise 2/K6 wo noise.csv",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 measurement wo noise 2/K7 wo noise.csv",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 measurement wo noise 2/K8 wo noise.csv",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 measurement wo noise 2/K9 wo noise.csv",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 measurement wo noise 2/K10 wo noise.csv",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 measurement wo noise 2/K11 wo noise.csv",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 measurement wo noise 2/K12 wo noise.csv"
    ]

'''
Give the charge states corresponding to data files. Must be in numerical order!
'''
cStates = [1,2,3,4,5,6,7,8,9,10,11,12]

color_idxs = np.linspace(0,1,len(cStates))
colors = [plt.cm.cool(x) for x in color_idxs]

fig, ax = plt.subplots()
for cState, fName, col in zip(cStates, dataFilePaths, colors):
    
    df = pd.read_csv(fName, names=["time", "current"])
    t = df["time"]
    i = df["current"]
        
    ax.plot(t,i, color = col)
    
plt.savefig("data_2020-06-09_wo_noise_2.png", dpi=300)