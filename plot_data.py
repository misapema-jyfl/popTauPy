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
        'size'   : 15}

matplotlib.rc('font', **font)


'''
Give the data file file paths. Must be in numerical order!
'''
dataFilePaths = [
"/home/miha/uni/research/ppp/code/popTauPy/data/data_2020-06-29_w_noise/K1+.csv",
"/home/miha/uni/research/ppp/code/popTauPy/data/data_2020-06-29_w_noise/K2+.csv",
"/home/miha/uni/research/ppp/code/popTauPy/data/data_2020-06-29_w_noise/K3+.csv",
"/home/miha/uni/research/ppp/code/popTauPy/data/data_2020-06-29_w_noise/K4+.csv",
"/home/miha/uni/research/ppp/code/popTauPy/data/data_2020-06-29_w_noise/K5+.csv",
"/home/miha/uni/research/ppp/code/popTauPy/data/data_2020-06-29_w_noise/K6+.csv",
"/home/miha/uni/research/ppp/code/popTauPy/data/data_2020-06-29_w_noise/K7+.csv",
"/home/miha/uni/research/ppp/code/popTauPy/data/data_2020-06-29_w_noise/K8+.csv",
"/home/miha/uni/research/ppp/code/popTauPy/data/data_2020-06-29_w_noise/K9+.csv",
"/home/miha/uni/research/ppp/code/popTauPy/data/data_2020-06-29_w_noise/K10+.csv",
"/home/miha/uni/research/ppp/code/popTauPy/data/data_2020-06-29_w_noise/K11+.csv",
"/home/miha/uni/research/ppp/code/popTauPy/data/data_2020-06-29_w_noise/K12+.csv"
]

'''
Give the charge states corresponding to data files. Must be in numerical order!
'''
cStates = [1,2,3,4,5,6,7,8,9,10,11,12]

color_idxs1 = np.linspace(0,1,len(cStates))

colors=[]
for color_id in color_idxs1:
    colors.append(plt.cm.spring(color_id))

# Pack the relevant information into a zipped iterable
ss = zip(cStates, dataFilePaths, colors)


# Make the figures
# fig, ax = plt.subplots()
# fig2, ax2 = plt.subplots()
# fig3, ax3 = plt.subplots()

# # Make the plots
# for s in ss:
    
#     # unpack the iterable
#     cState = s[0]
#     fName = s[1]
#     color = s[2]
    
#     # get the data from the data file
#     df = pd.read_csv(fName, names=["time", "current"])
#     t = df["time"]*1e3 # Convert to ms
#     i = df["current"]*1e3 # Convert to nA
    
    
#     # Create the data label
#     c = "{" + str(cState) + "+}"
#     lbl = r"K$^{}$".format(c)
    
#     if cState < 5:
        
#         # plot the data
#         ax.plot(t,i,label=lbl)
#         ax.set_xlim(left=-10, right=30)
#         ax.set_xlabel("Time (ms)")
#         ax.set_ylabel("Current (enA)")
#         ax.legend()
        
        
#     elif cState < 9:
    
#         # plot the data
#         ax2.plot(t,i,label=lbl)
        
#         ax2.set_xlim(left=-10, right=150)
#         ax2.set_xlabel("Time (ms)")
#         ax2.set_ylabel("Current (enA)")
        
#         ax2.legend()
        
#     else:
#         # plot the data
#         ax3.plot(t,i,label=lbl)
#         ax3.set_xlabel("Time (ms)")
#         ax3.set_ylabel("Current (enA)")
        
#         ax3.set_xlim(left=-10,right=360)
#         ax3.legend()
        


# fig.savefig("./results/data_upto4+.png",dpi=500)
# fig2.savefig("./results/data_upto8+.png",dpi=500)
# fig3.savefig("./results/data_upto12+.png",dpi=500)

# Make the figures
fig, ax = plt.subplots()
fig2, ax2 = plt.subplots()


# Pack the relevant information into a zipped iterable
ss = zip(cStates, dataFilePaths, colors)


# Make the plots
for s in ss:
    
    # unpack the iterable
    cState = s[0]
    fName = s[1]
    color = s[2]
    
    # get the data from the data file
    df = pd.read_csv(fName, names=["time", "current"])
    t = df["time"]*1e3 # Convert to ms
    i = df["current"]*1e3 # Convert to nA
    
    
    # Create the data label
    c = "{" + str(cState) + "+}"
    lbl = r"K$^{}$".format(c)
    
    if cState < 7:
        
        # plot the data
        ax.plot(t,i,label=lbl)
        ax.set_xlim(left=-10, right=50)
        ax.set_xlabel("Time (ms)")
        ax.set_ylabel("Current (enA)")
        ax.set_ylim(bottom=-5, top=70)
        ax.legend()
        
        
    else:
    
        # plot the data
        ax2.plot(t,i,label=lbl)
        ax2.set_xlabel("Time (ms)")
        ax2.set_ylabel("Current (enA)")        
        ax2.set_xlim(left=-10, right=350)
        ax2.set_ylim(bottom=-5, top=70)
        ax2.legend()
        

fig.savefig("./results/data_upto6+.png",dpi=500)
fig2.savefig("./results/data_upto12+.png",dpi=500)


