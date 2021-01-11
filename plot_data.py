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


from parameters import d

# Set font for plots
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 15}

matplotlib.rc('font', **font)


dataFilePaths = d["parsed_data_files"]
cStates = d["charge_states"]

color_idxs1 = np.linspace(0,1,len(cStates))

colors=[]
for color_id in color_idxs1:
    colors.append(plt.cm.spring(color_id))

# Pack the relevant information into a zipped iterable
ss = zip(cStates, dataFilePaths, colors)


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
        
        
# Plot the 1+ control signal
fName = d["1+_control_signals"][0]
print(fName)
df = pd.read_csv(fName, names=["t", "i"])
t,i = df["t"], df["i"]
t_off = t[ i > 0.5*max(i) ].values[0]
t = t-t_off
ax.plot(t*1e3,50*i/max(i),c="b",ls="--",label="1+ command")
ax.legend()
fig.savefig("./results/data_upto6+.eps", format="eps")
fig2.savefig("./results/data_upto12+.eps", format="eps")


