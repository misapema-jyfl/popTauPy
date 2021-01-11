#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 16:11:20 2020

@author: miha
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from parameters import parsing_parameters as pp


# Set the paths to files to use
onePlusFilePaths = pp["1+_control_signals"]
nPlusFilePaths = pp["n+_signals"]
cStates = pp["charge_states"]

# Pack filepaths to data frame
df = pd.DataFrame()
df["cState"] = cStates
df["onePlusFilePath"] = onePlusFilePaths
df["nPlusFilePath"] = nPlusFilePaths


def get_data_from_file(fileName):
    '''
    Given path to file, retrieve t,i data
    '''
    data = pd.read_csv( fileName, names = ["t", "i"] )
    
    t, i = data["t"], data["i"]
    
    return t, i
    

if not os.path.isdir(pp["save_to_path"]+"/figures/"):
    os.mkdir(pp["save_to_path"]+"/figures/")
    
    
# Parse all files
for charge_state in cStates:
    
    # Make figure to check results
    fig, ax = plt.subplots()
    
    # Get 1+ command signal and extracted n+ signal
    # corresponding to desired charge state
    
    c = df["cState"] == charge_state
    
    onePlusFileName = df[c]["onePlusFilePath"].values[0]
    nPlusFileName = df[c]["nPlusFilePath"].values[0]
    
    t_1, i_1 = get_data_from_file(onePlusFileName)
    t_n, i_n = get_data_from_file(nPlusFileName)
    
    
    # Remove time offset, i.e. set 1+ rise onset as t=0
    t_off = t_1[ i_1 > 0.1*max(i_1) ].values[0]
    t_1 = t_1 - t_off
    t_n = t_n - t_off
    
    # Remove background from n+ signal
    bg = np.average(i_n[t_n<0])
    i_n = i_n-bg
    
    
    # Remove background from 1+ signal
    bg = np.average(i_1[t_1<0])
    i_1 = i_1-bg
    
    # Normalize 1+ command signal
    
    i_1 = max(i_n)*i_1/max(i_1)
    
    # Output parsed files
    onePlusOutName = pp["save_to_path"] + "/1+_" + pp["elemental_symbol"] + str(charge_state) + "+.csv"
    nPlusOutName = pp["save_to_path"] + "/" + pp["elemental_symbol"] + str(charge_state) + "+.csv"
    
    df_1 = pd.DataFrame([t_1,i_1])
    df_n = pd.DataFrame([t_n,i_n])
    
    df_1.to_csv(onePlusOutName, index=None)
    df_n.to_csv(nPlusOutName, index=None)

    # Plot for sanity check
    ax.plot(t_1, i_1)
    ax.plot(t_n, i_n)
    
    
    fig.savefig(pp["save_to_path"] + "/figures/" + pp["elemental_symbol"] + str(charge_state) + "+.png")