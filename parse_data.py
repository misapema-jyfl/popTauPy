#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 16:11:20 2020

@author: miha
"""

# Comment

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


onePlusFilePaths=["/home/miha/uni/research/population_confinement_time/data/2020 06 29 Measurements/1+_K1+.txt",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 Measurements/1+_K2+.txt",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 Measurements/1+_K3+.txt",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 Measurements/1+_K4+.txt",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 Measurements/1+_K5+.txt",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 Measurements/1+_K6+.txt",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 Measurements/1+_K7+.txt",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 Measurements/1+_K8+.txt",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 Measurements/1+_K9+.txt",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 Measurements/1+_K10+.txt",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 Measurements/1+_K11+.txt",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 Measurements/1+_K12+.txt"]

nPlusFilePaths=["/home/miha/uni/research/population_confinement_time/data/2020 06 29 Measurements/K1+.txt",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 Measurements/K2+.txt",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 Measurements/K3+.txt",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 Measurements/K4+.txt",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 Measurements/K5+.txt",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 Measurements/K6+.txt",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 Measurements/K7+.txt",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 Measurements/K8+.txt",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 Measurements/K9+.txt",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 Measurements/K10+.txt",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 Measurements/K11+.txt",
"/home/miha/uni/research/population_confinement_time/data/2020 06 29 Measurements/K12+.txt"]

cStates = [1,2,3,4,5,6,7,8,9,10,11,12]


def denoise(fName1, fNameN, q):
    
    df_1 = pd.read_csv(fName1, sep="\t", skipfooter=12, names=["time","current"],engine="python")
    df_n = pd.read_csv(fNameN, sep="\t", skipfooter=12, names=["time","current"],engine="python")
    
    t_1, i_1 = df_1["time"],df_1["current"]
    t_n, i_n = df_n["time"],df_n["current"]
    
    #
    # Find injection onset time
    #
    
    # Normalize 1+ injection signal
    i_1 = i_1/max(i_1)
    
    # Find onset time
    t_onset = t_1[i_1 > 0.1*max(i_1)].values[0]
    
    #
    # Time shift s.t. t=0 <-> onset time
    #
    t_1 = t_1 - t_onset
    t_n = t_n - t_onset
    
    #
    # Remove background
    #
    bg = np.average(i_n[t_n<0])
    i_n = i_n - bg
    
    bg = np.average(i_1[t_1<0])
    i_1 = i_1 - bg
    
    #
    # Scale 1+ signal to current bump height
    #
    i_1 = i_1*max(i_n) / max(i_1)
    
    return t_1, i_1, t_n, i_n


def fourier(t, signal):
    
    n = len(signal)
    shat = np.fft.fft(signal, n)
    psd = shat * np.conj(shat) / n
    
    dt = t.values[1]-t.values[0]
    freq = (1/(dt*n))*np.arange(n)
    
    L = np.arange(1, np.floor(n/2),dtype="int")
    
    freq=freq[L]
    psd = psd[L]
    
    return freq, psd

q = 4
fName1 = "/home/miha/uni/research/population_confinement_time/data/2020 06 29 Measurements/1+_K4+.txt"
fNameN = "/home/miha/uni/research/population_confinement_time/data/2020 06 29 Measurements/K4+.txt"



fig, ax = plt.subplots()
t_1, i_1, t_n, i_n = denoise(fName1, fNameN, q)

noise = i_n[t_n<0]
current = i_n[t_n>0]
t = t_n[t_n < 0]
freq_noise, psd_noise = fourier(t, noise)
t = t_n[t_n > 0]
freq_current, psd_current = fourier(t_n, i_n)

ax.scatter(freq_current,psd_current,label="n+")
ax.scatter(freq_noise,psd_noise, label="noise")
ax.set_xscale("log")
ax.set_yscale("log")
ax.legend()

ihat = np.fft.fft(i_n)

indices = freq_current < 1e3
psd_current = psd_current*indices
freq_current = freq_current*indices
ifilt = np.fft.ifft(ihat)
fig, ax = plt.subplots()
ax.plot(i_n,label="n+")
ax.plot(ifilt,label="n+ clean")



# fig,ax = plt.subplots()
# for fName1, fNameN, q in zip(onePlusFilePaths, nPlusFilePaths, cStates):
    
#     t_1, i_1, t_n, i_n = denoise(fName1, fNameN, q)
    
#     ax.plot(t_1, i_1)
#     ax.plot(t_n, i_n)

    



    # Pack parsed data into df and save
    # s = "/home/miha/uni/research/population_confinement_time/data"
    # saveToPath = s + "/K" + str(q) + "+.csv"
    
    # d = {"time":t_n, "current": i_n}
    # df = pd.DataFrame(d)
    # df.to_csv(saveToPath, sep=",",index=None,header=None)


