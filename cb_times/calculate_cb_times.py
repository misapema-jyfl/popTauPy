#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 18 10:04:25 2021

Script for calculating the charge breeding times.
Use the parsed currents as input.

@author: miha
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

font = {"size":15}
matplotlib.rc("font", **font)

# Import the files
onePlusFiles = ["/home/miha/Work/codes/Consecutive_Transients_Analyzer/test/results/1+_K3+.csv",
"/home/miha/Work/codes/Consecutive_Transients_Analyzer/test/results/1+_K4+.csv",
"/home/miha/Work/codes/Consecutive_Transients_Analyzer/test/results/1+_K5+.csv",
"/home/miha/Work/codes/Consecutive_Transients_Analyzer/test/results/1+_K6+.csv",
"/home/miha/Work/codes/Consecutive_Transients_Analyzer/test/results/1+_K7+.csv",
"/home/miha/Work/codes/Consecutive_Transients_Analyzer/test/results/1+_K8+.csv",
"/home/miha/Work/codes/Consecutive_Transients_Analyzer/test/results/1+_K9+.csv",
"/home/miha/Work/codes/Consecutive_Transients_Analyzer/test/results/1+_K10+.csv",
"/home/miha/Work/codes/Consecutive_Transients_Analyzer/test/results/1+_K11+.csv",
"/home/miha/Work/codes/Consecutive_Transients_Analyzer/test/results/1+_K12+.csv"]
nPlusFiles = ["/home/miha/Work/codes/Consecutive_Transients_Analyzer/test/results/K3+.csv",
"/home/miha/Work/codes/Consecutive_Transients_Analyzer/test/results/K4+.csv",
"/home/miha/Work/codes/Consecutive_Transients_Analyzer/test/results/K5+.csv",
"/home/miha/Work/codes/Consecutive_Transients_Analyzer/test/results/K6+.csv",
"/home/miha/Work/codes/Consecutive_Transients_Analyzer/test/results/K7+.csv",
"/home/miha/Work/codes/Consecutive_Transients_Analyzer/test/results/K8+.csv",
"/home/miha/Work/codes/Consecutive_Transients_Analyzer/test/results/K9+.csv",
"/home/miha/Work/codes/Consecutive_Transients_Analyzer/test/results/K10+.csv",
"/home/miha/Work/codes/Consecutive_Transients_Analyzer/test/results/K11+.csv",
"/home/miha/Work/codes/Consecutive_Transients_Analyzer/test/results/K12+.csv"]
cStates = [3,4,5,6,7,8,9,10,11,12]


def calculate_total_number_extracted(t, i, q):
    '''Integrate the given signal t, i to find the total number of 
    charge state q particles extracted.
    
    Signal i must be in units of A, while t must be in s.'''
    e = 1.6021773e-19
    particles = i/(q*e)
    numOfExtractedParticles = np.trapz(y=particles, x=t)
    return numOfExtractedParticles

def determine_cb_time(t, i, q):
    '''Determines the charge breeding time - defined as the 
    time required to extract 90% of total particles extracted -
    when given the signal t, i and charge state q of ion 
    in question.'''
    # Integrate to find the number of extracted N+ particles
    numOfExtractedParticles = calculate_total_number_extracted(t=t,
                                                               i=i,
                                                               q=cState)
    
    # Determine the point in time when 90% of all particles have been extracted
    e = 1.6021773e-19
    particles = i/(q*e)
    condition = True
    j = 0
    while condition:
        t_tmp = t.values[0:j+1]
        particles_tmp = particles.values[0:j+1]
        numberExtracted = np.trapz(y=particles_tmp, x=t_tmp)
        
        if numberExtracted/numOfExtractedParticles > 0.90:
            condition=False
        else:
            j+=1
    
    t_cb = t.values[j+1] # Charge breeding time

    return t_cb

def select_extraction_current(t, i):
    '''Select only the non-zero extraction current.'''
    # Determine peak location
    c = (i == max(i))
    t_peak = t[c].values[0]
    
    # Choose portion of signal on RHS of peak
    iRHS = i[t > t_peak]
    tRHS = t[t > t_peak]
    
    # Find point in time when current 
    # has decayed back to zero
    c = iRHS < max(iNP)*0.01
    t_decay = tRHS[c].values[0]
    
    # Choose signal from t=0 to t=t_decay
    c = (t > 0)&(t<t_decay)
    t = t[c]
    i = i[c]
    
    return t, i
    
def determine_cb_efficiency(t, i, q, pulseWidth, onePlusIntensity):
    '''Give injection pulse width in units of s,
    and the intensity in units of A. 
    Similarly, the extracted transient signal, t in s and i in A.'''
    
    # Calculate total number of extracted particles
    numExtracted = calculate_total_number_extracted(t, i, q)
    
    # Number of injected particles
    e = 1.6021773e-19
    numInjected = pulseWidth*onePlusIntensity/e
    
    efficiency = 100*numExtracted/numInjected 
    
    return efficiency

cb_times = []
efficiencies = []
for i in range(len(cStates)):
    onePlusFile = onePlusFiles[i]
    nPlusFile = nPlusFiles[i]
    cState = cStates[i]
    
    # Select the N+ current
    nPlus = pd.read_csv(nPlusFile)
    windowSize = int(.1*len(nPlus)/100) # rolling average
    # nPlus = nPlus.rolling(window=windowSize, min_periods=1).mean() 
    tNP = nPlus["t"]
    iNP = nPlus["i"]
    
    # Select only the extraction transient
    # neglecting zero current
    t, i = select_extraction_current(tNP,iNP)
    
    # Determine the CB time
    t_cb = determine_cb_time(t=t, i=i*1e-6, q=cState)
    
    # Store cb time in dictionary
    cb_times.append(t_cb)
    
    # Conversion factor 
    F = (1/1e6 + 1/5.7e6)
    
    # Determine cb efficiency
    eff = determine_cb_efficiency(t=t,
                                  i=i*F,
                                  q=cState,
                                  pulseWidth=5*1e-3,
                                  onePlusIntensity=500e-9)
    efficiencies.append(eff)
    
    # Choose data for fill_between plotting
    c = (t>0)&(t<t_cb)
    t_shade = t[c]
    i_shade = i[c]
    
    # Create a figure for visualization
    fig, ax = plt.subplots()
    ax.plot(t*1e3, i, "firebrick")
    ax.fill_between(t_shade*1e3, 0, i_shade, color="crimson")
    
    # Adjust plot settings
    # Save
    # Close
    ax.set(xlim=[0, None], ylim=[0, None],
           xlabel="Time (ms)",
           ylabel=r"Current ($e\mu$A)")
    fig.tight_layout()
    plt.savefig("fig_transient_K{}+.eps".format(str(cState)), format="eps")
    plt.savefig("fig_transient_K{}+.png".format(str(cState)), dpi=300)
    plt.close()
    
    
    
# Write the cb times to file
with open("out_cb_times.csv", "w") as f:
    f.write("cState,tau_cb (s),efficiency (%)\n")
    for cState, t_cb, eff in zip(cStates, cb_times, efficiencies):
        s = (str(cState),str(t_cb),str(eff))
        line = ",".join(s) + "\n"
        f.write(line)
f.close()