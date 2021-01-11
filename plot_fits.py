#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 09:26:45 2021

@author: miha
"""

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import interpolate

from parameters import d
from parameters import p
from parameters import plotting_fits


# Set font for plots
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 15}

matplotlib.rc('font', **font)




# Dataframe for holding abc parameters
df_abc = pd.read_csv(p["abc_file_path"],index_col=0)
colNames = [ int(name) for name in df_abc.columns.values ] # Convert columns to int
df_abc.columns = colNames

# Dataframe for holding data file paths
df = pd.DataFrame()
df["cState"] = d["charge_states"]
df["fPath"] = d["parsed_data_files"]



def get_data(charge_state):
    '''
    Get measurement data corresponding to desired charge state.
    '''
    c = df["cState"] == charge_state
    data = pd.read_csv(df[c]["fPath"].values[0])
    
    ndf = pd.DataFrame(data)
    ndf.columns = ["t", "i"]
    
    t,i = ndf["t"],ndf["i"]
    
    return t,i



def get_abc(charge_state):
    '''
    Get the abc parameters for fit corresponding to desired charge state.
    '''
    a,b,c = df_abc[charge_state]["a"], df_abc[charge_state]["b"], df_abc[charge_state]["c"]

    return a,b,c



def get_interp_data(charge_state):
    '''
    Interpolate data corresponding to given charge state,
    and return an interpolate object for the data.
    '''
    t,i = get_data(charge_state)
    
    I = interpolate.interp1d(t,i,kind="cubic")
    
    return I


def make_interpolate_dataframe(charge_state):
    '''
    Make a dataframe out of the interpolate objects neighboring 
    a given charge state.
    '''
    df_interp = pd.DataFrame()
    interps = []
    for q in [charge_state-1, charge_state, charge_state+1]:
        I = get_interp_data(q)
        interps.append(I)
    df_interp["cState"] = [charge_state-1, charge_state, charge_state+1]
    df_interp["interp"] = interps
    
    return df_interp

def RHS(t,y,q,a,b,c):
    '''
    Defines the Right-Hand-Side of the differential equation, 
    which determines the time-evolution of the ion extraction current
    with respect to the neighboring charge state currents.
    '''
    a,b,c = get_abc(q)
    
    I_lo = df_interp[df_interp["cState"]==q-1]["interp"].values[0]
    I_hi = df_interp[df_interp["cState"]==q+1]["interp"].values[0]
    
    return a*I_lo(t) - b*y + c*I_hi(t)



def rk4(t0, y0, tf, f,**kwargs):
    '''
    Fourth order Runge-Kutta integrator for the arbitrary function
    f=f(t,y,**kwargs), where kwargs are the variables necessary to
    completely define the right-hand-side of the differential equation
    dy/dt = f(t,y; **kwargs)

    Parameters
    ----------
    t0 : float
        Moment in time to begin solution.
    y0 : float
        Initial value at t = t0.
    tf : float
        Final time.
    f : function
        Function which defines the DE to solve.
    **kwargs : dict
        Keyword arguments to completely define the right-hand-side of
        the DE.

    Returns
    -------
    T : array
        Time array.
    Y : array
        Function value array.

    '''
    
    h = plotting_fits["h"] # Time step (units of s) for the RK4 integrator. 
    Y = []
    T = []
    yn = y0
    tn = t0
    
    while tn <= tf-h:
        #
        # Update time T and solution Y
        #
        T.append(tn)
        Y.append(np.float(yn))
        #
        # Calculate new coefficients k
        #
        k1 = h*f(tn, yn, **kwargs)
        k2 = h*f(tn+0.5*h, yn + 0.5*k1, **kwargs)
        k3 = h*f(tn+0.5*h, yn + 0.5*k2, **kwargs)
        k4 = h*f(tn+h, yn + k3, **kwargs)
        #
        # Calculate next solution value and time 
        #
        yn += (1/6)*(k1 + 2*k2 + 2*k3 + k4)                  # Update yn
        tn += h                                              # Update tn
    
    return np.array(T),np.array(Y)





def determine_interpolation_limits(charge_state):
    '''
    Determines the interpolation limits for three consecutive transients,
    with t_i set to 0 s, always.
    '''    
    t_lo, i_lo = get_data(charge_state-1)
    t_hi, i_hi = get_data(charge_state+1)
    
    t_i = 0
    t_f = min( max(t_lo), max(t_hi) )
    
    return t_i, t_f






'''
Do the plotting. 

Set plotting parameters in the corresponding dictionary
in the file 'parameters.py'
'''
for q in plotting_fits["charge_states"]:
    
    df_interp = make_interpolate_dataframe(q)
    
    print("Plotting fit: " + str(q) + "+...")
    
    t_i, t_f = determine_interpolation_limits(q)
    a,b,c = get_abc(q)
    T,Y = rk4(t_i, 0, t_f, RHS, q=q, a=a,b=b,c=c)
    t,i = get_data(q)
    
    fig, ax = plt.subplots()
    
    ax.plot(t*1e3,i*1e3)
    ax.plot(T*1e3,Y*1e3)
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Current (nA)")
    
    outName = plotting_fits["output_file_name"]
    outType = plotting_fits["output_file_type"]
    
    fig.tight_layout()
    
    fig.savefig("./results/" + outName + "_q=" + str(q) + "+" + outType)

print("Finished!")







