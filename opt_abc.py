#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 12:43:01 2020

@author: miha
"""

import matplotlib
import pandas as pd
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import numpy as np
from scipy import interpolate


# Set font for plots
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 12}

matplotlib.rc('font', **font)


'''
Give the data file file paths. Must be in numerical order!
'''

dataFilePaths = [
"/home/miha/uni/research/ppp/data/data_2020-06-29_w_noise/K1+.csv",
"/home/miha/uni/research/ppp/data/data_2020-06-29_w_noise/K2+.csv",
"/home/miha/uni/research/ppp/data/data_2020-06-29_w_noise/K3+.csv",
"/home/miha/uni/research/ppp/data/data_2020-06-29_w_noise/K4+.csv",
"/home/miha/uni/research/ppp/data/data_2020-06-29_w_noise/K5+.csv",
"/home/miha/uni/research/ppp/data/data_2020-06-29_w_noise/K6+.csv",
"/home/miha/uni/research/ppp/data/data_2020-06-29_w_noise/K7+.csv",
"/home/miha/uni/research/ppp/data/data_2020-06-29_w_noise/K8+.csv",
"/home/miha/uni/research/ppp/data/data_2020-06-29_w_noise/K9+.csv",
"/home/miha/uni/research/ppp/data/data_2020-06-29_w_noise/K10+.csv",
"/home/miha/uni/research/ppp/data/data_2020-06-29_w_noise/K11+.csv",
"/home/miha/uni/research/ppp/data/data_2020-06-29_w_noise/K12+.csv"
]

'''
Give the charge states corresponding to data files. Must be in numerical order!
'''
chargeStates = [1,2,3,4,5,6,7,8,9,10,11,12]




def get_data_dict(chargeStates, dataFilePaths):
    '''
    Pack filepaths into a dictionary.

    Parameters
    ----------
    chargeStates : list
        List of charge states corresponding to data file paths.
    dataFilePaths : list
        List of data file paths.

    Returns
    -------
    dataDictionary : dict
        Data file paths organized into a dictionary with 
        the charge states as keys.

    '''
    #
    # Pack the data file paths to a dictionary for easy access later.
    #
    dataDictionary = {}
    for chargeState,dataFilePath in zip(chargeStates,dataFilePaths):
        dataDictionary[chargeState] = dataFilePath

    return dataDictionary



def get_interp_dict(dataDictionary):
    '''
    Read the data files in the data dictionary and interpolate the data 
    contained within. Then pack the interpolation objects into a 
    dictionary with the charge states as keys.

    Parameters
    ----------
    dataDictionary : dict
        Dictionary containing the data file paths.

    Returns
    -------
    interpDictionary : dict
        Dictionary containing the interpolate objects
        corresponding to data files given in dataDict.

    '''
    #
    # Interpolate the data files and pack them to a dictionary
    #
    interpDictionary = {}
    for key in dataDictionary:
        #
        # Get transient data and interpolate
        #
        df = pd.read_csv(dataDictionary[key], names=["time","current"])
        t = df["time"]
        i = df["current"]
        
        #
        # Get the interpolate object and pack to dictionary
        #
        I = interpolate.interp1d(t,i,kind="cubic")
        interpDictionary[key] = I

    return interpDictionary


    
def runge_kutta_DE(t,y,q,l,m,n,interpCurrents):
    '''
    The function, which defines the DE to be solved by the RK4 method.

    Parameters
    ----------
    t : float
        Time at which the DE is evaluated.
    y : float
        y(t).
    q : int
        Charge state around which the DE is solved.
    parameters : Dictionary
        The parameters l,m,n.
    interpCurrents : Dictionary
        The interpolated currents organized into a dictionary. 
        Keys must be the corresponding charge states (integer)! 

    Returns
    -------
    float
        Right hand side of the DE:
            d/dt I^q = l*I^(q-1) - b*I^q + c*I^(q+1).

    '''
    
    I_lo = interpCurrents[q-1](t)
    I_hi = interpCurrents[q+1](t)
    
    return l*I_lo - m*y + n*I_hi
    



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
    
    h = 1000e-6 # Time step (units of s) for the RK4 integrator. #TODO! Change this as per need.
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



def determineNoise(chargeState, dataDictionary):
    '''
    Determines the noise for the data corresponding to the key chargeState
    in the dataDictionary. The noise is determined by calculating the 
    standard deviation of the signal in the region t < 0 s.

    Parameters
    ----------
    chargeState : int
        Charge state corresponding to data of interest.
    dataDictionary : dict
        Dictionary of the data file paths.

    Returns
    -------
    noise : float
        Noise level estimation for the signal.

    '''
    
    q = chargeState
    
    #
    # Load the data
    #
    filename = dataDictionary[q]
    data=pd.read_csv(filename, names=["time","current"])
    t = data["time"]*1e3 # convert to ms
    i = data["current"]
    
    #
    # Determine the area where the noise is calculated
    #
    c = t < 0
    
    #
    # Calculate the noise level
    #
    noise = np.std(i[c])
    
    return noise


def determine_interp_limits(chargeState, dataDictionary):
    '''
    Determines the interpolation limits for the data files marked by the 
    data file paths in dataDictionary.

    Parameters
    ----------
    chargeState : int
        Charge state corresponding to data of interest.
    dataDictionary : dict
        Dictionary of the data file paths.


    Returns
    -------
    list
        List containing the lower and upper interpolation limit for the data.
        I.e. the largest time interval which is common to all signals.

    '''
    datas = {}
    for q in [chargeState-1,chargeState,chargeState+1]:
        datas[q] = dataDictionary[q]
        
    #
    # Find the largest common time instance for all signals
    # N.B. Fit must begin at t=0 (where also y0 = 0), and therefore t_lo = 0
    #
    t_his = []
    for key in datas:
        df = pd.read_csv(datas[key],names=["time","current"])
        time = df["time"]
        t_hi = time.values[-1]
        t_his.append(t_hi)

    t_hi = min(t_his)
    
    interpLimits = [0, t_hi]
    
    return interpLimits



def do_rk4(a,b,c,chargeState, interpDictionary, dataDictionary):
    '''
    Perform the Runge-Kutta integration for the given chargeState signal.

    Parameters
    ----------
    a : float
        Parameter a.
    b : float
        Parameter b.
    c : float
        Parameter b.
    chargeState : int
        Charge state corresponding to data of interest.
    dataDictionary : dict
        Dictionary of the data file paths.
    interpDictionary : dict
        Dictionary of the interpolate objects.
        
    Returns
    -------
    T : array
        Time array.
    Y : array
        Solution array.

    '''
    
    #
    # Get the necessary currents (q-1,q,q+1)
    #
    interpCurrents = {}
    for q in [chargeState-1, chargeState, chargeState+1]:
        interpCurrents[q] = interpDictionary[q]


    #
    # Determine interpolation limits
    # 
    interpLimits = determine_interp_limits(chargeState, dataDictionary)
    
    #
    # Get the initial values
    #
    t0=interpLimits[0]
    tf=interpLimits[-1]
    y0=0
    
    #
    # Do the RK4 integration
    #
    T,Y = rk4(t0=t0,y0=y0,tf=tf,f=runge_kutta_DE,q=chargeState,l=a,m=b,n=c,interpCurrents=interpCurrents)
    
    return T,Y





def fitting_function(t, a, b, c):
    '''
    Function to fit to the data. Defined as the solution to the differential
    equation of runge_kutta_DE.

    Parameters
    ----------
    t : float
        Moment in time.
    a : float
        Parameter a.
    b : float
        Parameter b.
    c : float
        Parameter c.

    Returns
    -------
    y : float
        Value of the fitting function at time t.

    '''
    global q, interpDict, dataDict
            
    T,Y = do_rk4(a,b,c,q,interpDict,dataDict)
    
    Y = interp1d(T, Y, kind="cubic",fill_value="extrapolate")
        
    y = Y(t)
    
    return y





#
# Run the optimization for desired charge states, and save results to file.
#

dataDict = get_data_dict(chargeStates, dataFilePaths)
interpDict = get_interp_dict(dataDict)

df_res = pd.DataFrame()

for q in chargeStates[1:-1]:
    
    print(q)
    
    #
    # Get the interpolation limits
    #
    interpLimits = determine_interp_limits(q, dataDict)
    lo, hi = interpLimits[0], interpLimits[1]
    
    #
    # Get the data to fit within limits
    #
    df = pd.read_csv(dataDict[q], names=["time", "current"])
    c = ( lo < df["time"])&(df["time"] < hi)
    
    df = df[c]

    # Limit fitting range
    hi = df["time"].values[-1] * 1.0
    c = df["time"] < hi    
    
    xdata = df["time"][c].values
    ydata = df["current"][c].values


    
    #
    # Get the noise to determine the uncertainty
    #
    noise = determineNoise(q,dataDict)
    noise = np.ones(len(xdata))*noise
    
    #
    # Do the fitting
    #
    lo_bnds = [0,0,0]
    hi_bnds = [np.inf,np.inf,np.inf]
    bnds = (lo_bnds, hi_bnds)
    popt,pcov = curve_fit(fitting_function, xdata, ydata, p0=[1000,1000,100],
                          sigma=noise, absolute_sigma=True,
                          method="trf", bounds=bnds)
    
    #
    # Calculate the chi2
    #
    chi2 = np.sum( np.square( (ydata - fitting_function(xdata, *popt))/noise ) )
    
    #
    # Calculate the (reduced) chi2
    #
    chi2_reduced = chi2 / (len(ydata)-3)
    
    #
    # Unpack the result
    #
    [a, b, c] = popt
    [da, db, dc] = np.sqrt(np.diag(pcov))
    
    df_res[q]=[a,b,c,da,db,dc,chi2_reduced]



df_res.index=["a", "b", "c", "da", "db", "dc", "chi2"]

print(df_res)

# df_res.to_csv("./results/solution_abc_2020-06-29_w_noise_h=10e-6.csv") #TODO! Change save to directory as you wish.








