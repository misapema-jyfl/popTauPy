#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 12:57:24 2021

@author: miha
"""

import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy import interpolate
# from parameters import d
import matplotlib.pyplot as plt
import matplotlib


# Set font for plots
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 15}

matplotlib.rc('font', **font)






def rk4(t0, y0, tf, f, h, **kwargs):
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
    h : float
        Time step of the method (s).
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
    
    # h = d["h"] # Time step (units of s) for the RK4 integrator. 
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



def RHS(t,y,a,b,c,I_lo,I_hi):
    '''
    Right-Hand-Side of the DE, which determines the time-evolution
    of the ion extraction current.
    
    Takes the a,b,c parameters and the interpolate objects
    for I_lo and I_hi (i.e. measurement data for I^q-1 & I^q+1) 
    as arguments.
    
    Returns the instantaneous value of dy/dt = RHS.
    
    Use with rk4 to solve for y(t).
    '''
    
    return a*I_lo(t) - b*y + c*I_hi(t)



class Fitter:
    
    def __init__(self, charge_state, params):
        
        self.params = params
        
        d = self.params["data"]
        o = self.params["optimizer"]
        p = self.params["parsing"]
        
        if not d["available_charge_states"] == False:
            cStates = d["available_charge_states"]
        else:
            cStates = p["available_charge_states"]
        
        if not d["parsed_data_path"] == False:
            fileLocDir = d["parsed_data_path"]
        else:
            fileLocDir = p["parsed_data_path"]
            
        species = d["species"].upper()
        
        self.h = o["rk_time_step"]

        self.charge_state = charge_state
        
        # Get paths to necessary files and charge state information.
        
        # Make the list of control signals
        # ctrl_files = []
        # for q in cStates:
        #     s = ("1+_{}".format(species),str(q),"+.csv")
        #     s = "".join(s)
        #     f = fileLocDir + s
        #     ctrl_files.append(f)
        
        # Make the list of N+ signals
        n_files = []
        for q in cStates:
            s = (species, str(q),"+.csv")
            s = "".join(s)
            f = fileLocDir + s
            n_files.append(f)
        self.dataFilePaths = n_files
        self.chargeStates = cStates
        
        #
        # Make a dataframe for holding data file paths
        #
        self.df = pd.DataFrame()
        self.df["cState"] = self.chargeStates
        self.df["fPath"] = self.dataFilePaths
        
        #
        # Make an interpolate dataframe out of three consecutive 
        # charge state currents around q
        #
        self.df_interp = self.make_interpolate_dataframe()
        
        
    def get_data(self, charge_state):
        '''
        Get measurement data corresponding to desired charge state.
        '''
        
        # charge_state = self.charge_state
        
        c = self.df["cState"] == charge_state
        data = pd.read_csv(self.df[c]["fPath"].values[0])
        
        ndf = pd.DataFrame(data) # Temporary dataframe for reading data to lists.
        ndf.columns = ["t", "i"]
        
        t,i = ndf["t"],ndf["i"]
        
        return t,i


    def get_interp_data(self, charge_state):
        '''
        Interpolate data corresponding to given charge state,
        and return an interpolate object for the data.
        '''
        
        # charge_state = self.charge_state
        
        t,i = self.get_data(charge_state)
        
        I = interpolate.interp1d(t,i,kind="cubic")
        
        return I
    
    
    def make_interpolate_dataframe(self):
        '''
        Make a dataframe out of the interpolate objects neighboring 
        a given charge state.
        '''
        
        charge_state = self.charge_state
        
        df_interp = pd.DataFrame()
        interps = []
        for q in [charge_state-1, charge_state, charge_state+1]:
            I = self.get_interp_data(q)
            interps.append(I)
        df_interp["cState"] = [charge_state-1, charge_state, charge_state+1]
        df_interp["interp"] = interps
        
        return df_interp


    def determine_interpolation_limits(self):
        '''
        Determines the interpolation limits for three consecutive transients,
        with t_i set to 0 s, always.
        '''    
        
        charge_state = self.charge_state
        
        t_lo, i_lo = self.get_data(charge_state-1)
        t_mid, i_mid = self.get_data(charge_state)
        t_hi, i_hi = self.get_data(charge_state+1)
        
        t_i = 0
        t_f = min( max(t_lo), max(t_hi), max(t_mid) )
    
        return t_i, t_f



    def fit_rk4(self, a, b, c):
        
        charge_state = self.charge_state
        
        #
        # Fetch the interpolate objects from the dataframe
        #
        I_lo = self.df_interp[ self.df_interp["cState"]==charge_state-1 ]["interp"].values[0]
        I_hi = self.df_interp[ self.df_interp["cState"]==charge_state+1 ]["interp"].values[0]
                
        t0, tf = self.determine_interpolation_limits()
        
        T,Y = rk4(t0=t0,y0=0,tf=tf,f=RHS,h=self.h,a=a,b=b,c=c,I_lo=I_lo,I_hi=I_hi)
        
        y = interp1d(T,Y,kind="cubic", fill_value="extrapolate")
        
        return y


    def fitting_function(self,t,a,b,c):
        
        y = self.fit_rk4(a,b,c)
        
        return y(t)

    

    def determine_noise(self):
        '''
        Determine the noise present in the measurement. 
        Calculated as the standard deviation of the signal
        in the region t<0, i.e. before onset of the 1+ injection pulse.
        '''
        t,i = self.get_data(self.charge_state)
        
        noise = np.std( i[t<0] )
        
        return noise



    def do_fit(self):
        
        #
        # Get the interpolation limits
        #
        t_i, t_f = self.determine_interpolation_limits()

        #
        # Get the data to fit 
        #
        t,i = self.get_data(self.charge_state)
        xdata, ydata = t[(t>t_i)&(t<t_f)], i[(t>t_i)&(t<t_f)]
        
        #
        # Determine the noise
        #
        noise = self.determine_noise()
        noise = np.ones(len(xdata))*noise

        lo_bnds = [0,0,0]
        hi_bnds = [np.inf,np.inf,np.inf]
        bnds = (lo_bnds, hi_bnds)
        p0 = [1000,1000,100]
        popt,pcov = curve_fit(f=self.fitting_function,
                              xdata=xdata,
                              ydata=ydata,
                              p0=p0,
                              sigma=noise,
                              absolute_sigma=True,
                              method="trf",
                              bounds=bnds)  
        
        # Calculate the chi2
        chi2 = np.sum( np.square( (ydata - self.fitting_function(xdata, *popt))/noise ) )
        
        # Calculate the (reduced) chi2
        chi2_reduced = chi2 / (len(ydata)-3)


        return popt, pcov, chi2_reduced




