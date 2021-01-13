#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  5 10:53:37 2020

@author: miha
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib import cm
from scipy.ndimage.filters import gaussian_filter
import matplotlib

from parameters import plotting_results as pr
from parameters import p



# Set font for plots
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 15}

matplotlib.rc('font', **font)






class Data:
    
   
    def __init__(self):
        
        self.fNames = pr["solution_set_files"]
        self.cStates= pr["charge_states"]  
        self.outDir = pr["output_directory"]
        
        # Set plotting limits
        self.Te_lo = pr["Te_lo"]
        self.Te_hi = pr["Te_hi"]
        self.ne_lo = pr["ne_lo"]
        self.ne_hi = pr["ne_hi"]
        self.F_hi = pr["F"]
        self.conf = pr["conf"]
        
        # 
        # Get solution set files into a dataframe,
        # organized according to the charge state.
        #
        self.dataFiles = pd.DataFrame({"q": self.cStates, "f": self.fNames})

    

    def set_limits(self, df):
        
        # Set the desired ne, Te, F limits on the data
        df = df[df["F"]<self.F_hi]
        df = df[(df["n"]<self.ne_hi)&(df["n"]>self.ne_lo)]
        df = df[(df["T"]<self.Te_hi)&(df["T"]>self.Te_lo)]
        
        return df
    
    

    def get_df(self,q):
            
            # Retrieve data from file corresponding to charge state q
            df = pd.read_csv(self.dataFiles["f"][self.dataFiles["q"]==q].values[0],index_col=0)
            
            # Make certain there are no negative solutions.        
            df = df[df["tau"] > 0]
            df = df[df["cx_rate"] > 0]
            df = df[df["inz_rate"] > 0]
            
            # Set upper and lower limits on the data
            df = self.set_limits(df)
            
            return df
        
    
    
    def plot_heatmap_solution_set(self, q, fig, ax):
        '''
        Plot a heatmap of solutions for charge state q
        '''
        
        df = self.get_df(q)
        
        x = df["T"]
        y = df["n"]
        
        # rng = [[10,10e3],[1e11,2.61e12]]
        rng = [ [p["Te_lo"],p["Te_hi"]], [10**p["ne_lo"],10**p["ne_hi"]]]        # TODO! Note how awful this is.
        
        heatmap, xedges, yedges = np.histogram2d(x, y, 
                                                 bins=1000, 
                                                 range=rng,
                                                 density=True)
        
        heatmap = gaussian_filter(heatmap, sigma=12)
        
        img = heatmap.T
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
        
        im = ax.imshow(img,origin='lower',
                       cmap=cm.hot,
                       extent=extent,
                       aspect="auto",
                       interpolation="none")
        
        cb = fig.colorbar(im)
        
        # Scatter plot the data points
        ax.scatter(x,y,s=.05,c="w")
        
        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_xlim(left=p["Te_lo"],right=p["Te_hi"])                             # TODO! This is similarly awful.
        ax.set_ylim(bottom=10**p["ne_lo"],top=10**p["ne_hi"])
        ax.set_xlabel(r"$T_e$ (eV)")
        ax.set_ylabel(r"$n_e$ (cm$^{-3}$)")
        
        
        
d = Data()
for q in d.cStates:
    fig, ax = plt.subplots()
    d.plot_heatmap_solution_set(q,fig,ax)    
    fig.tight_layout()
    plt.savefig( d.outDir + "solution_set_q={}+.eps".format(str(q)),format="eps")
    plt.close()
    
    
    
    
    
    
    