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

# Set font for plots
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 15}

matplotlib.rc('font', **font)



fNames = [
"/home/miha/uni/research/ppp/code/popTauPy/results/unc_lo=-60%_unc_hi=60%_MC_iters=1000_N=1000_q=5.csv",
"/home/miha/uni/research/ppp/code/popTauPy/results/unc_lo=-60%_unc_hi=60%_MC_iters=1000_N=1000_q=6.csv",
"/home/miha/uni/research/ppp/code/popTauPy/results/unc_lo=-60%_unc_hi=60%_MC_iters=1000_N=1000_q=7.csv",
"/home/miha/uni/research/ppp/code/popTauPy/results/unc_lo=-60%_unc_hi=60%_MC_iters=1000_N=1000_q=8.csv",
"/home/miha/uni/research/ppp/code/popTauPy/results/unc_lo=-60%_unc_hi=60%_MC_iters=1000_N=1000_q=9.csv",
"/home/miha/uni/research/ppp/code/popTauPy/results/unc_lo=-60%_unc_hi=60%_MC_iters=1000_N=1000_q=10.csv"
]


cStates=[5,6,7,8,9,10]



class Data:
    
    # Data limits 
    F_hi = 1e-4
    ne_lo = 0
    ne_hi = +np.inf
    Te_lo = 0
    Te_hi = +np.inf    
    
    def __init__(self, fNames, cStates):
            
        self.dataFiles = pd.DataFrame({"q": cStates, "f": fNames})

    

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
        
        rng = [[10,10e3],[1e11,2.61e12]]
        
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
        ax.set_xlim(left=10,right=10e3)
        ax.set_ylim(bottom=1e11,top=2.61e12)
        ax.set_xlabel(r"$T_e$ (eV)")
        ax.set_ylabel(r"$n_e$ (cm$^{-3}$)")
        fig.tight_layout()
        
        plt.savefig("./results/solution_set_q={}+.png".format(str(q)),dpi=500)
        plt.close()
        
        
d = Data(fNames, cStates)
# cStates=[5]
for q in cStates:
    fig, ax = plt.subplots()
    d.plot_heatmap_solution_set(q,fig,ax)    
    
    
    
    
    
    
    