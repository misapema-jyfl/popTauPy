#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 10:53:26 2020

@author: miha
"""

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


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
    
    # The plotting limits.
    F_hi = 1e-4
    ne_lo = 0
    ne_hi = +np.inf
    Te_lo = 0
    Te_hi = +np.inf
    
    # Confidence interval
    conf = 0.341
    
    
    def __init__(self, fNames, cStates):
        
        self.dataFiles = pd.DataFrame({"q": cStates, "f": fNames})
    
    
    
    def find_confidence_interval(self, list_of_values, condition_percentage):
        '''
        Seek the lower and upper limit in a list of values, which enclose 
        between themselves and the median a given percentage of all values 
        in the list. 
        
    
        Parameters
        ----------
        list_of_values : list
            List of the values, for which the analysis will be carried out.
        condition_percentage : float
            The percentage limit.
    
        Returns
        -------
        x_lo : float
            Lower limit.
        median : float
            Median value.
        x_hi : float
            Upper limit.
    
        '''
        xs = list_of_values
        p = condition_percentage
        
        median = np.median(xs)
        percentages = []
        x_his = []
        for x in xs:
            if x > median:
                # Select elements between current element and median
                interval = xs[(xs<=x)&(xs>=median)]
                # Calculate percentage of all values within current interval
                percentage = (len(interval))/len(xs)
                # If this interval satisfies the condition, add it to list
                if percentage >= p:
                    percentages.append(percentage)
                    x_his.append(x)
        # Find the minimum percentage satisfying the condition
        # along with the corresponding element.
        percentages=np.array(percentages)
        result = np.where(percentages==min(percentages))
        idx = result[0][0]
        x_hi = x_his[idx]
        
        # Find value x_lo, for which p fraction of all results are between 
        # it and the median
        percentages = []
        x_los = []
        for x in xs:
            if x < median:
                # Select elements between current element and median
                interval = xs[(xs>=x)&(xs<=median)]
                # Calculate percentage of all values within current interval
                percentage = (len(interval))/len(xs)
                # If this interval satisfies the condition, add it to list
                if percentage >= p:
                    percentages.append(percentage)
                    x_los.append(x)
        # Find the minimum percentage satisfying the condition
        # along with the corresponding element.
        percentages=np.array(percentages)
        result = np.where(percentages==min(percentages))
        idx = result[0][0]
        x_lo = x_los[idx]
        
        return x_lo, median, x_hi
    
    
    
    def get_df(self,q):
        
        # Retrieve data from file corresponding to charge state q
        df = pd.read_csv(self.dataFiles["f"][self.dataFiles["q"]==q].values[0],index_col=0)
        
        # Make certain there are no negative solutions.        
        df = df[df["tau"] > 0]
        df = df[df["cx_rate"] > 0]
        df = df[df["inz_rate"] > 0]
        
        return df
        


    def set_limits(self, df):
        
        # Set the desired ne, Te, F limits on the data
        df = df[df["F"]<self.F_hi]
        df = df[(df["n"]<self.ne_hi)&(df["n"]>self.ne_lo)]
        df = df[(df["T"]<self.Te_hi)&(df["T"]>self.Te_lo)]
        
        return df


    def plot_solution_set(self, q, fig, ax):
        
        df = self.get_df(q)
        
        df = self.set_limits(df)
        
        n, T = df["n"], df["T"]
        
        lbl_charge = "{"+str(q)+"+}"
        lbl = r"K$^{}$".format(lbl_charge)
        
        ax.scatter(T, n, label=lbl)
        
        if not self.ne_hi == np.inf:
            ax.set_ylim(self.ne_lo,self.ne_hi)
        if not self.Te_hi == np.inf:
            ax.set_xlim(self.Te_lo,self.Te_hi)
        ax.set_yscale("log")
        
        return fig, ax

    def plot_solution_sets(self, qs):
        
        fig, ax = plt.subplots()
        for q in qs:
            fig, ax = self.plot_solution_set(q, fig, ax)
        ax.set_ylabel(r"$n_e$ (1/cm$^{3}$)")
        ax.set_xlabel(r"$T_e$ (eV)")
        ax.legend(bbox_to_anchor=[.1,1.],ncol=3)
        fig.tight_layout()
        
        return fig, ax


    def plot_number_of_solutions(self, Fs, q):
        
        nSols=[]
        fig, ax = plt.subplots()
        for F in Fs:
            Data.F_hi = F
            df = self.get_df(q)
            df = self.set_limits(df)
            nSols.append(len(df))
        
        ax.scatter(Fs, nSols, marker="s", s=48, color="k")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$F$")
        ax.set_ylabel("Number of solutions")
        ax.set_xticks(Fs)
        
        return fig, ax
    
    
    def plot_time_against_F(self, Fs, q, key, marker, color):
        
        lo_errs=[]
        hi_errs=[]
        medians=[]
        fig, ax = plt.subplots()
        for F in Fs:
            
            Data.F_hi=F
            
            df = self.get_df(q)
            df = self.set_limits(df)
            
            
            if key == "inz_time":
                data=np.array(df["inz_rate"])
                data = data**(-1)
                data = 1e3*data
            elif key == "cx_time":
                data=np.array(df["cx_rate"])
                data = data**(-1)
                data = 1e3*data
            else:
                data=np.array(df[key])*1e3
            
            
            lo,median,hi = self.find_confidence_interval(data,self.conf)
            
            lo_err = median-lo
            hi_err = hi-median
            
            lo_errs.append(lo_err)
            hi_errs.append(hi_err)
            medians.append(median)
            
        ax.errorbar(x=np.array(Fs),y=medians,yerr=[lo_errs,hi_errs],
                    fmt="",
                    ls="",
                    lw=2,
                    capsize=8,
                    marker=marker,
                    markersize=8,
                    color=color)
        ax.set_xscale("log")
        ax.set_ylabel("Milliseconds")
        ax.set_xlabel(r"$F$")
        ax.set_xticks(Fs)
        ax.set_ylim(bottom=0)
        
        return fig, ax
    
    
    def plot_time_against_q(self, key, qs, marker, color):
        '''
        Make sure to set F_hi (and other limits) as you wish before 
        running this plotting!
        '''
        lo_errs=[]
        hi_errs=[]
        medians=[]
        fig, ax = plt.subplots()
        for q in qs:
            
            df = self.get_df(q)
            df = self.set_limits(df)
            
            
            if key == "inz_time":
                data=np.array(df["inz_rate"])
                data = data**(-1)
                data = 1e3*data
                lbl = r"$[\left\langle\sigma v\right\rangle^{\mathrm{inz}}_{q\to q+1} n_e]^{-1}$"
                
            elif key == "cx_time":
                data=np.array(df["cx_rate"])
                data = data**(-1)
                data = 1e3*data
                lbl = r"$[\left\langle\sigma v\right\rangle^{\mathrm{cx}}_{q\to q-1} n_0]^{-1}$"

            else:
                data=np.array(df[key])*1e3
                lbl = r"$\tau^q$"
                
            lo,median,hi = self.find_confidence_interval(data,self.conf)
            
            lo_err = median-lo
            hi_err = hi-median
            
            lo_errs.append(lo_err)
            hi_errs.append(hi_err)
            medians.append(median)
            
        ax.errorbar(x=np.array(qs),y=medians,yerr=[lo_errs,hi_errs],
                    fmt="",
                    ls="",
                    lw=2,
                    capsize=8,
                    marker=marker,
                    markersize=8,
                    color=color,
                    label=lbl)
        
        # ax.set_xscale("log")
        ax.set_ylabel("Milliseconds")
        ax.set_xlabel("Charge state")
        ax.set_xticks(qs)
        ax.set_ylim(bottom=0)
        ax.legend()
        fig.tight_layout()
        
        return fig, ax
    
    def output_time_results(self, qs, key):
        '''
        N.B. all times are output in units of milliseconds!
        '''
        lo_errs=[]
        hi_errs=[]
        medians=[]
        for q in qs:
            
            df = self.get_df(q)
            df = self.set_limits(df)
            
            if key == "inz_time":
                data=np.array(df["inz_rate"])
                data = data**(-1)
                data = 1e3*data
                
            elif key == "cx_time":
                data=np.array(df["cx_rate"])
                data = data**(-1)
                data = 1e3*data

            else:
                data=np.array(df[key])*1e3
                
            lo,median,hi = self.find_confidence_interval(data,self.conf)
            
            lo_err = median-lo
            hi_err = hi-median
            
            lo_errs.append(lo_err)
            hi_errs.append(hi_err)
            medians.append(median)
    
        # Output file
        df_out = pd.DataFrame()
        df_out["qs"]=qs
        df_out["lo_errs"]=lo_errs
        df_out["medians"]=medians
        df_out["hi_errs"]=hi_errs
        
        return df_out
    
    def plot_nT_against_q(self, qs, marker, color):
        '''
        Make sure to set F_hi (and other limits) as you wish before 
        running this plotting!
        '''
        lo_errs=[]
        hi_errs=[]
        medians=[]
        fig, ax = plt.subplots()
        for q in qs:
            
            df = self.get_df(q)
            df = self.set_limits(df)
            
            data = np.array(df["n"])*np.array(df["T"])
                
            lo,median,hi = self.find_confidence_interval(data,self.conf)
            
            lo_err = median-lo
            hi_err = hi-median
            
            lo_errs.append(lo_err)
            hi_errs.append(hi_err)
            medians.append(median)
            
        ax.errorbar(x=np.array(qs),y=medians,yerr=[lo_errs,hi_errs],
                    fmt="",
                    ls="",
                    lw=2,
                    capsize=8,
                    marker=marker,
                    markersize=8,
                    color=color,
                    )
        
        # ax.set_xscale("log")
        ax.set_ylabel(r"$n_e \left\langle E_e \right\rangle$ (eV/cm$^{3}$)")
        ax.set_xlabel("Charge state")
        ax.set_xticks(qs)
        ax.set_yscale("log")
        fig.tight_layout()
        
        return fig, ax
    
    def output_nT_results(self,qs):
        lo_errs=[]
        hi_errs=[]
        medians=[]
        for q in qs:
            
            df = self.get_df(q)
            df = self.set_limits(df)
            
            data = np.array(df["n"])*np.array(df["T"])*1.5
                
            lo,median,hi = self.find_confidence_interval(data,self.conf)
            
            lo_err = median-lo
            hi_err = hi-median
            
            lo_errs.append(lo_err)
            hi_errs.append(hi_err)
            medians.append(median)
            
        # Output file
        df_out = pd.DataFrame()
        df_out["qs"]=qs
        df_out["lo_errs"]=lo_errs
        df_out["medians"]=medians
        df_out["hi_errs"]=hi_errs
        
        return df_out
    
    def plot_nTtau_against_q(self, qs, marker, color):
        '''
        Make sure to set F_hi (and other limits) as you wish before 
        running this plotting!
        '''
        lo_errs=[]
        hi_errs=[]
        medians=[]
        fig, ax = plt.subplots()
        for q in qs:
            
            df = self.get_df(q)
            df = self.set_limits(df)
            
            data = np.array(df["n"])*np.array(df["T"])*np.array(df["tau"])*1.5
                
            lo,median,hi = self.find_confidence_interval(data,self.conf)
            
            lo_err = median-lo
            hi_err = hi-median
            
            lo_errs.append(lo_err)
            hi_errs.append(hi_err)
            medians.append(median)
            
        ax.errorbar(x=np.array(qs),y=medians,yerr=[lo_errs,hi_errs],
                    fmt="",
                    ls="",
                    lw=2,
                    capsize=8,
                    marker=marker,
                    markersize=8,
                    color=color,
                    )
        
        # ax.set_xscale("log")
        ax.set_ylabel(r"$n_e \left\langleE_e\right\rangle\tau^q$ (eVs/cm$^{3}$)")
        ax.set_xlabel("Charge state")
        ax.set_xticks(qs)
        ax.set_yscale("log")
        fig.tight_layout()
        
        return fig, ax
    
    def output_nTtau_results(self,qs):
        lo_errs=[]
        hi_errs=[]
        medians=[]
        for q in qs:
            
            df = self.get_df(q)
            df = self.set_limits(df)
            
            data = np.array(df["n"])*np.array(df["T"])*np.array(df["tau"])*1.5
                
            lo,median,hi = self.find_confidence_interval(data,self.conf)
            
            lo_err = median-lo
            hi_err = hi-median
            
            lo_errs.append(lo_err)
            hi_errs.append(hi_err)
            medians.append(median)
            
        # Output file
        df_out = pd.DataFrame()
        df_out["qs"]=qs
        df_out["lo_errs"]=lo_errs
        df_out["medians"]=medians
        df_out["hi_errs"]=hi_errs
        
        return df_out
    
    
d = Data(fNames, cStates)



'''
Plot characteristic time vs F.
Use the keys "tau", "cx_time", "inz_time" to plot each respectively.
Make sure to change figure name appropriately!
Remember to change the marker!
Remember to change the markercolor!
'''
for q in[5,6,7,8,9,10]:
    fig,ax=d.plot_time_against_F(Fs=[1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6],
                                  q=q,
                                  key="tau",
                                  marker=".",
                                  color="r")
    fig.tight_layout()
    plt.savefig("./results/time_{}_vs_F_q={}.eps".format("CONF",q),format="eps")
    plt.close(fig)


'''
Plot number of solutions vs F.
'''
# for q in [5,6,7,8,9,10]:
#     fig, ax = d.plot_number_of_solutions(Fs=[1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6], q=q)
#     ax.set_ylim(bottom=10,top=1e5)
#     fig.tight_layout()
#     plt.savefig("./results/number_of_solutions_vs_F_q={}.eps".format(q), format="eps")
#     plt.close(fig)


'''
Plot characteristic time vs charge state, with error bars.
Use the keys "tau", "cx_time", "inz_time" to plot each respectively.
Remember to change the plot name!
Remember to change the marker!
Remember to change the markercolor!
Remember to change the ylim!
Remember to change the yscale!
'''
# Data.F_hi = 1e-4
# fig, ax = d.plot_time_against_q("tau",
#                                 qs=[5,6,7,8,9,10],
#                                 marker=".", 
#                                 color="r")
# # ax.set_ylim(bottom=1,top=1e3)
# # ax.set_yscale("log")
# plt.savefig("./results/fig_{}_vs_q.eps".format("CONF"), format="eps")
# plt.close()

'''
Output characteristic times in a .csv file.
'''
# df_out=d.output_time_results(qs=[5,6,7,8,9,10], key="tau")
# df_out.to_csv("./results/output_time_{}_vs_q.csv".format("CONF"))

'''
Plot plasma energy contents against charge state.
'''
# Data.F_hi=1e-4
# fig, ax = d.plot_nT_against_q(qs=[5,6,7,8,9,10], marker="s", color="k")
# ax.set_ylim(bottom=1e13,top=1e16)
# plt.savefig("./results/fig_nT_vs_q.eps", format="eps")

'''
Output energy contents in a .csv file.
'''
# df_out=d.output_nT_results(qs=[5,6,7,8,9,10])
# df_out.to_csv("./results/output_nT_vs_q.csv")



'''
Plot plasma triple product against charge state.
'''
# Data.F_hi=1e-4
# fig, ax = d.plot_nTtau_against_q(qs=[5,6,7,8,9,10], marker="s", color="m")
# ax.set_ylim(bottom=1e11,top=1e14)
# plt.savefig("./results/fig_nTtau_vs_q.eps",format="eps")


'''
Output triple products in a .csv file.
'''
# df_out=d.output_nTtau_results(qs=[5,6,7,8,9,10])
# df_out.to_csv("./results/output_nTtau_vs_q.csv")










'''
Plot solution sets.
N.B! This is old! Don't use this! Use plot_solution_sets.py!
'''
# Data.F_hi=1e-4
# fig, ax = d.plot_solution_sets(qs=[5,6,7,8,9,10])
# ax.set_xscale("log")
# plt.savefig("./results/fig_solution_sets.eps",format="eps")

