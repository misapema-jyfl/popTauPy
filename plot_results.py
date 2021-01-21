#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 08:44:31 2021

@author: miha
"""

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from parameters import plotting_results as pr


# Set font for plots
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 15}

matplotlib.rc('font', **font)






class Plotter:
    
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
        
        if not xs.size > 0:
            print("Given list of values is empty!")
            return np.nan, np.nan, np.nan
        
        else:
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
            
            if not percentages.size > 0:
                print("There are no elements within the confidence interval!")
                x_hi = 0
            else:    
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
            if not percentages.size > 0:
                print("There are no elements within the confidence interval!")
                x_lo = 0
            else:    
                result = np.where(percentages==min(percentages))
                idx = result[0][0]
                x_lo = x_los[idx]
        
        return x_lo, median, x_hi
    
    
    
    
    def get_df(self, q):
        '''
        Read data from file corresponding to charge state q
        into a dataframe. 
        
        Take out any unphysical solutions if such exist.
        '''
        
        # Retrieve data from file corresponding to charge state q
        df = pd.read_csv(self.dataFiles["f"][self.dataFiles["q"]==q].values[0],index_col=0)
        
        # Make certain there are no negative solutions.        
        df = df[df["tau"] > 0]
        df = df[df["cx_rate"] > 0]
        df = df[df["inz_rate"] > 0]
        
        return df
    
    
    def set_limits(self, df):
        '''
        Prune out undesired data points,
        as determined by plotting limits.
        '''
        
        # Set the desired ne, Te, F limits on the data
        # df = df[df["F"]<self.F_hi]
        df = df[(df["n"]<self.ne_hi)&(df["n"]>self.ne_lo)]
        df = df[(df["T"]<self.Te_hi)&(df["T"]>self.Te_lo)]
        
        return df
    
    
    def set_F_upper(self,df):
        '''
        Set the upper limit of F from the parameters.py
        '''
        
        df = df[df["F"]<self.F_hi]
        
        return df
    
    
    def plot_number_of_solutions(self, Fs, q):
        '''
        Plot the number of solutions as a function of the 
        upper limit of the penalty function.
        '''
        
        nSols=[]
        fig, ax = plt.subplots()
        for F in Fs:
            
            df = self.get_df(q)
            df = self.set_limits(df)
            df = df[df["F"]<F]
            
            nSols.append(len(df))
        
        ax.scatter(Fs, nSols, marker="s", s=48, color="k")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$F$")
        ax.set_ylabel("Number of solutions")
        ax.set_xticks(Fs)
        
        return fig, ax
    



    def plot_times_against_F(self, Fs, q, key, marker, color):
        '''
        Plot the plasma characteristic time against the upper limit
        of the penalty function value.
        '''
        
        lo_errs=[]
        hi_errs=[]
        medians=[]
        
        fig, ax = plt.subplots()
        
        for F in Fs:
            
            
            df = self.get_df(q)
            df = self.set_limits(df)
            df = df[df["F"]<F]
            
            #
            # Select the characteristic time to plot
            #
            if key == "inz_time":
                data=np.array(df["inz_rate"])
                data = data**(-1)
                data = 1e3*data
            elif key == "cx_time":
                data=np.array(df["cx_rate"])
                data = data**(-1)
                data = 1e3*data
            elif key == "tau":
                data=np.array(df["tau"])*1e3
            else:
                print("Erroneous key!")
                break
            
            #
            # Calculate the error bars and median values
            #
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
        ax.set_xlabel(r"$F$")
        ax.set_xticks(Fs)
        ax.set_ylim(bottom=0)
        
        return fig, ax


    def plot_time_against_q(self, key, qs, marker, color):
        '''
        Plot the plasma characteristic time vs charge state.
        '''
        
        lo_errs=[]
        hi_errs=[]
        medians=[]
        fig, ax = plt.subplots()
        for q in qs:
            
            df = self.get_df(q)
            df = self.set_limits(df)
            df = self.set_F_upper(df)
            
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
                
        if key == "inz_time":
            
            ax.set_ylabel(r"$[\left\langle\sigma v\right\rangle^{\mathrm{inz}}_{q\to q+1} n_e]^{-1}$ (ms)")
            
        elif key == "cx_time":
            
            ax.set_ylabel(r"$[\left\langle\sigma v\right\rangle^{\mathrm{cx}}_{q\to q-1} n_0]^{-1}$ (ms)")

        else:
            ax.set_ylabel(r"$\tau^q$ (ms)")
                
        ax.set_xlabel("Charge state")
        ax.set_xticks(qs)
        ax.legend()
        
        return fig, ax
    
    
    
    
    def plot_Ee_against_q(self, qs, marker, color):
        '''
        Plot the plasma average energy content vs q.
        '''
        lo_errs=[]
        hi_errs=[]
        medians=[]
        fig, ax = plt.subplots()
        for q in qs:
            
            df = self.get_df(q)
            df = self.set_limits(df)
            
            data = np.array(df["n"])*np.array(df["T"]*1.5) # Convert nT to <Ee>
                
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
    
    
    def plot_triple_against_q(self, qs, marker, color):
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
        
        return fig, ax


    
    def output_triple_results(self,qs):
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
    
    
    
    def output_time_results(self, qs, key):
        '''
        Output characteristic times vs q as .csv 
        
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


    def output_Ee_results(self,qs):
        '''
        Output average energy content as a function of charge state as .csv
        '''
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




#
# Plot the number of solutions vs upper limit of penalty function value F
#
pf = pr["plotting_vs_F"]
if pf["plot_num_of_solutions_vs_F"] == True:
    P = Plotter()
    for charge_state in pr["charge_states"]:
        
        fig, ax = P.plot_number_of_solutions(Fs=pf["list_of_Fs"], q=charge_state)
        
        errFlag = 0
        if pf["y_lo"] <= 0:
            print("Can't set non-positive bottom limit on log-scale.")    
            errFlag = 1
            ax.set_ylim(bottom = 1, top = pf["y_hi"])
        else:
            ax.set_ylim(bottom = pf["y_lo"], top = pf["y_hi"])
            
        fig.tight_layout()
        fig.savefig(P.outDir + "fig_number_of_solutions_vs_F_q={}.eps".format(charge_state), format="eps")
        
    if errFlag == 1:
        print("Check y_lo in plotting_vs_F in parameters.py!")



#
# Plot the characteristic time vs upper limit of penalty function value F
#
pcf = pr["plotting_time_vs_F"]
if pcf["plot_or_not"] == True:
    P = Plotter()
    for charge_state in pr["charge_states"]:
        
        # Plot confinement times vs F
        if pcf["plot_tau"] == True:
            fig, ax = P.plot_times_against_F(Fs = pcf["list_of_Fs"],
                                             q = charge_state,
                                             key = "tau",
                                             marker = pcf["tau_marker"],
                                             color = pcf["tau_color"]
                                             )
            ax.set_ylabel("Confinement time (ms)")
            fig.tight_layout()
            fig.savefig(P.outDir + "fig_CONF_n_of_sols_vs_F_q={}.eps".format(charge_state), format="eps")

        # Plot ionization times vs F
        if pcf["plot_inz"] == True:
            fig, ax = P.plot_times_against_F(Fs = pcf["list_of_Fs"],
                                             q = charge_state,
                                             key = "inz_time",
                                             marker = pcf["inz_marker"],
                                             color = pcf["inz_color"]
                                             )
            ax.set_ylabel("Ionisation time (ms)")
            fig.tight_layout()
            fig.savefig(P.outDir + "fig_INZ_n_of_sols_vs_F_q={}.eps".format(charge_state), format="eps")
            
        # Plot charge exchange times vs F
        if pcf["plot_cx"] == True:
            fig, ax = P.plot_times_against_F(Fs = pcf["list_of_Fs"],
                                             q = charge_state,
                                             key = "cx_time",
                                             marker = pcf["cx_marker"],
                                             color = pcf["cx_color"]
                                             )
            ax.set_ylabel("Charge exchange time (ms)")
            fig.tight_layout()
            fig.savefig(P.outDir + "fig_CX_n_of_sols_vs_F_q={}.eps".format(charge_state), format="eps")





#
# Plot the characteristic time vs charge state
#
pcq = pr["plotting_time_vs_q"]
if pcq["plot_or_not"] == True:
    
    charge_states = pr["charge_states"]
    
    P = Plotter()    
    
    # Plot confinement times vs F
    if pcq["plot_tau"] == True:
        fig, ax = P.plot_time_against_q( qs = charge_states,
                                         key = "tau",
                                         marker = pcq["tau_marker"],
                                         color = pcq["tau_color"]
                                         )
        ax.set_ylabel("Confinement time (ms)")
        ax.set_yscale(pcq["tau_yscale"])
        
        if not pcq["tau_yscale"] == "log":
            ax.set_ylim(bottom=0)
        
        fig.tight_layout()
        fig.savefig(P.outDir + "fig_time_CONF_vs_q.eps", format="eps")

    # Plot ionization times vs F
    if pcq["plot_inz"] == True:
        fig, ax = P.plot_time_against_q( qs = charge_states,
                                         key = "inz_time",
                                         marker = pcq["inz_marker"],
                                         color = pcq["inz_color"]
                                         )
        ax.set_ylabel("Ionisation time (ms)")
        ax.set_yscale(pcq["inz_yscale"])
        if not pcq["tau_yscale"] == "log":
            ax.set_ylim(bottom=0)
        fig.tight_layout()
        fig.savefig(P.outDir + "fig_time_INZ_vs_q.eps", format="eps")
        
    # Plot charge exchange times vs F
    if pcq["plot_cx"] == True:
        fig, ax = P.plot_time_against_q(qs = charge_states,
                                         key = "cx_time",
                                         marker = pcq["cx_marker"],
                                         color = pcq["cx_color"]
                                         )
        ax.set_ylabel("Charge exchange time (ms)")
        ax.set_yscale(pcq["cx_yscale"])
        if not pcq["tau_yscale"] == "log":
            ax.set_ylim(bottom=0)
        fig.tight_layout()
        fig.savefig(P.outDir + "fig_time_CX_vs_q.eps", format="eps")
        
        


#
# Plot average energy content vs charge state
#
pce = pr["plotting_Ee_vs_q"]
if pce["plot_or_not"] == True:
    P = Plotter()
    fig, ax = P.plot_Ee_against_q(qs=pr["charge_states"],
                        marker=pce["marker"],
                        color=pce["color"])
    ax.set_ylim(bottom=pce["y_lo"], top=pce["y_hi"])
    fig.tight_layout()
    fig.savefig(P.outDir + "fig_Ee_vs_q.eps", format="eps")
    
    
    
#
# Plot triple product vs charge state
#
pct = pr["plotting_triple_vs_q"]
if pct["plot_or_not"] == True:
    P = Plotter()
    fig, ax = P.plot_triple_against_q(qs=pr["charge_states"],
                                      marker=pct["marker"],
                                      color=pct["color"])
    
    ax.set_ylim(bottom=pct["y_lo"], top=pct["y_hi"])
    fig.tight_layout()
    fig.savefig(P.outDir + "fig_triple_vs_q.eps", format="eps")
    
    
'''
Output characteristic times in a .csv file.
'''
if pr["output_tau_vs_q"] == True:
    P = Plotter()
    df_out=P.output_time_results(qs=pr["charge_states"], key="tau")
    df_out.to_csv(P.outDir + "output_time_CONF_vs_q.csv")

if pr["output_inz_vs_q"] == True:
    P = Plotter()
    df_out=P.output_time_results(qs=pr["charge_states"], key="inz_time")
    df_out.to_csv(P.outDir + "output_time_INZ_vs_q.csv")
    
if pr["output_cx_vs_q"] == True:
    P = Plotter()
    df_out=P.output_time_results(qs=pr["charge_states"], key="cx_time")
    df_out.to_csv(P.outDir + "output_time_CX_vs_q.csv")
    
    
'''
Output average energy content in a .csv file
'''
if pr["output_Ee_vs_q"] == True:
    P = Plotter()
    df_out=P.output_Ee_results(qs=pr["charge_states"])
    df_out.to_csv(P.outDir + "output_Ee_vs_q.csv")
    
    
    
'''
Output triple products in a .csv file.
'''
if pr["output_triple_vs_q"] == True:
    P = Plotter()
    df_out=P.output_triple_results(qs=pr["charge_states"])
    df_out.to_csv(P.outDir + "output_triple_vs_q.csv")
    
