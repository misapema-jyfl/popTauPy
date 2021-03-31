#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 08:44:31 2021

@author: miha
"""

import seaborn as sns
# sns.set()
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.ndimage.filters import gaussian_filter
import pandas as pd
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
# from parameters import plotting_results as pr
# from parameters import p


# Set font for plots
font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 15}

matplotlib.rc('font', **font)


def get_histogram(data, binCount, lo, hi, spacing="log"):
    '''
    Takes in an array of data and determines the histogram of the array. 
    
    Parameters:
        data: The data array
        
        binCount: The number of bins into which to distribute the data.
        
        lo: Lower limit of the data values.
        
        hi: Upper limit of the data values.
        
        spacing: One of ["linear", "log"]
    '''
    if spacing == "linear":
        bins = np.linspace(lo, hi, binCount)
    elif spacing == "log":
        bins = np.logspace(np.log10(lo), np.log10(hi), binCount)
    
    # Sort
    data = np.sort(data)

    histogramData = []
    for j in range(binCount):
        
        if j+1 == len(bins):
            histogramData.append(0)
            break
        
        rBin = bins[j+1]
        lBin = bins[j]

        binned = 0
        
        # TODO! 
        # It's unnecessary to iterate over the entire array
        # each time, since the data is sorted.
        # Make this smarter!
        for i,x in enumerate(data):
            if x <= rBin and x > lBin:
                binned+=1
        
        histogramData.append(binned)
    
    return bins, histogramData



class SolSetPlotter:
    
   
    def __init__(self, params):
        
        
        self.plotting = params["plotting"]
        self.optimizer = params["optimizer"]
        self.data = params["data"]
        self.g = params["general"]
        self.solSettings = self.plotting["plot_solution_sets"]
        self.parsing = params["parsing"]
        # Set the charge states for which data is available
        if self.plotting["available_charge_states"] == False:
            if self.data["available_charge_states"] == False:
                self.cStates = self.parsing["available_charge_states"][2:-2]
            else:
                self.cStates = self.data["available_charge_states"][2:-2]
        else:
            self.cStates = self.plotting["available_charge_states"]
        
        # Set the available solution set file names
        if self.plotting["solution_set_files"] == False:
            self.fNames = ["solset_MC_iters-{}_N-{}_q-{}.csv"\
                           .format(self.optimizer["number_of_MC"],
                                   self.optimizer["number_of_ne"],
                                   cState)\
                               for cState in self.cStates]
            self.fNames = [self.g["save_to_path"] + f for f in self.fNames]
        else:
            self.fNames = self.plotting["solution_set_files"]
        
        # Set the output directory
        self.outDir = params["general"]["save_to_path"]
        
        # Set plotting limits
        if not self.plotting["Ee_lo"] == False:
            self.Ee_lo = self.plotting["Ee_lo"]
        else: 
            self.Ee_lo = self.optimizer["Ee_lo"]
        if not self.plotting["Ee_hi"] == False:
            self.Ee_hi = self.plotting["Ee_hi"]
        else:
            self.Ee_hi = self.optimizer["Ee_hi"]
        if not self.plotting["ne_lo"] == False:
            self.ne_lo = self.plotting["ne_lo"]
        else:
            self.ne_lo = self.optimizer["ne_lo"]
        if not self.plotting["ne_hi"] == False:
            self.ne_hi = self.plotting["ne_hi"]
        else:
            self.ne_hi = self.optimizer["ne_hi"]
        
        
        self.F_hi = self.plotting["F_hi"]
        self.conf = self.plotting["confidence"]
        
        
        self.sigma = self.solSettings["sigma"]
        
        # Set solution set files into a dataframe,
        # organized according to the charge state.
        self.dataFiles = pd.DataFrame({"q": self.cStates, "f": self.fNames})

    

    def set_limits(self, df):
        
        # Set the desired ne, Te, F limits on the data
        df = df[df["F"]<self.F_hi]
        df = df[(df["n"]<self.ne_hi)&(df["n"]>self.ne_lo)]
        df = df[(df["E"]<self.Ee_hi)&(df["E"]>self.Ee_lo)]
        
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
        
    
    
    def plot_heatmap_solution_set(self, q):
        '''
        Plot a heatmap of solutions for charge state q
        '''
        
        df = self.get_df(q)
        x=df["E"]
        y=df["n"]
        
        # Create figure and subplots
        
        # =====================================================================
        #         Plot the heatmap
        # =====================================================================
        
        # Select plotting range
        rng = [ [self.Ee_lo,self.Ee_hi], [self.ne_lo,self.ne_hi]]
        
        # Distribute the solution set into bins
        heatmap, xedges, yedges = np.histogram2d(x=x, y=y, 
                                                  bins=1000, 
                                                  range=rng,
                                                  density=True)
        # Apply a gaussian filter
        heatmap = gaussian_filter(heatmap, sigma=self.sigma)
        
        img = heatmap.T
        
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
        
# =============================================================================
#         Generate the figure
# =============================================================================
        fig = plt.figure()

        # Settings for the subplots to add to the figure
        
        # Main heatmap
        mainx = 0.15
        mainy = 0.275
        mainw = 0.7
        mainh = 0.6
        
        # Cumulation plot (x-projection)
        xmargx = mainx
        xmargy = mainy + mainh
        xmargw = mainw
        xmargh = 0.1
        
        # Cumulation plot (y-projection)
        ymargx = mainx + mainw
        ymargy = mainy
        ymargw = 0.1
        ymargh = mainh
        
        # Colorbar
        cbaxx = mainx
        cbaxy = mainy - 0.15
        cbaxw = mainw
        cbaxh = 0.02
        
        # Plot the heatmap
        ax1 = fig.add_axes([mainx, mainy, mainw, mainh])
        im = plt.imshow(img, origin="lower",
                  cmap = cm.Blues,
                  extent = extent,
                  aspect = "auto",
                  interpolation = "none")

        # Plot x-projection histogram
        xmarg = fig.add_axes([xmargx, xmargy, xmargw, xmargh])
        xmarg.hist(x, 250)
        xmarg.set(xscale=self.solSettings["x_scale"],
        xlim=(self.Ee_lo, self.Ee_hi))
        
        xmarg.spines["left"].set_visible(False)
        xmarg.spines["right"].set_visible(False)
        xmarg.spines["top"].set_visible(False)
        xmarg.spines["bottom"].set_visible(False)
        
        # Plot y-projection histogram
        ymarg = fig.add_axes([ymargx, ymargy, ymargw, ymargh])
        ymarg.hist(y, 250, orientation="horizontal")
        ymarg.set(yscale=self.solSettings["y_scale"],
        ylim=(self.ne_lo, self.ne_hi))
        
        ymarg.spines["top"].set_visible(False)
        ymarg.spines["bottom"].set_visible(False)
        ymarg.spines["left"].set_visible(False)
        ymarg.spines["right"].set_visible(False)
        
        
        # Plot the colorbar
        cbax = fig.add_axes([cbaxx, cbaxy, cbaxw, cbaxh])
        plt.colorbar(im, cax=cbax, orientation="horizontal")
        
        # Set axis settings
        ax1.set(xlabel=r"$\left\langle E_e\right\rangle$ (eV)",
        ylabel=r"$n_e$ (cm$^{-3}$)",
        xscale=self.solSettings["x_scale"],
        yscale=self.solSettings["y_scale"],
        xlim=(self.Ee_lo, self.Ee_hi),
        ylim=(self.ne_lo, self.ne_hi))
        
        # Marginal axis labels are unnecessary
        xmarg.set_yticklabels([])
        xmarg.set_xticklabels([])
        xmarg.set_yticks([])
        ymarg.set_yticklabels([])
        ymarg.set_xticklabels([])
        ymarg.set_xticks([])
        
        # Save figure
        plt.savefig( self.outDir + "solution_set_q-{}+.eps".format(str(q)), format="eps")
        plt.savefig( self.outDir + "solution_set_q-{}+.png".format(str(q)), format="png", dpi=300)

        
      
    
    def do_solution_set_plots(self):
        # Run the solution set plotting.
        for q in self.cStates:
            self.plot_heatmap_solution_set(q)




class Plotter:
    
    def __init__(self, params):
        
        self.plotting = params["plotting"]
        self.optimizer = params["optimizer"]
        self.data = params["data"]
        self.g = params["general"]
        self.solSettings = self.plotting["plot_solution_sets"]
        parsing = params["parsing"]
        
        # Set the charge states for which data is available
        if self.plotting["available_charge_states"] == False:
            if not self.data["available_charge_states"] == False:
                self.cStates = self.data["available_charge_states"][2:-2]
            else:
                self.cStates = parsing["available_charge_states"][2:-2]
        else:
            self.cStates = self.plotting["available_charge_states"]
        
        # Set the available solution set file names
        if self.plotting["solution_set_files"] == False:
            self.fNames = ["solset_MC_iters-{}_N-{}_q-{}.csv"\
                           .format(self.optimizer["number_of_MC"],
                                   self.optimizer["number_of_ne"],
                                   cState)\
                               for cState in self.cStates]
            self.fNames = [self.g["save_to_path"] + f for f in self.fNames]
        else:
            self.fNames = self.plotting["solution_set_files"]
        
        # Set the output directory
        self.outDir = params["general"]["save_to_path"]
        
        # Set plotting limits
        if not self.plotting["Ee_lo"] == False:
            self.Ee_lo = self.plotting["Ee_lo"]
        else: 
            self.Ee_lo = self.optimizer["Ee_lo"]
        if not self.plotting["Ee_hi"] == False:
            self.Ee_hi = self.plotting["Ee_hi"]
        else:
            self.Ee_hi = self.optimizer["Ee_hi"]
        if not self.plotting["ne_lo"] == False:
            self.ne_lo = self.plotting["ne_lo"]
        else:
            self.ne_lo = self.optimizer["ne_lo"]
        if not self.plotting["ne_hi"] == False:
            self.ne_hi = self.plotting["ne_hi"]
        else:
            self.ne_hi = self.optimizer["ne_hi"]
        
        self.F_hi = self.plotting["F_hi"]
        self.conf = self.plotting["confidence"]
        
        # Set solution set files into a dataframe,
        # organized according to the charge state.
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
        
        # Set the desired ne, Ee, F limits on the data
        # df = df[df["F"]<self.F_hi]
        df = df[(df["n"]<self.ne_hi)&(df["n"]>self.ne_lo)]
        df = df[(df["E"]<self.Ee_hi)&(df["E"]>self.Ee_lo)]
        
        return df
    
    
    def set_F_upper(self,df):
        '''
        Set the upper limit of F from the parameters.py
        '''
        
        df = df[df["F"]<self.F_hi]
        
        return df
    
    
    def plot_nSol_v_F(self, Fs, q, ax):
        '''
        Plot the number of solutions as a function of the 
        upper limit of the penalty function 
        for a given charge state q.
        '''
        
        nSols=[]
        
        for F in Fs:
            
            df = self.get_df(q)
            df = self.set_limits(df)
            df = df[df["F"]<F]
            
            nSols.append(len(df))
        
        ax.scatter(Fs, nSols, s=48, label="{}+".format(str(q)))
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$F$")
        ax.set_ylabel("Number of solutions")
        ax.set_xticks(Fs)
        
        return ax
    
    def plot_nSol_v_F_all(self):
        '''
        Plot the number of solutions vs upper limit of penalty function value F.
        Plot each charge state data in one figure.
        '''
        p = self.plotting["plot_number_of_solutions_vs_F"]
        if p["do_plot"]:
            fig, ax = plt.subplots()
            for charge_state in self.cStates:
                
                ax = self.plot_nSol_v_F(Fs=p["list_of_F_values"], 
                                         q=charge_state,
                                         ax=ax)
                
                try:
                    if not p["y_lo"] == False:
                        ax.set_ylim(bottom=p["y_lo"])
                    if not p["y_hi"] == False:
                        ax.set_ylim(bottom=p["y_hi"])
                except:
                    print("Erroneous value encountered in y-limits\
                          when plotting number of solutions vs F.")
                          
            ax.legend()
            fig.tight_layout()
            fig.savefig(self.outDir + "fig_number_of_solutions_vs_F.eps", format="eps")
            
            

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
    
    
    
    def plot_time_v_q_all(self):
        '''
        Plot the characteristic time vs charge state
        '''
        pcq = self.plotting["plot_characteristic_time_vs_charge"]
        if pcq["do_plot"] == True:
            
            charge_states = self.cStates
            
            # Plot confinement times vs F
            fig, ax = self.plot_time_against_q( qs = charge_states,
                                              key = "tau",
                                              marker = pcq["conf_marker"],
                                              color = pcq["conf_color"]
                                              )
            ax.set_ylabel("Confinement time (ms)")
            ax.set_yscale(pcq["conf_yscale"])
            
            if not pcq["conf_yscale"] == "log":
                ax.set_ylim(bottom=0)
            
            fig.tight_layout()
            fig.savefig(self.outDir + "fig_time_CONF_vs_q.eps", format="eps")
            fig.savefig(self.outDir + "fig_time_CONF_vs_q.png", format="png", dpi=300)
    
            # Plot ionization times vs F
            fig, ax = self.plot_time_against_q( qs = charge_states,
                                              key = "inz_time",
                                              marker = pcq["inz_marker"],
                                              color = pcq["inz_color"]
                                              )
            ax.set_ylabel("Ionisation time (ms)")
            ax.set_yscale(pcq["inz_yscale"])
            if not pcq["inz_yscale"] == "log":
                ax.set_ylim(bottom=0)
            fig.tight_layout()
            fig.savefig(self.outDir + "fig_time_INZ_vs_q.eps", format="eps")
            fig.savefig(self.outDir + "fig_time_INZ_vs_q.png", format="png", dpi=300)
            
            # Plot charge exchange times vs F
        
            fig, ax = self.plot_time_against_q(qs = charge_states,
                                              key = "cx_time",
                                              marker = pcq["cx_marker"],
                                              color = pcq["cx_color"]
                                              )
            ax.set_ylabel("Charge exchange time (ms)")
            ax.set_yscale(pcq["cx_yscale"])
            if not pcq["cx_yscale"] == "log":
                ax.set_ylim(bottom=0)
            fig.tight_layout()
            fig.savefig(self.outDir + "fig_time_CX_vs_q.eps", format="eps")
            fig.savefig(self.outDir + "fig_time_CX_vs_q.png", format="png", dpi=300)
    
    
    
    def plot_Ee_against_q(self):
        '''
        Plot the plasma average energy content vs q.
        '''
        qs = self.cStates
        pce = self.plotting["plot_energy_content_vs_charge"]
        lo_errs=[]
        hi_errs=[]
        medians=[]
        
        fig, ax = plt.subplots()
        for q in qs:
            
            df = self.get_df(q)
            df = self.set_limits(df)
            df = self.set_F_upper(df)
            
            data = np.array(df["n"])*np.array(df["E"]) 
                
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
                    marker=pce["marker"],
                    markersize=8,
                    color=pce["color"],
                    )
        
        # ax.set_xscale("log")
        ax.set_ylabel(r"$n_e \left\langle E_e \right\rangle$ (eV/cm$^{3}$)")
        ax.set_xlabel("Charge state")
        ax.set_xticks(qs)
        ax.set_yscale("log")
        
        if not pce["y_lo"] == False:
            ax.set_ylim(bottom=pce["y_lo"])
        if not pce["y_hi"] == False:
            ax.set_ylim(top=pce["y_hi"])
        
        fig.tight_layout()
        
        fig.savefig(self.outDir + "fig_Ee_vs_q.eps", format="eps")
        fig.savefig(self.outDir + "fig_Ee_vs_q.png", format="png", dpi=300)
    
    
    
    def plot_triple_against_q(self):
        '''
        Make sure to set F_hi (and other limits) as you wish before 
        running this plotting!
        '''
        qs = self.cStates
        pce = self.plotting["plot_triple_product_vs_charge"]
        
        lo_errs=[]
        hi_errs=[]
        medians=[]
        fig, ax = plt.subplots()
        for q in qs:
            
            df = self.get_df(q)
            df = self.set_limits(df)
            df = self.set_F_upper(df)
            
            data = np.array(df["n"])*np.array(df["E"])*np.array(df["tau"])
                
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
                    marker=pce["marker"],
                    markersize=8,
                    color=pce["color"],
                    )
        
        ax.set_ylabel(r"$n_e \left\langleE_e\right\rangle\tau^q$ (eVs/cm$^{3}$)")
        ax.set_xlabel("Charge state")
        ax.set_xticks(qs)
        ax.set_yscale("log")
        
        if not pce["y_lo"]==False:
            ax.set_ylim(bottom=pce["y_lo"])
        if not pce["y_hi"]==False:
            ax.set_ylim(top=pce["y_hi"])
        
        fig.tight_layout()
        fig.savefig(self.outDir + "fig_triple_vs_q.eps", format="eps")
        fig.savefig(self.outDir + "fig_triple_vs_q.png", format="png", dpi=300)

        
    
    
    
    def output_triple_results(self):
        
        qs = self.cStates
        
        lo_errs=[]
        hi_errs=[]
        medians=[]
        for q in qs:
            
            df = self.get_df(q)
            df = self.set_limits(df)
            df = self.set_F_upper(df)
            
            data = np.array(df["n"])*np.array(df["E"])*np.array(df["tau"])
                
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
        
        df_out.to_csv(self.outDir + "output_triple.csv")    
    
    
    def output_time_results(self, key):
        '''
        Output characteristic times vs q as .csv 
        
        N.B. all times are output in units of milliseconds!
        '''
        qs = self.cStates
        
        lo_errs=[]
        hi_errs=[]
        medians=[]
        for q in qs:
            
            df = self.get_df(q)
            df = self.set_limits(df)
            df = self.set_F_upper(df)
            
            if key == "inz_time":
                data=np.array(df["inz_rate"])
                data = data**(-1)
                data = 1e3*data
                
            elif key == "cx_time":
                data=np.array(df["cx_rate"])
                data = data**(-1)
                data = 1e3*data

            elif key == "conf_time":
                data=np.array(df["tau"])*1e3
                
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
        
        if key == "conf_time":
            df_out.to_csv(self.outDir + "output_time_CONF_vs_q.csv")
        if key == "inz_time":
            df_out.to_csv(self.outDir + "output_time_INZ_vs_q.csv")
        if key == "cx_time":
            df_out.to_csv(self.outDir + "output_time_CX_vs_q.csv")
            
        
        

    
# '''
# Output characteristic times in a .csv file.
# '''
# if pr["output_tau_vs_q"] == True:
#     P = Plotter()
#     df_out=P.output_time_results(qs=pr["charge_states"], key="tau")
#     df_out.to_csv(P.outDir + "output_time_CONF_vs_q.csv")

# if pr["output_inz_vs_q"] == True:
#     P = Plotter()
#     df_out=P.output_time_results(qs=pr["charge_states"], key="inz_time")
#     df_out.to_csv(P.outDir + "output_time_INZ_vs_q.csv")
    
# if pr["output_cx_vs_q"] == True:
#     P = Plotter()
#     df_out=P.output_time_results(qs=pr["charge_states"], key="cx_time")
#     df_out.to_csv(P.outDir + "output_time_CX_vs_q.csv")
    
    
    def output_Ee_results(self):
        '''
        Output average energy content as a function of charge state as .csv
        '''
        qs = self.cStates
        
        lo_errs=[]
        hi_errs=[]
        medians=[]
        for q in qs:
            
            df = self.get_df(q)
            df = self.set_limits(df)
            df = self.set_F_upper(df)
            
            data = np.array(df["n"])*np.array(df["E"])
                
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
        
        df_out.to_csv(self.outDir + "output_Ee.csv")
        


    def do_result_plots(self):
        
        if self.plotting["plot_number_of_solutions_vs_F"]["do_plot"]:
            self.plot_nSol_v_F_all()
        if self.plotting["plot_characteristic_time_vs_charge"]["do_plot"]:
            self.plot_time_v_q_all()
        if self.plotting["plot_energy_content_vs_charge"]["do_plot"]:
            self.plot_Ee_against_q()
        if self.plotting["plot_triple_product_vs_charge"]["do_plot"]:
            self.plot_triple_against_q()


        if self.plotting["plot_characteristic_time_vs_charge"]["output_data"]:
            self.output_time_results(key="conf_time")
            self.output_time_results(key="inz_time")
            self.output_time_results(key="cx_time")
        if self.plotting["plot_energy_content_vs_charge"]["output_data"]:
            self.output_Ee_results()
        if self.plotting["plot_triple_product_vs_charge"]["output_data"]:
            self.output_triple_results()



# #
# # Plot the characteristic time vs upper limit of penalty function value F
# #
# pcf = pr["plotting_time_vs_F"]
# if pcf["plot_or_not"] == True:
#     P = Plotter()
#     for charge_state in pr["charge_states"]:
        
#         # Plot confinement times vs F
#         if pcf["plot_tau"] == True:
#             fig, ax = P.plot_times_against_F(Fs = pcf["list_of_Fs"],
#                                               q = charge_state,
#                                               key = "tau",
#                                               marker = pcf["tau_marker"],
#                                               color = pcf["tau_color"]
#                                               )
#             ax.set_ylabel("Confinement time (ms)")
#             fig.tight_layout()
#             fig.savefig(P.outDir + "fig_CONF_vs_F_q={}.eps".format(charge_state), format="eps")
#             fig.savefig(P.outDir + "fig_CONF_vs_F_q={}.png".format(charge_state), format="png", dpi=300)
#             plt.close(fig)
            
#         # Plot ionization times vs F
#         if pcf["plot_inz"] == True:
#             fig, ax = P.plot_times_against_F(Fs = pcf["list_of_Fs"],
#                                               q = charge_state,
#                                               key = "inz_time",
#                                               marker = pcf["inz_marker"],
#                                               color = pcf["inz_color"]
#                                               )
#             ax.set_ylabel("Ionisation time (ms)")
#             fig.tight_layout()
#             fig.savefig(P.outDir + "fig_INZ_vs_F_q={}.eps".format(charge_state), format="eps")
#             fig.savefig(P.outDir + "fig_INZ_vs_F_q={}.png".format(charge_state), format="png", dpi=300)
#             plt.close(fig)
            
#         # Plot charge exchange times vs F
#         if pcf["plot_cx"] == True:
#             fig, ax = P.plot_times_against_F(Fs = pcf["list_of_Fs"],
#                                               q = charge_state,
#                                               key = "cx_time",
#                                               marker = pcf["cx_marker"],
#                                               color = pcf["cx_color"]
#                                               )
#             ax.set_ylabel("Charge exchange time (ms)")
#             fig.tight_layout()
#             fig.savefig(P.outDir + "fig_CX_F_q={}.eps".format(charge_state), format="eps")
#             fig.savefig(P.outDir + "fig_CX_F_q={}.png".format(charge_state), format="png")
#             plt.close(fig)



        



    
    

    

    
# '''
# Output average energy content in a .csv file
# '''
# if pr["output_Ee_vs_q"] == True:
#     P = Plotter()
#     df_out=P.output_Ee_results(qs=pr["charge_states"])
#     df_out.to_csv(P.outDir + "output_Ee_vs_q.csv")
    
    
    
# '''
# Output triple products in a .csv file.
# '''
# if pr["output_triple_vs_q"] == True:
#     P = Plotter()
#     df_out=P.output_triple_results(qs=pr["charge_states"])
#     df_out.to_csv(P.outDir + "output_triple_vs_q.csv")
    
