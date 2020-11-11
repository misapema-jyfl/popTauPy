#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 09:52:53 2020

Just give the paths to files in the fNames list, and the charge states 
corresponding to said files in the list cStates. Arrange in numerical order! 

You may also change the penalty function upper limit 
through the variable penaltyLim.

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


penaltyLim = 1e-4 # The penalty function value upper limit


def find_confidence_interval(list_of_values, condition_percentage):
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





def plot_distribution(list_of_values, confidence_interval, bins,
                      xlabel, ylabel):
    
    x = list_of_values
    p = confidence_interval
    
    fig, ax = plt.subplots()
    
    x_lo, median, x_hi = find_confidence_interval(x, p)
    
    if x_lo < 0:
        print("Warning! Negative value for x_lo when plotting distribution of {:s}.".format(xlabel))
    
    ax.hist(x, bins=bins, density=False)
    ax.axvline(x_lo, color="r", ls="--", lw=2)
    ax.axvline(x_hi, color="r", ls="--", lw=2)
    ax.axvline(median, color="k", ls="--", lw=2)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    
    return fig, ax



'''
Run code below!
'''



# Choose the confidence interval for plots
# 0.341 corresponds to 1 sigma confidence.
conf = 0.341







'''
Output characteristic times, energy contents and triple products to .csv file.
'''

for q,fName in zip(cStates,fNames):
    
    df = pd.read_csv(fName)
    
    # Make sure that no unphysical solutions are counted
    c = df["tau"] > 0
    df = df[c]
    c = df["cx_rate"] > 0
    df = df[c]
    c = df["inz_rate"] > 0
    df = df[c]
    
    
    df = df[df["F"] < penaltyLim] # Set upper limit of the penalty function
    
    ns = df["n"]
    Ts = df["T"]
    
    taus = df["tau"]
    taus = np.array( [1e3*t for t in taus] ) # Convert to ms
    
    inz_rates = df["inz_rate"]
    inz_times = np.array( [1e3/s for s in inz_rates] ) # Conversion from rate to time (ms)
    
    cx_rates = df["cx_rate"]
    cx_times = np.array( [1e3/s for s in cx_rates] ) # Conversion from rate to time (ms)
    
    eCs = df["energy_content"]
    eCs = np.array([e*1.5 for e in eCs]) # Convert to average energy
    
    triprods = np.array([e*t*1e-3 for e, t in zip(eCs, taus)]) # N.B. the conversion to seconds!
    
    ntaus = np.array([n*t*1e-3 for n,t in zip(ns, taus)]) # N.B. conversion to s
    
    
    
    # Set up the output file data frame
    df_out = pd.DataFrame(index=["lo_err", "median", "hi_err"])
    
    
    
    # taus
    x_lo, median, x_hi = find_confidence_interval(taus, conf)
    lo_err = median - x_lo
    hi_err = x_hi - median
    df_out["taus"] = [lo_err, median, hi_err]
    
    
    
    # inz_times
    x_lo, median, x_hi = find_confidence_interval(inz_times, conf)
    lo_err = median - x_lo
    hi_err = x_hi - median
    
    df_out["inz_time"] = [lo_err, median, hi_err]
    
    
    
    # cx_times
    x_lo, median, x_hi = find_confidence_interval(cx_times, conf)
    lo_err = median - x_lo
    hi_err = x_hi - median
    
    df_out["cx_time"] = [lo_err, median, hi_err]
    
    
    
    # energy contents
    x_lo, median, x_hi = find_confidence_interval(eCs, conf)
    lo_err = median - x_lo
    hi_err = x_hi - median
    
    df_out["energy_content"] = [lo_err, median, hi_err]
    
    
    
    # triple products
    x_lo, median, x_hi = find_confidence_interval(triprods, conf)
    lo_err = median - x_lo
    hi_err = x_hi - median
    
    df_out["triple_product"] = [lo_err, median, hi_err]
    
    
    
    # ne*tau values
    x_lo, median, x_hi = find_confidence_interval(ntaus, conf)
    lo_err = median - x_lo
    hi_err = x_hi - median
    
    df_out["neTaus"] = [lo_err, median, hi_err]
    
    
    
    df_out.to_csv("./results/outfile_results_q={}.csv".format(q))







'''
Solution set plot
'''
fig, ax = plt.subplots()
for q, fName in zip(cStates, fNames):
    
    df = pd.read_csv(fName)
    
    # Make sure that no unphysical solutions are counted
    c = df["tau"] > 0
    df = df[c]
    c = df["cx_rate"] > 0
    df = df[c]
    c = df["inz_rate"] > 0
    df = df[c]
    
    df = df[df["F"] < penaltyLim] # Set upper limit of the penalty function
    
    ns = df["n"]
    Ts = df["T"]
    
    lbl = "{" + str(q) + "+}"
    
    ax.scatter(Ts, ns, label="K$^{}$".format(lbl))
ax.set_xlabel(r"$T_e$ (eV)")
ax.set_ylabel(r"$n_e$ (cm$^{-3}$)")
ax.set_yscale("log")
ax.set_xscale("log")
ax.legend()
plt.tight_layout()
plt.savefig("./results/fig_solution_set.png", dpi=300)
plt.close()
    
    
    
    
    

'''
Plot of resultant characteristic time distributions
'''    

for q,fName in zip(cStates,fNames):
    
    df = pd.read_csv(fName)
    
    # Make sure that no unphysical solutions are counted
    c = df["tau"] > 0
    df = df[c]
    c = df["cx_rate"] > 0
    df = df[c]
    c = df["inz_rate"] > 0
    df = df[c]
    
    df = df[df["F"] < penaltyLim] # Set upper limit of the penalty function
    
    ns = df["n"]
    Ts = df["T"]
    
    taus = df["tau"]
    taus = np.array( [1e3*t for t in taus] ) # Convert to ms
    
    inz_rates = df["inz_rate"]
    inz_times = np.array( [1e3/s for s in inz_rates] ) 
    
    cx_rates = df["cx_rate"]
    cx_times = np.array( [1e3/s for s in cx_rates] ) 
    
    
    u = np.floor(max(taus))+1
    l = np.floor(min(taus))
    d = int(np.floor(u-l) + 1)
    bns = np.linspace(l,u,d)
    fig, ax = plot_distribution(taus, conf, bins=bns, 
                                xlabel=r"$\tau^q$ (ms)", ylabel="Number of solutions")
    ax.set_xlim(0,100)
    plt.tight_layout()
    plt.savefig("./results/fig_distribution_conf-time_q={}.png".format(q), dpi=300)
    plt.close(fig)
    
    
    u = np.floor(max(inz_times))+1
    l = np.floor(min(inz_times))
    d = int(u-l)
    bns = np.linspace(l,u,d)
    fig, ax = plot_distribution(inz_times, conf, bins=bns, 
                                xlabel=r"$1/n_e\left\langle\sigma v\right\rangle_{q\to q+1}^{{inz}}$", 
                                ylabel="Number of solutions")
    ax.set_xlim(0,40)
    plt.tight_layout()
    plt.savefig("./results/fig_distribution_inz-time_q={}.png".format(q), dpi=300)
    plt.close(fig)
    
    
    u = np.floor(max(cx_times))+1
    l = np.floor(min(cx_times))
    d = int(u-l)
    bns = np.linspace(l,u,d)
    fig, ax = plot_distribution(cx_times, conf, bins=bns, 
                                xlabel=r"$1/n_0\left\langle\sigma v\right\rangle_{q\to q-1}^{{cx}}$", 
                                ylabel="Number of solutions")
    ax.set_xlim(right=40)
    plt.tight_layout()
    plt.savefig("./results/fig_distribution_cx-time_q={}.png".format(q), dpi=300)
    plt.close(fig)
    
    


    


    
    
'''
Plasma characteristic times plot.

Plot all in one figure, and then each separately.
'''
 
M = 8 # <- markersize
   
#
# CX-time
#
err_los = []
err_his = []
medians = []
for q, fName in zip(cStates, fNames):
    
    df = pd.read_csv(fName)
    
    # Make sure that no unphysical solutions are counted
    c = df["tau"] > 0
    df = df[c]
    c = df["cx_rate"] > 0
    df = df[c]
    c = df["inz_rate"] > 0
    df = df[c]
    
    df = df[df["F"] < penaltyLim] # Set upper limit of the penalty function
    
    cx_rates = df["cx_rate"]
    cx_times = np.array( [1e3/s for s in cx_rates] ) 
    
    x_lo, median, x_hi = find_confidence_interval(cx_times, conf)
    
    err_los.append(median - x_lo)
    err_his.append(x_hi - median)
    medians.append(median)



# CX times in their own figure.
fig2,ax2=plt.subplots()
x = np.array(cStates)
ax2.errorbar(x, medians, yerr=[err_los, err_his], fmt="^", barsabove=True,
            capsize=7, color = "b", markersize=M,
            label = r"$1/n_0\left\langle\sigma v\right\rangle^{cx}_{q\to q-1}$")
ax2.set_ylim(bottom=1, top=1e3)
ax2.set_xticks(cStates)
ax2.set_xlabel("Charge state")
ax2.set_ylabel("Milliseconds")
ax2.legend()
ax2.set_yscale("log")
fig2.tight_layout()
fig2.savefig("./results/fig_times_CX.png", dpi=300)
plt.close(fig2)


# Figure object where to plat all plots in one.
fig,ax=plt.subplots()

#
# Offset cStates by -0.1 to avoid overlap in plot
#
x = np.array(cStates) - 0.05

ax.errorbar(x, medians, yerr=[err_los, err_his], fmt="^", barsabove=True,
            capsize=7, color = "b", markersize=M,
            label = r"$1/n_0\left\langle\sigma v\right\rangle^{cx}_{q\to q-1}$")

#
# INZ-time
#
err_los = []
err_his = []
medians = []
for q, fName in zip(cStates, fNames):
    
    df = pd.read_csv(fName)
    
    df = df[df["F"] < penaltyLim] # Set upper limit of the penalty function
    
    inz_rates = df["inz_rate"]
    inz_times = np.array( [1e3/s for s in inz_rates] ) 
    
    x_lo, median, x_hi = find_confidence_interval(inz_times, conf)
    
    err_los.append(median - x_lo)
    err_his.append(x_hi - median)
    medians.append(median)

#
# Offset cStates by 0.1 to avoid overlap in plot
#
x = np.array(cStates) + 0.05

ax.errorbar(x, medians, yerr=[err_los, err_his], fmt="s", barsabove=True,
            capsize=7, color = "k", markersize=M,
            label=r"$1/n_e\left\langle\sigma v\right\rangle^{inz}_{q\to q+1}$")


# INZ times in their own figure.
fig2,ax2=plt.subplots()
x = np.array(cStates)
ax2.errorbar(x, medians, yerr=[err_los, err_his], fmt="s", barsabove=True,
            capsize=7, color = "k", markersize=M,
            label = r"$1/n_e\left\langle\sigma v\right\rangle^{inz}_{q\to q+1}$")
ax2.set_ylim(bottom=0)
ax2.set_xticks(cStates)
ax2.set_xlabel("Charge state")
ax2.set_ylabel("Milliseconds")
ax2.legend()
fig2.tight_layout()
fig2.savefig("./results/fig_times_INZ.png", dpi=300)
plt.close(fig2)


#
# Confinement time
#
err_los = []
err_his = []
medians = []
for q, fName in zip(cStates, fNames):
    
    df = pd.read_csv(fName)
    
    # Make sure that no unphysical solutions are counted
    c = df["tau"] > 0
    df = df[c]
    c = df["cx_rate"] > 0
    df = df[c]
    c = df["inz_rate"] > 0
    df = df[c]
    
    df = df[df["F"] < penaltyLim] # Set upper limit of the penalty function
    
    taus = df["tau"]
    taus = np.array( [1e3*t for t in taus] )
    
    x_lo, median, x_hi = find_confidence_interval(taus, conf)
    
    err_los.append(median - x_lo)
    err_his.append(x_hi - median)
    medians.append(median)
    
ax.errorbar(cStates, medians, yerr=[err_los, err_his], fmt="o", barsabove=True,
            capsize=7, color = "r", markersize=M,
            label=r"$\tau^q$")

ax.set_ylim(bottom=-1)
ax.set_xticks(cStates)
# ax.set_xlim(left=7.5, right=10.5)                                             # x-limits for the plot
ax.set_xlabel("Charge state")
ax.set_ylabel("Milliseconds")
ax.legend()
plt.tight_layout()
plt.savefig("./results/fig_times_characteristic.png", dpi=300)
plt.close(fig)


# Conf. times in their own figure.
fig2,ax2=plt.subplots()
x = np.array(cStates)
ax2.errorbar(x, medians, yerr=[err_los, err_his], fmt="o", barsabove=True,
            capsize=7, color = "r", markersize=M,
            label = r"$\tau^q$")
ax2.set_ylim(bottom=0)
ax2.set_xticks(cStates)
ax2.set_xlabel("Charge state")
ax2.set_ylabel("Milliseconds")
ax2.legend()
fig2.tight_layout()
fig2.savefig("./results/fig_times_CONF.png", dpi=300)
plt.close(fig2)







'''
Energy content plot
'''



err_los = []
err_his = []
medians = []
for q, fName in zip(cStates, fNames):
    
    df = pd.read_csv(fName)
    
    # Make sure that no unphysical solutions are counted
    c = df["tau"] > 0
    df = df[c]
    c = df["cx_rate"] > 0
    df = df[c]
    c = df["inz_rate"] > 0
    df = df[c]
    
    df = df[df["F"] < penaltyLim] # Set upper limit of the penalty function
    
    eCs = df["energy_content"]
    eCs = np.array([e*1.5 for e in eCs]) # Convert to average energy
    
    x_lo, median, x_hi = find_confidence_interval(eCs, conf)
    
    err_los.append(median - x_lo)
    err_his.append(x_hi - median)
    medians.append(median)


M = 8 # <- markersize
fig,ax=plt.subplots()
ax.errorbar(cStates, medians, yerr=[err_los, err_his], fmt="s", barsabove=True,
            capsize=3, color = "k", markersize=M)
ax.set_xlabel("Charge state")
ax.set_ylabel(r"$n_e\left\langleE_e\right\rangle$ (eV/cm$^3$)")
ax.set_yscale("log")
ax.set_ylim(bottom=1e12)
plt.tight_layout()
plt.savefig("./results/fig_energy-contents.png", dpi=300)
plt.close()






'''
Triple product plot
'''

err_los = []
err_his = []
medians = []
for q, fName in zip(cStates, fNames):
    
    df = pd.read_csv(fName)
    
    # Make sure that no unphysical solutions are counted
    c = df["tau"] > 0
    df = df[c]
    c = df["cx_rate"] > 0
    df = df[c]
    c = df["inz_rate"] > 0
    df = df[c]
    
    df = df[df["F"] < penaltyLim] # Set upper limit of the penalty function
    
    eCs = df["energy_content"]
    eCs = np.array([e*1.5 for e in eCs]) # Convert to average energy

    taus = df["tau"]
    triprods = np.array([e*t for e, t in zip(eCs, taus)])
    
    x_lo, median, x_hi = find_confidence_interval(triprods, conf)
    
    err_los.append(median - x_lo)
    err_his.append(x_hi - median)
    medians.append(median)


M = 8 # <- markersize
fig,ax=plt.subplots()
ax.errorbar(cStates, medians, yerr=[err_los, err_his], fmt="s", barsabove=True,
            capsize=3, color = "m", markersize=M)
ax.set_xlabel("Charge state")
ax.set_ylabel(r"$n_e\left\langle E_e\right\rangle\tau^q$ (eVs/cm$^3$)")
ax.set_yscale("log")
plt.tight_layout()
plt.savefig("./results/fig_triple-products.png", dpi=300)
plt.close()